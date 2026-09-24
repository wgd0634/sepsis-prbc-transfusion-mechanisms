# -*- coding: utf-8 -*-
"""
竞争风险敏感性分析：7天AKI（死亡为竞争事件）
加权 Aalen-Johansen CIF 差异（IPTW 权重来自原 PS 模型）+ bootstrap CI
队列: MIMIC 重建队列 (mimic_cr.csv) + eICU 验证队列 (eicu_cr.csv)
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(sys.executable).parent.parent.parent))
import numpy as np
import pandas as pd
from scipy import stats as _st

rng = np.random.default_rng(42)
B = 500

def weighted_cif_rd(df, t_eval=7.0):
    """每组: 加权 Aalen-Johansen CIF1(AKI); 返回 RD% (PRBC - NoPRBC)"""
    out = {}
    for g in (0, 1):
        d = df[df.prbc_24h_flag == g]
        w = d.w.values.astype(float)
        # 事件时间(天): AKI onset 或死亡; 无事件 -> 7 天删失
        t_aki = d.t_aki.values.astype(float)
        t_death = d.t_death.values.astype(float)
        t = np.where(np.isnan(t_aki), t_death, np.where(np.isnan(t_death), t_aki,
                     np.minimum(t_aki, t_death)))
        t = np.where(np.isnan(t), t_eval, t)
        # 事件类型: 0 删失, 1 AKI, 2 死亡; 先发生者优先
        ev = np.zeros(len(d), dtype=int)
        ev[~np.isnan(t_aki) & (t_aki <= t_eval) & (np.isnan(t_death) | (t_aki <= t_death))] = 1
        ev[~np.isnan(t_death) & (t_death <= t_eval) & (np.isnan(t_aki) | (t_death < t_aki))] = 2
        t = np.clip(t, 0, t_eval)
        order = np.argsort(t, kind='mergesort')
        t_s, ev_s, w_s = t[order], ev[order], w[order]
        # 加权 Nelson-Aalen 各原因累积风险 + 加权 KM 生存
        uniq = np.unique(t_s)
        H1 = H2 = 0.0
        S, prev_S = 1.0, 1.0
        cif1 = 0.0
        for u in uniq:
            at = (t_s == u)
            n_risk = w_s[t_s >= u].sum()
            if n_risk <= 0:
                continue
            d1 = w_s[at & (ev_s == 1)].sum()
            d2 = w_s[at & (ev_s == 2)].sum()
            dc = w_s[at & (ev_s == 0)].sum()
            d_all = d1 + d2
            cif1 += prev_S * d1 / n_risk          # Aalen-Johansen
            H1 += d1 / n_risk
            H2 += d2 / n_risk
            if d_all > 0:
                S *= (1 - d_all / n_risk)
            prev_S = S
        out[g] = cif1 * 100
    return out[1] - out[0]

def load_mimic():
    d = pd.read_csv('mimic_cr.csv')
    d['t_aki'] = (pd.to_datetime(d.aki_onset_t) - pd.to_datetime(d.baseline_t)).dt.total_seconds() / 86400
    d['t_death'] = (pd.to_datetime(d.death_t) - pd.to_datetime(d.baseline_t)).dt.total_seconds() / 86400
    return d

def load_eicu():
    d = pd.read_csv('eicu_cr.csv')
    d['t_aki'] = (d.aki_onset_min - d.baseline_min) / 1440.0
    d['t_death'] = (d.death_min - d.baseline_min) / 1440.0
    return d

def add_weights(df, conf_map, cohort):
    """PS 模型权重（与原分析一致）"""
    need = ['stay_id', 'prbc_24h_flag'] + [c for c in conf_map if c != 'male']
    if 'male' in conf_map and 'male' not in pd.read_csv(cohort, nrows=0).columns:
        need.append('gender')
    d = df.merge(pd.read_csv(cohort)[need],
                 on='stay_id', how='inner', suffixes=('', '_y'))
    assert d.prbc_24h_flag_x.equals(d.prbc_24h_flag_y) if 'prbc_24h_flag_x' in d else True
    if 'prbc_24h_flag_x' in d.columns:
        d = d.drop(columns=['prbc_24h_flag_y']).rename(columns={'prbc_24h_flag_x': 'prbc_24h_flag'})
    try:
        import statsmodels.api as sm
    except ImportError:
        from statsmodels.discrete.discrete_model import Logit
        from statsmodels.tools.tools import add_constant
        import types
        sm = types.SimpleNamespace(Logit=Logit, add_constant=add_constant)
    if 'gender' in conf_map:
        d['male'] = (d.gender == 'M').astype(int)
        conf_map = ['male' if c == 'gender' else c for c in conf_map]
    if 'male' in conf_map and 'male' not in d.columns:
        d['male'] = (d.gender == 'Male').astype(int)
    d = d.dropna(subset=conf_map + ['prbc_24h_flag']).copy()
    X = sm.add_constant(d[conf_map])
    m = sm.Logit(d.prbc_24h_flag, X).fit(disp=0)
    ps = m.predict(X)
    w = np.where(d.prbc_24h_flag == 1, 1 / ps, 1 / (1 - ps))
    lo, hi = np.quantile(w, [0.01, 0.99])
    d['w'] = np.clip(w, lo, hi)
    return d

def run(name, df, conf):
    d = add_weights(df, conf, {'MIMIC': 'mimic_cohort.csv', 'eICU': 'eicu_cohort.csv'}[name])
    rd = weighted_cif_rd(d)
    n = len(d)
    boots = []
    for _ in range(B):
        b = d.iloc[rng.integers(0, n, n)]
        try:
            boots.append(weighted_cif_rd(b))
        except Exception:
            continue
    ci = np.percentile(boots, [2.5, 97.5])
    # 未加权 CIF（描述）
    d0 = d.assign(w=1.0)
    rd_unw = weighted_cif_rd(d0)
    print(f"[{name}] n={n}, 输血={int(d.prbc_24h_flag.sum())} "
          f"加权CIF RD(AKI,7d) = {rd:+.2f}% (95%CI {ci[0]:+.2f} ~ {ci[1]:+.2f}); 未加权 {rd_unw:+.2f}%")
    return rd, ci

if __name__ == '__main__':
    print("== MIMIC 重建队列 ==")
    run('MIMIC', load_mimic(), ['age', 'hb_min_24h', 'sofa', 'gender', 'charlson', 'ckd', 'diabetes', 'mech_vent_24h', 'vasopressor_24h', 'rrt_24h'])
    print("== eICU 验证队列 ==")
    run('eICU', load_eicu(), ['age', 'hb_min_24h', 'apachescore', 'male', 'ckd', 'diabetes', 'mech_vent_24h', 'vasopressor_24h', 'rrt_24h'])
