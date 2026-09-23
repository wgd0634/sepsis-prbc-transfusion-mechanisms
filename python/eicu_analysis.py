# -*- coding: utf-8 -*-
"""
eICU 外部验证分析流水线
用法: python eICU分析代码.py <提取结果.csv>
输入列名见交接说明第五节。分析流程与 MIMIC 主分析完全一致：
IPTW -> SMD平衡 -> 加权结局比较 -> 三路径中介分析(含多重插补) -> EF亚组交互
依赖: pandas numpy statsmodels scipy scikit-learn
"""
import sys
import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats

CSV = sys.argv[1] if len(sys.argv) > 1 else 'eicu_cohort.csv'
B = 500          # bootstrap 次数
rng = np.random.default_rng(42)

d = pd.read_csv(CSV)
d['male'] = (d['gender'] == 'Male').astype(int)
print(f"队列 n={len(d)}, 输血组 n={int(d.prbc_24h_flag.sum())} ({d.prbc_24h_flag.mean()*100:.1f}%)")

CONF = ['age', 'hb_min_24h', 'apachescore', 'male',
        'ckd', 'diabetes', 'mech_vent_24h', 'vasopressor_24h', 'rrt_24h']

# ---------- 1. 倾向评分 + IPTW ----------
d_ps = d.dropna(subset=CONF + ['prbc_24h_flag']).copy()
X = sm.add_constant(d_ps[CONF])
psm = sm.Logit(d_ps.prbc_24h_flag, X).fit(disp=0)
d_ps['ps'] = psm.predict(X)
w = np.where(d_ps.prbc_24h_flag == 1, 1/d_ps.ps, 1/(1-d_ps.ps))
lo, hi = np.quantile(w, [0.01, 0.99])
d_ps['w'] = w.clip(lo, hi)
print(f"\nPS范围 {d_ps.ps.min():.3f}-{d_ps.ps.max():.3f}, 权重截断 {lo:.2f}-{hi:.2f}")

def smd(col, weighted=False):
    t_ = d_ps[d_ps.prbc_24h_flag == 1]; c_ = d_ps[d_ps.prbc_24h_flag == 0]
    if not weighted:
        return (t_[col].mean() - c_[col].mean()) / np.sqrt((t_[col].var() + c_[col].var()) / 2)
    mt = np.average(t_[col], weights=t_.w); mc = np.average(c_[col], weights=c_.w)
    vt = np.average((t_[col]-mt)**2, weights=t_.w); vc = np.average((c_[col]-mc)**2, weights=c_.w)
    return (mt - mc) / np.sqrt((vt + vc) / 2)

print("\n协变量平衡 (加权前 -> 加权后 SMD):")
for c in CONF:
    print(f"  {c:18s} {smd(c):+.3f} -> {smd(c, True):+.3f}")

# ---------- 2. 加权结局比较 + bootstrap ----------
def weighted_effects(data):
    res = {}
    for out in ['death_inhosp', 'aki_7d']:
        t_ = data[data.prbc_24h_flag == 1]; c_ = data[data.prbc_24h_flag == 0]
        res[out] = (np.average(t_[out], weights=t_.w) - np.average(c_[out], weights=c_.w)) * 100
    return res

base = weighted_effects(d_ps)
boot = {k: [] for k in base}
n = len(d_ps)
for _ in range(B):
    b = d_ps.iloc[rng.integers(0, n, n)]
    try:
        m = sm.Logit(b.prbc_24h_flag, sm.add_constant(b[CONF])).fit(disp=0)
        ps = m.predict(sm.add_constant(b[CONF]))
        b = b.assign(w=np.where(b.prbc_24h_flag == 1, 1/ps, 1/(1-ps)).clip(lo, hi))
        r = weighted_effects(b)
        for k in boot: boot[k].append(r[k])
    except Exception:
        continue
print("\nIPTW加权效应 (输血-不输血, %):")
for k, v in boot.items():
    ci = np.percentile(v, [2.5, 97.5])
    print(f"  {k}: {base[k]:+.2f} (95%CI {ci[0]:+.2f} ~ {ci[1]:+.2f})")

# ---------- 3. 中介分析 ----------
def mediation(data, med, outcome):
    w_, T = data.w.values, data.prbc_24h_flag.values.astype(float)
    Mv, Y = data[med].values.astype(float), data[outcome].values.astype(float)
    C = sm.add_constant(data[CONF].values); Cm = C.mean(axis=0)
    a = sm.WLS(Mv, np.column_stack([np.ones(len(T)), T, C]), weights=w_).fit().params
    Xy = np.column_stack([np.ones(len(T)), T, Mv, C])
    b = sm.GLM(Y, Xy, family=sm.families.Binomial(), freq_weights=w_).fit().params
    M0 = a[0] + a[2:] @ Cm; M1 = M0 + a[1]
    risk = lambda Tv, M_: 1/(1 + np.exp(-(b[0] + b[1]*Tv + b[2]*M_ + b[3:] @ Cm)))
    return risk(1,M1)-risk(0,M0), risk(1,M0)-risk(0,M0), risk(1,M1)-risk(1,M0)  # TE, NDE, NIE

def run_mediation(sub, med, out, label, n_boot=300):
    TE, NDE, NIE = mediation(sub, med, out)
    nies = []
    nn = len(sub)
    for _ in range(n_boot):
        bb = sub.iloc[rng.integers(0, nn, nn)]
        try:
            m = sm.Logit(bb.prbc_24h_flag, sm.add_constant(bb[CONF])).fit(disp=0)
            ps = m.predict(sm.add_constant(bb[CONF]))
            bb = bb.assign(w=np.where(bb.prbc_24h_flag == 1, 1/ps, 1/(1-ps)).clip(lo, hi))
            nies.append(mediation(bb, med, out)[2])
        except Exception:
            continue
    ci = np.percentile(nies, [2.5, 97.5])
    print(f"\n中介[{label}] -> {out}: TE={TE*100:+.2f}% NDE={NDE*100:+.2f}% "
          f"NIE={NIE*100:+.2f}% (95%CI {ci[0]*100:+.2f}~{ci[1]*100:+.2f}) 中介比例={NIE/TE*100:.1f}%")

d_ps['lac0'] = d_ps.lactate_0h.clip(upper=20); d_ps['lac24'] = d_ps.lactate_24h.clip(upper=20)
d_ps['lac_clear'] = (d_ps.lac0 - d_ps.lac24) / d_ps.lac0
d_ps['ne0'] = d_ps.norepi_0h.clip(upper=1.0); d_ps['ne24'] = d_ps.norepi_24h.clip(upper=1.0)
d_ps['ne_delta'] = d_ps.ne24 - d_ps.ne0

sub_fluid = d_ps.dropna(subset=['fluid_balance_ml_24h'])
sub_fluid = sub_fluid.assign(fluid_L=sub_fluid.fluid_balance_ml_24h/1000)
run_mediation(sub_fluid, 'fluid_L', 'death_inhosp', '24h液体平衡')
sub_lac = d_ps.dropna(subset=['lac_clear'])
print(f"\n(乳酸完整案例 n={len(sub_lac)})")
if len(sub_lac) > 500: run_mediation(sub_lac, 'lac_clear', 'death_inhosp', '乳酸清除率-完整案例')
sub_ne = d_ps.dropna(subset=['ne_delta'])
print(f"(升压药完整案例 n={len(sub_ne)})")
if len(sub_ne) > 300: run_mediation(sub_ne, 'ne_delta', 'death_inhosp', '升压药减量')

# 乳酸中介 - 多重插补版
from sklearn.experimental import enable_iterative_imputer  # noqa
from sklearn.impute import IterativeImputer
from sklearn.linear_model import BayesianRidge
imp_cols = CONF + ['prbc_24h_flag', 'death_inhosp', 'lac0', 'lac24']
TEs, NIEs = [], []
for m_ in range(5):
    imp = IterativeImputer(estimator=BayesianRidge(), max_iter=10,
                           random_state=100+m_, min_value=0, sample_posterior=True)
    Xi = imp.fit_transform(d_ps[imp_cols])
    t = d_ps.copy(); t[imp_cols] = Xi
    t['lac0'] = t.lac0.clip(0.3, 20); t['lac24'] = t.lac24.clip(0.3, 20)
    t['lac_clear'] = (t.lac0 - t.lac24)/t.lac0
    t['ps'] = psm.predict(sm.add_constant(t[CONF]))
    t['w'] = np.where(t.prbc_24h_flag == 1, 1/t.ps, 1/(1-t.ps)).clip(lo, hi)
    r = mediation(t, 'lac_clear', 'death_inhosp')
    TEs.append(r[0]); NIEs.append(r[2])
print(f"\n乳酸中介[多重插补合并] -> death_inhosp: TE={np.mean(TEs)*100:+.2f}% "
      f"NIE={np.mean(NIEs)*100:+.2f}% 中介比例={np.mean(NIEs)/np.mean(TEs)*100:.1f}%")

# ---------- 4. EF亚组交互 ----------
def ef_grp(r):
    if r.sys_hf == 1: return 'HFrEF'
    if r.dia_hf == 1: return 'HFpEF'
    if r.any_hf == 1: return 'OtherHF'
    return 'NoHF'
d_ps['ef_group'] = d_ps.apply(ef_grp, axis=1)
print("\nEF亚组加权风险差 (death_inhosp):")
for g in ['HFrEF', 'HFpEF', 'OtherHF', 'NoHF']:
    s = d_ps[d_ps.ef_group == g]
    if len(s) < 50:
        print(f"  {g}: n={len(s)} (样本过少, 跳过)"); continue
    t_ = s[s.prbc_24h_flag == 1]; c_ = s[s.prbc_24h_flag == 0]
    rd = (np.average(t_.death_inhosp, weights=t_.w) - np.average(c_.death_inhosp, weights=c_.w)) * 100
    print(f"  {g:8s} n={len(s):5d}  RD={rd:+.2f}%")
