# -*- coding: utf-8 -*-
"""
MIMIC 脓毒症输血队列 · EF 亚组交互再分析（超声 LVEF 版）
输入: mimic_cohort.csv（本地重建队列 + ICD 心衰标志 + 超声 LVEF）
输出: IPTW 平衡 / 主效应 / ICD 亚组(复现论文) / LVEF 亚组(新) / 交互检验
与论文一致: 加权 logistic 交互模型（处理×亚组乘积项）, 500 次 bootstrap
"""
import sys
import numpy as np
import pandas as pd
try:
    import statsmodels.api as sm
except ImportError:
    from statsmodels.tools.tools import add_constant
    from statsmodels.discrete.discrete_model import Logit
    from statsmodels.genmod.generalized_linear_model import GLM
    from statsmodels.genmod import families
    import types
    sm = types.SimpleNamespace(Logit=Logit, GLM=GLM, families=families, add_constant=add_constant)

CSV = sys.argv[1] if len(sys.argv) > 1 else 'mimic_cohort.csv'
B = 500
rng = np.random.default_rng(42)

d = pd.read_csv(CSV)
d['male'] = (d['gender'] == 'M').astype(int)
print(f"队列 n={len(d)}, 输血组 n={int(d.prbc_24h_flag.sum())} ({d.prbc_24h_flag.mean()*100:.1f}%)")
print(f"有超声 LVEF 者 n={d.lvef.notna().sum()} ({d.lvef.notna().mean()*100:.1f}%)")

# 与论文 Table1 对照的 ICD 亚组分布
def ef_icd(r):
    if r.sys_hf == 1: return 'HFrEF(ICD)'
    if r.dia_hf == 1: return 'HFpEF(ICD)'
    if r.any_hf == 1: return 'OtherHF(ICD)'
    return 'NoHF'
d['ef_icd'] = d.apply(ef_icd, axis=1)
print("\nICD 心衰亚组分布:", d.ef_icd.value_counts().to_dict())

CONF = ['age', 'hb_min_24h', 'sofa', 'male', 'charlson',
        'ckd', 'diabetes', 'mech_vent_24h', 'vasopressor_24h', 'rrt_24h']

# ---------- 1. IPTW ----------
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

# ---------- 2. 主效应 + bootstrap ----------
def weighted_rd(data, out='death_28d'):
    t_ = data[data.prbc_24h_flag == 1]; c_ = data[data.prbc_24h_flag == 0]
    return (np.average(t_[out], weights=t_.w) - np.average(c_[out], weights=c_.w)) * 100

base = weighted_rd(d_ps)
boots = []
n = len(d_ps)
for _ in range(B):
    b = d_ps.iloc[rng.integers(0, n, n)]
    try:
        m = sm.Logit(b.prbc_24h_flag, sm.add_constant(b[CONF])).fit(disp=0)
        ps = m.predict(sm.add_constant(b[CONF]))
        b = b.assign(w=np.where(b.prbc_24h_flag == 1, 1/ps, 1/(1-ps)).clip(lo, hi))
        boots.append(weighted_rd(b))
    except Exception:
        continue
ci = np.percentile(boots, [2.5, 97.5])
print(f"\n[主效应] 28天死亡 RD = {base:+.2f}% (95%CI {ci[0]:+.2f} ~ {ci[1]:+.2f})  [论文: -5.67% (-7.24~-4.02)]")

# ---------- 3. 亚组 RD（通用函数） ----------
def subgroup_rd(data, grp_col, groups, out='death_28d'):
    res = {}
    for g in groups:
        s = data[data[grp_col] == g]
        if len(s) < 30 or s.prbc_24h_flag.nunique() < 2:
            res[g] = (len(s), np.nan, (np.nan, np.nan)); continue
        rd = weighted_rd(s, out)
        bs = []
        for _ in range(200):
            b = s.iloc[rng.integers(0, len(s), len(s))]
            if b.prbc_24h_flag.nunique() < 2: continue
            try:
                bs.append(weighted_rd(b, out))
            except Exception:
                continue
        ci_g = np.percentile(bs, [2.5, 97.5]) if len(bs) > 50 else (np.nan, np.nan)
        res[g] = (len(s), rd, ci_g)
    return res

def print_sub(res, title):
    print(f"\n[{title}] 加权 28 天死亡 RD (%):")
    for g, (nn, rd, ci_g) in res.items():
        if np.isnan(rd):
            print(f"  {g:14s} n={nn:5d}  (无法估计)")
        else:
            print(f"  {g:14s} n={nn:5d}  RD={rd:+.2f}% (95%CI {ci_g[0]:+.2f} ~ {ci_g[1]:+.2f})")

# 3a. ICD 亚组（复现论文 Table/Figure 2）
icd_groups = ['HFrEF(ICD)', 'HFpEF(ICD)', 'OtherHF(ICD)', 'NoHF']
print_sub(subgroup_rd(d_ps, 'ef_icd', icd_groups), "ICD 编码亚组（复现论文）")

# 3b. 超声 LVEF 亚组（核心新分析）
d_ps['ef_echo'] = np.where(d_ps.lvef.notna(),
                    np.where(d_ps.lvef < 40, 'LVEF<40(HFrEF)', 'LVEF>=40'),
                    'NoEcho')
echo_groups = ['LVEF<40(HFrEF)', 'LVEF>=40', 'NoEcho']
print_sub(subgroup_rd(d_ps, 'ef_echo', echo_groups), "超声 LVEF 亚组（新）")

# ---------- 4. 交互检验（论文方法：加权 outcome 模型加乘积项） ----------
def interaction_p(data, grp_col, ref, out='death_28d'):
    """加权 logistic: y ~ T + G + T:G，G 为二分类(参照=ref)，返回 T×G 的 P 值"""
    s = data[data[grp_col].isin([ref]) | (data[grp_col] != ref)].copy()
    s['G'] = (s[grp_col] != ref).astype(float)
    T = s.prbc_24h_flag.values.astype(float)
    X = np.column_stack([np.ones(len(s)), T, s.G.values, T * s.G.values] + [s[c].values for c in CONF])
    try:
        m = sm.GLM(s[out], X, family=sm.families.Binomial(), freq_weights=s.w.values).fit()
        return m.pvalues.iloc[3], m.params.iloc[3]
    except Exception as e:
        return np.nan, np.nan

# ICD HFrEF vs NoHF（复现论文 P=0.515）
s = d_ps[d_ps.ef_icd.isin(['HFrEF(ICD)', 'NoHF'])].copy()
s['G'] = (s.ef_icd == 'HFrEF(ICD)').astype(float)
T = s.prbc_24h_flag.values.astype(float)
X = np.column_stack([np.ones(len(s)), T, s.G.values, T * s.G.values] + [s[c].values for c in CONF])
m = sm.GLM(s.death_28d, X, family=sm.families.Binomial(), freq_weights=s.w.values).fit()
print(f"\n[交互检验-ICD] HFrEF vs NoHF: OR_interaction={np.exp(m.params.iloc[3]):.3f}, P={m.pvalues.iloc[3]:.3f}  [论文: P=0.515]")

# LVEF<40 vs >=40（超声定义，限有超声者）
s2 = d_ps[d_ps.ef_echo.isin(['LVEF<40(HFrEF)', 'LVEF>=40'])].copy()
s2['G'] = (s2.ef_echo == 'LVEF<40(HFrEF)').astype(float)
T2 = s2.prbc_24h_flag.values.astype(float)
X2 = np.column_stack([np.ones(len(s2)), T2, s2.G.values, T2 * s2.G.values] + [s2[c].values for c in CONF])
m2 = sm.GLM(s2.death_28d, X2, family=sm.families.Binomial(), freq_weights=s2.w.values).fit()
print(f"[交互检验-超声] LVEF<40 vs LVEF>=40: OR_interaction={np.exp(m2.params.iloc[3]):.3f}, P={m2.pvalues.iloc[3]:.3f} (n={len(s2)}, 其中LVEF<40: {(s2.G==1).sum()})")

# LVEF<40 vs NoEcho（纳入全部队列）
s3 = d_ps[d_ps.ef_echo.isin(['LVEF<40(HFrEF)', 'NoEcho'])].copy()
s3['G'] = (s3.ef_echo == 'LVEF<40(HFrEF)').astype(float)
T3 = s3.prbc_24h_flag.values.astype(float)
X3 = np.column_stack([np.ones(len(s3)), T3, s3.G.values, T3 * s3.G.values] + [s3[c].values for c in CONF])
m3 = sm.GLM(s3.death_28d, X3, family=sm.families.Binomial(), freq_weights=s3.w.values).fit()
print(f"[交互检验-超声] LVEF<40 vs 无超声: OR_interaction={np.exp(m3.params.iloc[3]):.3f}, P={m3.pvalues.iloc[3]:.3f}")

# ---------- 5. LVEF 亚组特征描述 ----------
print("\n[超声亚组基线特征]:")
for g in echo_groups:
    s = d_ps[d_ps.ef_echo == g]
    if len(s) == 0: continue
    print(f"  {g:16s} n={len(s):5d} 输血率={s.prbc_24h_flag.mean()*100:.1f}% "
          f"28天死亡={s.death_28d.mean()*100:.1f}% LVEF中位数={s.lvef.median():.0f} "
          f"SOFA中位数={s.sofa.median():.0f} Hb中位数={s.hb_min_24h.median():.1f}")
