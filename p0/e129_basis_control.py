# -*- coding: utf-8 -*-
r"""实验 129：e128 的 G2（3 参数指数式）过了 [G6] 判据，但在写进 board 之前必须先排除两种平凡解释：
  (i) 增益来自"地板项"而不是"指数"（γ_fit 中位 0.0256、大量行贴下界 ⇒ 数据偏向 γ→0 即双曲）；
  (ii) 增益只是"7 个点配 3 个旋钮"的过拟合 ⇒ 用留出（holdout）检验，而不是只看外推。

同一个池（seed 20260930，rank 1–3 各 12）、同一条曲线，基函数逐个换：
  H2  K/ΔI                          （§74/75 的原律，1 参）
  HF  δ + K/ΔI                      （原律 + 自由地板，2 参）
  P2  δ + K·ΔI                      （2 参多项式）
  P3  δ + K1·ΔI + K2·ΔI²            （3 参多项式——与 G2 同旋钮数的**假基**）
  E2  K/(e^{γΔI}−1)                 （2 参纯指数，= e128 的 G1）
  EF3 δ + K/(e^{γΔI}−1)             （3 参指数 + 地板，= e128 的 G2）
  PL3 δ + K/ΔI^ν                    （3 参幂律 + 地板）
两条评测：
  [X1] 外推：拟合窗 ΔI∈[0.30,0.60]，定价 ΔI=1.831461（与 e128 同口径，便于对照）；
  [X2] 留出：仍只在 [0.30,0.60] 拟合，去预测**从未参与拟合**的 [0.80,1.00] 上的点（逐点相对误差中位）。
       ⇒ 过拟合会在 [X2] 上现形；真基函数在 [X2] 上也该站得住。
 [X3] 判据（跑前写死）：若 P3（同 knob 数的假基）在两条评测上都不劣于 EF3，则 e128 的 G2
       **不配**被称为律，只配称为"三参数插值"；此时 §75-C 的限制句可软化但不许替换成指数式。
       若 HF（2 参）就达到 EF3 的精度，则结论应写成"错在地板，不在指数"。
"""
import sys
import numpy as np
from scipy.linalg import solve_discrete_are, eigh
from scipy.optimize import curve_fit
sys.stdout.reconfigure(encoding='utf-8')

_src = open('p0/exp_c_audit.py', encoding='utf-8').read().split("print(r'== E53")[0]
_ns = {'__name__': 'p'}
exec(compile(_src, 'p0/exp_c_audit.py[preamble]', 'exec'), _ns)
A, W, TH, JC, n = _ns['A'], _ns['W'], _ns['TH'], _ns['JC'], _ns['n']
sym = _ns['sym']
R_EXP = 1.168539
TARGET = 3.0
DIT = TARGET - R_EXP
SGRID = 10.0 ** np.linspace(-6.0, 2.5, 90)
FLO, FI = 0.30, 0.60      # 拟合带
HLO, HI = 0.80, 1.00      # 留出带


def curve(Z, s):
    r = Z.shape[1]
    Ir = np.eye(r)
    C = np.sqrt(s) * Z.T
    Pm = sym(solve_discrete_are(A.T, C.T, W, Ir))
    Sm = C @ Pm @ C.T + Ir
    Lk = Pm @ C.T @ np.linalg.inv(Sm)
    Pp = sym(Pm - Lk @ Sm @ Lk.T)
    I = 0.5 * np.log(np.linalg.det(Ir + C @ Pm @ C.T)) / np.log(2.0)
    return I, JC + np.trace(TH @ Pp)


def true_at(Z, target=TARGET):
    lo, hi = -7.0, 5.0
    if not (curve(Z, 10 ** lo)[0] < target < curve(Z, 10 ** hi)[0]):
        return None
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if curve(Z, 10 ** mid)[0] < target:
            lo = mid
        else:
            hi = mid
    return curve(Z, 10 ** (0.5 * (lo + hi)))[1]


FORMS = {
    'H2':  (lambda di, K: K / di,                              [1.0],                    1,
            [0.0],            [np.inf]),
    'HF':  (lambda di, K, d: d + K / di,                       [1.0, 0.1],               2,
            [0.0, 0.0],       [np.inf, np.inf]),
    'P2':  (lambda di, K, d: d + K * di,                       [1.0, 0.1],               2,
            [-np.inf, -np.inf], [np.inf, np.inf]),
    'P3':  (lambda di, K, d, c: d + K * di + c * di ** 2,      [1.0, 0.1, 0.1],          3,
            [-np.inf, -np.inf, -np.inf], [np.inf, np.inf, np.inf]),
    'E2':  (lambda di, K, g: K / (np.exp(g * di) - 1.0),       [1.0, 1.3863],            2,
            [1e-6, 1e-3],     [np.inf, 20.0]),
    'EF3': (lambda di, K, g, d: d + K / (np.exp(g * di) - 1.0), [1.0, 1.3863, 0.1],       3,
            [1e-6, 1e-3, 0.0], [np.inf, 20.0, np.inf]),
    'PL3': (lambda di, K, nu, d: d + K / di ** nu,             [1.0, 1.0, 0.1],          3,
            [0.0, 0.05, 0.0],  [np.inf, 5.0, np.inf]),
}

_, Vth = eigh(TH)
rng = np.random.default_rng(20260930)
pool = []
for rk in (1, 2, 3):
    for c in range(10):
        Z, _ = np.linalg.qr(rng.standard_normal((n, rk)))
        pool.append((f'r{rk}-{c:02d}', rk, Z))
    pool.append((f'THmin-{rk}', rk, Vth[:, :rk]))
    pool.append((f'THmax-{rk}', rk, Vth[:, n - rk:]))

acc = {k: {'ext': [], 'ho': [], 'rk': [], 'g': []} for k in FORMS}
nuse = []
for name, rk, Z in pool:
    Dtrue = true_at(Z)
    if Dtrue is None:
        continue
    Is, Ds = [], []
    for s in SGRID:
        try:
            I, D = curve(Z, s)
        except Exception:
            continue
        if np.isfinite(I) and np.isfinite(D):
            Is.append(I); Ds.append(D)
    Is = np.array(Is); Ds = np.array(Ds)
    di = Is - R_EXP
    fitm = (di >= FLO) & (di <= FI)
    hom = (di >= HLO) & (di <= HI)
    if fitm.sum() < 3 or hom.sum() < 2:
        print(f'  {name}: 拟合带 {fitm.sum()} 点 / 留出带 {hom.sum()} 点 —— 带不足，剔除')
        continue
    nuse.append((name, rk, fitm.sum(), hom.sum()))
    daf, Daf = di[fitm], Ds[fitm] - JC
    dah, Dah = di[hom], Ds[hom] - JC
    for k, (f, p0, npar, lo, ub) in FORMS.items():
        try:
            p, _ = curve_fit(f, daf, Daf, p0=p0, maxfev=60000, bounds=(lo, ub))
            ext = (JC + f(DIT, *p)) / Dtrue - 1.0
            pred_h = f(dah, *p)
            ho = float(np.median(np.abs(pred_h / Dah - 1.0)))
            acc[k]['ext'].append(ext); acc[k]['ho'].append(ho); acc[k]['rk'].append(rk)
            if k in ('E2', 'EF3', 'PL3'):
                acc[k]['g'].append(p[1])
        except Exception as e:
            print(f'  {name} {k} fit fail: {type(e).__name__}')

print(f'\n== 有效设计 {len(nuse)} 个；拟合带 [%.2f,%.2f] 3~%d 点，留出带 [%.2f,%.2f] 2~%d 点 =='
      % (FLO, FI, max(x[2] for x in nuse), HLO, HI, max(x[3] for x in nuse)))
print(f'{"形式":<6}{"knob":>5}   {"外推 med|e|":>12}{"最差":>9}   {"留出 med|e|":>12}{"最差":>9}   外推有符号med')
for k in FORMS:
    a = np.array(acc[k]['ext']); b = np.array(acc[k]['ho'])
    if len(a) == 0:
        print(f'  {k}: 无有效拟合')
        continue
    print(f'  {k:<5}{FORMS[k][2]:>4}   {100*np.median(np.abs(a)):>11.2f}%{100*a.max():>8.1f}%'
          f'   {100*np.median(b):>11.2f}%{100*b.max():>8.1f}%   {100*np.median(a):>+10.2f}%')

print('\n== [X3] 的两条判定 ==')
def mm(k, which):
    a = np.array(acc[k][which])
    return np.median(np.abs(a)) if len(a) else np.nan
eE, hE = mm('EF3', 'ext'), mm('EF3', 'ho')
eP, hP = mm('P3', 'ext'), mm('P3', 'ho')
eH, hH = mm('HF', 'ext'), mm('HF', 'ho')
eL, hL = mm('PL3', 'ext'), mm('PL3', 'ho')
print(f'  P3（同 3 旋钮、符号自由的假基）：外推 {100*eP:.2f}% / 留出 {100*hP:.2f}%')
print(f'  EF3（=e128 的 G2，3 参指数+地板）：外推 {100*eE:.2f}% / 留出 {100*hE:.2f}%')
print(f'  PL3（3 参幂律+地板）：外推 {100*eL:.2f}% / 留出 {100*hL:.2f}%')
print(f'  HF（2 参，原律+自由地板）：外推 {100*eH:.2f}% / 留出 {100*hH:.2f}%')
print(f'  ⇒ 判定一（G2 是否只是"三参数插值"）：'
      f'{"P3 不劣 ⇒ 是，不配叫律" if eP <= eE else "P3 显著更差 ⇒ 不是；旋钮数不是解释（外推差 %.1f×）"}'
      % ((eP / eE) if eE > 0 else float("nan")))
print(f'     留出带上同结论？{"是" if hP > hL else "否（留出上 P3 反而不差，须警惕）"}'
      f'  P3 留出 {100*hP:.2f}% vs 留出最优 {100*min(x for x in [hH,hL,hE] if x==x):.2f}%')
print(f'  ⇒ 判定二（增益来自地板还是指数）：HF(2 参含地板) {100*eH:.2f}% vs PL3(3 参含地板+自由ν) {100*eL:.2f}%'
      f' vs EF3(3 参含地板+自由γ) {100*eE:.2f}%')

gE = np.array(acc['EF3']['g']); gP3 = np.array(acc['PL3']['g']); g2 = np.array(acc['E2']['g'])
print('\n== 指数/幂指数的可识别性 ==')
if len(g2):
    print(f'  E2  γ_fit: med={np.median(g2):.4f}  全距=[{g2.min():.4f},{g2.max():.4f}]  贴下界占比={np.mean(g2<0.01):.2f}')
if len(gE):
    print(f'  EF3 γ_fit: med={np.median(gE):.4f}  IQR=[{np.percentile(gE,25):.4f},{np.percentile(gE,75):.4f}]'
          f'  全距=[{gE.min():.4f},{gE.max():.4f}]  贴下界(<0.01)占比={np.mean(gE<0.01):.2f}  '
          f'落在[1.0,2.0](含理论 2ln2=1.3863)占比={np.mean((gE>1.0)&(gE<2.0)):.2f}')
if len(gP3):
    print(f'  PL3 ν_fit: med={np.median(gP3):.4f}  IQR=[{np.percentile(gP3,25):.4f},{np.percentile(gP3,75):.4f}]'
          f'  全距=[{gP3.min():.4f},{gP3.max():.4f}]  （对照 e127 F3 近地板 ν med=1.0044）')

print('\n== rank 分层（留出误差的符号由谁决定）==')
rka = np.array(acc['EF3']['rk'])
for k in ('H2', 'HF', 'P3', 'EF3'):
    if not acc[k]['ext']:
        continue
    line = []
    for rk in (1, 2, 3):
        a = np.array([v for v, r in zip(acc[k]['ext'], acc[k]['rk']) if r == rk])
        line.append(f'rank{rk}:{100*np.median(a):+7.2f}%(n={len(a)})')
    print(f'  {k:<5}' + '  '.join(line))
