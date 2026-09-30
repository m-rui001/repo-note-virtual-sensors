# -*- coding: utf-8 -*-
r"""R63-D 实验 127：把 §75-C 的"K 只能近地板排序、不能定价"变成可修的定量限制。

前沿律的真身是两条腿：$\Delta I=I-R_{\exp}=a\,s+O(s^2)$，$D-J_c=b/s+O(1)$
$\Rightarrow D-J_c=K/\Delta I+O(1)$，$K=ab$。所以"定价失效"的**原因**是那个 $O(1)$：
$P(\Delta I):=(D-J_c)\,\Delta I$ 沿前沿并不恒等于 $K$，它随 $\Delta I$ 漂。
本节量的就是漂移的形状，并回答"多花一次 DARE 解能买回多少精度"：

 [F1] 单参数（§75 的规则）：$D_{\rm pred}=J_c+K/\Delta I_{\rm tgt}$，$K$ 只取自 $\Delta I\le0.3$。
 [F2] 两参数线性：$P(\Delta I)=K+K_2\,\Delta I$ 在 $\Delta I\le0.3$ 上最小二乘，再回代。
 [F3] 幂律带指数：$\log(D-J_c)=\log \tilde K-\nu_{\rm fit}\log\Delta I$（$\nu$ 逐设计拟合，不写死 1）。
 三者都**只准用近地板数据拟合**，在 $\Delta I_{\rm tgt}=3-R_{\exp}=1.831461$ 处对真值比误差。
 [F4] 真值由二分（80 次解）给出，不许用同一批点既拟合又当答案。
 [F5] 判据：$\mathrm{med}|e|$ 从 F1 到 F2 若下降 $<2\times$，则"两次解换精度"不划算，
      §75-C 那句限制原样保留；若下降 $\ge3\times$ 且符号不再随 rank 翻，就替换成新的可写句。
"""
import sys
import numpy as np
from scipy.linalg import solve_discrete_are, eigh
sys.stdout.reconfigure(encoding='utf-8')

_src = open('p0/exp_c_audit.py', encoding='utf-8').read().split("print(r'== E53")[0]
_ns = {'__name__': 'p'}
exec(compile(_src, 'p0/exp_c_audit.py[preamble]', 'exec'), _ns)
A, W, TH, JC, n = _ns['A'], _ns['W'], _ns['TH'], _ns['JC'], _ns['n']
sym = _ns['sym']
R_EXP = 1.168539
TARGET = 3.0
DIT = TARGET - R_EXP
SGRID = 10.0 ** np.linspace(-6.0, 1.5, 46)
FITMAX = 0.30


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


_, Vth = eigh(TH)
rng = np.random.default_rng(20260930)
pool = []
for rk in (1, 2, 3):
    for c in range(10):
        Z, _ = np.linalg.qr(rng.standard_normal((n, rk)))
        pool.append((f'r{rk}-{c:02d}', rk, Z))
    pool.append((f'THmin-{rk}', rk, Vth[:, :rk]))
    pool.append((f'THmax-{rk}', rk, Vth[:, n - rk:]))

print('== F1/F2/F3：只用 $\\Delta I\\le%.2f$ 的近地板点拟合，在 $\\Delta I=%.6f$（$I=3$）定价 =='
      % (FITMAX, DIT))
print('  design    rank  n_fit  D_true     e(F1)%    e(F2)%    e(F3)%    nu_fit   P漂(相对)')
res = []
for name, rk, Z in pool:
    Dtrue = true_at(Z)
    if Dtrue is None:
        print(f'  {name:<9}{rk:>4}   (I=3 不在扫描范围)')
        continue
    Is, Ds = [], []
    for s in SGRID:
        try:
            I, D = curve(Z, s)
        except Exception:
            continue
        if not (np.isfinite(I) and np.isfinite(D)):
            continue
        Is.append(I); Ds.append(D)
    Is = np.array(Is); Ds = np.array(Ds)
    di = Is - R_EXP
    ok = (di > 1e-4) & (di <= FITMAX) & (Ds > JC)
    if ok.sum() < 4:
        print(f'  {name:<9}{rk:>4}   近地板点不足（{ok.sum()}）')
        continue
    dif, Df = di[ok], Ds[ok] - JC
    P = Df * dif
    K = np.median(P)
    e1 = (JC + K / DIT) / Dtrue - 1.0
    cf = np.polyfit(dif, P, 1)
    e2 = (JC + np.polyval(cf, DIT) / DIT) / Dtrue - 1.0
    lf = np.polyfit(np.log(dif), np.log(Df), 1)
    nu = -lf[0]
    e3 = (JC + np.exp(lf[1]) * DIT ** (-nu)) / Dtrue - 1.0
    drift = (P.max() - P.min()) / K
    res.append((name, rk, e1, e2, e3, nu, drift, ok.sum()))
    print(f'  {name:<9}{rk:>4}  {ok.sum():>4}   {Dtrue:>9.4f}  {100*e1:>8.2f}  {100*e2:>8.2f}  '
          f'{100*e3:>8.2f}   {nu:>6.4f}   {drift:>7.3f}')

R = np.array([[r[2], r[3], r[4]] for r in res])
nuu = np.array([r[5] for r in res])
dr = np.array([r[6] for r in res])
print(f'\n== 汇总（{len(R)} 个设计）==')
for j, lab in enumerate(('F1 单参数 K/ΔI', 'F2 两参数 K+K2·ΔI', 'F3 逐设计幂律 ν')):
    a = R[:, j]
    print(f'  {lab:<20} med|e|={100*np.median(np.abs(a)):>7.2f}%  最差={100*np.max(np.abs(a)):>8.2f}%  '
          f'med有符号={100*np.median(a):>+8.2f}%  |e|<5%占比={np.mean(np.abs(a)<0.05):.2f}')
print(f'  $\\nu$ 拟合：med={np.median(nuu):.4f}  IQR=[{np.percentile(nuu,25):.4f},{np.percentile(nuu,75):.4f}]  '
      f'全距=[{nuu.min():.4f},{nuu.max():.4f}]   $\\nu\\equiv1.0001$ 的偏离是否随 $\\Delta I$ 变宽：见逐行')
print(f'  $P$ 的相对漂移（$\\Delta I\\le0.3$ 窗口内）：med={np.median(dr):.3f}  max={dr.max():.3f}')

print('\n== 符号随 rank 是否翻（§75-C 的那条硬约束）==')
for rk in (1, 2, 3):
    for j, lab in enumerate(('F1', 'F2', 'F3')):
        a = np.array([r[2 + j] for r in res if r[1] == rk])
        if len(a) == 0:
            continue
        print(f'  rank={rk} {lab}: n={len(a)}  med={100*np.median(a):>+7.2f}%  '
              f'负/正={int(np.sum(a<0))}/{int(np.sum(a>0))}')

print('\n== 成本口径：F1 用 1 次 DARE 解；F2/F3 用 %d~%d 次近地板解 + 1 次拟合；真值用 80 次二分 =='
      % (min(r[7] for r in res), max(r[7] for r in res)))
print('判据 [F5]：F2 的 med|e| 相对 F1 若降不到 2 倍，限制句保留原样。')
