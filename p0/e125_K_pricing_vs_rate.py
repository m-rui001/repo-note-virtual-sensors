# -*- coding: utf-8 -*-
"""R62-D 实验 125：K 定价误差随目标率 I 的走向——判"K 只能排序不能定价"是全局失败还是渐近段失败。
每个设计先用 s=1e-6 定一次 K，再对 I in {1.25,1.35,1.5,2,3,4} 二分解真值，报 (D_pred-JC)/(D_true-JC)。
"""
import numpy as np
from scipy.linalg import solve_discrete_are, eigh

_src = open('p0/exp_c_audit.py', encoding='utf-8').read().split("print(r'== E53")[0]
_ns = {'__name__': 'p'}
exec(compile(_src, 'p0/exp_c_audit.py[preamble]', 'exec'), _ns)
A, W, TH, JC, n = _ns['A'], _ns['W'], _ns['TH'], _ns['JC'], _ns['n']
sym = _ns['sym']
R_EXP, S = 1.168539, 1e-6
TG = (1.25, 1.35, 1.50, 2.00, 3.00, 4.00)
out = []


def p(s):
    out.append(str(s)); print(s, flush=True)


def curve(Z, s):
    r = Z.shape[1]
    Ir = np.eye(r)
    C = np.sqrt(s) * Z.T
    Pm = sym(solve_discrete_are(A.T, C.T, W, Ir))
    Sm = C @ Pm @ C.T + Ir
    Lk = Pm @ C.T @ np.linalg.inv(Sm)
    Pp = sym(Pm - Lk @ Sm @ Lk.T)
    I = 0.5 * np.log(np.linalg.det(Ir + C @ Pm @ C.T)) / np.log(2.0)
    return I, JC + np.trace(TH @ Pp), np.max(np.abs(np.linalg.eigvals((np.eye(n) - Lk @ C) @ A)))


def at_rate(Z, target):
    lo, hi = -7.0, 5.0
    if not (curve(Z, 10 ** lo)[0] < target < curve(Z, 10 ** hi)[0]):
        return None
    for _ in range(70):
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
    for c in range(30):
        Z, _ = np.linalg.qr(rng.standard_normal((n, rk)))
        pool.append((f'r{rk}-{c:02d}', rk, Z))
    pool.append((f'THmin-{rk}', rk, Vth[:, :rk]))
    pool.append((f'THmax-{rk}', rk, Vth[:, n - rk:]))

Ks = {}
for nm, rk, Z in pool:
    I0, D0, rho0 = curve(Z, S)
    if rho0 < 1.0 and I0 > R_EXP:
        Ks[nm] = (rk, (I0 - R_EXP) * (D0 - JC), Z)

p('== K 定价误差 vs 目标率（中位 / 90 分位 / 最差，单位 %%，对 (D_pred-JC)/(D_true-JC)）==')
p('  rank')
for rk in (1, 2, 3):
    p('   rank %d   I 目标   中位    90分位   最差     n' % rk)
    for t in TG:
        e = []
        for nm, (r, K, Z) in Ks.items():
            if r != rk:
                continue
            d3 = at_rate(Z, t)
            if d3 is None:
                continue
            e.append(((JC + K / (t - R_EXP)) - JC) / (d3 - JC))
        e = np.array(e)
        p('            %5.2f  %+7.2f%%  %6.2f%%  %7.2f%%  %2d'
          % (t, (np.median(e) - 1) * 100, np.percentile(np.abs(e - 1), 90) * 100,
             np.max(np.abs(e - 1)) * 100, len(e)))
open('p0/e125_out.txt', 'w', encoding='utf-8').write('\n'.join(out) + '\n')
