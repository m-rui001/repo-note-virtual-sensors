# -*- coding: utf-8 -*-
"""R62-D 实验 124：K 不只是排序量——能不能给 I=3 的代价【定价】。
前沿律给 D@I=3 的预测值 D_pred = JC + K/(3 - R_exp)，K 只用 s=1e-6 的一次 DARE 解。
报同 rank 的相对误差分布，以及 1 次解 vs 80 次二分解的成本比。
"""
import numpy as np
from scipy.linalg import solve_discrete_are, eigh

_src = open('p0/exp_c_audit.py', encoding='utf-8').read().split("print(r'== E53")[0]
_ns = {'__name__': 'p'}
exec(compile(_src, 'p0/exp_c_audit.py[preamble]', 'exec'), _ns)
A, W, TH, JC, n = _ns['A'], _ns['W'], _ns['TH'], _ns['JC'], _ns['n']
sym = _ns['sym']
R_EXP = 1.168539
TARGET = 3.0
S = 1e-6
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


def at_rate(Z, target=TARGET):
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
    for c in range(30):
        Z, _ = np.linalg.qr(rng.standard_normal((n, rk)))
        pool.append((f'r{rk}-{c:02d}', rk, Z))
for rk in (1, 2, 3):
    pool.append((f'THmin-{rk}', rk, Vth[:, :rk]))
    pool.append((f'THmax-{rk}', rk, Vth[:, n - rk:]))

p('== D_pred(3) = JC + K/(3-R_exp)，K 来自 s=1e-6 一次解；对比 80 次二分真解出的 D(3) ==')
p('  rank   n   中位相对误差   90 分位   最差   预测最便宜=真最便宜?   THmax- 预测/真值')
for rk in (1, 2, 3):
    errs, pairs = [], []
    for nm, r, Z in pool:
        if r != rk:
            continue
        I0, D0, rho0 = curve(Z, S)
        if rho0 >= 1.0 or I0 <= R_EXP:
            continue
        K = (I0 - R_EXP) * (D0 - JC)
        d3 = at_rate(Z)
        if d3 is None:
            continue
        dpred = JC + K / (TARGET - R_EXP)
        errs.append((dpred - JC) / (d3 - JC))
        pairs.append((nm, dpred, d3))
    e = np.array(errs)
    pred_best = min(pairs, key=lambda q: q[1])[0]
    true_best = min(pairs, key=lambda q: q[2])[0]
    truemax = [q for q in pairs if q[0] == f'THmax-{rk}'][0]
    p('  %d    %2d      %+.3f%%      %+.3f%%   %+.3f%%      %s (%s / %s)   %.4f / %.4f'
      % (rk, len(e), (np.median(e) - 1) * 100, np.percentile(np.abs(e - 1), 90) * 100,
         np.max(np.abs(e - 1)) * 100,
         '是 ✓' if pred_best == true_best else '否 ✗', pred_best, true_best,
         truemax[1], truemax[2]))
p('\n注：误差是 (D_pred-JC)/(D_true-JC)，即**超出 LQR 地板那一段**的比例误差（JC=%.6f 本身精确）✓' % JC)
p('成本：定价一个设计 = 1 次 DARE 解；真解 I=3 = 80 次 ⇒ %.0f× 节省。' % 80)
open('p0/e124_out.txt', 'w', encoding='utf-8').write('\n'.join(out) + '\n')
