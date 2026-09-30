r"""E64b = 把 E64 的等地板配对做到可引用强度：
 (1) 收紧容差到 0.1%/0.01%，看率散布是否仍然存在；
 (2) 交叉对计数：地板更差(Phi0 更大)却率更省(J 更小)的方向对——单调性(率只通过地板)的直接反例；
 (3) 同一地板水平线上取一圈方向，看 J(I) 沿圈的变化幅度。
"""
import sys, time
sys.stdout.reconfigure(encoding='utf-8')
import numpy as np
from scipy.linalg import solve_discrete_are as sdare
_src = open('p0/exp_c_audit.py', encoding='utf-8').read().split("print(r'== E53")[0]
_ns = {'__name__': 'p'}; exec(compile(_src, 'p0/exp_c_audit.py[preamble]', 'exec'), _ns)
A, W, TH, JC, n, sym = _ns['A'], _ns['W'], _ns['TH'], _ns['JC'], _ns['n'], _ns['sym']
print('== E64b：等地板配对的可引用强度 ==')

def JIfast(S):
    ev, EV = np.linalg.eigh(sym(S)); keep = ev > 1e-11 * ev.max()
    if not keep.any(): return np.inf, np.inf
    F = (EV[:, keep] * np.sqrt(ev[keep])).T; k = F.shape[0]
    try: Pt = sym(sdare(A.T, F.T, sym(W), np.eye(k)))
    except Exception: return np.inf, np.inf
    if not np.all(np.isfinite(Pt)): return np.inf, np.inf
    FPt = F @ Pt
    try: M = np.linalg.solve(FPt @ F.T + np.eye(k), FPt)
    except Exception: return np.inf, np.inf
    P = sym(Pt - FPt.T @ M)
    cl = A - A @ Pt @ F.T @ M
    if np.max(np.abs(np.linalg.eigvals(cl))) >= 1.0 - 1e-9: return np.inf, np.inf
    ld = np.linalg.slogdet(Pt)[1] - np.linalg.slogdet(P)[1]
    return JC + float(np.trace(TH @ P)), 0.5 * ld / np.log(2)

def J_hit1(z, It, lo=-4.0, hi=14.0, steps=26):
    S1 = np.outer(z, z)
    if not np.isfinite(JIfast(10.0 ** hi * S1))[0]: return np.nan
    for _ in range(steps):
        m = 0.5 * (lo + hi); j, i_ = JIfast(10.0 ** m * S1)
        if np.isfinite(i_) and i_ >= It: hi = m
        else: lo = m
    return JIfast(10.0 ** hi * S1)[0]

rng = np.random.default_rng(5)
Z = rng.normal(0, 1, (500, n)); Z /= np.linalg.norm(Z, axis=1, keepdims=True)
t0 = time.time()
FL = np.array([JIfast(10.0 ** 14 * np.outer(z, z))[0] for z in Z])
J2 = np.array([J_hit1(z, 2.0) for z in Z])
J5 = np.array([J_hit1(z, 5.0) for z in Z])
ok = np.isfinite(FL) & np.isfinite(J2) & np.isfinite(J5)
FL, J2, J5 = FL[ok], J2[ok], J5[ok]
print('  可用方向 %d/%d (%.0f s)' % (ok.sum(), len(Z), time.time() - t0))

for tol in (0.005, 0.001, 0.0001):
    sp = []; npair = 0
    for i in range(len(FL)):
        d = np.abs(FL - FL[i]) / np.maximum(FL, FL[i])
        j = np.where((d < tol) & (np.arange(len(FL)) > i))[0]
        for k in j:
            npair += 1
            sp.append(abs(J2[i] - J2[k]) / min(J2[i], J2[k]))
            sp.append(abs(J5[i] - J5[k]) / min(J5[i], J5[k]))
    if not sp:
        print('  地板相对差 <%.4f: 无配对' % tol); continue
    sp = np.array(sp)
    print('  地板相对差 <%.4f: 配对 %d 组，桶内率相对极差 中位 %.1f%% / 最大 %.1f%%'
          % (tol, npair, 100 * np.median(sp), 100 * sp.max()))

for nm, JV in (('I=2', J2), ('I=5', J5)):
    cross = 0; worst = 0.0; ex = None
    for i in range(len(FL)):
        worse = (FL > FL[i]) & (JV < JV[i])
        cross += int(worse.sum())
        if worse.any():
            k = np.where(worse)[0][np.argmin(JV[worse] / JV[i])]
            rel = (FL[k] - FL[i]) / FL[i]
            if rel > worst: worst, ex = rel, (i, k)
    print('  %s: 交叉对（地板更差但率更省）共 %d / %d 个方向参与；'
          '最大的一例地板差 %+.1f%% 而 J 由 %.3f 降到 %.3f'
          % (nm, cross, len(FL), 100 * worst, JV[ex[0]], JV[ex[1]]) if ex else '(无)')

print('\n[D] 一圈方向，地板固定在 Phi0=100 附近：沿圈 J(I=2) 的变化')
target = 100.0 + JC
# 从一个已知方向出发，在等值面上随机走（用二维子空间内旋转），只报告实际落到的地板散布
base = Z[np.argmin(np.abs(FL - target))]
u, v = np.linalg.qr(np.column_stack([base, rng.normal(0, 1, n)]))[:2]
uu, vv = u[:, 0], v
vv = vv - np.dot(vv, uu) * uu; vv /= np.linalg.norm(vv)
pts = []
for th in np.linspace(0, 2 * np.pi, 25):
    z = np.cos(th) * uu + np.sin(th) * vv
    f = JIfast(10.0 ** 14 * np.outer(z, z))[0]
    if np.isfinite(f): pts.append((f, J_hit1(z, 2.0)))
P = np.array(pts)
print('  同一大圆上：地板 [%.3f, %.3f]（跨度 %.1f%%），J(I=2) [%.3f, %.3f]（跨度 %.1f%%）'
      % (P[:, 0].min(), P[:, 0].max(), 100 * (P[:, 0].max() / P[:, 0].min() - 1),
         np.nanmin(P[:, 1]), np.nanmax(P[:, 1]),
         100 * (np.nanmax(P[:, 1]) / np.nanmin(P[:, 1]) - 1)))
print('  相关系数 corr(地板, J) = %.4f —— 若 <1，地板不是率的充分统计量'
      % np.corrcoef(P[:, 0], np.nan_to_num(P[:, 1], nan=np.nanmin(P[:, 1])))[0, 1])
