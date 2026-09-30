r"""E64 = 直接测 C 的 §12.4/§12.6 那句话："子空间只通过地板起作用，率对 range S 不敏感(<=0.5%)"。

他的证据是"同一 rank 族里换子空间，率差 <=0.5%"（Schur 尾 vs Theta 投影，rank-1 嵌入后）。
那句话的**正确量词**是：rate frontier 是 range(S) 的函数，且该函数只通过 Phi_0(range S) 依赖 range(S)。
可判定的反例形状：找两个 range（r=1 就是两个方向）满足 Phi_0 相同（相对差 <0.5%）
但同一目标率上的最优 J 显著不同（>5%）。若存在，命题为假；若扫遍球面都找不到，那是一个
真实的经验陈述（值得写下来，但要说清是"扫过的 400 个方向里没找到"）。

r=1 是 §12.6 那条的最小情形，也是我已知的 Grassmannian 有双稳态（10.169 vs 86.666）的维数，
所以地板谱是散布的、不是单值的，正好用来做等地板配对。
"""
import sys, time
sys.stdout.reconfigure(encoding='utf-8')
import numpy as np
from scipy.linalg import schur, solve_discrete_are as sdare

np.set_printoptions(precision=4, suppress=True, linewidth=170)
_src = open('p0/exp_c_audit.py', encoding='utf-8').read().split("print(r'== E53")[0]
_ns = {'__name__': 'p'}
exec(compile(_src, 'p0/exp_c_audit.py[preamble]', 'exec'), _ns)
A, W, TH, JC, n, sym = _ns['A'], _ns['W'], _ns['TH'], _ns['JC'], _ns['n'], _ns['sym']
print('== E64：等地板配对检验"率只通过地板依赖子空间" ==')
In_ = np.eye(n)


def JIfast(S):
    ev, EV = np.linalg.eigh(sym(S))
    keep = ev > max(1e-11 * ev.max(), 1e-300)
    if not keep.any():
        return np.inf, np.inf
    F = (EV[:, keep] * np.sqrt(ev[keep])).T
    k = F.shape[0]
    try:
        Pt = sym(sdare(A.T, F.T, sym(W), np.eye(k)))
    except Exception:
        return np.inf, np.inf
    if not np.all(np.isfinite(Pt)):
        return np.inf, np.inf
    FPt = F @ Pt
    try:
        P = sym(Pt - FPt.T @ np.linalg.solve(FPt @ F.T + np.eye(k), FPt))
    except Exception:
        return np.inf, np.inf
    cl = A - A @ Pt @ F.T @ np.linalg.solve(FPt @ F.T + np.eye(k), F)
    if np.max(np.abs(np.linalg.eigvals(cl))) >= 1.0 - 1e-9:
        return np.inf, np.inf
    ld = np.linalg.slogdet(Pt)[1] - np.linalg.slogdet(P)[1]
    return JC + float(np.trace(TH @ P)), 0.5 * ld / np.log(2)


def J_hit1(z, It, lo=-4.0, hi=14.0, steps=28):
    S1 = np.outer(z, z)
    if not np.isfinite(JIfast(10.0 ** hi * S1))[0]:
        return np.nan
    for _ in range(steps):
        m = 0.5 * (lo + hi)
        j, i_ = JIfast(10.0 ** m * S1)
        if np.isfinite(i_) and i_ >= It:
            hi = m
        else:
            lo = m
    return JIfast(10.0 ** hi * S1)[0]


rng = np.random.default_rng(5)
Z = rng.normal(0, 1, (500, n))
Z /= np.linalg.norm(Z, axis=1, keepdims=True)
rows = []
t0 = time.time()
for z in Z:
    fl = JIfast(10.0 ** 14 * np.outer(z, z))[0]
    if not np.isfinite(fl):
        continue
    j2 = J_hit1(z, 2.0); j3 = J_hit1(z, 3.0)
    rows.append((fl, j2, j3, float(np.dot(z, np.linalg.eigh(TH)[1][:, -1]))))
R = np.array(rows)
print('  可检测方向 %d/%d，%.0f s；地板范围 [%.4f, %.4f] = j_c+Phi0, Phi0 in [%.4f, %.4f]'
      % (len(R), len(Z), time.time() - t0, R[:, 0].min(), R[:, 0].max(), R[:, 0].min() - JC, R[:, 0].max() - JC))
print('  有限 J(I=2) 的比例 %.2f；J(I=2) 范围 [%.3f, %.3f]'
      % (np.isfinite(R[:, 1]).mean(), np.nanmin(R[:, 1]), np.nanmax(R[:, 1])))

print('\n[A] 等地板配对：地板相对差 <0.5%% 的方向对，率差多少？（命题为假 <=> 存在 >5%% 的对）')
best_pair = None
cnt = 0
for i in range(len(R)):
    for j in range(i + 1, len(R)):
        f1, f2 = R[i, 0], R[j, 0]
        if abs(f1 - f2) / max(f1, f2) > 0.005:
            continue
        for c, nm in ((1, 'I=2'), (2, 'I=3')):
            a, b = R[i, c], R[j, c]
            if not (np.isfinite(a) and np.isfinite(b)):
                continue
            cnt += 1
            rel = abs(a - b) / min(a, b)
            if best_pair is None or rel > best_pair[0]:
                best_pair = (rel, nm, f1, f2, a, b, i, j)
print('  配对比较次数 %d；最大相对率差 %.2f%%' % (cnt, 100 * best_pair[0]))
print('  取到最大者：%s，地板 %.4f vs %.4f（相对差 %.3f%%），J=%.4f vs %.4f'
      % (best_pair[1], best_pair[2], best_pair[3], 100 * abs(best_pair[2] - best_pair[3]) / best_pair[2],
         best_pair[4], best_pair[5]))

print('\n[B] 分桶：把地板切 10 桶，看桶内 J(I=2) 的散布（桶内相对极差）')
q = np.quantile(R[:, 0], np.linspace(0, 1, 11))
for b in range(10):
    sel = (R[:, 0] >= q[b]) & (R[:, 0] <= q[b + 1]) & np.isfinite(R[:, 1])
    if sel.sum() < 3:
        continue
    v = R[sel, 1]
    print('  桶%2d 地板[%.3f,%.3f] n=%3d: J(I=2) [%.3f, %.3f] 桶内极差 %.2f%%'
          % (b, q[b], q[b + 1], sel.sum(), v.min(), v.max(), 100 * (v.max() - v.min()) / v.min()))

print('\n[C] 地板相同(同一个 Phi0)但方向不同——最干净的对照：地板取分位数后，在每桶里取极值方向与 Theta 头的夹角')
for b in (2, 5, 8):
    sel = (R[:, 0] >= q[b]) & (R[:, 0] <= q[b + 1]) & np.isfinite(R[:, 1])
    if sel.sum() < 3:
        continue
    idx = np.where(sel)[0]
    lo_i = idx[np.argmin(R[sel, 1])]; hi_i = idx[np.argmax(R[sel, 1])]
    print('  桶%d: 最好方向 cos(z,theta_1)=%.3f 地板 %.4f J=%.4f | 最差 cos=%.3f 地板 %.4f J=%.4f | 率差 %.2f%%'
          % (b, abs(R[lo_i, 3]), R[lo_i, 0], R[lo_i, 1], abs(R[hi_i, 3]), R[hi_i, 0], R[hi_i, 1],
             100 * (R[hi_i, 1] - R[lo_i, 1]) / R[lo_i, 1]))
print('== E64 done ==')
