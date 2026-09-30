"""E19：两件事。

(1) 交出我 p≈0.53 的原始 (D, Δ) 数列（C3 §3 要求对齐数据而不是对齐结论）。
    E13 的 F2 网格是 D_min+0.3, D_min+2.0, 50, 65, 80 —— 前两点离地板 0.3 与 2.0，
    我当时用**最左两点**做对数斜率，得到 −0.53。这里把逐对本地斜率全算出来。

(2) 回答 C3 §5 第二个公开问题：range(S*(65)) 的受约束 Φ(+2.4%) 为什么会贴住 Φ_min(2)(+2.3%)。
    猜想：这是 **D 大的极限重合**，不是结构巧合。
    检验：对一列 D，取无约束最优 S*(D) 的 leading r 维跨度，与 C 存下来的 Φ_min(r) 极小元
    (.work3/c13_V_r{1,2,3}.npy，只读) 做主角/对齐度，并算两边 Φ。
    若 D 增大时 range(S*(D)) 的 leading 方向逐列转向 Φ_min(r) 的极小子空间（也即 Θ 的主方向），
    那 C 的"率侧 r*(D) 与地板侧合成一个定理"就有了正确的形式：**极限命题，不是等式命题**。
"""
import sys
import numpy as np
sys.stdout.reconfigure(encoding='utf-8')
from p0.p0_replicate_letter import A, W, n, sym, ln2, ctrl, unconstrained_sdp, init_gamma
from p0.exp_authoritative import noiseless

np.set_printoptions(precision=4, suppress=True)
Pc, K, Th = ctrl(np.eye(n))
unc = float(np.trace(W @ Pc))
wTh, VTh = np.linalg.eigh(sym(Th))            # Θ 的特征方向（升序）


def Phi(Y, iters=60000, tol=1e-14):
    Y = np.atleast_2d(Y)
    if Y.shape[0] != n and Y.shape[1] == n:      # 容错：给的是 r×n 就转置
        Y = Y.T
    if Y.shape[1] >= n:
        return 0.0
    Qb = np.linalg.qr(Y)[0]                       # n×r 正交基
    Pt = sym(W.copy()); P = None
    for _ in range(iters):
        P = sym(Pt - Pt @ Qb @ np.linalg.solve(Qb.T @ Pt @ Qb, Qb.T @ Pt))
        Pn = sym(A @ P @ A.T + W)
        if np.max(np.abs(Pn - Pt)) < tol * max(1.0, np.max(np.abs(Pn))):
            Pt = Pn; break
        Pt = Pn
    if P is None or not np.all(np.isfinite(P)) or np.max(np.abs(P)) > 1e12:
        return np.inf
    return float(np.trace(Th @ P))


def range_onb(S, tol=1e-9):
    w, V = np.linalg.eigh(sym(S))
    keep = [int(i) for i in range(n - 1, -1, -1) if w[i] > tol * w[n - 1]]
    return np.ascontiguousarray(V[:, keep])


def align(Y1, Y2):
    """Y1(dim r1) 的每个方向落在 Y2(dim r2, r1<=r2) 上的最小投影长度 = 最大主角的余弦。"""
    if Y1.shape[1] > Y2.shape[1]:
        Y1, Y2 = Y2, Y1
    s = np.linalg.svd(Y2.T @ Y1, compute_uv=False)
    return float(s[-1]) if len(s) else np.nan


E13 = [(46.423, 3.8180), (48.123, 1.4012), (50.000, 0.7761), (65.000, 0.1083), (80.000, 0.0320)]
DMIN = 46.1231

if __name__ == '__main__':
    print('==== (1) 我的 E13 原始数列与逐对本地指数 ====')
    print('%9s %9s %10s %10s' % ('D', 'D-D_min', 'Δ', '本地 p（与上一点）'))
    xs, ys = [], []
    for D, d in E13:
        x = D - DMIN
        xs.append(np.log(x)); ys.append(np.log(d))
        s = ('%.3f' % ((ys[-2] - ys[-1]) / (xs[-2] - xs[-1]))) if len(xs) > 1 else '—'
        print('%9.3f %9.3f %10.4f %10s' % (D, x, d, s))
    xs, ys = np.array(xs), np.array(ys)
    g = np.polyfit(xs, ys, 1)
    print('  全局 OLS 对数斜率 = %.3f （|残差| 最大 %.2f）' % (g[0], np.max(np.abs(np.polyval(g, xs) - ys))))
    print('  只用最左两点的斜率 = %.3f   ← 我 #5/#6 里的 p≈0.53 就是这么来的'
          % ((ys[0] - ys[1]) / (xs[0] - xs[1])))
    print('  ⇒ C3 §3 成立：不存在单一 p，我那条幂律作废（我认）。')

    print('\n==== (2) range(S*(D)) 与 Φ_min(r) 极小元：D 大时是否转向同一处 ====')
    Vr = {}
    for r in (1, 2, 3):
        try:
            Y = np.load('.work3/c13_V_r%d.npy' % r)
            if Y.shape[0] != n:
                Y = Y.T
            Vr[r] = np.ascontiguousarray(Y)            # 统一成 n×r
        except Exception as e:
            print('  读 c13_V_r%d.npy 失败: %s' % (r, type(e).__name__))
    for r, Y in Vr.items():
        print('  C 的 Φ_min(%d) 极小元：Φ=%.4f (+%.1f%%)，与 Θ 前 %d 主方向的对齐度=%.4f'
              % (r, Phi(Y), 100 * Phi(Y) / unc, r, align(VTh[:, n - r:], Y)))

    print('\n%9s %5s %10s %12s %12s %10s %10s'
          % ('D', 'r*', 'I_unc', 'Φ(range S*)', 'Φ_min(r*)', '对齐度', 'Δ相对差'))
    for D in (50.0, 65.0, 80.0, 120.0, 200.0, 400.0, 1000.0):
        u = unconstrained_sdp(D, np.eye(n), Th, Pc)
        S = np.linalg.inv(u['P']) - np.linalg.inv(u['Pt'])
        Y = range_onb(S)
        r = Y.shape[1]
        ph = Phi(Y)
        pm = Phi(Vr[r]) if r in Vr else np.nan
        al = align(Y, Vr[r]) if r in Vr else np.nan
        print('%9.1f %5d %10.4f %12s %12s %10s %10s'
              % (D, r, u['I_true'], '%.4f' % ph,
                 '%.4f' % pm if np.isfinite(pm) else '—',
                 '%.3f' % al if np.isfinite(al) else '—',
                 '%+.1f%%' % (100 * (ph - pm) / pm) if (np.isfinite(pm) and pm > 1e-6) else '—'))
