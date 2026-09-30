"""E32：E29 的修复值必须独立复核——解析点 $V^\\ast$ 到底可不可行。

E29 报"解析点在 $x\\le0.05$ 与原数值链重合到 $10^{-4}$ bit，深处低 $0.6\\sim1.0$ bit"。这条话的重量全压在
一件事上：$J(V^\\ast)=D_{\\min}(\\mathcal K_F)+x$ 是真的，不是我的定点迭代自说自话。三条独立路径同时算 $(I,J)$：
  (1) 我自己的信息形式定点迭代（E29/E31 用的那个）；
  (2) `p0_replicate_letter.rate_cost`：scipy 的 `solve_discrete_are` 解预报协方差，完全不同的实现；
  (3) `face_direct` 的锥上 SDP：把 $S=F^\\top\\Gamma F$ 摆进 L-CSS 的锥坐标，用 cvxpy/CLARABEL 求解，
      看它给的最优值是否 $\\le$ 解析点（若 SDP 也找不到更低值，而我的旧 SLSQP 链更高，那旧链就是垃圾）。
另外单独测：把解析点 $\\Gamma^\\ast$ 原样喂给 `solve_cone`，它在 $x=0.005$ 处交回 $10.58$（比解析点差 $0.92$），
说明它连自己的可行起点都没守住——查清是 `pack/unpack` 的精度损失还是 $10^{-9}$ 的可行性容差，这决定我后面
还要不要用这条链。
"""
import sys
import numpy as np
sys.stdout.reconfigure(encoding='utf-8')
from p0.p0_replicate_letter import A, W, n, sym, ln2, ctrl, rate_cost, unconstrained_sdp
from p0.exp_authoritative import noiseless
from p0.exp_nearfloor_law import Fs, floor_of, Th, Pc, solve_cone, pack, unpack

UNC = float(np.trace(W @ Pc))
P0 = [None]


def fp_solve(F, V, iters=60000, tol=1e-16):
    P = P0[0][0].copy()
    for k in range(iters):
        Pt = sym(A @ P @ A.T + W)
        Pn = sym(Pt - Pt @ F.T @ np.linalg.solve(F @ Pt @ F.T + V, F @ Pt))
        if k > 300 and np.max(np.abs(Pn - P)) < tol:
            P = Pn
            break
        P = Pn
    return P, sym(A @ P @ A.T + W)


def calib(F, Ci, x_target, fl):
    lo, hi = 1e-16, 1e2
    for _ in range(70):
        mid = np.sqrt(lo * hi)
        P, _ = fp_solve(F, mid * Ci)
        if float(np.trace(Th @ (P - P0[0][0]))) > x_target:
            hi = mid
        else:
            lo = mid
    a = np.sqrt(lo * hi)
    P, Pt = fp_solve(F, a * Ci)
    V = sym(a * Ci)
    I_fp = 0.5 * np.log(np.linalg.det(F @ Pt @ F.T + V) / np.linalg.det(V)) / ln2
    return V, P, Pt, I_fp, float(np.trace(Th @ P)) + UNC


if __name__ == '__main__':
    F = Fs['F2 尾部2(封闭, r=2)']
    r = F.shape[0]
    fl, phi, _ = floor_of(F)
    Pn, res, it = noiseless(F, A, W, iters=80000, tol=1e-15)
    P0[0] = (Pn, sym(A @ Pn @ A.T + W))
    C = np.load('p0/fig/Cmat_F2.npy')
    Ci = np.linalg.inv(C)
    print('F2 地板 %.4f，检验解析点可行性与三条独立路径\n' % fl)
    print('%8s %12s %12s %12s %12s %12s %12s'
          % ('x', 'J(定点)', 'J(DARE)', 'J−地板−x', 'I(定点)', 'I(DARE)', 'solve_cone'))
    for x in (0.05, 0.02, 0.01, 0.005, 0.002, 0.001):
        V, P, Pt, I_fp, J_fp = calib(F, Ci, x, fl)
        G = np.linalg.inv(V)
        I_dare, J_dare, _, _ = rate_cost(F, sym(G), Th, Pc)
        sc = solve_cone(F, fl + x, fl, [pack(sym(G))])
        print('%8.4f %12.5f %12.5f %12.1e %12.5f %12.5f %12s'
              % (x, J_fp, J_dare, J_fp - fl - x, I_fp, I_dare,
                 '%+.4f' % (sc[0] - I_fp) if sc else '失败'))

    print('\n精度诊断：把 $\\Gamma^\\ast$ 走一遍 pack→unpack，看参数化损失多少相对误差')
    V, P, Pt, I_fp, J_fp = calib(F, Ci, 0.005, fl)
    G = np.linalg.inv(V)
    G2 = unpack(pack(G))
    print('  $\\kappa(\\Gamma^\\ast)$=%.3e，$\\|G-G_2\\|_F/\\|G\\|_F$=%.3e'
          % (np.linalg.cond(G), np.linalg.norm(G - G2) / np.linalg.norm(G)))
    _, J2, _, _ = rate_cost(F, G2, Th, Pc)
    print('  重建后的 $J$ 比预算超 %+.3e（相对 %+.3e）；`solve_cone` 的容差是 $J\\le D(1+10^{-9})$'
          % (J2 - (fl + 0.005), (J2 - (fl + 0.005)) / (fl + 0.005)))
