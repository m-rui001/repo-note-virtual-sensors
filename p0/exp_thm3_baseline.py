"""E7：Theorem 3 的"精确归约"到底归约到哪个基线？泄漏项是谁？

自包含（不改 p0 的全局 A,B,W）：所有函数接受 (A,B,W) 参数。
(a) reduced 基线：任务子系统 (A_TT,B_T,W_TT) 上的 Tanaka 无约束 SDP —— Thm.3 声称 TRV 问题等于它；
(b) full 基线：原 4 维系统、Q=FᵀF 的无约束 SDP —— 原文 Eq.(13) 的 Δ 定义用它；
(c) 泄漏实验：把 W 的交叉块 W_TR 置零 / 把 A_RT 置零，看 Δ(full 基线) 是否塌回 0。
"""
import numpy as np
import cvxpy as cp
from scipy.linalg import solve_discrete_are, schur
from scipy.optimize import minimize

ln2 = np.log(2.0)
B4 = np.array([[0.66, -0.58, 0.03, -0.20],
               [2.61, -0.91, 0.87, -0.07],
               [-0.64, -1.12, -0.19, 0.61],
               [0.93, 0.58, -1.18, -1.21]])


def sym(M):
    return 0.5 * (M + M.T)


def ctrl(A, B, Q, R):
    Pc = solve_discrete_are(A, B, sym(Q), sym(R))
    K = np.linalg.solve(sym(R + B.T @ Pc @ B), B.T @ Pc @ A)
    return Pc, K, sym(K.T @ (R + B.T @ Pc @ B) @ K)


def ratecost(A, W, F, Gam, Theta, Pc, R):
    V = np.linalg.inv(Gam)
    Pt = sym(solve_discrete_are(A.T, F.T, W, sym(V)))
    P = sym(np.linalg.inv(np.linalg.inv(Pt) + F.T @ Gam @ F))
    I = 0.5 * np.log(np.linalg.det(Pt) / np.linalg.det(P)) / ln2
    J = float(np.trace(Theta @ P) + np.trace(W @ Pc))
    return I, J, P, Pt


def sdp(A, B, W, Q, R, D):
    """Tanaka 式(18)，任意维；返回 (I_true, S, P, Pc, Theta)"""
    nn = A.shape[0]
    Pc, K, Theta = ctrl(A, B, Q, R)
    P = cp.Variable((nn, nn), symmetric=True)
    Pi = cp.Variable((nn, nn), symmetric=True)
    c0 = 0.5 * np.log(np.linalg.det(W)) / ln2
    c2 = float(np.trace(W @ Pc))
    prob = cp.Problem(cp.Minimize(-0.5 * cp.log_det(Pi) / ln2 + c0),
                      [Pi >> 0, P >> 0, cp.trace(Theta @ P) + c2 <= D,
                       P << A @ P @ A.T + W,
                       cp.bmat([[P - Pi, P @ A.T], [A @ P, A @ P @ A.T + W]]) >> 0])
    prob.solve(solver=cp.CLARABEL, gp=False)
    if P.value is None:
        return None
    Pv, Piv = sym(P.value), sym(Pi.value)
    Pt = A @ Pv @ A.T + W
    S = np.linalg.inv(Pv) - np.linalg.inv(Pt)
    return 0.5 * np.log(np.linalg.det(Pt) / np.linalg.det(Pv)) / ln2, S, Pv, Pc, Theta


def cone_opt(A, B, W, Q, R, D, F, S0, starts=5, budget=1500):
    """min I s.t. J<=D，限制 S=F'ΓF, Γ⪰0（罚函数两段，与 p0 版一致）"""
    Pc, K, Theta = ctrl(A, B, Q, R)
    r = F.shape[0]
    rng = np.random.default_rng(0)

    def unpack(p):
        L = np.zeros((r, r)); k = 0
        for i in range(r):
            for j in range(i + 1):
                L[i, j] = p[k]; k += 1
        L[np.diag_indices(r)] = np.abs(L[np.diag_indices(r)]) + 1e-8
        return sym(L @ L.T)

    def ij(p):
        try:
            I, J, _, _ = ratecost(A, W, F, unpack(p), Theta, Pc, R)
        except Exception:
            return 1e6, 1e6
        return (I, J) if np.isfinite(I) and np.isfinite(J) else (1e6, 1e6)

    def pen(p):
        I, J = ij(p)
        v = max(0.0, (J - D) / max(abs(D), 1e-9))
        return I + 2000.0 * v * v + 200.0 * v

    FFt = np.linalg.inv(F @ F.T)
    G0 = sym(FFt @ (F @ S0 @ F.T) @ FFt) + 1e-6 * np.eye(r)
    cands = [G0, np.eye(r) * 1.0, np.eye(r) * 20.0]
    for _ in range(starts - 3):
        L = rng.standard_normal((r, r)) * 0.7
        cands.append(sym(L @ L.T) + 0.05 * np.eye(r))
    best = None
    for G in cands:
        p = np.array([np.linalg.cholesky(sym(G) + 1e-9 * np.eye(r))[i, j]
                      for i in range(r) for j in range(i + 1)])
        res = minimize(pen, p, method='Nelder-Mead', options=dict(maxfev=budget, fatol=1e-12, xatol=1e-8))
        res = minimize(lambda q: ij(q)[0], res.x, method='SLSQP',
                       constraints=[dict(type='ineq', fun=lambda q: D - ij(q)[1])],
                       options=dict(maxiter=200, ftol=1e-11))
        try:
            I, J, _, _ = ratecost(A, W, F, unpack(res.x), Theta, Pc, R)
        except Exception:
            continue
        if np.isfinite(I) and J <= D * 1.0001 and (best is None or I < best[0]):
            best = (I, J, unpack(res.x))
    return best


def tail_F(A, k):
    T, U, sdim = schur(A, output='real', sort=lambda ar, ai: abs(complex(ar, ai)) < 1.0)
    idx = list(range(A.shape[0] - k, A.shape[0]))
    return U[:, idx].T, U, idx, T


if __name__ == '__main__':
    A4 = np.array([[0.12, 0.63, -0.52, 0.33], [0.26, -1.28, 1.57, 1.13],
                   [-1.77, -0.30, 0.77, 0.25], [-0.16, 0.20, -0.58, 0.56]])
    W4 = sym(np.array([[4.94, -0.10, 1.29, 0.35], [-0.10, 5.55, 2.07, 0.31],
                       [1.29, 2.07, 2.02, 1.43], [0.35, 0.31, 1.43, 3.10]]))
    I4 = np.eye(4)
    k = 2
    F, U, idx, T = tail_F(A4, k)
    rest = [i for i in range(4) if i not in idx]
    Ut = U[:, idx]; Ur = U[:, rest]
    A_TR = T[np.ix_(idx, rest)]; A_RT = T[np.ix_(rest, idx)]
    W_TR = Ut.T @ W4 @ Ur
    print('k=%d 尾部坐标：max|A_TR|=%.2e  max|A_RT|=%.2e  ‖W_TR‖_max=%.2e' %
          (k, np.abs(A_TR).max(), np.abs(A_RT).max(), np.abs(W_TR).max()))
    Qf = sym(F.T @ F)
    Ds = [40.0, 65.0, 80.0]
    print('\n%6s %10s %10s %10s %10s' % ('D', 'I_full-unc', 'I_red-unc', 'I_TRV', 'Δ_full'))
    for D in Ds:
        full = sdp(A4, B4, W4, Qf, np.eye(4), D)
        red = sdp(Ut.T @ A4 @ Ut, Ut.T @ B4, sym(Ut.T @ W4 @ Ut), np.eye(k), np.eye(4), D)
        trv = cone_opt(A4, B4, W4, Qf, np.eye(4), D, F, full[1])
        print('%6.1f %10.4f %10.4f %10s %10s' % (
            D, full[0], red[0], ('%.4f' % trv[0]) if trv else '不可行',
            ('%+.4f' % (trv[0] - full[0])) if trv else '—'))

    print('\n(c) 泄漏实验（D=65，Q=FᵀF）：在 Schur 坐标里逐个置零耦合块，看 Δ_full 是否塌回 0')
    T_W = U.T @ W4 @ U          # 噪声在 Schur 坐标下的分块
    T_A = T.copy()
    def back(TmatA, TmatW):
        return U @ TmatA @ U.T, U @ sym(TmatW) @ U.T
    variants = [('原样', T_A, T_W)]
    a1 = T_A.copy(); a1[np.ix_(rest, idx)] = 0.0                      # A_RT=0
    variants.append(('A_RT=0', a1, T_W))
    w1 = T_W.copy(); w1[np.ix_(idx, rest)] = 0.0; w1[np.ix_(rest, idx)] = 0.0   # W_TR=0
    variants.append(('W_TR=0', T_A, w1))
    variants.append(('A_RT=0 且 W_TR=0', a1, w1))
    D = 65.0
    for tag, Ta, Tw in variants:
        Aa, Ww = back(Ta, Tw)
        Qa = sym(F.T @ F)
        full = sdp(Aa, B4, Ww, Qa, np.eye(4), D)
        trv = cone_opt(Aa, B4, Ww, Qa, np.eye(4), D, F, full[1])
        red = sdp(Ut.T @ Aa @ Ut, Ut.T @ B4, sym(Ut.T @ Ww @ Ut), np.eye(k), np.eye(4), D)
        print('  %-16s I_full=%7.4f  I_red=%7.4f  I_TRV=%s  Δ_full=%s  Δ_vs_red=%s' % (
            tag, full[0], red[0],
            ('%7.4f' % trv[0]) if trv else '  不可行',
            ('%+7.4f' % (trv[0] - full[0])) if trv else '   —  ',
            ('%+7.4f' % (trv[0] - red[0])) if trv else '   —  '))
