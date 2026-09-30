"""复现 Tanaka TAC2018 (arXiv:1510.04214) 的最小有向信息 LQG 曲线，
并与 IEEE L-CSS 2026 (Task-Restricted Sensing) 的表示惩罚 Δ(D) 对齐。

两个独立求解器：
  A) unconstrained_sdp(D, Q)  —— Tanaka 式(18) 的 max-det SDP（cvxpy/CLARABEL）
  B) face_direct(D, F, Q)     —— 在 L-CSS 的锥 K_F={F'ΓF} 上做定点迭代 + 约束优化（全局多起点）
校验锚点（原文 [page 6]）：D=33→6.133, D=40→3.266, D=80→1.602 bits/sample；
稳定化下界 Σlog|λ_u| = 1.169 bits/sample（本文件验算 1.1690）。
"""
import numpy as np
from scipy.linalg import solve_discrete_are, logm
import cvxpy as cp

ln2 = np.log(2.0)

A = np.array([[0.12, 0.63, -0.52, 0.33],
              [0.26, -1.28, 1.57, 1.13],
              [-1.77, -0.30, 0.77, 0.25],
              [-0.16, 0.20, -0.58, 0.56]])
B = np.array([[0.66, -0.58, 0.03, -0.20],
              [2.61, -0.91, 0.87, -0.07],
              [-0.64, -1.12, -0.19, 0.61],
              [0.93, 0.58, -1.18, -1.21]])
W = np.array([[4.94, -0.10, 1.29, 0.35],
              [-0.10, 5.55, 2.07, 0.31],
              [1.29, 2.07, 2.02, 1.43],
              [0.35, 0.31, 1.43, 3.10]])
W = 0.5 * (W + W.T)
n = 4


def sym(M):
    return 0.5 * (M + M.T)


def plant_facts():
    ev = np.linalg.eigvals(A)
    u = np.abs(ev[np.abs(ev) >= 1.0])
    return ev, float(u.size), float(np.sum(np.log(u)) / ln2)


def ctrl(Qmat, Rmat=np.eye(n)):
    """Riccati 解 Pc、增益 K、Θ = K'(R + B'Pc B)K   （L-CSS Eq.16 / Tanaka Eq.11d）"""
    Pc = solve_discrete_are(A, B, sym(Qmat), sym(Rmat))
    K = solve_(Rmat + B.T @ Pc @ B, B.T @ Pc @ A)
    Theta = K.T @ (Rmat + B.T @ Pc @ B) @ K
    return Pc, K, sym(Theta)


def solve_(M, rhs):
    return np.linalg.solve(sym(M) + 1e-12 * np.eye(M.shape[0]), rhs)


def posterior(Pt, F, Gam):
    """给定先验 Pt 与精度 S=F'ΓF，返回后验 P=(Pt^{-1}+S)^{-1}"""
    S = F.T @ Gam @ F
    return np.linalg.inv(np.linalg.inv(Pt) + S)


def ss_filter(F, Gam):
    """稳态信息滤波。scipy 的 solve_discrete_are(a,b,q,r) 解 X=a'Xa−a'Xb(b'Xb+r)^{-1}b'Xa+q，
    代入 (A',F',W,V) 得到的 X 是**预报**协方差 P̃ 而非后验 P，需再做一次 information-form 更新。
    （早期版本把 X 当成 P，Γ 从 1e2 增至 1e8 时 I_dir 只从 3.61 变到 3.67 bits —— 已用
    20000 步定点迭代交叉校验修正。）"""
    V = np.linalg.inv(Gam)
    Pt = sym(solve_discrete_are(A.T, F.T, W, sym(V)))
    P = sym(np.linalg.inv(np.linalg.inv(Pt) + F.T @ Gam @ F))
    return P, Pt


def rate_cost(F, Gam, Theta, Pc):
    """L-CSS Eq.(10)(11)：有向信息率与代价（bits/sample）"""
    P, Pt = ss_filter(F, Gam)
    S = F.T @ Gam @ F
    I_dir = 0.5 * np.log(np.linalg.det(Pt) / np.linalg.det(P)) / ln2
    J = float(np.trace(Theta @ P) + np.trace(W @ Pc))
    return I_dir, J, P, Pt


def chol_to_psd(p, r):
    L = np.zeros((r, r))
    idx = 0
    for i in range(r):
        for j in range(i + 1):
            L[i, j] = p[idx]
            idx += 1
    L[np.diag_indices(r)] = np.abs(L[np.diag_indices(r)]) + 1e-6
    return L @ L.T


def init_gamma(F, S_ref):
    """把无约束最优信息矩阵 S_ref 压到锥 K_F 上，作为 TRV 求解的初值：
    Γ = (FF')^-1 · F S_ref F' · (FF')^-T   （即 range(Fᵀ) 上的分块）"""
    FFt_inv = np.linalg.inv(F @ F.T)
    G = FFt_inv @ (F @ S_ref @ F.T) @ FFt_inv
    return sym(G) + 1e-6 * np.eye(F.shape[0])


def face_direct(D, F, Theta, Pc, gam0=None, n_starts=2, seed=0, warm=None, budget=1200):
    """在锥 K_F 上求 min I_dir s.t. J<=D。
    两段式：先罚函数 Nelder-Mead 保证走到可行域（相对违背量做尺度无关的罚），
    再用 SLSQP 硬约束精修。初值优先 warm（沿 D 网格热启动）与 gam0（锥投影初值）。"""
    from scipy.optimize import minimize
    r = F.shape[0]

    def unpack(p):
        L = np.zeros((r, r))
        idx = 0
        for i in range(r):
            for j in range(i + 1):
                L[i, j] = p[idx]
                idx += 1
        L[np.diag_indices(r)] = np.abs(L[np.diag_indices(r)]) + 1e-8
        return sym(L @ L.T)

    def pack(G):
        G = sym(G) + 1e-9 * np.eye(r)
        L = np.linalg.cholesky(G)
        return np.array([L[i, j] for i in range(r) for j in range(i + 1)])

    def ij(p):
        try:
            I, J, _, _ = rate_cost(F, unpack(p), Theta, Pc)
        except Exception:
            return 1e6, 1e6
        if not np.isfinite(I) or not np.isfinite(J):
            return 1e6, 1e6
        return I, J

    def pen(p):
        I, J = ij(p)
        v = max(0.0, (J - D) / max(abs(D), 1e-9))
        return I + 2000.0 * v * v + 200.0 * v

    starts = [G for G in (warm, gam0) if G is not None]
    starts += [np.eye(r) * s for s in (1.0, 20.0)]
    rng = np.random.default_rng(seed)
    for _ in range(n_starts):
        L = rng.standard_normal((r, r)) * 0.7
        starts.append(sym(L @ L.T) + 0.05 * np.eye(r))

    best = None
    for G0 in starts:
        try:
            p = pack(G0)
        except np.linalg.LinAlgError:
            continue
        res = minimize(pen, p, method='Nelder-Mead',
                       options=dict(maxfev=budget, xatol=1e-7, fatol=1e-11))
        res = minimize(lambda q: ij(q)[0], res.x, method='SLSQP',
                       constraints=[dict(type='ineq', fun=lambda q: D - ij(q)[1])],
                       options=dict(maxiter=150, ftol=1e-11))
        Gam = unpack(res.x)
        try:
            I, J, _, _ = rate_cost(F, Gam, Theta, Pc)
        except Exception:
            continue
        if np.isfinite(I) and np.isfinite(J) and J <= D * 1.0001 and (best is None or I < best[0]):
            best = (I, J, Gam)
    return best


def unconstrained_sdp(D, Qmat, Theta, Pc, solver=cp.CLARABEL):
    """Tanaka 式(18)：变量 (P, Pi)，返回 SDP 目标值、½logdet(P̃/P)、 recovered S"""
    P = cp.Variable((n, n), symmetric=True)
    Pi = cp.Variable((n, n), symmetric=True)
    const = 0.5 * np.log(np.linalg.det(W)) / ln2
    obj = cp.Minimize(-0.5 * cp.log_det(Pi) / ln2 + const)
    cons = [Pi >> 0, P >> 0,
            cp.trace(Theta @ P) + const2(Qmat, Pc) <= D,
            P << A @ P @ A.T + W,
            cp.bmat([[P - Pi, P @ A.T], [A @ P, A @ P @ A.T + W]]) >> 0]
    prob = cp.Problem(obj, cons)
    prob.solve(solver=solver, gp=False)
    if P.value is None:
        return None
    Pv = sym(P.value)
    Piv = sym(Pi.value)
    Pt = A @ Pv @ A.T + W
    S = np.linalg.inv(Pv) - np.linalg.inv(Pt)
    I_true = 0.5 * np.log(np.linalg.det(Pt) / np.linalg.det(Pv)) / ln2
    return dict(sdp_obj=float(prob.value), I_true=I_true, P=Pv, Pt=Pt,
                rank_S=int(np.sum(np.linalg.eigvalsh(sym(S)) > 1e-6 * np.abs(np.linalg.eigvalsh(sym(S))).max())),
                feasible=prob.status == 'optimal')


def const2(Qmat, Pc):
    return float(np.trace(W @ Pc))


def schur_F(nu, keep=()):
    """L-CSS §V：A=UTU' 实 Schur，不稳定模态排在末尾；F 选取对应 Schur 坐标"""
    T, U = np.linalg.schur(A, 'real')
    ev = np.linalg.eigvals(T)
    order = np.argsort(np.abs(ev))          # 稳定在前，不稳定在后
    T = T[np.ix_(order, order)]
    U = U[:, order]
    rows = []
    idx = list(range(n - nu, n)) + [order.index(k) for k in keep] if keep else list(range(n - nu, n))
    F = U[:, idx].T
    return sym(U.T @ A @ U), F, U


if __name__ == '__main__':
    ev, nu, ib = plant_facts()
    print('eig(A)=', np.round(ev, 4))
    print('unstable modes ν=%d,  stabilization lower bound = %.4f bits/sample (原文 1.169)' % (nu, ib))
    Pc, K, Theta = ctrl(np.eye(n))
    print('Tr(W Pc)=%.4f (垂直渐近线 D_min)' % float(np.trace(W @ Pc)))
    for D, anchor in ((33, 6.133), (40, 3.266), (80, 1.602)):
        out = unconstrained_sdp(D, np.eye(n), Theta, Pc)
        if out:
            print('D=%3d  SDP obj=%.4f  I_true=%.4f  rank(S)=%d  | 原文 %.3f' %
                  (D, out['sdp_obj'], out['I_true'], out['rank_S'], anchor))
