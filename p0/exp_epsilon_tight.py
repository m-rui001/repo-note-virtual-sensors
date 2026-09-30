"""E15：ε-扫描的高精度版。E14 的 Δ 只报到 4 位小数，而 face_direct 接受 J ≤ 1.0001·D 的松弛，
   在 D=80 上这相当于 ~2e-5 bits 的乐观偏差——所以"Δ=+0.0000"并不等于 Δ=0。
   这里把约束收紧到 J ≤ D(1+1e-9)，逐 ε 报 Δ 的科学计数位、Φ、以及独立于我的局部求解器的诊断量 ρ。

   ρ 的定义：无约束（Tanaka 式18 SDP）最优 S*(D) 在不可测方向 ker F 上的相对能量
       ρ = ‖ZᵀS*Z‖_F / ‖S*‖_F。
   Cor.2 语言：ρ=0 ⟺ 无约束最优可被 F 实现 ⟺ Δ=0。ρ 由 cvxpy+CLARABEL 给出，
   与我的 SLSQP 完全无关，所以它是这条"Δ 与 Φ 谁先动"结论的独立仲裁者。
"""
import numpy as np
from scipy.linalg import schur
from scipy.optimize import minimize
from p0.p0_replicate_letter import (A, B, W, n, sym, rate_cost, unconstrained_sdp,
                                    init_gamma)
from p0.exp_authoritative import noiseless, ctrl_at
from p0.exp_thm3_qhypothesis import F2, Z2, align_family

np.set_printoptions(precision=4, suppress=True)
r2 = F2.shape[0]
D0 = 80.0


def pack3(G):
    G = sym(G) + 1e-10 * np.eye(r2)
    L = np.linalg.cholesky(G)
    return np.array([L[0, 0], L[1, 0], L[1, 1]])


def unpack3(p):
    L = np.array([[p[0], 0.0], [p[1], p[2]]])
    return sym(L @ L.T)


def ij(F, p, Th, Pc):
    try:
        I, J, _, _ = rate_cost(F, unpack3(p), Th, Pc)
    except Exception:
        return 1e6, 1e6
    if not (np.isfinite(I) and np.isfinite(J)):
        return 1e6, 1e6
    return I, J


def tight_face(F, Th, Pc, D, G0):
    """硬约束 min I s.t. J ≤ D(1+1e-9)，多起点 Cholesky 参数化 SLSQP。"""
    best = None
    cands = [G0, np.eye(r2), np.eye(r2) * 20.0, np.diag([5.0, 0.5])]
    rng = np.random.default_rng(0)
    cands += [sym(v @ v.T) + 0.05 * np.eye(r2) for v in rng.standard_normal((6, r2, r2))]
    for G in cands:
        try:
            p = pack3(G)
        except np.linalg.LinAlgError:
            continue
        res = minimize(lambda q: ij(F, q, Th, Pc)[0], p, method='SLSQP',
                       constraints=[dict(type='ineq', fun=lambda q: D - ij(F, q, Th, Pc)[1])],
                       options=dict(maxiter=400, ftol=1e-13))
        Gam = unpack3(res.x)
        I, J = ij(F, res.x, Th, Pc)
        if J <= D * (1 + 1e-9) and (best is None or I < best[0]):
            best = (I, J)
    return best


if __name__ == '__main__':
    Pc0, K0, _ = ctrl_at(A, B, np.eye(n), np.eye(B.shape[1]))
    print('D=%.1f  下界步长参考：dI/dD ≈ %.2e bits/单位代价 ⟹ Δ 的可分辨下限' % (D0, (1.3419 - 1.3069) / 15))
    print('%-8s %10s %11s %11s %11s %12s %12s' %
          ('ε', 'Φ', 'max|ΘZ|', 'ρ(SDP)', 'I_unc', 'I_TRV(紧)', 'Δ'))
    prev = None
    for eps in (0.0, 1e-4, 1e-2, 3e-2, 0.1, 0.3, 1.0):
        Q = align_family(F2, Z2, eps)
        Pc, K, Th = ctrl_at(A, B, Q, np.eye(B.shape[1]))
        floor = float(np.trace(W @ Pc))
        P, res, it = noiseless(F2, A, W)
        phi = float(np.trace(Th @ P))
        u = unconstrained_sdp(D0, Q, Th, Pc)
        S = np.linalg.inv(u['P']) - np.linalg.inv(u['Pt'])
        rho = float(np.linalg.norm(Z2.T @ S @ Z2) / np.linalg.norm(S))
        f = tight_face(F2, Th, Pc, D0, init_gamma(F2, S))
        if f is None:
            print('%-8g %+10.5f %11.3e %11.3e %11.4f  紧约束求解失败' % (eps, phi, rho, np.max(np.abs(Th @ Z2)), u['I_true']))
            continue
        print('%-8g %+10.5f %11.3e %11.3e %11.4f %12.4f %+12.3e' %
              (eps, phi, np.max(np.abs(Th @ Z2)), rho, u['I_true'], f[0], f[0] - u['I_true']))
