"""E8（2026-09-28 第 3 轮，B 侧）：Φ(F) 的几何刻画 + 秩 1 情形的精确 Δ(D) 曲线。

主张 1（可判定预测）：Φ(F) 只依赖 (A,F) 的最大不可观测块 V* = ∩ ker(F A^k)，
    与 F 的秩、与 F 在其他方向上的取向都无关。
    闭式：取 V* 的正交基 Z，X = dlyap(A|V*, W|V*)，P_∞ = Z X Zᵀ，Φ = tr(Θ P_∞)。
    推论：存在秩 1 的 F 使 (A,F) 可观测 ⟹ Φ=0 ⟹ L-CSS 的 Prop.3（满秩才无损）远非必要条件。

主张 2（秩 1 锥问题精确可解）：rank F = 1 时 Γ 只有一个标量自由度，
    I_TRV(D) = I_dir(γ*(D))，γ* 由 J(γ)=D 唯一决定（J 随 γ 单调降、I_dir 随 γ 单调升）。
    于是 Δ(D) 不需要 Algorithm 1 的 SCA，可以逐点精确算，直接看 D→地板 端的行为。
"""
import numpy as np
from scipy.linalg import solve_discrete_lyapunov
from p0.p0_replicate_letter import (A, B, W, n, sym, ln2, ctrl, rate_cost,
                                    unconstrained_sdp, ss_filter)

np.set_printoptions(precision=4, suppress=True)
Pc, K, Theta = ctrl(np.eye(n))
unc_floor = float(np.trace(W @ Pc))
TH = sym(Theta)
THsq = sym(np.linalg.eigh(TH)[1] @ np.diag(np.sqrt(np.maximum(np.linalg.eigh(TH)[0], 0)))
           @ np.linalg.eigh(TH)[1].T)


def blind_spot(F, tol=1e-8):
    """V* = 最大的、含于 ker F 的 A-不变子空间：可观测性矩阵的零空间"""
    r = F.shape[0]
    O = np.vstack([F @ np.linalg.matrix_power(A, k) for k in range(n)])
    _, s, Vt = np.linalg.svd(O, full_matrices=True)
    k_obs = int((s > tol * s[0]).sum())
    Z = Vt[k_obs:].T                      # n × (n-k_obs)，V* 的正交基
    return Z, n - k_obs, np.linalg.matrix_rank(O) >= n


def Phi_geom(F):
    Z, dim, obs = blind_spot(F)
    if dim == 0:
        return 0.0, dim, obs
    Av = Z.T @ A @ Z
    Wv = sym(Z.T @ W @ Z)
    try:
        X = sym(solve_discrete_lyapunov(Av, Wv))
    except Exception:
        return np.inf, dim, obs
    return float(np.trace(TH @ (Z @ X @ Z.T))), dim, obs


def Phi_num(F, prec=1e12):
    _, J, _, _ = rate_cost(F, prec * np.eye(F.shape[0]), Theta, Pc)
    return J - unc_floor


def curve_rank1(F, gmin=1.0, gmax=1e14, npts=400):
    """秩 1 锥问题的精确前沿：γ 网格上 (J(γ), I_dir(γ))"""
    gs = np.logspace(np.log10(gmin), np.log10(gmax), npts)
    pts = []
    for g in gs:
        try:
            I, J = rate_cost(F, np.array([[g]]), Theta, Pc)[:2]
        except Exception:
            continue
        if np.isfinite(I) and np.isfinite(J):
            pts.append((g, J, I))
    return np.array(pts)


def I_trv_rank1(pts, D):
    """J ≤ D 的最小 I_dir：单调性下即 J=D 的交点"""
    Js, Is = pts[:, 1], pts[:, 2]
    ok = Js <= D
    if not ok.any():
        return np.nan
    return float(Is[ok].min())


def rand_F(r, seed):
    rng = np.random.default_rng(seed)
    F = rng.standard_normal((r, n))
    return F - F @ blind_spot(np.eye(n))[0] @ blind_spot(np.eye(n))[0].T if False else F


def show_F(tag, F, Ds):
    Z, dim, obs = blind_spot(F)
    pg, pd, po = Phi_geom(F)
    pn = Phi_num(F) if F.shape[0] < n else 0.0
    print('%-26s rank=%d  dim V*=%d  可观测=%s  Φ_几何=%+9.4f  Φ_数值=%+9.4f' %
          (tag, np.linalg.matrix_rank(F), dim, po, pg, pn))
    if F.shape[0] == 1:
        pts = curve_rank1(F)
        floor = unc_floor + pg
        print('    锥地板 D_min(K_F)=%.4f；γ 网格 J 最小值=%.4f' % (floor, pts[:, 1].min()))
        for D in Ds:
            u = unconstrained_sdp(D, np.eye(n), Theta, Pc)
            it = I_trv_rank1(pts, D)
            d = it - u['I_true'] if np.isfinite(it) else np.nan
            print('    D=%8.4f  I_unc=%8.4f  I_TRV=%9.4f  Δ=%+9.4f' % (D, u['I_true'], it, d))


if __name__ == '__main__':
    print('D_min^unc = tr(W Pc) = %.4f\n' % unc_floor)

    print('---- 主张 1：Φ 只依赖 V*（几何式 vs Γ→∞ 数值式对账）')
    from p0.exp_trv_penalty import schur_task
    F2, _ = schur_task(2)
    F3, c3 = schur_task(3)
    cases = [('Schur 尾部 k=2', F2), ('Schur 尾部 k=3', F3),
             ('单坐标 e1', np.array([[1., 0, 0, 0]])),
             ('单坐标 e3', np.array([[0., 0, 1., 0]])),
             ('单坐标 e4', np.array([[0., 0, 0, 1.]]))]
    for i, s in enumerate([1.3, -0.7, 2.1, 0.4, -1.9]):
        cases.append(('随机秩1 #%d' % i, np.array([[1., s, 0.3 * s, s * s - 1]])))
    for tag, F in cases:
        Z, dim, obs = blind_spot(F)
        pg, _, _ = Phi_geom(F)
        print('  %-14s dim V*=%d 可观测=%-5s Φ_几何=%+10.4f  Φ_数值=%+10.4f  差=%.2e'
              % (tag, dim, obs, pg, Phi_num(F), abs(pg - Phi_num(F))))

    print('\n---- 同一 V* 的不同 F 给同一 Φ（秩 2 例）')
    Fs = [np.array([[1., 0., 0., 0.], [0., 1., 0., 0.]]),
          np.array([[1., 0., 0., 0.], [0., 1., 0.2, 0.]]),
          np.array([[1., 1., 0., 0.], [0., 1., 0., 0.]])]
    for F in Fs:
        Z, dim, obs = blind_spot(F)
        print('  dim V*=%d 可观测=%-5s Φ_几何=%+9.4f Φ_数值=%+9.4f' % (dim, obs, Phi_geom(F)[0], Phi_num(F)))

    print('\n---- 主张 2：秩 1 可观测 F 的精确 Δ(D)（Φ=0 而 Δ≠0 的角）')
    for tag, F in [('e1(不可观测块非平凡)', np.array([[1., 0, 0, 0]])),
                   ('随机可观测秩1', np.array([[1., 1.3, 0.39, 0.49]]))]:
        Z, dim, obs = blind_spot(F)
        floor = unc_floor + Phi_geom(F)[0]
        Ds = [floor + e for e in (1e-3, 1e-2, 0.1, 0.5, 1.0, 2.0, 5.0)] + [40.0, 80.0]
        show_F('%s 可观测=%s' % (tag, obs), F, Ds)
