"""E9（2026-09-28 第 3 轮，B 侧）：Φ(F) 的真身 = 无噪声（V=0）滤波 DARE 的定点，
   不是我在 entry #3 里用的 "Γ=1e10 代理"。E8 证伪了我的几何猜想，这里把量修正。

事实链：
  (a) 量测噪声 V→0 时，后验必然满足 range(P) ⊆ ker F（每步更新后 Fe=0），故 rank P_∞ ≤ n−r；
      即便 (A,F) 可观测，x_t 也只由 (Fx_t, Fx_{t+1},…) 的**未来**组合决定，一步延迟不可消除
      ⟹ Φ(F)=0 ⟸ 只能 r=n（与 L-CSS Prop.3 一致；E8 的"秩1可观测⟹零地板"是错的）。
  (b) 精确算法：在预测协方差上迭代 P̃ ← A(I−KF)P̃Aᵀ + W，K = P̃Fᵀ(FP̃Fᵀ)⁻¹，收敛到定点。
      与 Γ=1e10 的直接求逆相比无条件数问题（后者会给出负的 J，已被 E8 抓到）。
  (c) 对账：Φ_定点 vs Φ_代理(Γ=1e4/1e6/1e8/1e10)，看代理是否单调收敛到定点；不收敛则 entry #3 的数作废。
"""
import numpy as np
from p0.p0_replicate_letter import (A, B, W, n, sym, ln2, ctrl, rate_cost, ss_filter,
                                    unconstrained_sdp)
from p0.exp_trv_penalty import schur_task, cone_resid

np.set_printoptions(precision=4, suppress=True)
Pc, K, Theta = ctrl(np.eye(n))
unc_floor = float(np.trace(W @ Pc))


def noiseless_posterior(F, iters=200000, tol=1e-14):
    """V=0 稳态 Kalman：迭代预测协方差的投影递归，返回 (后验 P, 定点残差, 实际步数)"""
    FF = F @ F.T
    Pt = sym(W.copy())
    res = np.inf
    for t in range(iters):
        P = sym(Pt - Pt @ F.T @ np.linalg.solve(F @ Pt @ F.T + 1e-300 * np.eye(F.shape[0]), F @ Pt))
        Pt_new = sym(A @ P @ A.T + W)
        res = np.max(np.abs(Pt_new - Pt))
        Pt = Pt_new
        if res < tol:
            break
    P = sym(Pt - Pt @ F.T @ np.linalg.solve(F @ Pt @ F.T, F @ Pt))
    return P, res, t + 1


def phi_exact(F):
    P, res, t = noiseless_posterior(F)
    return float(np.trace(Theta @ P)), P, res, t


def phi_proxy(F, precs=(1e4, 1e6, 1e8, 1e10)):
    out = []
    for p in precs:
        try:
            _, J, P, Pt = rate_cost(F, p * np.eye(F.shape[0]), Theta, Pc)
        except Exception:
            J = np.nan
        out.append(J - unc_floor if np.isfinite(J) else np.nan)
    return out


def report(tag, F):
    if F.shape[0] >= n:
        print('%-24s 满秩：Φ=0（Prop.3）' % tag)
        return
    ph, P, res, t = phi_exact(F)
    Z, sV = np.linalg.svd(P, hermitian=True)[:2]
    rk = int((sV > 1e-9 * sV[0]).sum())
    # range(P) ⊆ ker F ?
    FFt = F @ F.T
    proj_err = np.max(np.abs(F @ (Z[:, :rk] * sV[:rk]) @ Z[:, :rk].T)) if rk else 0.0
    prox = phi_proxy(F)
    print('%-24s rankF=%d  Φ_定点=%+10.4f  Φ_代理(1e4..1e10)=%s' % (tag, np.linalg.matrix_rank(F), ph,
          ' '.join('%+.3f' % x for x in prox)))
    print('      定点残差=%.1e@%d步  rank P_∞=%d(=n−r=%d)  max|F P_∞|=%.2e  tr(WP_c)=%.4f' %
          (res, t, rk, n - np.linalg.matrix_rank(F), proj_err, unc_floor))
    return ph


if __name__ == '__main__':
    print('D_min^unc = %.4f\n' % unc_floor)
    F2, _ = schur_task(2)
    F3, _ = schur_task(3)
    print('---- (b)(c) 定点 vs 代理：entry #3 报过的所有 F')
    rows = []
    for tag, F in [('Schur 尾部 k=2', F2), ('Schur 尾部 k=3', F3),
                   ('e1,e2', np.array([[1., 0, 0, 0], [0, 1., 0, 0]])),
                   ('e1,e3', np.array([[1., 0, 0, 0], [0, 0, 1., 0]])),
                   ('e1,e4', np.array([[1., 0, 0, 0], [0, 0, 0, 1.]])),
                   ('单坐标 e1', np.array([[1., 0, 0, 0]])),
                   ('单坐标 e3', np.array([[0., 0, 1., 0]])),
                   ('随机秩2 #1', np.array([[0.8, -0.3, 0.5, 0.2], [0.1, 0.9, -0.4, 0.3]])),
                   ('随机秩2 #2', np.array([[1., 2., -1., 0.5], [0.2, 0.1, 1., -0.7]])),
                   ('随机秩3', np.random.default_rng(7).standard_normal((3, n)))]:
        rows.append((tag, F, report(tag, F)))
    print('\n  全状态 I：', end='')
    report('I (r=n)', np.eye(n))

    print('\n---- 同一 V* / 同一 ker F 的两个 F 是否同 Φ（E8 猜想的正确版本：Φ 只依赖 ker F？）')
    # ker 相同：F 与 R·F（R 可逆）ker 完全相同
    F = np.array([[1., 0., 0., 0.], [0, 1., 0., 0.]])
    R = np.array([[2., 1.], [-0.5, 3.]])
    print('  F 与 R·F（同一行空间/同一 ker）: Φ=%.6f vs %.6f' %
          (phi_exact(F)[0], phi_exact(R @ F)[0]))
    Fb = np.array([[1., 0., 0., 0.], [0.3, 1., 0., 0.]])
    print('  同行空间第三种基: Φ=%.6f' % phi_exact(Fb)[0])

    print('\n---- 结论性检查：Φ_定点 是否 ≥ 0 且严格递增地随 r 下降而增大（同 ker 嵌套）')
    for r in (3, 2, 1):
        Fr = np.eye(n)[:r]
        print('  前 %d 个坐标：Φ=%+.5f' % (r, phi_exact(Fr)[0]))
