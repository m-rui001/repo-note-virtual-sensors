"""E11：机制判决——Δ(率/表示损失) 由 Schur 块 A_TR 控制，Φ(代价地板) 由 A_RT 控制。

全部在实 Schur 坐标里做：被控对象 T，噪声 UᵀWU，量测 F=[0 I_ν]（尾部 ν 个坐标）。
四个变体只动非对角块：
  (1) 原始 T：A_TR=0, A_RT≠0  → 预测 Δ≡0（原文 Thm.3）、Φ>0。
  (2) 置零 A_RT：两个非对角块都 0（子系统解耦）→ 预测 Φ=0。
  (3) 人为加出 A_TR（保持 A_RT≠0）→ 预测 Φ 基本不变（地板由 A_RT 决定），而 Δ>0。
  (4) 只置零 A_RT、保留人造 A_TR → 预测 Φ=0 而 Δ>0。
四个角填满 ⟹ Δ 与 Φ 是两个独立开关，分别由 T 的两个非对角块控制。
"""
import numpy as np
from scipy.linalg import schur, solve_discrete_are
from p0.p0_replicate_letter import A, B, W, n, sym, solve_

np.set_printoptions(precision=4, suppress=True)
T, U, sdim = schur(A, output='real', sort=lambda a: np.abs(a) < 1.0)
nu = n - sdim
iT, iR = list(range(n - nu, n)), list(range(n - nu))
F = np.zeros((nu, n)); F[:, iR:] = np.eye(nu)
Ws, Bs = sym(U.T @ W @ U), U.T @ B
print('Schur 排序后 |λ| =', np.round(np.abs(np.linalg.eigvals(T)), 4))
print('max|A_TR| = %.4f   max|A_RT| = %.4f' % (np.max(np.abs(T[np.ix_(iT, iR)])),
                                                np.max(np.abs(T[np.ix_(iR, iT)]))))


def ctrl_at(Am, Bm, Qm, Rm=None):
    Rm = np.eye(Bm.shape[1]) if Rm is None else Rm
    Pc = solve_discrete_are(Am, Bm, sym(Qm), sym(Rm))
    K = solve_(Rm + Bm.T @ Pc @ Bm, Bm.T @ Pc @ Am)
    return Pc, K, sym(K.T @ (Rm + Bm.T @ Pc @ Bm) @ K)


def noiseless(Am, Wm):
    Pt = sym(Wm.copy()); res = np.inf
    for k in range(200000):
        M = F @ Pt @ F.T + 1e-300 * np.eye(nu)
        P = sym(Pt - Pt @ F.T @ np.linalg.solve(M, F @ Pt))
        Pn = sym(Am @ P @ Am.T + Wm)
        res = np.max(np.abs(Pn - Pt)); Pt = Pn
        if res < 1e-14:
            break
    M = F @ Pt @ F.T
    return sym(Pt - Pt @ F.T @ np.linalg.solve(M + 1e-300 * np.eye(nu), F @ Pt)), res, k + 1


def report(tag, Tm, Qm):
    try:
        Pc, K, Th = ctrl_at(Tm, Bs, Qm)
    except Exception as e:
        print('  %-40s LQR 不可解: %s' % (tag, e)); return
    unc = float(np.trace(Ws @ Pc))
    P, res, it = noiseless(Tm, Ws)
    phi = float(np.trace(Th @ P))
    rk = np.linalg.matrix_rank(P, tol=1e-9 * max(np.linalg.eigvalsh(P).max(), 1e-30))
    print('  %-40s Φ=%+9.4f  相对=%+7.1f%%  D_min(K_F)=%9.4f  rankP∞=%d  res=%.0e@%d' %
          (tag, phi, 100 * phi / unc, unc + phi, rk, res, it))


if __name__ == '__main__':
    Qalign = sym(F.T @ F)
    TRv = T[np.ix_(iT, iR)].copy(); RTv = T[np.ix_(iR, iT)].copy()
    T_noRT = T.copy(); T_noRT[np.ix_(iR, iT)] = 0.0
    T_addTR = T.copy(); T_addTR[np.ix_(iT, iR)] = 0.7 * np.ones((nu, n - nu))
    T_both = T_addTR.copy(); T_both[np.ix_(iR, iT)] = 0.0

    for Qtag, Qm in (('Q = FᵀF（原文 Eq.(4) 的 natural choice）', Qalign), ('Q = I', np.eye(n))):
        print('\n---- %s' % Qtag)
        report('(1) 原始：A_TR=0, A_RT≠0', T, Qm)
        report('(2) 置零 A_RT（两非对角块全 0）', T_noRT, Qm)
        report('(3) 人造 A_TR=0.7（A_RT 保留）', T_addTR, Qm)
        report('(4) 人造 A_TR=0.7 且 A_RT=0', T_both, Qm)

    print('\n---- 机制读法：Φ 只随 A_RT 变（(1)vs(2)、(3)vs(4)），Δ 只随 A_TR 变')
    print('  A_TR 影响的是"任务子系统能否自我复现"（率损失），A_RT 影响的是"看不见的坐标在驱动任务"（代价地板）')
