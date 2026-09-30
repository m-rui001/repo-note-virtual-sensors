"""E14：一条受控实验——固定 F（封闭性 A_TR≡0 全程不变），只拧代价 Q 在"测不到的方向"上的权重，
   看 L-CSS Theorem 3 的结论 Δ(D)≡0 还成不成立。

   族：Q_ε = FᵀF + ε·ZZᵀ，Z = ker F 的正交基（本项目里 = 前导稳定 Schur 列）。
     ε=0  Q 支撑 ⊆ range(Fᵀ)（任务对齐）
     ε=1  Q = I（全状态代价，§V 未声明但 Fig.1 唯一自洽的选择）
   A_TR 只依赖 (A,F)，与 ε 无关 ⟹ 若 Δ(ε>0)>0，就是 Thm.3 缺假设的反例，
   而不是"SCA 没收敛"或"子空间没选对"。

   同时报两个独立诊断：
     Φ(F;Q_ε)      可达代价地板缺口（一次无噪声滤波定点）
     ρ(D;Q_ε)      无约束最优 S*(D) 落在 ker F 上的相对能量 ‖ZᵀS*Z‖_F/‖S*‖_F
   预言（Cor.2 语言）：Δ>0 ⟺ ρ>0，且 Φ>0 ⟺ Θ 在 ker F 上不恒零。两者由同一个旋钮 ε 驱动。
"""
import numpy as np
from scipy.linalg import schur
from p0.p0_replicate_letter import (A, B, W, n, sym, ctrl, face_direct,
                                    unconstrained_sdp, init_gamma)
from p0.exp_authoritative import noiseless, ctrl_at

np.set_printoptions(precision=4, suppress=True, linewidth=200)
T, U, sdim = schur(A, output='real', sort=lambda a: np.abs(a) < 1.0)
nu = n - sdim
iR, iT = list(range(n - nu)), list(range(n - nu, n))
F2 = U[:, iT].T                     # 封闭（A_TR=0），原文 §V 的 ν=2 曲线
Z2 = U[:, iR]                       # ker F2
F3 = U[:, [1, 2, 3]].T              # 不封闭（把 2×2 复块劈开），Fig.1 的 (ii) 曲线
Z3 = np.linalg.svd(F3)[2][3:].T

EPS = (0.0, 1e-4, 1e-2, 0.1, 0.5, 1.0)
DS = (50.0, 65.0, 80.0, 110.0, 200.0)


def row_ortho(F):
    return np.max(np.abs(F @ F.T - np.eye(F.shape[0])))


def align_family(F, Z, eps):
    return sym(F.T @ F + eps * (Z @ Z.T))


def report(tag, F, Z, eps, closed):
    Q = align_family(F, Z, eps)
    m = B.shape[1]
    Pc, K, Th = ctrl_at(A, B, Q, np.eye(m))
    floor = float(np.trace(W @ Pc))
    P, res, it = noiseless(F, A, W)
    phi = float(np.trace(Th @ P)) if np.isfinite(P).all() else np.inf
    tz = float(np.max(np.abs(Th @ Z)))
    line = '%-22s ε=%-7g %s  Φ=%+9.4f  D_min(K_F)=%9.4f  max|ΘZ|=%.2e' % (
        tag, eps, '封闭' if closed else '不封闭', phi, floor + phi, tz)
    print(line)
    print('     %8s %10s %10s %10s %12s' % ('D', 'I_unc', 'I_TRV', 'Δ', 'ρ=‖ZᵀS*Z‖/‖S*‖'))
    for D in DS:
        if D <= floor + phi + 1e-6:
            print('     %8.2f  低于地板 %.4f，跳过' % (D, floor + phi))
            continue
        u = unconstrained_sdp(D, Q, Th, Pc)
        S = np.linalg.inv(u['P']) - np.linalg.inv(u['Pt'])
        r = float(np.linalg.norm(Z.T @ S @ Z) / np.linalg.norm(S))
        f = face_direct(D, F, Th, Pc, gam0=init_gamma(F, S), n_starts=2)
        if not f:
            print('     %8.2f %10.4f %10s  锥上求解失败' % (D, u['I_true'], '—'))
            continue
        print('     %8.2f %10.4f %10.4f %+10.4f %12.3e' % (D, u['I_true'], f[0], f[0] - u['I_true'], r))


if __name__ == '__main__':
    print('F2 行正交性 max|FFᵀ−I|=%.2e ; F3=%.2e' % (row_ortho(F2), row_ortho(F3)))
    print('A_TR: F2=%.3e  F3=%.3e   A_RT: F2=%.3e' % (
        np.max(np.abs(T[np.ix_(iT, iR)])), np.max(np.abs(T[np.ix_([1, 2, 3], [0])])),
        np.max(np.abs(T[np.ix_(iR, iT)]))))
    print('\n======== 主实验：F 固定为封闭的 ν=2，只拧 ε ========')
    for eps in EPS:
        report('F2(ν=2)', F2, Z2, eps, True)
    print('\n======== 对照：把 ε=0 装到不封闭的 F3 上（只有封闭性破了）========')
    for eps in (0.0, 1.0):
        report('F3(ν=3)', F3, Z3, eps, False)
