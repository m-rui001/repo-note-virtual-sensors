"""判伪实验：L-CSS 2026 的表示惩罚 Δ(D) 到底在度量什么？

E1 多对一：同一 range(S) 的不同表示（F 与 c·F、F 与任意可逆行变换）给出完全相同的 (I,J)
    → 在线性高斯 LQG 内，"表示结构"被压缩成一个矩阵 S=F'ΓF，下游 (I,J) 看不见其他信息。
E2 复现 Fig.1：Q=I（全状态代价）下三种 TRV 架构（ν=2 不稳定 Schur 坐标 / +1 稳定坐标 / 全状态）
E3 关键检验：把代价换成任务对齐的 Q=F'F（原文 Eq.(4) 自称的 natural choice），
    任务子系统若动态封闭（Thm 3 前提），Δ(D) 应恒为 0。数值验证之。
"""
import numpy as np
from scipy.linalg import schur
from p0.p0_replicate_letter import (A, B, W, n, sym, ln2, ctrl, rate_cost, face_direct,
                                    unconstrained_sdp, plant_facts, init_gamma)

np.set_printoptions(precision=4, suppress=True)


CLOSURE = {}          # schur_task(ν) 的 (max|A_TR|, max|A_RT|) 诊断表


def schur_task(nu):
    """原文 §V 的 F=[0 I_ν]Ūᵀ：实 Schur，稳定块在前、不稳定块在尾，取尾部 ν 行。

    2026-09-28 修正（见 community.md 第 4 条第 0 节）：旧实现先取未排序的 schur(A,'real')，
    再用朴素下标置换把不稳定块搬到末尾。置换本身仍是正交相似，但会破坏拟三角性，
    于是 (a) max|A_TR| 出现 1.42 这类伪影，(b) 更致命的是当稳定模态是一对复共轭（2×2 块）时，
    按列搬运后尾部 ν 列张成的子空间并不对应"任务坐标 dynamics 封闭"那一块，
    Φ 与 Δ 全被带偏（同一个 F 的主角度与真不变子空间差 35.39°）。
    现在直接用 sort= 让 LAPACK 做正交块重排，并显式检查封闭性。
    """
    T, U, sdim = schur(A, output='real', sort=lambda a: abs(a) < 1.0)
    idx = list(range(n - nu, n))
    comp = [k for k in range(n) if k not in idx]
    F = U[:, idx].T
    A_TR = float(np.max(np.abs(T[np.ix_(idx, comp)])))     # 封闭性：ν=2 应为 0；ν=3 劈开复块 ⟹ 非零
    A_RT = float(np.max(np.abs(T[np.ix_(comp, idx)])))
    CLOSURE[(nu)] = (A_TR, A_RT)                            # 诊断用，不参与计算
    return F, A_TR


def full_state_F():
    return np.eye(n)


def cone_resid(S, F):
    """S 到锥 K_F={S⪰0: range(S)⊆range(Fᵀ)} 的违约量：正交补块 + 交叉块。
    对应原文 Corollary 2 的 range(S*)⊆range(Fᵀ) 条件。"""
    Q, _ = np.linalg.qr(F.T)
    M = np.eye(n) - Q @ Q.T
    return (float(np.linalg.norm(M @ S @ M) / np.linalg.norm(S)),
            float(np.linalg.norm(M @ S @ Q) / np.linalg.norm(S)))


def run(Ds, F, Qmat):
    """沿 D 从宽到窄 warm-start；同时记录无约束最优 S* 的锥违约量"""
    Pc, K, Theta = ctrl(Qmat)
    warm, rows = None, []
    for D in sorted(Ds, reverse=True):
        unc = unconstrained_sdp(D, Qmat, Theta, Pc)
        S_unc = np.linalg.inv(unc['P']) - np.linalg.inv(unc['Pt']) if unc else None
        g0 = init_gamma(F, S_unc) if S_unc is not None else None
        fac = face_direct(D, F, Theta, Pc, gam0=g0, warm=warm, n_starts=1)
        if unc and fac:
            warm = fac[2]
            I_trv, J_trv = rate_cost(F, fac[2], Theta, Pc)[:2]
            a, b = cone_resid(S_unc, F)
            rows.append((D, unc['I_true'], I_trv, I_trv - unc['I_true'], J_trv, a, b))
        else:
            rows.append((D, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan))
    return sorted(rows, key=lambda r: r[0])


if __name__ == '__main__':
    ev, nu, ib = plant_facts()
    print('ν=%d, 下界=%.4f bits/sample\n' % (nu, ib))

    F2, closed2 = schur_task(2)
    F3, closed3 = schur_task(3)
    print('F(ν=2) 行空间封闭残差 max|A_TR| = %.2e' % closed2)
    print('F(ν=3) 行空间封闭残差 max|A_TR| = %.2e\n' % closed3)

    # ---- E1 多对一：F 与 2F、F 与 QF（Q 可逆）应给出同一 (I,J)
    Pc, K, Theta = ctrl(np.eye(n))
    Gam = np.eye(2)
    i1, j1 = rate_cost(F2, Gam, Theta, Pc)[:2]
    i2, j2 = rate_cost(2 * F2, 0.25 * Gam, Theta, Pc)[:2]
    R = np.array([[0.9, 0.1], [0.2, 0.8]])
    i3, j3 = rate_cost(R @ F2, np.linalg.inv(R @ R.T), Theta, Pc)[:2]
    print('E1  F,Γ=(I):        I=%.6f J=%.6f' % (i1, j1))
    print('E1  2F,Γ/4:         I=%.6f J=%.6f' % (i2, j2))
    print('E1  RF,(RRᵀ)^-1:    I=%.6f J=%.6f' % (i3, j3))
    print('E1  结论：三者 S=FᵀΓF 相同 ⟹ 下游 (I,J) 相同：',
          np.allclose([i1, j1], [i2, j2], rtol=1e-6) and np.allclose([i1, j1], [i3, j3], rtol=1e-4), '\n')

    Ds = [34.0, 40.0, 50.0, 65.0, 80.0]
    hdr = '%6s %9s %9s %9s %9s %9s %9s' % ('D', 'I_unc', 'I_TRV', 'Δ', 'J_trv', '‖MS*M‖/‖S‖', '‖MS*Q‖/‖S‖')
    print('---- E2  Q = I（全状态代价，对应原文 Fig.1 的情形）')
    for tag, F in (('ν=2 不稳定', F2), ('ν=3 +1稳定', F3)):
        print(' [%s]' % tag)
        print(hdr)
        for D, iu, it, d, jt, a, b in run(Ds, F, np.eye(n)):
            print('  %6.1f %9.4f %9.4f %9.4f %9.3f %9.4f %9.4f' % (D, iu, it, d, jt, a, b))

    print('\n---- E3  Q = FᵀF（任务对齐代价，原文 Eq.(4) 的 natural choice）')
    for tag, F in (('ν=2 不稳定(封闭残差见上)', F2), ('ν=3 +1稳定', F3)):
        Qm = sym(F.T @ F)
        print(' [%s]  Q=FᵀF' % tag)
        print(hdr)
        for D, iu, it, d, jt, a, b in run(Ds, F, Qm):
            print('  %6.1f %9.4f %9.4f %9.4f %9.3f %9.4f %9.4f' % (D, iu, it, d, jt, a, b))
