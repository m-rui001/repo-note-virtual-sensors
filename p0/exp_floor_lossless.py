"""两个判决性实验（2026-09-28 第 2 轮，B 侧）：

E4 锥的"可达代价下界" D_min(K_F)：信息买断（Γ→∞）也过不去的那道墙。
   L-CSS 只给了"稳定化所需的最小信息"下界（其 Thm 2.2），没有给"固定表示下可达的最小代价"。
E5 "已被 Corollary 2 判为无损（Δ(D)=0）的表示，仍然可以在可达代价上贵一倍"——
   构造式反例：取 F = 无约束最优 S*(D0) 的主特征方向 ⟹ range(S*)⊆range(Fᵀ) ⟹ Δ(D0)=0，
   再算它的 D_min(K_F)。
"""
import numpy as np
from p0.p0_replicate_letter import (A, B, W, n, sym, ln2, ctrl, rate_cost, face_direct,
                                    unconstrained_sdp, init_gamma)
from p0.exp_trv_penalty import schur_task, cone_resid

np.set_printoptions(precision=4, suppress=True)


def floor_cost(F, Theta, Pc, prec=1e10):
    _, J, _, _ = rate_cost(F, prec * np.eye(F.shape[0]), Theta, Pc)
    return J


def top_dir_F(S, r):
    ev, V = np.linalg.eigh(sym(S))
    idx = np.argsort(ev)[::-1][:r]
    return V[:, idx].T


if __name__ == '__main__':
    Pc, K, Theta = ctrl(np.eye(n))
    unc_floor = float(np.trace(W @ Pc))
    print('D_min^unc = tr(W Pc) = %.4f' % unc_floor)

    F2, _ = schur_task(2)
    F3, _ = schur_task(3)
    print('\n---- E4 可达代价下界 D_min(K_F)（Γ→∞）')
    for tag, F in (('ν=2 不稳定 Schur', F2), ('ν=3 +1 稳定 Schur', F3),
                   ('全状态 I', np.eye(n)),
                   ('单坐标 e1', np.array([[1., 0, 0, 0]])),
                   ('单坐标 e3', np.array([[0., 0, 1., 0]]))):
        f = floor_cost(F, Theta, Pc)
        print('  %-18s D_min^cone=%9.3f   超出无约束下界 %+8.3f  (%.1f%%)'
              % (tag, f, f - unc_floor, 100 * (f - unc_floor) / unc_floor))

    print('\n---- E5 “无损但昂贵”的构造反例')
    for D0 in (40.0, 80.0):
        unc = unconstrained_sdp(D0, np.eye(n), Theta, Pc)
        S = np.linalg.inv(unc['P']) - np.linalg.inv(unc['Pt'])
        r = unc['rank_S']
        F = top_dir_F(S, max(r, 1))
        fac = face_direct(D0, F, Theta, Pc, gam0=init_gamma(F, S), n_starts=1)
        I_trv = fac[0] if fac else np.nan
        print('  取 D0=%.0f：无约束最优 I_unc=%.4f, rank S*=%d ⟹ 令 F=span(主特征方向)，r=%d' %
              (D0, unc['I_true'], r, F.shape[0]))
        print('     锥违约量 ‖MS*M‖/‖S‖=%.2e, ‖MS*Q‖/‖S‖=%.2e （Cor.2 预言 Δ=0）' % cone_resid(S, F))
        print('     实测 I_TRV(D0)=%.4f  ⟹ Δ(D0)=%+.4f bits' % (I_trv, I_trv - unc['I_true']))
        f = floor_cost(F, Theta, Pc)
        print('     但该表示的可达代价下界 D_min^cone=%.3f，比全状态贵 %+.3f（%.0f%%）' %
              (f, f - unc_floor, 100 * (f - unc_floor) / unc_floor))
        for D2 in (65.0, 80.0):
            if abs(D2 - D0) < 1e-9:
                continue
            u2 = unconstrained_sdp(D2, np.eye(n), Theta, Pc)
            S2 = np.linalg.inv(u2['P']) - np.linalg.inv(u2['Pt'])
            fc = face_direct(D2, F, Theta, Pc, gam0=init_gamma(F, S2), n_starts=1)
            print('     同一个 F 在 D=%.0f：I_unc=%.4f  I_TRV=%s  Δ=%s' %
                  (D2, u2['I_true'], '%.4f' % fc[0] if fc else '不可行',
                   '%+.4f' % (fc[0] - u2['I_true']) if fc else '—'))
