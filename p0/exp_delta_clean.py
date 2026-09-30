"""E13：在排序 Schur 的正确 F 上重跑 Δ(D)（entry #3 与 E9/E10 的两个版本都不作数）。
   与 E12 的地板 D_min(K_F) 对照：D 必须 > 地板才可行；顺便检验 SLSQP 是否在地板附近伪造可行点。
"""
import numpy as np
from scipy.linalg import schur
from p0.p0_replicate_letter import (A, B, W, n, sym, ln2, ctrl, rate_cost, face_direct,
                                    unconstrained_sdp, init_gamma)
from p0.exp_authoritative import noiseless, ctrl_at

np.set_printoptions(precision=4, suppress=True)
T, U, sdim = schur(A, output='real', sort=lambda a: np.abs(a) < 1.0)
nu = n - sdim
iT = list(range(n - nu, n))
Fs = {'F2 尾部2(封闭)': U[:, iT].T,
      'F3 尾部3(不封闭)': U[:, [1, 2, 3]].T}
Pc, K, Theta = ctrl(np.eye(n))
unc_floor = float(np.trace(W @ Pc))

if __name__ == '__main__':
    print('D_min^unc=%.4f  (Q=I)' % unc_floor)
    for tag, F in Fs.items():
        P, res, it = noiseless(F, A, W)
        phi = float(np.trace(Theta @ P))
        dm = unc_floor + phi
        print('\n[%s]  rank F=%d  Φ=%+.4f  D_min(K_F)=%.4f' %
              (tag, np.linalg.matrix_rank(F), phi, dm))
        print('%8s %10s %10s %10s %10s %10s' % ('D', 'I_unc', 'I_TRV', 'Δ', 'J_trv', '锥违约'))
        for D in sorted({dm + 0.3, dm + 2.0, 50.0, 65.0, 80.0}):
            if D <= dm + 1e-6:
                print('%8.3f  低于地板，跳过' % D)
                continue
            u = unconstrained_sdp(D, np.eye(n), Theta, Pc)
            S = np.linalg.inv(u['P']) - np.linalg.inv(u['Pt'])
            f = face_direct(D, F, Theta, Pc, gam0=init_gamma(F, S), n_starts=2)
            if not f:
                print('%8.3f %10.4f %10s  锥上求解失败' % (D, u['I_true'], '—'))
                continue
            I_trv, J_trv = f[0], f[1]
            Gam = f[2]
            Sc = F.T @ Gam @ F
            viol = np.max(np.abs(Sc - Sc @ F.T @ np.linalg.solve(F @ F.T, F)) if F.shape[0] == n
                          else Sc - np.linalg.qr(F.T)[0] @ (np.linalg.qr(F.T)[0].T @ Sc))
            print('%8.3f %10.4f %10.4f %+10.4f %10.3f %10.2e' %
                  (D, u['I_true'], I_trv, I_trv - u['I_true'], J_trv, viol))
