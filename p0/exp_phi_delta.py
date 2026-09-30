"""E6：Δ(D)（信息代价）与 Φ(F)（可达代价下界间距）是否互相独立？
若能把 (Δ,Φ) 的四个角都填上，就证明“只看 Δ(D) 的无损性判据不完备”。
Φ(F) := D_min(K_F) − D_min^unc = tr(Θ P_∞(F))，一次 DARE 即可算，无需 SCA。
"""
import numpy as np
from p0.p0_replicate_letter import (W, n, sym, ln2, ctrl, rate_cost, face_direct,
                                    unconstrained_sdp, init_gamma)
from p0.exp_trv_penalty import schur_task, cone_resid

Pc, K, Theta = ctrl(np.eye(n))
unc_floor = float(np.trace(W @ Pc))


def Phi(F, prec=1e10):
    _, J, _, _ = rate_cost(F, prec * np.eye(F.shape[0]), Theta, Pc)
    return J - unc_floor


def delta(F, Ds=(40.0, 80.0)):
    out = []
    for D in Ds:
        u = unconstrained_sdp(D, np.eye(n), Theta, Pc)
        S = np.linalg.inv(u['P']) - np.linalg.inv(u['Pt'])
        f = face_direct(D, F, Theta, Pc, gam0=init_gamma(F, S), n_starts=1)
        out.append((D, u['I_true'], f[0] if f else np.nan,
                    (f[0] - u['I_true']) if f else np.inf, np.linalg.matrix_rank(S)))
    return out


def show(tag, F):
    ph = Phi(F)
    ds = delta(F)
    print('%-22s Φ=%+9.3f (%+6.1f%%) |' % (tag, ph, 100 * ph / unc_floor), end=' ')
    for D, iu, it, d, r in ds:
        print(' D=%.0f:Δ=%s' % (D, ('%+.3f' % d) if np.isfinite(d) else '不可行'), end=' ')
    print()


if __name__ == '__main__':
    ev, u = np.linalg.eigh(sym(Theta))
    Vth = u[:, np.argsort(ev)[::-1][:2]]           # range(Θ)（2 维）
    print('rank(Θ)=%d, D_min^unc=%.4f' % (int(np.sum(ev > 1e-9)), unc_floor))
    print('%-22s %10s %s' % ('表示 F', 'Φ(代价下界差)', 'Δ(D) 率惩罚'))
    show('Schur 不稳定 ν=2', schur_task(2)[0])
    show('Schur ν=3', schur_task(3)[0])
    show('全状态 I', np.eye(n))
    show('span(range Θ)', Vth.T)
    show('坐标 e1,e2', np.array([[1., 0, 0, 0], [0, 1., 0, 0]]))
    show('坐标 e1,e3', np.array([[1., 0, 0, 0], [0, 0, 1., 0]]))
    show('坐标 e1,e4', np.array([[1., 0, 0, 0], [0, 0, 0, 1.]]))
    rng = np.random.default_rng(3)
    for k in range(4):
        M = rng.standard_normal((2, n))
        show('随机 rank-2 #%d' % k, M)
    # 用无约束最优方向造的“Cor.2 无损”F（E5 的角）
    u80 = unconstrained_sdp(80.0, np.eye(n), Theta, Pc)
    S80 = np.linalg.inv(u80['P']) - np.linalg.inv(u80['Pt'])
    e, Vv = np.linalg.eigh(sym(S80))
    show('span(S* 主方向) r=1', Vv[:, np.argsort(e)[::-1][:1]].T)
