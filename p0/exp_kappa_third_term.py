"""E31：把残差的 $O(x)$ 系数 $\\kappa$ 也推出来并验证，让渐近式成为**三项全显式**。

沿最优射线 $V=tC^{-1}$（$C\\succ0$ 保证这是条合法射线）一维展开：
  $x(t)=\\operatorname{tr}(\\Theta P(tC^{-1}))-\\operatorname{tr}(\\Theta P_0)=rt+b\\,t^2+O(t^3)$
    （$rt=\\operatorname{tr}(C\\cdot tC^{-1})$ 就是步(b)的线性项，$b$ 是二阶项，直接拟合可测）
  $I(t)=-\\tfrac r2\\log_2 t+\\tfrac12\\log_2\\bigl(\\det C\\cdot\\det M_0\\bigr)+\\tfrac{t}{2\\ln2}\\operatorname{tr}(M_0^{-1}C^{-1})+O(t^2)$
    （$M_0=F\\tilde P_0F^\\top$；把 $\\det(F\\tilde PF^\\top+V)$ 在 $V=tC^{-1}$ 上展开即得那个迹）
反解 $t=(x/r)(1-bx/r^2+\\dots)$ 代回，再扣掉 $I_{unc}(地板+x)=I_{unc}(地板)+xI'_{unc}$：
$$\\Delta(x)=\\tfrac r2\\log_2\\tfrac1x+b_{pred}+\\kappa x+O(x^2),\\qquad
\\boxed{\\ \\kappa=\\frac{b+\\operatorname{tr}(M_0^{-1}C^{-1})}{2r\\ln2}+\\bigl|I'_{unc}(D_{\\min}(\\mathcal K_F))\\bigr|\\ }$$
右边三个量全可测：$b$ 由 $x(t)$ 的二次拟合，迹由 $M_0,C$ 的矩阵运算，$I'_{unc}$ 由无约束 SDP 的中心差分。
于是这是一个"用低阶数据预言高阶项"的闭环检验：若四个 $F$ 的 $\\kappa$ 都命中，Prop E 就升级成三项展开。
"""
import sys
import numpy as np
sys.stdout.reconfigure(encoding='utf-8')
from scipy.linalg import schur
from p0.p0_replicate_letter import A, W, n, sym, ln2, ctrl, unconstrained_sdp
from p0.exp_authoritative import noiseless
from p0.exp_nearfloor_law import Fs, floor_of, Th, Pc

UNC = float(np.trace(W @ Pc))
P0 = [None]


def filter_V(F, V, iters=40000, tol=1e-16):
    P = P0[0][0].copy()
    for k in range(iters):
        Pt = sym(A @ P @ A.T + W)
        Pn = sym(Pt - Pt @ F.T @ np.linalg.solve(F @ Pt @ F.T + V, F @ Pt))
        if k > 300 and np.max(np.abs(Pn - P)) < tol:
            P = Pn
            break
        P = Pn
    return P, sym(A @ P @ A.T + W)


def I_of(F, V, Pt=None):
    if Pt is None:
        _, Pt = filter_V(F, V)
    return 0.5 * np.log(np.linalg.det(F @ Pt @ F.T + V) / np.linalg.det(V)) / ln2


if __name__ == '__main__':
    import p0.exp_oos_asymptote as e30
    from p0.p0_replicate_letter import n as _n
    rng = np.random.default_rng(11)
    Qr = np.linalg.qr(rng.standard_normal((_n, 2)))[0].T
    def coord(i, j):
        M = np.zeros((2, _n)); M[0, i] = 1.; M[1, j] = 1.
        return M
    e1 = np.zeros((1, _n)); e1[0, 0] = 1.
    cases = [(t, F) for t, F in Fs.items()] + [
        ('F6 坐标{x1,x2}', coord(0, 1)), ('F7 坐标{x2,x4}', coord(1, 3)),
        ('F8 随机正交行对', Qr), ('F9 单通道 y=x1', e1)]
    print('%-22s %3s %9s %10s %9s %10s %10s %10s'
          % ('F', 'r', 'b(二阶)', 'tr(M0⁻¹C⁻¹)', 'κ_预言', 'κ_实测', "I'_unc", 'b_锚定差'))
    for tag, F in cases:
        r = F.shape[0]
        fl, phi, _ = floor_of(F)
        Pn, res, it = noiseless(F, A, W, iters=80000, tol=1e-15)
        P0[0] = (Pn, sym(A @ Pn @ A.T + W))
        M0 = sym(F @ P0[0][1] @ F.T)
        fp = 'p0/fig/Cmat_%s.npy' % tag.split()[0]
        try:
            C = np.load(fp)
        except Exception:
            e30.e28.P0[0] = (Pn, sym(A @ Pn @ A.T + W))   # measure_C 用的是 e28 的全局
            C = e30.measure_C(F)
            np.save(fp, C)
        Ci = np.linalg.inv(C)

        ts = np.logspace(-5, -1.5, 14)
        xs = np.array([float(np.trace(Th @ (filter_V(F, t * Ci)[0] - Pn))) for t in ts])
        # x = r t + b t²：固定线性项为 r，只回归 b
        b = float(np.dot(ts ** 2, xs - r * ts) / np.dot(ts ** 2, ts ** 2))
        trm = float(np.trace(np.linalg.solve(M0, Ci)))
        dIu = 0.5e-3
        sl = (unconstrained_sdp(fl + dIu, np.eye(n), Th, Pc)['I_true']
              - unconstrained_sdp(fl - dIu, np.eye(n), Th, Pc)['I_true']) / (2 * dIu)
        k_pred = (b + trm) / (2 * r * ln2) + abs(sl)

        def Delt(t):
            P, Pt = filter_V(F, t * Ci)
            xr = float(np.trace(Th @ (P - Pn)))
            Iu = unconstrained_sdp(fl + xr, np.eye(n), Th, Pc)['I_true']
            return xr, I_of(F, t * Ci, Pt) - Iu
        pts = np.array([Delt(t) for t in (1e-4, 3e-4, 1e-3, 3e-3, 1e-2)])
        asym = (r / 2.) * np.log2(10.)
        b_theory = (r / 2) * np.log2(r) + .5 * np.log2(np.linalg.det(C) * np.linalg.det(M0)) - \
            unconstrained_sdp(fl, np.eye(n), Th, Pc)['I_true']
        b_anchor = pts[0, 1] - asym * np.log10(1 / pts[0, 0])
        k_emp = np.polyfit(pts[:, 0], pts[:, 1] - (asym * np.log10(1 / pts[:, 0]) + b_anchor), 1)[0]
        print('%-22s %3d %9.4f %10.4f %9.5f %10.5f %10.5f %10.5f'
              % (tag, r, b, trm, k_pred, k_emp, sl, b_theory - b_anchor))
