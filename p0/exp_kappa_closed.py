"""E31b：补上 E31 漏掉的那一块，把 $\\kappa$ 的预言闭合。

E31 的失败项很清楚：我写了 $\\tfrac{d}{dt}\\big|_0\\det(F\\tilde P(tC^{-1})F^\\top+tC^{-1})$，
然后把 $\\tilde P$ 对 $V$ 的方向导数**错换成了 $C^{-1}$ 本身**。正确的对象是
$$D_{\\mathbb M}:=\\frac{d}{dt}\\Big|_{0}F\\big(\\tilde P(tC^{-1})-\\tilde P_0\\big)F^\\top\\quad(\\text{秩 }r\\times r)$$
可测：$[F\\tilde P(tC^{-1})F^\\top-F\\tilde P_0F^\\top]/t$ 在三个 $t$ 上做 Richardson。于是
$$\\kappa=\\frac{b+\\operatorname{tr}\\!\\big(M_0^{-1}(D_{\\mathbb M}+C^{-1})\\big)}{2r\\ln2}+\\bigl|I'_{unc}(D_{\\min}(\\mathcal K_F))\\bigr|$$
（$b$ 是 $x(t)=rt+bt^2$ 的二阶系数，来自"标定 $t$ 让真实 $x$ 命中目标"这一步的偏移。）
对照口径也改深：拟合 $\\kappa_{\\text{实测}}$ 只用 $x\\le5\\times10^{-3}$ 的点，避免把 $O(x^2)$ 曲率算进斜率。
"""
import sys
import numpy as np
sys.stdout.reconfigure(encoding='utf-8')
from p0.p0_replicate_letter import A, W, n, sym, ln2, ctrl, unconstrained_sdp
from p0.exp_authoritative import noiseless
from p0.exp_nearfloor_law import Fs, floor_of, Th, Pc
from p0.exp_kappa_third_term import filter_V, I_of
import p0.exp_oos_asymptote as e30
import p0.exp_kappa_third_term as e31

UNC = float(np.trace(W @ Pc))
rng = np.random.default_rng(11)
Qr = np.linalg.qr(rng.standard_normal((n, 2)))[0].T


def coord(i, j):
    M = np.zeros((2, n)); M[0, i] = 1.; M[1, j] = 1.
    return M


if __name__ == '__main__':
    e1 = np.zeros((1, n)); e1[0, 0] = 1.
    cases = [(t, F) for t, F in Fs.items()] + [
        ('F6 坐标{x1,x2}', coord(0, 1)), ('F7 坐标{x2,x4}', coord(1, 3)),
        ('F8 随机正交行对', Qr), ('F9 单通道 y=x1', e1)]
    print('%-22s %2s %9s %11s %11s %9s %9s %9s'
          % ('F', 'r', 'b', 'tr(M0⁻¹C⁻¹)', 'tr(M0⁻¹D_M)', 'κ_预言', 'κ_实测', '比值'))
    for tag, F in cases:
        r = F.shape[0]
        fl, phi, _ = floor_of(F)
        Pn, res, it = noiseless(F, A, W, iters=80000, tol=1e-15)
        e31.P0[0] = (Pn, sym(A @ Pn @ A.T + W))
        M0 = sym(F @ e31.P0[0][1] @ F.T)
        fp = 'p0/fig/Cmat_%s.npy' % tag.split()[0]
        try:
            C = np.load(fp)
        except Exception:
            e30.e28.P0[0] = (Pn, sym(A @ Pn @ A.T + W))
            C = e30.measure_C(F)
            np.save(fp, C)
        Ci = np.linalg.inv(C)

        ts = np.logspace(-5, -1.5, 14)
        xs = np.array([float(np.trace(Th @ (filter_V(F, t * Ci)[0] - Pn))) for t in ts])
        b = float(np.dot(ts ** 2, xs - r * ts) / np.dot(ts ** 2, ts ** 2))
        trm = float(np.trace(np.linalg.solve(M0, Ci)))
        # D_M：Richardson 外推到 t=0
        dm_est = []
        for t in (3e-5, 1e-5, 3e-6):
            _, Pt = filter_V(F, t * Ci)
            dm_est.append((sym(F @ Pt @ F.T - M0) / t, t))
        D_M = dm_est[-1][0]
        trD = float(np.trace(np.linalg.solve(M0, D_M)))
        dIu = 0.5e-3
        sl = (unconstrained_sdp(fl + dIu, np.eye(n), Th, Pc)['I_true']
              - unconstrained_sdp(fl - dIu, np.eye(n), Th, Pc)['I_true']) / (2 * dIu)
        k_pred = (b + trm + trD) / (2 * r * ln2) + abs(sl)

        pts = []
        for t in (1e-5, 3e-5, 1e-4, 3e-4, 1e-3):
            P, Pt = filter_V(F, t * Ci)
            xr = float(np.trace(Th @ (P - Pn)))
            if xr > 5e-3:
                continue
            Iu = unconstrained_sdp(fl + xr, np.eye(n), Th, Pc)['I_true']
            pts.append((xr, I_of(F, t * Ci, Pt) - Iu))
        pts = np.array(pts)
        asym = (r / 2.) * np.log2(10.)
        b_anchor = pts[0, 1] - asym * np.log10(1 / pts[0, 0])
        k_emp = np.polyfit(pts[:, 0], pts[:, 1] - (asym * np.log10(1 / pts[:, 0]) + b_anchor), 1)[0]
        print('%-22s %2d %9.4f %11.4f %11.4f %9.5f %9.5f %9.3f'
              % (tag, r, b, trm, trD, k_pred, k_emp, k_pred / k_emp if k_emp else np.nan))
        print('    拟合窗 $x\\in[%.1e,%.1e]$（%d 点），$b_{pred}-b_{锚定}=%+.5f$   $t$=$%s$ 处 $D_{\\mathbb M}$ 三分量 %s'
              % (pts[0, 0], pts[-1, 0], len(pts),
                 (r / 2) * np.log2(r) + .5 * np.log2(np.linalg.det(C) * np.linalg.det(M0))
                 - unconstrained_sdp(fl, np.eye(n), Th, Pc)['I_true'] - b_anchor,
                 ' '.join('%.0e' % d[1] for d in dm_est),
                 np.array2string(np.array([np.trace(np.linalg.solve(M0, d[0])) for d in dm_est]),
                                 precision=4)))
