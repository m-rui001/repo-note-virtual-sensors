r"""E34b：$\kappa$ 的最终表——把 $b=q(B)$ 换成"$\delta t+bt^2$ 双参拟合"，并全程避开 SDP。

E34 的判定：沿 $V=tC^{-1}$，$g(t):=x(t)-rt$ 的 $q̂(t)=g/t^2$ 在深端发散，但**不是**亚二次物理项。
拟合 $q̂(t)=b+\delta/t$（$\delta=-8.05\times10^{-6}$，$b=-0.840$）把 $t=10^{-7}\dots5\times10^{-4}$ 的 16 个点
全部打到 3 位有效数字以内。$\delta$ 的来源就是我测的 $C$ 有 $4\times10^{-6}$ 的相对误差
（E28 的可加性核对给的就是这个量级）：$x=t\,\mathrm{tr}(C_{真}C^{-1})=t(r+\delta)$。
DARE 残差 $\sim2\times10^{-14}$、"从 DARE 解再迭代 20 万步"的位移 $\sim3\times10^{-14}$，
而深端 $|g|$ 在 $t=10^{-7}$ 仍有 $8\times10^{-13}$ —— 信噪比 30 倍，不是求根噪声。

于是 E31b 那两个错都清楚了：
 (i) $b$ 的拟合窗跨到 $t=0.032$，三阶项被吸进二阶系数（F2 给 $-0.8145$，真值 $-0.840$）；
 (ii) 就算窗对了，用 $b=\lim g/t^2$ 的朴素投影仍会被 $\delta/t$ 拖走。
现在 $b$ 用 $(\delta,b)$ 双参在深窗拟合，$\kappa$ 的**率侧部分** $\kappa'=\kappa-|I'_{unc}|$
直接从 $I_{dir}$ 对 $x$ 拟合（完全不含 SDP）：
$$\kappa'_{pred}=\frac{b+\operatorname{tr}(M_0^{-1}(D_{\mathbb M}+C^{-1}))}{2r\ln2}$$
对照实测窗 $x\in[5\times10^{-4},3\times10^{-2}]$、模型 $I_{dir}=-\tfrac r2\log_2x+B+\kappa'x+\zeta x^2$。
"""
import sys
import numpy as np
sys.stdout.reconfigure(encoding='utf-8')
from scipy.linalg import solve_discrete_are
from p0.p0_replicate_letter import A, W, n, sym, ln2, ctrl, unconstrained_sdp
from p0.exp_authoritative import noiseless
from p0.exp_nearfloor_law import Fs, floor_of, Th, Pc

UNC = float(np.trace(W @ Pc))


def xI(F, V, Pn):
    Pt = sym(solve_discrete_are(A.T, F.T, W, V))
    S = F @ Pt @ F.T + V
    P = sym(Pt - Pt @ F.T @ np.linalg.solve(S, F @ Pt))
    x = float(np.trace(Th @ (P - Pn)))
    _, l1 = np.linalg.slogdet(S)
    _, l2 = np.linalg.slogdet(V)
    return x, .5 * (l1 - l2) / ln2, Pt


def coord(i, j):
    M = np.zeros((2, n)); M[0, i] = 1.; M[1, j] = 1.
    return M


def fit_db(T, G):
    M = np.vstack([T, T ** 2]).T
    cf = np.linalg.lstsq(M, G, rcond=None)[0]
    return cf[0], cf[1]


if __name__ == '__main__':
    rng = np.random.default_rng(11)
    Qr = np.linalg.qr(rng.standard_normal((n, 2)))[0].T
    e1 = np.zeros((1, n)); e1[0, 0] = 1.
    cases = [(t, F) for t, F in Fs.items()] + [
        ('F6 坐标{x1,x2}', coord(0, 1)), ('F7 坐标{x2,x4}', coord(1, 3)),
        ('F8 随机正交行对', Qr), ('F9 单通道 y=x1', e1)]
    if __import__('os').environ.get('E34_ONE'):
        cases = cases[:1]
    print('%-22s %1s %10s %10s %9s %9s %9s %9s %7s'
          % ('F', 'r', '$\\delta$', '$b$', 'trm', 'trD', "$\\kappa'_{测}$", "$\\kappa'_{预}$", '比值'))
    rows = []
    for tag, F in cases:
        r = F.shape[0]
        fl, phi, _ = floor_of(F)
        Pn, res0, it0 = noiseless(F, A, W, iters=200000, tol=1e-17)
        Pt0 = sym(A @ Pn @ A.T + W)
        M0 = sym(F @ Pt0 @ F.T)
        C = np.load('p0/fig/Cmat_%s.npy' % tag.split()[0])
        Ci = np.linalg.inv(C)
        trm = float(np.trace(np.linalg.solve(M0, Ci)))
        # 深窗 $(\delta,b)$
        T1 = np.logspace(-6.3, -3.4, 9)
        G1 = np.array([xI(F, t * Ci, Pn)[0] - r * t for t in T1])
        d1, b1 = fit_db(T1, G1)
        T2 = np.logspace(-4.2, -2.6, 9)
        G2 = np.array([xI(F, t * Ci, Pn)[0] - r * t for t in T2])
        d2, b2 = fit_db(T2, G2)
        # 朴素单参投影（E31b 口径）对照
        b_naive = float(np.dot(T2 ** 2, G2) / np.dot(T2 ** 2, T2 ** 2))
        # $D_{\mathbb M}$
        dm = [float(np.trace(np.linalg.solve(M0, sym(F @ xI(F, t * Ci, Pn)[2] @ F.T - M0) / t)))
              for t in (1e-4, 1e-5, 1e-6)]
        trD = dm[-1]
        # 率侧斜率：$I_{dir}$ 对 $x$
        T3 = np.logspace(-4.3, -1.6, 9)
        X, IV = [], []
        for t in T3:
            x, I, _ = xI(F, t * Ci, Pn)
            X.append(x); IV.append(I)
        X, IV = np.array(X), np.array(IV)
        L = IV + (r / 2) * np.log2(X)
        cf = np.polyfit(X, L, 2)
        k_meas = cf[-2]
        k_pred = (b1 + trm + trD) / (2 * r * ln2)
        # 对数系数放开，独立验证 $r/2$
        Mfree = np.vstack([np.log2(1 / X), np.ones_like(X), X]).T
        cf3 = np.linalg.lstsq(Mfree, IV, rcond=None)[0]
        print('%-22s %1d %10.2e %10.5f %9.5f %9.5f %9.5f %9.5f %7.4f'
              % (tag, r, d1, b1, trm, trD, k_meas, k_pred, k_pred / k_meas))
        print('    浅窗 $\\delta$=%.1e $b$=%.5f（朴素投影 $b$=%.5f，E31b 用它）  $D_M$ 三点 %s  '
              '$\\zeta$=%.2e  拟合残差max=%.1e  自由对数系数=%.6f（理论 %.6f）'
              % (d2, b2, b_naive, np.array2string(np.array(dm), precision=5),
                 cf[0], np.abs(L - np.polyval(cf, X)).max(), cf3[0], r / 2.))
        rows.append((tag, r, b1, k_meas, k_pred, cf3[0]))
    print('\n汇总：率侧 $\kappa\'$ 比值 %s' % '  '.join('%s=%.4f' % (t.split()[0], p / m)
                                                       for t, r, b, m, p, lf in rows))
    print('对数系数偏差 %s' % '  '.join('%s=%+.1e' % (t.split()[0], lf - rr / 2.)
                                        for t, rr, b, m, p, lf in rows))
