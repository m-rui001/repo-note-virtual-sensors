r"""E33：$\kappa$ 的 5% 残差到底是谁的锅——重测 $b=q(B)$，并在深水区直接测率侧斜率。

推导先自查一遍（二阶 KKT 我重新算过）：沿 $V=tB,\ B=C^{-1}$，
  $x(t)=rt+q(B)t^2+O(t^3)$，
  $I_{dir}(t)=\tfrac12\log_2(\det C\det M_0)-\tfrac r2\log_2 t+\tfrac{t}{2\ln2}\operatorname{tr}(M_0^{-1}(D_{\mathbb M}+C^{-1}))+O(t^2)$。
把约束反解 $t=(x/r)(1-b_{eff}x/r^2)$ 代回时，"非对数项梯度 $G$"与"代价二次型梯度 $\mathcal Q$"产生的两处
$s^2$ 修正**正好互相抵消**，留下的确实是 $b=q(B)$。所以公式没有漏项，漏的是**测量**：
E31b 用定点迭代取 $x(t)$，而迭代残差在 $10^{-12}$ 量级——$q̂(t)=(x-rt)/t^2$ 的噪声地板是
$10^{-12}/t^2$，在 $t=10^{-4}$ 已达 $10^{-4}$、$t=10^{-6}$ 达 $10^{0}$，把整个拟合窗（$t$ 跨到 $0.032$）
往偏的方向拉。E31b 的 $b_{宽窗}=-0.8326$ 与 Richardson 的 $-0.8407$ 之差正是这个量级。

这轮全部换成精确路径：$\tilde P(V)$ 用 `solve_discrete_are`$(A^\top,F^\top,W,V)$ 直接解
（$10^{-15}$ 相对残差），$P=\tilde P-\tilde PF^\top(F\tilde PF^\top+V)^{-1}F\tilde P$，
$x=\operatorname{tr}(\Theta(P-P_0))$。于是
 (1) $q̂(t)$ 在 $t\in[3\times10^{-4},10^{-2}]$ 应当干净地趋于常数；用 $q̂=b+ct+dt^2$ 外推取 $b$；
 (2) 率侧斜率 $\kappa'$ 从 $I_{dir}$ 对 $x$ 的拟合取，窗选在 $x\ge5\times10^{-4}$：
     深水点虽然 $I$ 准，但 $x$ 的绝对误差 $\epsilon$ 会被对数项放大成 $(r/2)\epsilon/(x\ln2)$，
     $x<10^{-4}$ 处那一项比 $\kappa'x$ 还大——这是 E31b 口径的第二个隐患；
 (3) 预测 $\kappa'_{pred}=[b+\operatorname{tr}(M_0^{-1}(D_{\mathbb M}+C^{-1}))]/(2r\ln2)$，
     $\kappa=\kappa'+|I'_{unc}|$。$D_{\mathbb M}$ 也用 DARE 重新外推。
"""
import sys
import numpy as np
sys.stdout.reconfigure(encoding='utf-8')
from scipy.linalg import solve_discrete_are
from p0.p0_replicate_letter import A, W, n, sym, ln2, ctrl, unconstrained_sdp
from p0.exp_authoritative import noiseless
from p0.exp_nearfloor_law import Fs, floor_of, Th, Pc

UNC = float(np.trace(W @ Pc))


def Pt_of(F, V):
    """精确预测协方差：DARE 解 $\tilde P=A\tilde PA^\top-(A\tilde PF^\top)(V+F\tilde PF^\top)^{-1}(F\tilde PA^\top)+W$。"""
    return sym(solve_discrete_are(A.T, F.T, W, V))


def xI_of(F, V, Pn, Pt0):
    Pt = Pt_of(F, V)
    P = sym(Pt - Pt @ F.T @ np.linalg.solve(F @ Pt @ F.T + V, F @ Pt))
    x = float(np.trace(Th @ (P - Pn)))
    s1, l1 = np.linalg.slogdet(F @ Pt @ F.T + V)
    s2, l2 = np.linalg.slogdet(V)
    return x, .5 * (l1 - l2) / ln2, Pt


def coord(i, j):
    M = np.zeros((2, n)); M[0, i] = 1.; M[1, j] = 1.
    return M


if __name__ == '__main__':
    rng = np.random.default_rng(11)
    Qr = np.linalg.qr(rng.standard_normal((n, 2)))[0].T
    e1 = np.zeros((1, n)); e1[0, 0] = 1.
    cases = [(t, F) for t, F in Fs.items()] + [
        ('F6 坐标{x1,x2}', coord(0, 1)), ('F7 坐标{x2,x4}', coord(1, 3)),
        ('F8 随机正交行对', Qr), ('F9 单通道 y=x1', e1)]
    if __import__('os').environ.get('E33_ONE'):
        cases = cases[:1]
    TS = np.logspace(-5, -2, 7)
    print('自检：DARE 的 $V\to$ 定点解 vs `noiseless` 的迭代解（$t=10^{-3}$）')
    for tag, F in cases:
        r = F.shape[0]
        fl, phi, _ = floor_of(F)
        Pn, res0, it0 = noiseless(F, A, W, iters=80000, tol=1e-15)
        Pt0 = sym(A @ Pn @ A.T + W)
        M0 = sym(F @ Pt0 @ F.T)
        C = np.load('p0/fig/Cmat_%s.npy' % tag.split()[0])
        Ci = np.linalg.inv(C)
        trm = float(np.trace(np.linalg.solve(M0, Ci)))
        Iu_fl = unconstrained_sdp(fl, np.eye(n), Th, Pc)['I_true']
        b_pred = (r / 2) * np.log2(r) + .5 * np.log2(np.linalg.det(C) * np.linalg.det(M0)) - Iu_fl
        Vd = 1e-3 * Ci
        x1, _, Pt1 = xI_of(F, Vd, Pn, Pt0)
        print('\n===== %s  r=%d  地板=%.6f  $\\det C$=%.4f  $b_{pred}$=%.6f  定点残差=%.1e ====='
              % (tag, r, fl, np.linalg.det(C), b_pred, res0))

        # ---------- (1) $q̂(t)$：DARE 精确路径 ----------
        xs, qs, Is = [], [], []
        for t in TS:
            x, I, _ = xI_of(F, t * Ci, Pn, Pt0)
            xs.append(x); Is.append(I); qs.append((x - r * t) / t ** 2)
        print('  %10s %13s %14s %14s' % ('t', 'x(t)', 'q̂(t)', 'log残差 $I+\\frac r2\log_2x$'))
        for t, x, q, I in zip(TS, xs, qs, Is):
            print('  %10.2e %13.6e %14.6f %14.7f' % (t, x, q, I + (r / 2) * np.log2(x)))
        qq = np.array(qs)
        # $q̂=b+ct+dt^2$ 外推（用整窗，但窗已避开噪声地板）
        cb = np.polyfit(TS, qq, 2)
        b_deep = cb[-1]
        b_wide = float(np.dot(TS ** 2, np.array(xs) - r * TS) / np.dot(TS ** 2, TS ** 2))
        print('  $b=q(B)$ 外推=%.5f   E31b 式宽窗投影=%.5f   差=%+.5f   ($ct$ 系数=%.3f)'
              % (b_deep, b_wide, b_deep - b_wide, cb[1]))

        # ---------- (2) $D_{\mathbb M}$ 外推 ----------
        dm = []
        for t in (1e-3, 3e-4, 1e-4, 3e-5, 1e-5):
            _, _, Pt = xI_of(F, t * Ci, Pn, Pt0)
            dm.append(float(np.trace(np.linalg.solve(M0, sym(F @ Pt @ F.T - M0) / t))))
        print('  $\\operatorname{tr}(M_0^{-1}D_{\\mathbb M})$ @ $t=10^{-3}\dots10^{-5}$：%s  取末值=%.5f'
              % (np.array2string(np.array(dm), precision=5), dm[-1]))
        trD = dm[-1]

        # ---------- (3) 率侧斜率 ----------
        X = np.array(xs); IV = np.array(Is)
        asym = r / 2.
        L = IV + asym * np.log2(X)          # $=B+\kappa'x+\zeta x^2+\dots$
        for lo, hi, deg in ((5e-4, 3e-2, 1), (5e-4, 3e-2, 2), (1e-3, 2e-2, 1)):
            sel = (X >= lo) & (X <= hi)
            if sel.sum() < deg + 2:
                continue
            cf = np.polyfit(X[sel], L[sel], deg)
            k_emp = cf[-2]
            pred = (b_deep + trm + trD) / (2 * r * ln2)
            res = L[sel] - np.polyval(cf, X[sel])
            print('  窗 $x\in[%.0e,%.0e]$ 多项式阶 %d：$\\kappa\'_{实测}$=%.5f  $\\kappa\'_{预言}(b_{深})$=%.5f  比值=%.4f'
                  '  $(b_{宽}+\dots)$=%.5f 比值=%.4f  残差max=%.1e  $B$=%.5f'
                  % (lo, hi, deg, k_emp, pred, pred / k_emp,
                     (b_wide + trm + trD) / (2 * r * ln2),
                     (b_wide + trm + trD) / (2 * r * ln2) / k_emp, np.abs(res).max(), cf[-1]))
        # 对数项系数本身也放开，检验 $r/2$ 是否被数据支持
        M = np.vstack([np.log2(1 / X), np.ones_like(X), X]).T
        cf3, rs3 = np.linalg.lstsq(M, IV, rcond=None)[:2]
        sel = (X >= 5e-4)
        cov = np.linalg.inv(M[sel].T @ M[sel]) * (rs3[0] / (sel.sum() - 3))
        print('  自由三参拟合：对数系数=%.6f$\\pm$%.6f（理论 %.6f）  $\\kappa\'$=%.5f$\\pm$%.5f  $B$=%.5f'
              % (cf3[0], np.sqrt(cov[0, 0]), asym, cf3[2], np.sqrt(cov[2, 2]), cf3[1]))
