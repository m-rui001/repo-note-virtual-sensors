r"""E34：$x(t)$ 深端的亚二次项是真的，还是我的数值证书不够？

E33 的梯子（DARE 精确路径 + 独立定点迭代，两路在 $t=10^{-5}$ 处对 $x-rt$  agreeing 到 $2\times10^{-14}$）
显示 $q̂(t)=(x-rt)/t^2$ 在 $t\lesssim3\times10^{-4}$ 重新发散：
$g(t):=x-rt$ 的相邻十倍比值给 $\Delta D$ 出的局部指数是 $2.00,2.00,1.99,1.96,1.85,1.65$——
**越深越小**。若真实存在 $g=bt^2+ct^p\ (1<p<2)$，那么
 $\Delta(x)=\tfrac r2\log_2\tfrac1x+b_{pred}+\kappa x+O(x^p)$，三项形式仍然成立（$x^p=o(x)$ 不成立！
 $x^p$ 在 $p<2$ 时比 $x^2$ 大但比 $x$ 小，仍 $\to0$ 于 $x$ 之比），
 可是"$\kappa$ 的实测窗"里混的就是 $x^{p-1}$ 阶污染：$x^{p-2}$ 在 $x=10^{-3},p=1.65$ 已有 11 倍。
所以必须先把两件事分开：
 (1) **精度证书**：DARE 残差 $\max\bigl|\tilde P-(A\tilde PA^\top-A\tilde PF^\top(V+F\tilde PF^\top)^{-1}F\tilde PA^\top+W)\bigr|$、
     闭环本征值（判条件数）、以及"从 DARE 解再迭代 20 万步"的位移。若位移 $\ll |g|$ 的量级差，
     发散就不是求根误差。再用两个不同的 $\Phi_0$ 参照（$10^5$ 与 $10^6$ 步定点）看共模偏差。
 (2) **指数**：$t\in[10^{-7},3\times10^{-4}]$ 上拟合 $g=\beta t^{\rho}$（$\beta,\rho$ 自由），
     以及 $g=bt^2+ct^p$ 四参拟合；报 $\rho$ 的窗稳定性。$\rho\to2$ 就是我看错了，$\rho\approx1.5$ 就是
     退化 Riccati 的 Kato 半幂律——那将直接改写 Prop E 的第二项误差阶。
"""
import sys
import numpy as np
sys.stdout.reconfigure(encoding='utf-8')
from scipy.linalg import solve_discrete_are
from p0.p0_replicate_letter import A, W, n, sym, ln2, ctrl, unconstrained_sdp
from p0.exp_authoritative import noiseless
from p0.exp_nearfloor_law import Fs, floor_of, Th, Pc

UNC = float(np.trace(W @ Pc))


def Pt_dare(F, V):
    return sym(solve_discrete_are(A.T, F.T, W, V))


def dare_res(F, V, Pt):
    """Riccati 残差（绝对值，量纲同 $P$）。"""
    S = F @ Pt @ F.T + V
    rhs = A @ Pt @ A.T - A @ Pt @ F.T @ np.linalg.solve(S, F @ Pt @ A.T) + W
    return float(np.max(np.abs(sym(Pt - rhs))))


def P_from(F, V, Pt):
    return sym(Pt - Pt @ F.T @ np.linalg.solve(F @ Pt @ F.T + V, F @ Pt))


def refine(F, V, P, iters=200000):
    """从给定 $P$ 再迭代信息形态，看它移动多少（求根误差的独立证书）。"""
    for k in range(iters):
        Pt = sym(A @ P @ A.T + W)
        Pn = sym(Pt - Pt @ F.T @ np.linalg.solve(F @ Pt @ F.T + V, F @ Pt))
        if k > 2000 and np.max(np.abs(Pn - P)) < 1e-18:
            P = Pn
            break
        P = Pn
    return P


if __name__ == '__main__':
    tag = sys.argv[1] if len(sys.argv) > 1 else list(Fs)[0]
    F = Fs[tag]
    r = F.shape[0]
    fl, phi, _ = floor_of(F)
    Pn, res0, it0 = noiseless(F, A, W, iters=200000, tol=1e-17)
    Phi_ref = float(np.trace(Th @ Pn))
    Pn2, res1, it1 = noiseless(F, A, W, iters=2000000, tol=1e-18)
    print('===== %s  r=%d  地板=%.6f =====' % (tag, r, fl))
    print('  $\Phi_0$（20 万步）=%.14e   （200 万步）=%.14e   差=%.2e'
          % (Phi_ref, float(np.trace(Th @ Pn2)), float(np.trace(Th @ (Pn2 - Pn)))))
    Ci = np.linalg.inv(np.load('p0/fig/Cmat_%s.npy' % tag.split()[0]))
    Pt0 = sym(A @ Pn @ A.T + W)
    M0 = sym(F @ Pt0 @ F.T)
    print('  $\operatorname{cond}(M_0)$=%.3e   $\operatorname{cond}(C^{-1})$=%.3e'
          % (np.linalg.cond(M0), np.linalg.cond(Ci)))

    TS = np.logspace(-7, -3.3, 16)
    print('%11s %14s %12s %13s %12s %13s %11s'
          % ('t', 'g=x-rt', 'q̂=g/t²', 'DARE残差', '再迭代位移', '$|\\lambda_{cl}|_{max}$', '$\\rho_{loc}$'))
    gs, ts = [], []
    prev = None
    for t in TS:
        V = t * Ci
        Pt = Pt_dare(F, V)
        P = P_from(F, V, Pt)
        x = float(np.trace(Th @ (P - Pn)))
        g = x - r * t
        Pr = refine(F, V, P)
        shift = float(np.trace(Th @ (Pr - P)))
        Kl = np.linalg.eigvals(A - A @ Pt @ F.T @ np.linalg.solve(F @ Pt @ F.T + V, F))
        rho = np.nan if prev is None else np.log(gs[-1] / g) / np.log(TS[len(gs) - 1] / t) if g * gs[-1] > 0 else np.nan
        ts.append(t); gs.append(g)
        print('%11.3e %14.4e %12.3e %13.2e %12.2e %13.6f %11.3f'
              % (t, g, g / t ** 2, dare_res(F, V, Pt), shift, np.max(np.abs(Kl)), rho))
    g = np.array(gs); T = np.array(ts)
    # 幂律拟合（自由指数），分窗
    for lo, hi in ((-7, -4.5), (-5.5, -3.3), (-7, -3.3)):
        sel = (np.log10(T) >= lo) & (np.log10(T) <= hi)
        if sel.sum() < 3:
            continue
        cf = np.polyfit(np.log10(T[sel]), np.log10(np.abs(g[sel])), 1)
        print('  窗 $10^{%.1f}\dots10^{%.1f}$（%d 点）：$\rho$=%.4f  $\beta$=%.3e  残差rms=%.1e'
              % (lo, hi, sel.sum(), cf[0], 10 ** cf[1], np.std(np.polyval(cf, np.log10(T[sel]))
                                                             - np.log10(np.abs(g[sel])))))
    # 双项模型 $g=\beta_2 t^2+\beta_p t^p$：扫 $p$，线性最小二乘解 $(\beta_2,\beta_p)$
    best = None
    for p in np.linspace(1.05, 2.6, 311):
        M = np.vstack([T ** 2, T ** p]).T
        cf, rs = np.linalg.lstsq(M, g, rcond=None)[:2]
        rel = float(np.sqrt(np.sum((M @ cf - g) ** 2)) / np.abs(g).sum())
        if best is None or rel < best[0]:
            best = (rel, p, cf[0], cf[1])
    print('  全窗双项 $g=b t^2+c t^p$：$p$=%.3f  $b$=%.5f  $c$=%.3e  相对残差=%.2e'
          % (best[1], best[2], best[3], best[0]))
    # 只看深端一半，检验 $b$ 是否可辨识
    for cut in (-7, -6.3):
        sel = np.log10(T) <= cut
        best = None
        for p in np.linspace(1.05, 2.6, 311):
            M = np.vstack([T[sel] ** 2, T[sel] ** p]).T
            cf, rs = np.linalg.lstsq(M, g[sel], rcond=None)[:2]
            rel = float(np.sqrt(np.sum((M @ cf - g[sel]) ** 2)) / np.abs(g[sel]).sum())
            if best is None or rel < best[0]:
                best = (rel, p, cf[0], cf[1])
        print('  仅 $t\le10^{%.1f}$：%d 点  $p$=%.3f  $b$=%.5f  $c$=%.3e  相对残差=%.2e'
              % (cut, sel.sum(), best[1], best[2], best[3], best[0]))
