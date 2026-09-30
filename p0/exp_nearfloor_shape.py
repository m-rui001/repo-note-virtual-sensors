"""E22：近地板处 Δ(D) 的发散形状——用单参数精度族做**无优化**的参数化测量。

动机（对 #9 的自我更正）：C 的加密网格 c16 最左一对本地斜率 −0.516，我的粗网格同一 x 区间给
−0.528，两套独立实现差 2%。所以 p≈0.5 不是我的网格伪影。但"局部斜率在最左一对 ≈0.5"
与"极限指数 =0.5"是两件事，需要看 x→0 的**形状**而不是再拟合一个 p。

方法：绕开约束优化。取 $\Gamma=\gamma I_r$（$\gamma$ 从 $10^0$ 到 $10^{14}$），
用 `rate_cost` 直接算该精度下的 $(J,I)$（$J$ 随 $\gamma$ 单调下降到地板
$D_{\min}=D_{\min}^{unc}+\Phi$），再用 SDP 算 $I_{unc}(J)$。于是
    $x(\gamma)=J(\gamma)-D_{\min}$，$\Delta(\gamma)=I(\gamma)-I_{unc}(J(\gamma))$
给出一条**参数化的** $\Delta$–$x$ 曲线，没有任何求解器噪声。
若 $x\to0$ 时 $\Delta\sim x^{-p}$，则 $\log\Delta$ 对 $\log x$ 的局部斜率应趋于常数；
若 $\Delta\sim\log(1/x)$（我预计的另一支：$\gamma\to\infty$ 时 $I\sim r\log\gamma$ 而
$J-D_{\min}\sim c/\gamma$，则 $\Delta\sim\log(1/x)$），局部斜率应趋于 0。
"""
import sys
import numpy as np
sys.stdout.reconfigure(encoding='utf-8')
from p0.p0_replicate_letter import A, W, n, sym, ln2, ctrl, rate_cost, unconstrained_sdp
from p0.exp_authoritative import noiseless
from scipy.linalg import schur

Pc, K, Th = ctrl(np.eye(n))
unc = float(np.trace(W @ Pc))
T, U, sdim = schur(A, output='real', sort=lambda a: abs(a) < 1.0)
nu = n - sdim
iT = list(range(n - nu, n))
F2 = U[:, iT].T
r2 = F2.shape[0]
Pinf, _, _ = noiseless(F2, A, W, iters=80000, tol=1e-15)
phi = float(np.trace(Th @ Pinf))
floor = unc + phi

if __name__ == '__main__':
    print('F2：D_min^unc=%.4f  Φ=%.5f  地板=%.5f  rank F=%d' % (unc, phi, floor, r2))
    gams = np.logspace(0, 14, 29)
    xs, ds, js, iss = [], [], [], []
    for g in gams:
        Gam = g * np.eye(r2)
        try:
            I, J, _, _ = rate_cost(F2, Gam, Th, Pc)
        except Exception:
            continue
        if not (np.isfinite(I) and np.isfinite(J)) or J <= floor * (1 + 1e-12):
            continue
        u = unconstrained_sdp(J, np.eye(n), Th, Pc)
        if u is None:
            continue
        xs.append(J - floor); ds.append(I - u['I_true']); js.append(J); iss.append(I)
    xs = np.array(xs); ds = np.array(ds)
    print('%12s %12s %12s %12s %10s' % ('γ', 'x=J−地板', 'Δ', '本地斜率', 'log10 x'))
    prev = None
    for g, x, d in zip(gams[:len(xs)], xs, ds):
        sl = '' if prev is None else '%10.3f' % ((np.log(d) - np.log(prev[1])) / (np.log(x) - np.log(prev[0])))
        print('%12.3e %12.5f %12.5f %10s %10.3f' % (g, x, d, sl, np.log10(x)))
        prev = (x, d)
    m = xs < np.percentile(xs, 40)
    if m.sum() >= 3:
        b = np.polyfit(np.log(xs[m]), np.log(ds[m]), 1)
        print('\n最左 40%% 点（x∈[%.3g, %.3g]）的拟合指数 = %.3f' % (xs[m].min(), xs[m].max(), -b[0]))
    m2 = xs > np.percentile(xs, 60)
    if m2.sum() >= 3:
        b2 = np.polyfit(np.log(xs[m2]), np.log(ds[m2]), 1)
        print('最右 40%% 点（x∈[%.3g, %.3g]）的拟合指数 = %.3f' % (xs[m2].min(), xs[m2].max(), -b2[0]))
