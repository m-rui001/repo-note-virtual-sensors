"""E20：从可行域那一侧测 D_min(V)，独立于 noiseless 定点。

先前设计错了：本例 $J_{ol}=+\infty$（A 不稳定、$S=0$ 时 $\mathrm{tr}(\Theta P)$ 发散），
所以 $J\le D$ 在每个有限 $D$ 上都吃紧（E19 尾部实测 $J_{trv}=D$ 到 8000），
"大 $D$ 时代价收敛到地板"这个测试是不可能成立的。

改成正确的方向：**从下往上撞**。对 Prop D 预言的地板
    F2（封闭）：$D_{\min}^{unc}+\Phi=31.4833+14.6398=46.1231$
    F3（不封闭）：$31.4833+5.1199=36.6032$
在 $D<$ 地板处问 $\mathcal K_F$ 上的求解器：它能否找到 $J\le D$ 的点？
若 Prop D 对，则找不到，且它所能达到的最小 $J$ 会从上方贴住地板。
求解器只用到 $\Theta,\ \widetilde P$ 与锥约束，完全不走我的 $V\to0$ 定点，所以这是独立的一路。
"""
import sys
import numpy as np
sys.stdout.reconfigure(encoding='utf-8')
from scipy.optimize import minimize
from p0.p0_replicate_letter import A, W, n, sym, ln2, ctrl, rate_cost, init_gamma
from scipy.linalg import schur

Pc, K, Th = ctrl(np.eye(n))
unc = float(np.trace(W @ Pc))
T, U, sdim = schur(A, output='real', sort=lambda a: abs(a) < 1.0)
nu = n - sdim
iT = list(range(n - nu, n))
Fs = [('F2 尾部2(封闭)', U[:, iT].T, 46.1231),
      ('F3 尾部3(不封闭)', U[:, [1, 2, 3]].T, 36.6032)]


def minJ(F, gam0, D, budget=60000):
    """在 K_F 上最小化 J（不看 D），带硬可行性报告。返回 (J_min, I, 是否 J≤D)."""
    r = F.shape[0]

    def pack(G):
        G = sym(G) + 1e-10 * np.eye(r)
        L = np.linalg.cholesky(G)
        return np.array([L[i, j] for i in range(r) for j in range(i + 1)])

    def unpack(p):
        L = np.zeros((r, r))
        k = 0
        for i in range(r):
            for j in range(i + 1):
                L[i, j] = p[k]; k += 1
        return sym(L @ L.T)

    def ij(p):
        try:
            I, J, _, _ = rate_cost(F, unpack(p), Th, Pc)
        except Exception:
            return 1e6, 1e6
        if not (np.isfinite(I) and np.isfinite(J)):
            return 1e6, 1e6
        return I, J

    best = None
    rng = np.random.default_rng(0)
    starts = [gam0] + [np.eye(r) * s for s in (1.0, 10.0, 100.0, 1e4)]
    for _ in range(8):
        L = rng.standard_normal((r, r)) * 0.7
        starts.append(sym(L @ L.T) + 0.05 * np.eye(r))
    for G0 in starts:
        try:
            p0 = pack(G0)
        except np.linalg.LinAlgError:
            continue
        res = minimize(lambda p: ij(p)[1], p0, method='Nelder-Mead',
                       options=dict(maxfev=budget, xatol=1e-9, fatol=1e-12))
        v = ij(res.x)[1]
        if best is None or v < best:
            best = v
    I_at = ij(np.zeros(0)) if best is None else None
    return best


if __name__ == '__main__':
    print('D_min^unc = %.4f' % unc)
    for tag, F, floor in Fs:
        print('\n[%s]  Prop D 预言地板 = %.4f   (Φ = %.4f)' % (tag, floor, floor - unc))
        Gam0 = init_gamma(F, np.eye(n) * 3.0)
        jm = minJ(F, Gam0, None)
        print('  K_F 上能达到的最小 J（13 个起点，只看 J）= %.4f   与地板之差 = %+.4f' % (jm, jm - floor))
        print('  相对 Φ 的误差 = %+.2f%%' % (100 * (jm - floor) / (floor - unc)))
        for eps in (0.5, 0.2, 0.05):
            D = floor - eps
            print('  问：D=%.4f（地板下方 %.2f）→ 可达最小 J=%.4f ⟹ 可行？%s'
                  % (D, eps, jm, '是（Prop D 死）' if jm <= D else '否 ✓'))
        for eps in (0.05, 0.5):
            D = floor + eps
            print('  问：D=%.4f（地板上方 %.2f）→ 同上最小 J=%.4f ⟹ 可行？%s'
                  % (D, eps, jm, '是 ✓' if jm <= D else '否（求解器没够到，非否证）'))
