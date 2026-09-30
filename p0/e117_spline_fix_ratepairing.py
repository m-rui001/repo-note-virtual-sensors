import numpy as np
from scipy.optimize import minimize
from scipy.stats import pearsonr

lam = 0.1
N = 1201
xs = np.linspace(-6.0, 6.0, N)
dxs = xs[1] - xs[0]
px = np.exp(-0.5 * xs ** 2) / np.sqrt(2 * np.pi)
out = []
p = lambda s: out.append(str(s))


def power_match(phi):
    v = np.trapezoid(px * phi ** 2, xs) - np.trapezoid(px * phi, xs) ** 2
    return phi * np.sqrt(1.0 / max(v, 1e-12))


def stage(phi, g, V, nz=1001):
    phi = np.asarray(phi, float)
    zg = np.linspace(phi.min() - 6 * np.sqrt(V), phi.max() + 6 * np.sqrt(V), nz)
    K = np.exp(-0.5 * (zg[None, :] - phi[:, None]) ** 2 / V) / np.sqrt(2 * np.pi * V)
    m = (px * dxs)[:, None]
    pz = (K * m).sum(0)
    pz = np.maximum(pz, 1e-300)
    post = K * m / pz[None, :]
    mg = post.T @ g
    mg2 = post.T @ (g * g)
    cost = float(np.trapezoid(pz * (mg2 - mg * mg / (1 + lam)), zg))
    Hz = -np.trapezoid(pz * np.log(pz), zg)
    return cost, float((Hz - 0.5 * np.log(2 * np.pi * np.e * V)) / np.log(2))


def cubic(q, c):
    return power_match(c * (xs + q * xs ** 3))


def spline(theta, K=12):
    sp = np.log1p(np.exp(np.clip(theta, -30, 30)))
    base = np.linspace(xs.min(), xs.max(), K + 1)
    vals = np.concatenate([[0.0], np.cumsum(sp)])
    vals = vals - np.interp(0.0, base, vals)
    return power_match(np.interp(xs, base, vals))


def rate_at(F, C):
    o = np.argsort(F[:, 0])
    Cc, Rr = F[o, 0], F[o, 1]
    if C < Cc.min() or C > Cc.max():
        return None
    return float(np.interp(C, Cc, Rr))


th0 = np.zeros(12)
s0 = spline(th0)
b = np.linspace(xs.min(), xs.max(), 13)
sp = np.log1p(np.exp(np.clip(th0, -30, 30)))
vals = np.concatenate([[0.0], np.cumsum(sp)])
p('== (0) theta=0 的"spline" 是什么形状 ==')
p('   归一化后与设计 x 的相关系数 r = %.10f   对 x 的三次项拟合残差 = %.3e'
  % (pearsonr(xs, s0)[0], np.abs(s0 - (np.polyval(np.polyfit(xs, s0, 1), xs))).max()))
p('   => 增量全相等（softplus(0)=ln2 每段），分段线性 + 等距节点 = 线性函数；power_match 后仍是线性传感器')

for eps in (0.0, 0.3):
    g = xs + eps * xs ** 3
    Vg = 10.0 ** np.linspace(-1.6, 0.4, 6)
    F1 = min([np.array([stage(cubic(q, c), g, V) for V in Vg])
              for q in np.linspace(-0.5, 0.5, 7) for c in (0.6, 0.9, 1.3, 1.8)],
             key=lambda A: A[:, 0].min())
    F0 = np.array([stage(power_match(xs), g, V) for V in Vg])
    C0 = F1[:, 0].max() * 0.5 + F1[:, 0].min() * 0.5
    rng = np.random.default_rng(7)
    bug, fix = [], []
    for V in Vg:
        best, bx = None, None
        for s in range(4):
            th_init = np.zeros(12) if s == 0 else rng.normal(0, 0.8, 12)
            r = minimize(lambda th: stage(spline(th), g, V)[0], th_init, method='Nelder-Mead',
                         options={'maxiter': 400, 'xatol': 1e-3, 'fatol': 1e-6})
            if best is None or r.fun < best:
                best, bx = r.fun, r.x
        cost_opt = best
        rate_opt = stage(spline(bx), g, V)[1]
        rate_zer = stage(spline(np.zeros(12)), g, V)[1]
        bug.append((cost_opt, rate_zer))
        fix.append((cost_opt, rate_opt))
    F2b, F2f = np.array(bug), np.array(fix)
    r0, r1 = rate_at(F0, C0), rate_at(F1, C0)
    r2b, r2f = rate_at(F2b, C0), rate_at(F2f, C0)
    p('')
    p('== eps=%.1f  匹配代价 C0=%.4f ==' % (eps, C0))
    p('   R_lin=%.4f   R_cubic=%.4f' % (r0, r1))
    p('   [复现 C 的写法：rate 用 theta=0 的设计]   R_spline = %s  => 样条-三次 = %+.4f bit'
      % ('%.4f' % r2b, r2b - r1))
    p('   [修正配对：rate 用被优化的那个 theta]      R_spline = %s  => 样条-三次 = %+.4f bit'
      % ('%.4f' % r2f, r2f - r1))
    p('   率节省(相对线性): 三次 %.4f | 样条(修正) %.4f | 常数 F=%.6f'
      % (r0 - r1, r0 - r2f, 0.474891))
open('p0/e117_out.txt', 'w', encoding='utf-8').write('\n'.join(out) + '\n')
print('done')
