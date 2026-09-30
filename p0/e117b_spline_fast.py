# -*- coding: utf-8 -*-
"""R59-D 快速版：只跑 eps=0.3，把 c68 的样条行按【修正配对】重算。
对照：buggy 配对 (cost=spline(bx), rate=spline(theta=0)) vs corrected (cost,rate 同 bx)。
"""
import sys, time, numpy as np
from scipy.optimize import minimize
sys.stdout.reconfigure(encoding='utf-8')
lam = 0.1
N = 1201
xs = np.linspace(-6.0, 6.0, N); dxs = xs[1] - xs[0]
px = np.exp(-0.5 * xs**2) / np.sqrt(2*np.pi)
out = []


def p(s):
    print(s); out.append(str(s))


def power_match(phi):
    v = np.trapezoid(px*phi**2, xs) - np.trapezoid(px*phi, xs)**2
    return phi*np.sqrt(1.0/max(v, 1e-12))


def stage(phi, g, V, nz=1001):
    phi = np.asarray(phi, float)
    zg = np.linspace(phi.min()-6*np.sqrt(V), phi.max()+6*np.sqrt(V), nz)
    K = np.exp(-0.5*(zg[None, :]-phi[:, None])**2/V)/np.sqrt(2*np.pi*V)
    m = (px*dxs)[:, None]
    pz = (K*m).sum(0); pz = np.maximum(pz, 1e-300)
    post = K*m/pz[None, :]
    mg = post.T@g; mg2 = post.T@(g*g)
    cost = float(np.trapezoid(pz*(mg2 - mg*mg/(1+lam)), zg))
    Hz = -np.trapezoid(pz*np.log(pz), zg)
    rate = float((Hz - 0.5*np.log(2*np.pi*np.e*V))/np.log(2))
    return cost, rate


def spline(theta):
    sp = np.log1p(np.exp(np.clip(theta, -30, 30)))
    base = np.linspace(xs.min(), xs.max(), 13)
    vals = np.concatenate([[0.0], np.cumsum(sp)])
    vals = vals - np.interp(0.0, base, vals)
    return power_match(np.interp(xs, base, vals))


def cubic(q, c):
    return power_match(c*(xs + q*xs**3))


def rate_at(F, C):
    o = np.argsort(F[:, 0]); Cc = F[o, 0]; Rr = F[o, 1]
    if C < Cc.min() or C > Cc.max():
        return None
    return float(np.interp(C, Cc, Rr))


Vg = 10.0**np.linspace(-1.6, 0.4, 6)
eps = 0.3
g = xs + eps*xs**3

s0 = spline(np.zeros(12))
pm = power_match(xs)
p('== (0) theta=0 样条是否就是线性传感器 ==')
p('  12 个 softplus 增量 = %s' % np.round(np.log1p(np.exp(np.zeros(12))), 16)[:3])
p('  pearson r(xs, s0)          = %.14f' % np.corrcoef(xs, s0)[0, 1])
p('  max|s0 - pm(x)|            = %.3e' % np.max(np.abs(s0-pm)))

F0 = np.array([stage(pm, g, V) for V in Vg])
F1 = min([np.array([stage(cubic(q, c), g, V) for V in Vg])
          for q in np.linspace(-0.5, 0.5, 7) for c in (0.6, 0.9, 1.3, 1.8)],
         key=lambda A: A[:, 0].min())

rng = np.random.default_rng(7)
starts = [np.zeros(12)] + [rng.normal(0, 0.8, 12) for _ in range(3)]
bug, cor = [], []
t0 = time.time()
for i, V in enumerate(Vg):
    bx, bc = None, np.inf
    for s in starts:
        r = minimize(lambda th: stage(spline(th), g, V)[0], s, method='Nelder-Mead',
                     options={'maxiter': 400, 'xatol': 1e-3, 'fatol': 1e-6})
        if r.fun < bc:
            bc, bx = r.fun, r.x
    copt, ropt = stage(spline(bx), g, V)
    rz = stage(s0, g, V)[1]
    bug.append((copt, rz)); cor.append((copt, ropt))
    p('  V=%.5f  cost=%.5f | rate(bx)=%.4f  rate(theta=0)=%.4f  |dr|=%+.4f   [%.0fs]'
      % (V, copt, ropt, rz, ropt-rz, time.time()-t0))
bug = np.array(bug); cor = np.array(cor)

C0 = F1[:, 0].max()*0.5 + F1[:, 0].min()*0.5
r_lin, r_cub = rate_at(F0, C0), rate_at(F1, C0)
r_b, r_c = rate_at(bug, C0), rate_at(cor, C0)
p('\n== (1) 匹配代价 C0=%.5f 处的率读数 ==' % C0)
p('  线性传感器 R_lin = %s' % ('none' if r_lin is None else '%.4f' % r_lin))
p('  三次族     R_cub = %s' % ('none' if r_cub is None else '%.4f' % r_cub))
p('  样条(buggy 配对, 复现 c68) R_spl_bug = %s' % ('none' if r_b is None else '%.4f' % r_b))
p('  样条(修正配对)             R_spl_cor = %s' % ('none' if r_c is None else '%.4f' % r_c))
if None not in (r_lin, r_cub, r_b, r_c):
    p('  样条-三次  buggy = %+.4f | corrected = %+.4f' % (r_b-r_cub, r_c-r_cub))
    p('  样条-线性  buggy = %+.4f (=> 与线性同值即 bug 证据) | corrected = %+.4f' % (r_b-r_lin, r_c-r_lin))
open('p0/e117b_out.txt', 'w', encoding='utf-8').write('\n'.join(out)+'\n')
p('done %.0fs' % (time.time()-t0))
