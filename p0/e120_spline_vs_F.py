# -*- coding: utf-8 -*-
"""R59-D 实验 120（§73-D 的判决实验）：把 12 参数单调样条族按【修正配对】加进 C 的 c69 口径，
看 C=0.5 处的增益能否越过 F=0.474891。
判据：gain_spl > 0.474891 => §54 的"包络"与 §71 的收窄版一起作废，只能写"三次族的上包络"；
      否则 => "包络"在 eps=0.3 保留，但必须加族限定（样条族在该点没越过）。
口径全部照 c69：N=1401 xs in [-6.5,6.5], nz=1201, Vg=10^linspace(-2.2,0.6,16),
R = min over 该族各设计的 interp(代价->率) @ C，增益 = R_lin - R_nl。
"""
import sys, time, numpy as np
from scipy.optimize import minimize
sys.stdout.reconfigure(encoding='utf-8')
lam = 0.1; eps = 0.3; FC = 0.474891
N = 1401
xs = np.linspace(-6.5, 6.5, N); dxs = xs[1]-xs[0]
px = np.exp(-0.5*xs**2)/np.sqrt(2*np.pi)
g = xs + eps*xs**3
C = 0.5
Vg = 10.0**np.linspace(-2.2, 0.6, 16)


def p(s):
    print(s, flush=True)


def stage(phi, V, nz=1201):
    phi = np.asarray(phi, float)
    zg = np.linspace(phi.min()-6*np.sqrt(V), phi.max()+6*np.sqrt(V), nz)
    K = np.exp(-0.5*(zg[None, :]-phi[:, None])**2/V)/np.sqrt(2*np.pi*V)
    m = (px*dxs)[:, None]
    pz = (K*m).sum(0); pz = np.maximum(pz, 1e-300)
    post = K*m/pz[None, :]
    mg = post.T@g; mg2 = post.T@(g*g)
    cost = float(np.trapezoid(pz*(mg2-mg*mg/(1+lam)), zg))
    Hz = -np.trapezoid(pz*np.log(pz), zg)
    rate = float((Hz-0.5*np.log(2*np.pi*np.e*V))/np.log(2))
    return cost, rate


def curve(phi):
    return np.array([stage(phi, V) for V in Vg])


def rate_at(A, Cc):
    o = np.argsort(A[:, 0]); cc = A[o, 0]; rr = A[o, 1]
    if Cc < cc.min() or Cc > cc.max():
        return None
    return float(np.interp(Cc, cc, rr))


def pmatch(phi):
    v = np.trapezoid(px*phi**2, xs) - np.trapezoid(px*phi, xs)**2
    return phi*np.sqrt(1.0/max(v, 1e-12))


def spline(theta):
    sp = np.log1p(np.exp(np.clip(theta, -30, 30)))
    base = np.linspace(xs.min(), xs.max(), 13)
    vals = np.concatenate([[0.0], np.cumsum(sp)])
    vals = vals - np.interp(0.0, base, vals)
    return pmatch(np.interp(xs, base, vals))


lin = {(0.0, c): curve(0.0*(xs) + c*xs) for c in (0.6, 0.8, 1.0, 1.3, 1.7, 2.2)}
nlc = {(q, c): curve(c*(xs + q*xs**3)) for q in (0.1, 0.25, 0.45, 0.7)
       for c in (0.6, 0.8, 1.0, 1.3, 1.7, 2.2)}
Rl = [r for r in (rate_at(A, C) for A in lin.values()) if r is not None]
Rn = [r for r in (rate_at(A, C) for A in nlc.values()) if r is not None]
gain_cub = min(Rl) - min(Rn)
p('复现 c69 @ C=0.5 : R_lin=%.4f  R_nl(三次)=%.4f  增益=%.4f  相对 F=%.3f x'
  % (min(Rl), min(Rn), gain_cub, gain_cub/FC))

t0 = time.time()
rng = np.random.default_rng(11)
OPT_V = [10.0**-1.4, 10.0**-1.0, 10.0**-0.6, 10.0**-0.2]
thetas = []
for V in OPT_V:
    bx, bc = None, np.inf
    for s in range(2):
        th0 = np.zeros(12) if s == 0 else rng.normal(0, 0.8, 12)
        r = minimize(lambda th: stage(spline(th), V)[0], th0, method='Nelder-Mead',
                     options={'maxiter': 260, 'xatol': 1e-3, 'fatol': 1e-6})
        if r.fun < bc:
            bc, bx = r.fun, r.x
    thetas.append(bx)
    p('  优化 V=%.4f -> cost=%.5f  非零增量段 %d/12   [%.0fs]'
      % (V, bc, int(np.sum(np.abs(bx) > 1e-3)), time.time()-t0))

best = None
for i, th in enumerate(thetas):
    A = curve(spline(th))
    r = rate_at(A, C)
    if r is None:
        p('  样条设计 %d（V_opt）在 C=0.5 不可读，代价区间 [%.4f,%.4f]' % (i, A[:, 0].min(), A[:, 0].max()))
        continue
    p('  样条设计 %d : R@C=0.5 = %.4f  增益 = %.4f  相对 F = %.3f x'
      % (i, r, min(Rl)-r, (min(Rl)-r)/FC))
    if best is None or r < best:
        best = r
if best is not None:
    gs = min(Rl) - best
    p('\n== 判决 ==')
    p('  三次族增益 %.4f (%.3f x F)  |  加入样条族后 %.4f (%.3f x F)' % (gain_cub, gain_cub/FC, gs, gs/FC))
    p('  越过 F=0.474891 ? %s   => %s' % ('是 ✗ 包络口径作废' if gs > FC else '否 ✓ 包络在 eps=0.3 保留（需加族限定）',
                                        '只能写"三次族的上包络"' if gs > FC else '§54/§71 的"包络"字样可保留，但须写明含 12 维族'))
p('done %.0fs' % (time.time()-t0))
