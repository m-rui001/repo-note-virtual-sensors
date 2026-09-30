# -*- coding: utf-8 -*-
"""R59-D 实验 118：增益/F 比值的【误差棒】实测——把 C 的 6 点 V 网格加密到 21 点，
看同一 (eps, C) 上的增益移动多少。用于判 §56 的 1.044 / 1.004 是否显著大于 1，
也用于判我自己 §71 里"F 是包络"这句需要多强的 eps 限定。
"""
import sys, time, numpy as np
sys.stdout.reconfigure(encoding='utf-8')
lam = 0.1
N = 1201
xs = np.linspace(-6.0, 6.0, N); dxs = xs[1] - xs[0]
px = np.exp(-0.5 * xs**2) / np.sqrt(2*np.pi)
out = []


def p(s):
    print(s, flush=True); out.append(str(s))


def pm(phi):
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


def cubic(q, c):
    return pm(c*(xs + q*xs**3))


DESIGNS = [(q, c) for q in np.linspace(-0.5, 0.5, 7) for c in (0.6, 0.9, 1.3, 1.8)]


def front(fam, Vg):
    return np.array([stage(fam(V), g, V) for V in Vg])


def rate_at(F, C):
    o = np.argsort(F[:, 0]); Cc = F[o, 0]; Rr = F[o, 1]
    if C < Cc.min() or C > Cc.max():
        return None
    return float(np.interp(C, Cc, Rr))


def best_front(Vg):
    return min([front(lambda V, q=q, c=c: cubic(q, c), Vg) for q, c in DESIGNS],
               key=lambda A: A[:, 0].min())


F_DISP = {0.1: 0.127998, 0.2: 0.312347, 0.3: 0.474891, 0.45: 0.669372}
CASES = [(0.3, 0.500, 'C 的 §54 表首行 = 0.4649'), (0.1, 0.19982, 'C 的 §56 eps=0.1 = 0.13361'),
         (0.2, 0.31966, 'C 的 §56 eps=0.2 = 0.31373')]
GRIDS = {'6pt(c68/c69 原样)': 10.0**np.linspace(-1.6, 0.4, 6),
         '21pt(加密)': 10.0**np.linspace(-1.6, 0.4, 21)}
t0 = time.time()
p('== 增益对 V 网格的敏感性（同一族、同一设计集，只改网格密度）==')
for eps, C, tag in CASES:
    globals()['g'] = xs + eps*xs**3
    row = {}
    for nm, Vg in GRIDS.items():
        r_lin = rate_at(front(lambda V: pm(xs), Vg), C)
        r_nl = rate_at(best_front(Vg), C)
        row[nm] = (r_lin, r_nl, None if r_lin is None or r_nl is None else r_lin-r_nl)
    a, b = row['6pt(c68/c69 原样)'], row['21pt(加密)']
    d = '' if None not in (a[2], b[2]) else ' [不可读]'
    ratio_a = '' if (b[2] is None or eps not in F_DISP) else ' | 比值 %.3f->%.3f' % (a[2]/F_DISP[eps], b[2]/F_DISP[eps])
    p(' eps=%.2f (%s) C=%.4f : 增益 %.4f -> %.4f  漂移 %+.4f bit%s'
      % (eps, tag, C, a[2] if a[2] is not None else float('nan'),
         b[2] if b[2] is not None else float('nan'),
         (b[2]-a[2]) if (a[2] is not None and b[2] is not None) else float('nan'), ratio_a))
open('p0/e118_out.txt', 'w', encoding='utf-8').write('\n'.join(out)+'\n')
p('done %.0fs' % (time.time()-t0))
