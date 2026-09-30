import numpy as np

lam = 0.1
eps = 0.3
F = 0.474891
N = 1401
xs = np.linspace(-6.5, 6.5, N)
dxs = xs[1] - xs[0]
px = np.exp(-0.5 * xs ** 2) / np.sqrt(2 * np.pi)
g = xs + eps * xs ** 3
out = []
p = lambda s: out.append(str(s))


def build(q, c, V, nz=1201):
    ph = c * (xs + q * xs ** 3)
    zg = np.linspace(ph.min() - 6 * np.sqrt(V), ph.max() + 6 * np.sqrt(V), nz)
    K = np.exp(-0.5 * (zg[None, :] - ph[:, None]) ** 2 / V) / np.sqrt(2 * np.pi * V)
    m = (px * dxs)[:, None]
    pz = (K * m).sum(0)
    pz = np.maximum(pz, 1e-300)
    post = K * m / pz[None, :]
    mg = post.T @ g
    mg2 = post.T @ (g * g)
    cost = float(np.trapezoid(pz * (mg2 - mg * mg / (1 + lam)), zg))
    Hz = -np.trapezoid(pz * np.log(pz), zg)
    rate = float((Hz - 0.5 * np.log(2 * np.pi * np.e * V)) / np.log(2))
    return cost, rate


CS = (0.5, 0.8, 1.2, 1.8, 2.5, 3.5)
CL = (0.6, 0.8, 1.0, 1.3, 1.7, 2.2)
Vg = 10.0 ** np.linspace(-2.2, 0.6, 16)
QS_BASE = (0.1, 0.25, 0.45, 0.7)
QS_ADD = (0.2, 0.3, 0.35, 0.4)


def rate_at(Fc, C):
    o = np.argsort(Fc[:, 0])
    Cc = Fc[o, 0]
    Rr = Fc[o, 1]
    if C < Cc.min() or C > Cc.max():
        return None
    return float(np.interp(C, Cc, Rr))


lin = {(0.0, c): np.array([build(0.0, c, V) for V in Vg]) for c in CL}
nlb = {(q, c): np.array([build(q, c, V) for V in Vg]) for q in QS_BASE for c in CL}
nla = {(q, c): np.array([build(q, c, V) for V in Vg]) for q in QS_ADD for c in CL}
p('== q-lattice test on C\'s execution stage (g=x+0.3x^3, lam=0.1, F=0.474891) ==')
p('   QS_BASE = C\'s q grid {0.1,0.25,0.45,0.7}   |   +QS_ADD = {0.2,0.3,0.4,0.35}  (0.3 = eps is NOT in his grid)')
p('')
p('   C      R_lin     R_nl(base) gain(base)  rel     R_nl(base+add) gain(all)  rel     best q   q=0.3 gain')
for C in CS:
    rl = [rate_at(Fc, C) for Fc in lin.values()]
    rl = [r for r in rl if r is not None]
    if not rl:
        p('   %.1f   (lin envelope does not cover C)' % C)
        continue
    Rl = min(rl)
    rb = [rate_at(Fc, C) for Fc in nlb.values()]
    ra = [rate_at(Fc, C) for Fc in {**nlb, **nla}.values()]
    rb = [r for r in rb if r is not None]
    ra = [r for r in ra if r is not None]
    if not rb:
        p('   %.1f   (nl envelope does not cover C)' % C)
        continue
    gb = Rl - min(rb)
    ga = Rl - min(ra)
    per_q = {}
    for dic in (nlb, nla):
        for (qk, ck), Fc in dic.items():
            r = rate_at(Fc, C)
            if r is None:
                continue
            gv = Rl - r
            if qk not in per_q or gv > per_q[qk]:
                per_q[qk] = gv
    bq = max(per_q, key=lambda k: per_q[k])
    p('   %5.2f  %8.4f  %8.4f  %+7.4f  %6.3fx   %8.4f  %+7.4f  %6.3fx   q=%.2f(+%.4f)  q0.3=%s'
      % (C, Rl, min(rb), gb, gb / F, min(ra), ga, ga / F, bq, per_q[bq],
         ('%+.4f' % per_q[0.3]) if 0.3 in per_q else 'n/a'))
p('')
p('== per-q gain at each C (base grid marked *) ==')
tags = sorted({q for q, _ in {**nlb, **nla}.keys()})
p('   C    ' + '  '.join('%s%.2f' % ('*' if q in QS_BASE else ' ', q) for q in tags))
for C in CS:
    rl = [rate_at(Fc, C) for Fc in lin.values()]
    rl = [r for r in rl if r is not None]
    if not rl:
        continue
    Rl = min(rl)
    row = []
    for q in tags:
        cand = [rate_at(Fc, C) for (qq, cc), Fc in {**nlb, **nla}.items() if qq == q]
        cand = [x for x in cand if x is not None]
        r = min(cand) if cand else None
        row.append('     n/a' if r is None else '%+8.4f' % (Rl - r))
    p('   %5.2f ' % C + ' '.join(row))
open('p0/e114_out.txt', 'w', encoding='utf-8').write('\n'.join(out) + '\n')
print('done')
