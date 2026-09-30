# -*- coding: utf-8 -*-
"""R59-D 实验 119：c71 的"最紧可用代价点"是 C = Clo + 0.05*(Chi-Clo) —— 这个 0.05 是自由参数。
固定 C 的设计集与网格，只改 f，看 增益/F(eps) 比值摆动多少；再在 eps=0.3 上把 Vg 从 14 加密到 28 看网格漂移。
用途：判 §56 的 "0.98-1.04" 是精度陈述还是读点选择。
"""
import sys, time, numpy as np
sys.stdout.reconfigure(encoding='utf-8')
lam = 0.1
N = 1401
xs = np.linspace(-6.5, 6.5, N); dxs = xs[1]-xs[0]
px = np.exp(-0.5*xs**2)/np.sqrt(2*np.pi)
_q = np.linspace(-7, 7, 4001); _p = np.exp(-0.5*_q**2)/np.sqrt(2*np.pi)


def p(s):
    print(s, flush=True)


def F_of(eps):
    gp = 1 + 3*eps*_q**2
    return 0.5*np.log2(np.trapezoid(_p*gp*gp, _q)/np.exp(np.trapezoid(_p*np.log(gp*gp), _q)))


def build(eps, q, c, V, nz=1201):
    g = xs + eps*xs**3
    ph = c*(xs + q*xs**3)
    zg = np.linspace(ph.min()-6*np.sqrt(V), ph.max()+6*np.sqrt(V), nz)
    K = np.exp(-0.5*(zg[None, :]-ph[:, None])**2/V)/np.sqrt(2*np.pi*V)
    m = (px*dxs)[:, None]
    pz = (K*m).sum(0); pz = np.maximum(pz, 1e-300)
    post = K*m/pz[None, :]
    mg = post.T@g; mg2 = post.T@(g*g)
    cost = float(np.trapezoid(pz*(mg2-mg*mg/(1+lam)), zg))
    Hz = -np.trapezoid(pz*np.log(pz), zg)
    rate = float((Hz-0.5*np.log(2*np.pi*np.e*V))/np.log(2))
    return cost, rate


def rate_at(A, C):
    o = np.argsort(A[:, 0]); Cc = A[o, 0]; Rr = A[o, 1]
    if C < Cc.min() or C > Cc.max():
        return None
    return float(np.interp(C, Cc, Rr))


CL = (0.6, 0.8, 1.0, 1.3, 1.7)
QS = (0.1, 0.25, 0.45, 0.7)
FS = (0.02, 0.05, 0.10, 0.20, 0.35)
t0 = time.time()
for eps in (0.10, 0.20, 0.30):
    F = F_of(eps)
    for nV, tag in ((14, 'Vg14(c71 原样)'),):
        Vg = 10.0**np.linspace(-2.2, 0.6, nV)
        lin = [np.array([build(eps, 0.0, c, V) for V in Vg]) for c in CL]
        nl = [np.array([build(eps, q, c, V) for V in Vg]) for q in QS for c in CL]
        allF = lin + nl
        Clo = max(A[:, 0].min() for A in allF); Chi = min(A[:, 0].max() for A in allF)
        p('\n eps=%.2f  F(eps)=%.6f  共同代价区间 [%.5f, %.5f]  (%s, %.0fs)'
          % (eps, F, Clo, Chi, tag, time.time()-t0))
        p('    f(读点)     C          R_lin    R_nl     增益      增益/F')
        for f in FS:
            C = Clo + f*(Chi-Clo)
            Rl = [r for r in (rate_at(A, C) for A in lin) if r is not None]
            Rn = [r for r in (rate_at(A, C) for A in nl) if r is not None]
            if not Rl or not Rn:
                p('    %5.2f   %9.5f     (无可比读数)' % (f, C)); continue
            gain = min(Rl)-min(Rn)
            p('    %5.2f   %9.5f   %7.4f  %7.4f  %8.5f   %7.3f %s'
              % (f, C, min(Rl), min(Rn), gain, gain/F, '<- c71 用的就是这个' if f == 0.05 else ''))
    if eps == 0.30:
        Vg = 10.0**np.linspace(-2.2, 0.6, 28)
        lin = [np.array([build(eps, 0.0, c, V) for V in Vg]) for c in CL]
        nl = [np.array([build(eps, q, c, V) for V in Vg]) for q in QS for c in CL]
        allF = lin + nl
        Clo = max(A[:, 0].min() for A in allF); Chi = min(A[:, 0].max() for A in allF)
        p('\n eps=0.30  网格加密对照 (Vg28, %.0fs)' % (time.time()-t0))
        p('    f(读点)     C          R_lin    R_nl     增益      增益/F')
        for f in FS:
            C = Clo + f*(Chi-Clo)
            Rl = [r for r in (rate_at(A, C) for A in lin) if r is not None]
            Rn = [r for r in (rate_at(A, C) for A in nl) if r is not None]
            if not Rl or not Rn:
                continue
            gain = min(Rl)-min(Rn)
            p('    %5.2f   %9.5f   %7.4f  %7.4f  %8.5f   %7.3f' % (f, C, min(Rl), min(Rn), gain, gain/F))
p('done %.0fs' % (time.time()-t0))
