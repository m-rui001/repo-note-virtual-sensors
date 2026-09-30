import numpy as np
from scipy.linalg import solve_discrete_are, eigh

_src = open('p0/exp_c_audit.py', encoding='utf-8').read().split("print(r'== E53")[0]
_ns = {'__name__': 'p'}
exec(compile(_src, 'p0/exp_c_audit.py[preamble]', 'exec'), _ns)
A, W, TH, JC, n = _ns['A'], _ns['W'], _ns['TH'], _ns['JC'], _ns['n']
sym = _ns['sym']

out = []
p = lambda s: out.append(str(s))


def curve(Z, s):
    r = Z.shape[1]
    Ir = np.eye(r)
    C = np.sqrt(s) * Z.T
    Pm = sym(solve_discrete_are(A.T, C.T, W, Ir))
    Sm = C @ Pm @ C.T + Ir
    Lk = Pm @ C.T @ np.linalg.inv(Sm)
    Pp = sym(Pm - Lk @ Sm @ Lk.T)
    I = 0.5 * np.log(np.linalg.det(Ir + C @ Pm @ C.T)) / np.log(2.0)
    return I, JC + np.trace(TH @ Pp)


def at_rate(Z, target):
    lo, hi = -9.0, 6.0
    if not (curve(Z, 10 ** lo)[0] < target < curve(Z, 10 ** hi)[0]):
        return None, None
    for _ in range(90):
        mid = 0.5 * (lo + hi)
        if curve(Z, 10 ** mid)[0] < target:
            lo = mid
        else:
            hi = mid
    s = 10 ** (0.5 * (lo + hi))
    return s, curve(Z, s)[1]


wth, Vth = eigh(TH)
R_EXP = 1.168539
dI = 2e-3
t = R_EXP + dI
p('== b vs K anatomy at the matched near-floor rate  (dI = %.4f bit, JC = %.6f) ==' % (dI, JC))
p('   b   = s*(D-JC) at s=1e-6      (the quantity my [G1] ranked on)')
p('   K   = (D-JC)*dI  at I=R_exp+dI (the quantity that orders D at matched rate)')
p('   a   = K/b                      (the rate-side slope factor I ignored)')
p('')
p(' rank  leg       s(min)      D(min)       b(min)     K(min)      a(min)   |  s(max)      D(max)       b(max)     K(max)      a(max)')
for rk in (1, 2, 3):
    Zmin, Zmax = Vth[:, :rk], Vth[:, n - rk:]
    rec = {}
    for nm, Z in (('min', Zmin), ('max', Zmax)):
        s, D = at_rate(Z, t)
        I6, D6 = curve(Z, 1e-6)
        b = 1e-6 * (D6 - JC)
        K = (D - JC) * dI
        rec[nm] = (s, D, b, K, K / b)
    a = rec['min']
    x = rec['max']
    p(' %d      min %11.4e %12.4f %11.4f %11.4f %11.4f   |  max %11.4e %12.4f %11.4f %11.4f %11.4f'
      % (rk, a[0], a[1], a[2], a[3], a[4], x[0], x[1], x[2], x[3], x[4]))
    p('       ratios:  b_min/b_max=%.4f   a_min/a_max=%.4f   K_min/K_max=%.4f   D_min/D_max=%.4f'
      % (a[2] / x[2], a[4] / x[4], a[3] / x[3], a[1] / x[1]))
p('')
p('== frontier-law cross-check:  D - JC  ==  b / s  at the matched s ==')
for rk in (1, 2, 3):
    for nm in ('min', 'max'):
        Z = Vth[:, :rk] if nm == 'min' else Vth[:, n - rk:]
        s, D = at_rate(Z, t)
        I6, D6 = curve(Z, 1e-6)
        b = 1e-6 * (D6 - JC)
        pred = b / s
        p('  rank %d  %s:  b/s = %.6e   D-JC(obs) = %.6e   rel = %+.4f%%'
          % (rk, nm, pred, D - JC, 100.0 * (pred / (D - JC) - 1.0)))
open('p0/e112_out.txt', 'w', encoding='utf-8').write('\n'.join(out) + '\n')
print('done')
