import sys
import numpy as np
sys.stdout.reconfigure(encoding='utf-8')
import p0.exp_plants_oos as E35
import p0.exp_ladder as X37
import p0.exp_C_exact as X36

Pl = E35.Pl
for nm, pl, F in E35.plants():
    if 'T-F2' not in nm and 'PL' not in nm: continue
    C = X37.C_of(pl, F)
    rho, elam, U = X36.rank_C(C)
    Ci = U[:, :rho] @ np.diag(1.0/elam[:rho]) @ U[:, :rho].T
    fms, ee, uu, _ = X37.fams_of(pl, F, C, 14.639794977004 if 'T-F2' in nm else 0.037934348335)
    f0 = [f for f in fms if f['a'] == rho][0]
    print('\n===', nm, 'lam(C)=', elam)
    for t in (1e-5, 1e-7, 1e-9):
        V = t * np.diag(1.0/f0['w'])
        Fv = f0['F']
        gl, dev = X37.grads(pl, Fv, V)
        phi0 = 0.0
        h = 1e-4
        row = []
        for k, ((i, j), dI, dx) in enumerate(gl):
            if i != j: continue
            d = np.zeros_like(V); d[i, i] = h * V[i, i]
            I1 = X37.pv(pl, Fv, V + d, phi0)[1]; I2 = X37.pv(pl, Fv, V - d, phi0)[1]
            x1 = X37.pv(pl, Fv, V + d, phi0)[0]; x2 = X37.pv(pl, Fv, V - d, phi0)[0]
            fdI = (I1 - I2) / (2 * d[i, i]); fdx = (x1 - x2) / (2 * d[i, i])
            row.append('i=%d: dI %.4e vs fd %.4e | dx %.4e vs fd %.4e' % (i, dI, fdI, dx, fdx))
        print(' t=%.0e  dev=%.1e' % (t, dev)); [print('   ', s) for s in row]
