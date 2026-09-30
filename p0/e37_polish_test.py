import sys
import numpy as np
sys.stdout.reconfigure(encoding='utf-8')
import p0.exp_plants_oos as E35
import p0.exp_ladder as X37
import p0.exp_C_exact as X36

for nm, pl, F in E35.plants():
    if 'T-F2' not in nm and '随机正交' not in nm: continue
    C = X37.C_of(pl, F)
    rho, elam, U = X36.rank_C(C)
    fms, ee, uu, _ = X37.fams_of(pl, F, C, 0.0)
    f0 = [f for f in fms if f['a'] == rho][0]
    Fv = f0['F']
    Vref = None
    print('\n===', nm, 'lam(C)=', np.array2string(elam, precision=4))
    for pol in (600, 6000, 60000):
        line = []
        for t in (1e-5, 1e-7, 1e-9):
            V = t * np.diag(1.0 / f0['w'])
            x, I, Pp, Pt_, res = X37.pv(pl, Fv, V, 0.0, polish=pol)
            M = np.linalg.inv(Pt_)
            Pi = M + Fv.T @ np.linalg.inv(V) @ Fv
            dev = np.linalg.norm(np.linalg.inv(Pi) - Pp) / np.linalg.norm(Pp)
            # 参考解：直接 scipy DARE 的后验
            Pref = pl.P(Fv, V)
            dp = np.linalg.norm(Pp - Pref) / np.linalg.norm(Pref)
            line.append('t=%.0e res=%.1e dev=%.1e |P-scipy|=%.1e' % (t, res, dev, dp))
        print('  polish=%-6d ' % pol + '  |  '.join(line))
