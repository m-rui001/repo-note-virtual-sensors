import numpy as np
from scipy.linalg import solve_discrete_are, schur, eigh

_src = open('p0/exp_c_audit.py', encoding='utf-8').read().split("print(r'== E53")[0]
_ns = {'__name__': 'p'}
exec(compile(_src, 'p0/exp_c_audit.py[preamble]', 'exec'), _ns)
A, W, TH, JC, n = _ns['A'], _ns['W'], _ns['TH'], _ns['JC'], _ns['n']
sym = _ns['sym']

out = []


def p(s):
    out.append(str(s))


def curve(Z, s):
    r = Z.shape[1]
    Ir = np.eye(r)
    C = np.sqrt(s) * Z.T
    Pm = sym(solve_discrete_are(A.T, C.T, W, Ir))
    Sm = C @ Pm @ C.T + Ir
    Lk = Pm @ C.T @ np.linalg.inv(Sm)
    Pp = sym(Pm - Lk @ Sm @ Lk.T)
    I = 0.5 * np.log(np.linalg.det(Ir + C @ Pm @ C.T)) / np.log(2.0)
    D = JC + np.trace(TH @ Pp)
    return I, D


def at_rate(Z, target):
    lo, hi = -9.0, 6.0
    flo = curve(Z, 10 ** lo)[0]
    fhi = curve(Z, 10 ** hi)[0]
    if not (flo < target < fhi):
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

p('== (0) 目标：Theta-min 与 Theta-max 两个设计的代价次序在哪个率上翻转 ==')
p(f'  JC={JC:.6f}  R_exp={R_EXP}  lambda(TH) 升序={np.round(wth, 4)}  （两个设计都是 rank=r 的 r 根本征方向）')

GRID = [R_EXP + 2e-3, 1.175, 1.18, 1.20, 1.30, 1.50, 1.75, 2.00, 2.50, 3.00, 4.00, 5.00]
p('\n== (1) [C1] 逐率读数：D(Theta-min) 对 D(Theta-max)，比值 <1 表示 min 更便宜（近地板那侧）==')
cross = {}
for rk in (1, 2, 3):
    Zmin = Vth[:, :rk]
    Zmax = Vth[:, n - rk:]
    p(f'\n rank {rk}:')
    p('    I 目标      s(min)      D(min)      s(max)      D(max)     Dmin/Dmax   谁便宜')
    seq = []
    for t in GRID:
        sm, dm = at_rate(Zmin, t)
        sx, dx = at_rate(Zmax, t)
        if dm is None or dx is None:
            p(f'    {t:8.4f}   不可达')
            continue
        seq.append((t, dm, dx))
        ratio = dm / dx
        p(f'    {t:8.4f}  {sm:10.3e}  {dm:11.4f}  {sx:10.3e}  {dx:11.4f}   {ratio:8.4f}   '
          f'{"min" if ratio < 1 else "max"}')
    ts = np.array([a for a, _, _ in seq])
    gd = np.array([b / c for _, b, c in seq]) - 1.0
    Ix = None
    for k in range(len(ts) - 1):
        if gd[k] * gd[k + 1] < 0:
            lo_, hi_ = ts[k], ts[k + 1]
            g_ = lambda t: (at_rate(Zmin, t)[1] or np.nan) / (at_rate(Zmax, t)[1] or np.nan) - 1.0
            for _ in range(45):
                mid = 0.5 * (lo_ + hi_)
                if g_(lo_) * g_(mid) <= 0:
                    hi_ = mid
                else:
                    lo_ = mid
            Ix = 0.5 * (lo_ + hi_)
            break
    cross[rk] = Ix
    if Ix is None:
        p(f'   符号变化：无（在网格 {ts[0]:.4f}…{ts[-1]:.4f} 内一侧到底）')
    else:
        p(f'   [C1] 翻转率 I*({rk}) = {Ix:.6f} bit   （离地板 R_exp 的距离 = {Ix - R_EXP:.6f} bit）'
          f'   相对地板位置 = {(Ix - R_EXP) / R_EXP * 100:.1f}%')

p('\n== (2) [C2] 同一批随机设计：Theta-min/Theta-max 在每个率上的同 rank 分位（0=最便宜） ==')
_, U = schur(A, output='real', sort=lambda a: abs(a) < 1.0)[:2]
rng = np.random.default_rng(20260930)
pool = {}
for rk in (1, 2, 3):
    zs = []
    for c in range(30):
        G = rng.standard_normal((n, rk))
        Z, _ = np.linalg.qr(G)
        zs.append((f'r{rk}-{c:02d}', Z))
    pool[rk] = zs

for rk in (1, 2, 3):
    p(f'\n rank {rk}:')
    p('    I 目标     D(THmin)   分位%    D(THmax)  分位%    随机最好(设计)     随机最好 D')
    for t in [R_EXP + 2e-3, 1.20, 1.50, 2.00, 2.50, 3.00, 4.00]:
        vals = []
        for tag, Z in pool[rk]:
            s, d = at_rate(Z, t)
            if d is not None and np.isfinite(d):
                vals.append((tag, d))
        if not vals:
            continue
        _, dmin = at_rate(Vth[:, :rk], t)
        _, dmax = at_rate(Vth[:, n - rk:], t)
        vals.sort(key=lambda z: z[1])
        ds = [d for _, d in vals]

        def pct(x):
            return 100.0 * sum(1 for y in ds if y <= x) / len(ds)
        best = vals[0]
        p(f'    {t:8.4f}  {dmin:10.4f}  {pct(dmin):6.1f}   {dmax:10.4f}  {pct(dmax):6.1f}   '
          f'{best[0]:<10} {best[1]:10.4f}')

p('\n== (3) [C3] 近地板前因子 b 与翻转率的一致性核对（b 小者应当是近地板端便宜的那个） ==')
for rk in (1, 2, 3):
    row = []
    for nm, Z in [('min', Vth[:, :rk]), ('max', Vth[:, n - rk:])]:
        s, d0 = at_rate(Z, R_EXP + 2e-3)
        Iend, Dend = curve(Z, 1e-6)
        row.append((nm, 1e-6 * (Dend - JC), Iend - R_EXP, Dend))
    p(f'  rank {rk}: ' + '   '.join(f'{nm}: b={bb:.4f}, I(1e-6)-R_exp={di:.2e}' for nm, bb, di, _ in row)
      + f'   翻转 I*={cross[rk]:.6f}' if cross[rk] else f'  rank {rk}: 无翻转')

open('p0/e110_out.txt', 'w', encoding='utf-8').write('\n'.join(out) + '\n')
p('done')
