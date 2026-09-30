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
    rho = np.max(np.abs(np.linalg.eigvals((np.eye(n) - Lk @ C) @ A)))
    return I, D, rho


def at_rate(Z, target=3.0):
    lo, hi = -7.0, 5.0
    flo, fhi = curve(Z, 10 ** lo)[0], curve(Z, 10 ** hi)[0]
    if not (flo < target < fhi):
        return None
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if curve(Z, 10 ** mid)[0] < target:
            lo = mid
        else:
            hi = mid
    s = 10 ** (0.5 * (lo + hi))
    I, D, rho = curve(Z, s)
    return s, I, D, rho


def prefactor(Z):
    s = 1e-6
    I, D, rho = curve(Z, s)
    if rho >= 1.0 or I <= 1.168539:
        return np.nan
    return s * (D - JC)


_, U = schur(A, output='real', sort=lambda a: abs(a) < 1.0)[:2]
designs = []
for k, nm in [(4, 'tail-4'), (3, 'tail-3'), (2, 'tail-2')]:
    designs.append((nm, k, U[:, n - k:]))

rng = np.random.default_rng(20260930)
for rk in (1, 2, 3):
    for c in range(30):
        G = rng.standard_normal((n, rk))
        Z, _ = np.linalg.qr(G)
        designs.append((f'r{rk}-{c:02d}', rk, Z))

wth, Vth = eigh(TH)
for rk in (1, 2, 3):
    designs.append((f'THmin-{rk}', rk, Vth[:, :rk]))
    designs.append((f'THmax-{rk}', rk, Vth[:, n - rk:]))

R_EXP = 1.168539
rows = []
for tag, rk, Z in designs:
    try:
        got = at_rate(Z, 3.0)
        b = prefactor(Z)
    except Exception as e:
        got, b = None, np.nan
        p(f'  {tag} 失败 {type(e).__name__}: {e}')
        continue
    if got is None:
        rows.append((tag, rk, np.nan, np.nan, np.nan, np.nan, b))
        continue
    s, I, D, rho = got
    rows.append((tag, rk, s, I, D, rho, b))

p('== (0) 锚点 ==')
p(f'  JC={JC:.6f}  R_exp={R_EXP}  lambda(TH) 升序={np.round(wth, 4)}')
p('  口径：同一批 93 个设计（seed 20260930 + 坐标尾平面），把 s 解到 I=3.000 bit，读 D；')
p('  再与近地板前因子 b(V)=s(D-JC)|_{s=1e-6} 对照。')

p('\n== (1) 逐设计：D 在 I=3 的读数 vs 近地板 b ==')
p('  设计        rank   s(I=3)     D(I=3)     rho      b(1e-6)    判定')
for tag, rk, s, I, D, rho, b in rows:
    if np.isnan(D):
        p(f'  {tag:<11} {rk}    不可达 I=3（该设计的 I 上限不足 3 bit）')
        continue
    p(f'  {tag:<11} {rk}    {s:9.3e}  {D:10.4f}  {rho:6.4f}  {b:10.4f}  ')

p('\n== (2) [Q1] 同 rank 内按 D(I=3) 排名：Theta-规则落在哪 ==')
for rk in (1, 2, 3, 4):
    grp = [r for r in rows if r[1] == rk and not np.isnan(r[4])]
    if not grp:
        continue
    ds = sorted(grp, key=lambda r: r[4])
    p(f'\n rank {rk}: 可算 {len(grp)} 个，D 区间 {ds[0][4]:.4f} … {ds[-1][4]:.4f}（倍差 {ds[-1][4]/ds[0][4]:.3f}x）')
    for pos, r in enumerate(ds[:4]):
        p(f'   最小 {pos+1}: {r[0]:<10} D={r[4]:.4f}  b={r[6]:.4f}')
    for r in grp:
        if r[0].startswith('TH'):
            q = 100.0 * sum(1 for x in ds if x[4] <= r[4]) / len(ds)
            p(f'   {r[0]:<9} D={r[4]:.4f} 落在同 rank 的 {q:.1f}% 分位（100%=最贵）')

p('\n== (3) [Q2] 近地板前因子 b 能否预测 I=3 的 D（同 rank 的 Spearman） ==')


def spearman(x, y):
    m = [(a, b) for a, b in zip(x, y) if np.isfinite(a) and np.isfinite(b)]
    if len(m) < 6:
        return np.nan, len(m)
    ax = np.argsort(np.argsort([a for a, _ in m]))
    ay = np.argsort(np.argsort([b for _, b in m]))
    return float(np.corrcoef(ax, ay)[0, 1]), len(m)


for rk in (1, 2, 3):
    grp = [r for r in rows if r[1] == rk and not np.isnan(r[4])]
    rho_s, m = spearman([r[6] for r in grp], [r[4] for r in grp])
    p(f'  rank {rk}: rho(b, D@I=3) = {rho_s:+.3f}（n={m}）')
grp = [r for r in rows if not np.isnan(r[4])]
rho_s, m = spearman([r[6] for r in grp], [r[4] for r in grp])
p(f'  全样本（混 rank）: rho = {rho_s:+.3f}（n={m}）—— 混 rank 只作参考，rank 本身强烈改变 D')

p('\n== (4) [Q3] 与已发表口径/本车道已知最好值的对照 ==')
p('  42.0402 = 无约束精确 converse（[4] 的 SDP，I=3）；45.4537 = tail-3 平面 Powell（本车道已知最好平面设计）；')
p('  50.6023 = tail-2 平面 Powell。')
for rk in (2, 3, 4):
    grp = [r for r in rows if r[1] == rk and not np.isnan(r[4])]
    if not grp:
        continue
    best = min(grp, key=lambda r: r[4])
    p(f'  rank {rk}: 随机/坐标里最好 {best[0]} D={best[4]:.4f}；比 45.4537 {"便宜" if best[4] < 45.4537 else "贵"} {abs(best[4]-45.4537):.4f}')
open('p0/e109_out.txt', 'w', encoding='utf-8').write('\n'.join(out) + '\n')
print('rows', len(rows))
