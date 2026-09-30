# -*- coding: utf-8 -*-
"""R62-D 实验 123：打分器稳健性——K=a*b 的近地板读数换一个 s 还成不成立。
s 取 1e-6 / 1e-5 / 1e-4（都在 ν=1.0001 律的渐近段内、且不踩 s_break≈1e-7）。
若三个 s 的 rho 都高 ⇒ "K 是排序量"不是挑出来的读数；若只有 1e-6 高 ⇒ 降级为单点巧合。
"""
import numpy as np
from scipy.linalg import solve_discrete_are, eigh
from scipy.stats import spearmanr

_src = open('p0/exp_c_audit.py', encoding='utf-8').read().split("print(r'== E53")[0]
_ns = {'__name__': 'p'}
exec(compile(_src, 'p0/exp_c_audit.py[preamble]', 'exec'), _ns)
A, W, TH, JC, n = _ns['A'], _ns['W'], _ns['TH'], _ns['JC'], _ns['n']
sym = _ns['sym']
R_EXP = 1.168539
SS = (1e-6, 1e-5, 1e-4)
out = []


def p(s):
    out.append(str(s)); print(s, flush=True)


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
    if not (curve(Z, 10 ** lo)[0] < target < curve(Z, 10 ** hi)[0]):
        return None
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if curve(Z, 10 ** mid)[0] < target:
            lo = mid
        else:
            hi = mid
    return curve(Z, 10 ** (0.5 * (lo + hi)))[1]


_, Vth = eigh(TH)
rng = np.random.default_rng(20260930)
pool = []
for rk in (1, 2, 3):
    for c in range(30):
        Z, _ = np.linalg.qr(rng.standard_normal((n, rk)))
        pool.append((f'r{rk}-{c:02d}', rk, Z))
for rk in (1, 2, 3):
    pool.append((f'THmin-{rk}', rk, Vth[:, :rk]))
    pool.append((f'THmax-{rk}', rk, Vth[:, n - rk:]))

rec = {}
for nm, rk, Z in pool:
    legs = {}
    ok = True
    for s in SS:
        I0, D0, rho0 = curve(Z, s)
        if rho0 >= 1.0 or I0 <= R_EXP:
            ok = False; break
        legs[s] = (s * (D0 - JC), (I0 - R_EXP) / s)
    if not ok:
        continue
    d3 = at_rate(Z, 3.0)
    if d3 is None:
        continue
    rec[nm] = (rk, d3, legs)

p('== rho(打分器, D@I=3) 随近地板读数 s 的变化（同 rank 组内，n=%d/组）==' % 32)
p('  rank   s=1e-6        s=1e-5        s=1e-4       （每格：rho(b) / rho(a) / rho(K)）')
for rk in (1, 2, 3):
    sel = [v for k, v in rec.items() if v[0] == rk]
    dd = np.array([v[1] for v in sel])
    cells = []
    for s in SS:
        bb = np.array([v[2][s][0] for v in sel]); aa = np.array([v[2][s][1] for v in sel])
        kk = bb * aa
        cells.append('%+.3f/%+.3f/%+.3f' % (spearmanr(bb, dd).correlation,
                                            spearmanr(aa, dd).correlation,
                                            spearmanr(kk, dd).correlation))
    p('  %d     %s   %s   %s' % (rk, cells[0], cells[1], cells[2]))

p('\n== min-K 命中 min-D 的情况（按三个 s 各判一次）==')
for rk in (1, 2, 3):
    sel = [(k, v) for k, v in rec.items() if v[0] == rk]
    dd = np.array([v[1] for _, v in sel]); tags = [k for k, _ in sel]
    best_d = tags[int(np.argmin(dd))]
    line = []
    for s in SS:
        kk = np.array([v[2][s][0] * v[2][s][1] for _, v in sel])
        line.append('%s:%s' % ('{:.0e}'.format(s), '是 ✓' if tags[int(np.argmin(kk))] == best_d else '否(' + tags[int(np.argmin(kk))] + ')'))
    p('  rank %d  min-D=%-9s  %s' % (rk, best_d, '  '.join(line)))

p('\n== a、b、K 的跨设计倍差（s=1e-6）==')
sel = list(rec.values())
bb = np.array([v[2][1e-6][0] for v in sel]); aa = np.array([v[2][1e-6][1] for v in sel])
p('  b %.1fx   a %.1fx   K %.1fx' % (bb.max() / bb.min(), aa.max() / aa.min(), (bb * aa).max() / (bb * aa).min()))
open('p0/e123_out.txt', 'w', encoding='utf-8').write('\n'.join(out) + '\n')
