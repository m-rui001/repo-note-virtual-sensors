# -*- coding: utf-8 -*-
"""R62-D 实验 122：把 §69 的 K=a*b 做成【可交付的打分器】——同一批 96 个设计，
只用近地板单点 (s=1e-6) 的读数构造 b、a、K，再与 I=3.0 处精确二分解出的 D 比排序一致性。
对照：§68-C 报过 rho(b,D@3)=+0.860/+0.682/+0.158。本节判 rho(K,D@3) 是否把它救回来。
这是 §69-G-② "a,b,K 三列对表"的我侧交付件（C 报 D，我报两腿 + K）。
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
    return s, D


_, Vth = eigh(TH)
designs = [('tail-2', 2, None)]
rng = np.random.default_rng(20260930)
for rk in (1, 2, 3):
    for c in range(30):
        G = rng.standard_normal((n, rk))
        Z, _ = np.linalg.qr(G)
        designs.append((f'r{rk}-{c:02d}', rk, Z))
for rk in (1, 2, 3):
    designs.append((f'THmin-{rk}', rk, Vth[:, :rk]))
    designs.append((f'THmax-{rk}', rk, Vth[:, n - rk:]))

rows = []
for nm, rk, Z in designs:
    if Z is None:
        continue
    I0, D0, rho0 = curve(Z, 1e-6)
    if rho0 >= 1.0 or I0 <= R_EXP:
        rows.append((nm, rk, np.nan, np.nan, np.nan, None)); continue
    a = (I0 - R_EXP) / 1e-6
    b = 1e-6 * (D0 - JC)
    res = at_rate(Z, 3.0)
    rows.append((nm, rk, b, a, a * b, res[1] if res else None))

p('== 96 个设计：近地板单点 (s=1e-6) 的两腿 a、b 与 K=a*b，对 I=3.0 处精确 D 的排序一致性 ==')
p('   rho 是同 rank 组内的 Spearman(打分器, D@I=3)；§68-C 只报过 b 那一列')
p('')
p('  rank   n   rho(b,D@3)   rho(a,D@3)   rho(K,D@3)   min-K 设计 == min-D 设计 ?')
for rk in (1, 2, 3):
    sel = [r for r in rows if r[1] == rk and np.isfinite(r[2]) and r[5] is not None]
    bb = np.array([r[2] for r in sel]); aa = np.array([r[3] for r in sel])
    kk = np.array([r[4] for r in sel]); dd = np.array([r[5] for r in sel])
    rb = spearmanr(bb, dd).correlation; ra = spearmanr(aa, dd).correlation
    rk_ = spearmanr(kk, dd).correlation
    same = sel[int(np.argmin(kk))][0] == sel[int(np.argmin(dd))][0]
    p('  %d     %2d    %+.3f        %+.3f        %+.3f          %s  (%s / %s)'
      % (rk, len(sel), rb, ra, rk_, '是 ✓' if same else '否 ✗', sel[int(np.argmin(kk))][0], sel[int(np.argmin(dd))][0]))

p('\n== Theta 极值设计在两种打分器下的名次（同 rank 内，1=最便宜；K 只用 s=1e-6 一个读数）==')
for rk in (1, 2, 3):
    sel = [r for r in rows if r[1] == rk and np.isfinite(r[2]) and r[5] is not None]
    tags = [r[0] for r in sel]
    pos_d = {tags[i]: 1 + p_ for p_, i in enumerate(sorted(range(len(sel)), key=lambda i: sel[i][5]))}
    pos_k = {tags[i]: 1 + p_ for p_, i in enumerate(sorted(range(len(sel)), key=lambda i: sel[i][4]))}
    pos_b = {tags[i]: 1 + p_ for p_, i in enumerate(sorted(range(len(sel)), key=lambda i: sel[i][2]))}
    for tag in (f'THmax-{rk}', f'THmin-{rk}'):
        p('  rank %d  %-9s : D@3 第 %2d 名 | K 第 %2d 名 | b 第 %2d 名   (K 名次差 %+d, b 名次差 %+d)'
          % (rk, tag, pos_d[tag], pos_k[tag], pos_b[tag], pos_k[tag] - pos_d[tag], pos_b[tag] - pos_d[tag]))

p('\n== K 的跨度与 b 的跨度（同一个 s=1e-6 读数）==')
kk = np.array([r[4] for r in rows if np.isfinite(r[4])]); bb = np.array([r[2] for r in rows if np.isfinite(r[2])])
aa = np.array([r[3] for r in rows if np.isfinite(r[3])])
p('  b : %.4g .. %.4g  (%.1f x)' % (bb.min(), bb.max(), bb.max() / bb.min()))
p('  a : %.4g .. %.4g  (%.1f x)' % (aa.min(), aa.max(), aa.max() / aa.min()))
p('  K : %.4g .. %.4g  (%.1f x)' % (kk.min(), kk.max(), kk.max() / kk.min()))
p('\n判据：若 rho(K,D@3) 在 rank-3 从 +0.158 抬到接近 +1 ⇒ 前沿律的排序量确实是 K，b 单独用是错的腿；')
p('        若 K 与 b 的 rho 接近 ⇒ K 的说法只是重述，rank-3 的低一致另有原因。')
open('p0/e122_out.txt', 'w', encoding='utf-8').write('\n'.join(out) + '\n')
