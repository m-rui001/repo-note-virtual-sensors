# -*- coding: utf-8 -*-
r"""e167：秩-$r$ 最优支撑**是否唯一**的直接测法（等值带扫，不是多起点优化器）。

动机（两块互相矛盾的读数，都在我车道）：
 * `p0/e162b_out.txt` [B]：$D{=}80$ 的自由秩-1 最优射线离开 $\Theta$ 首轴 $6.02^\circ$；
 * `p0/e165_out.txt`：同格给 $62.5^\circ$，两处的 $\min I$ 都 $=$ SDP 值 $1.603105$。
 $\Rightarrow$ 要么"同一值的两个不同支撑"（$=\mathrm{Opt}(80)$ 的秩-1 面非平凡 $-$ $-$ 正是 C 板 73-A 的命门），
    要么其中一轮欠收敛/口径不同。多起点优化器分不清这两者 $-$ $-$ 它给"收敛到的点"，不给"值函数的等值带"。
本脚本改成**方向扫**：给定方向 $u$ 直接解增益使 $J=D_{\rm tgt}$（e162 的 `read_pt`），于是每个方向有确定的 $I(u)$；
等值带 $B_\epsilon:=\{u:\ I(u)\le\min I+\epsilon\}$ 的**角散布**就是"支撑不唯一"的可测形式。

判据（跑前写死；#46：每条判决行都带有效读数条数）：
 [Q0] 有效读数 $:=|D-D_{\rm tgt}|\le$ e162 的 RES 门；逐格打印扫向数/有效数/回退数。
 [Q1] 三个容差 $\epsilon\in\{10^{-6},10^{-4},10^{-2}\}$ bit，各报带内方向数与带内最大两两角距（度）。
 [Q2] 判决：取带内两两角距 $>5^\circ$ 的前 6 条方向，**各自独立细化**（Powell fine$=1500$）。
      细化后若仍存在两条方向角距 $\ge5^\circ$ 且值差 $\le10^{-6}$ bit $\\Rightarrow$ 判"**秩-$r$ 支撑在 $10^{-6}$ bit 内不唯一（数值级）**"；
      若细化读数两两 $<1^\circ$ $\\Rightarrow$ 判"**未见不唯一**"；其余 $\\Rightarrow$ "**未决**"，把具体角距与值差打出来。
 [Q3] 口径上限：本脚本不证"精确最优"；$\min I\ge I_{\rm unc}$（我的可行集 $\subseteq$ 凸集），所以 $B_\epsilon$ 是**近最优**集。
网格：rank-1 用 $\Theta$ 四轴 $+$ 600 随机单位向量；rank-2 用 $\Theta$ 前 2 列 $+$ 300 个随机正交 2-标架。
格：$(D,r)\in\{(80,1),(56.66,1),(34.41,2)\}$。
"""
import sys
import time
import numpy as np
sys.stdout.reconfigure(encoding='utf-8')

out = []
T0 = time.time()


def p(*a):
    s = ' '.join(str(x) for x in a)
    out.append(s)
    print(s)
    sys.stdout.flush()


P = {'__name__': 'p'}
exec(compile(open('p0/e162_ray_closure.py', encoding='utf-8').read().split('CRAB = {')[0],
             'p0/e162_ray_closure.py[prefix]', 'exec'), P)
AUD = P['_ns']

SDP = {34.41: 4.893329, 56.66: 2.053565, 80.00: 1.603105}
A, B, W = AUD['A'].copy(), AUD['B'].copy(), AUD['W'].copy()
Pc, K, TH, JC = P['ctrl_full'](A, B, W, np.eye(4), np.eye(4))
P['A'], P['B'], P['W'], P['TH'], P['JC'] = A, B, W, TH, JC
wth, U0 = np.linalg.eigh(TH)
Uth = U0[:, np.argsort(-wth)]
P['Uth'] = Uth
n = 4
rng = np.random.default_rng(7788)
EPS = [1e-6, 1e-4, 1e-2]


def gap(Pa, Pb):
    """两 $r$-平面的最大主角度（rank-1 时退化为射线夹角，用 $|\\cos|$ 所以对 $u\\to-u$ 不变）。"""
    sv = np.linalg.svd(Pa.T @ Pb, compute_uv=False)
    return float(np.degrees(np.arccos(np.clip(sv.min(), -1.0, 1.0))))


def pa_vs_th(G, r):
    return np.degrees(np.arccos(np.clip(np.linalg.svd(Uth[:, :r].T @ G, compute_uv=False), -1.0, 1.0)))


p('== 锚点 $J_C=%.6f$；$I_{\rm unc}$：%s ==' % (JC, SDP))
p('   $\Theta$ 特征值（降序）：%s' % np.array2string(wth[::-1], precision=4))
p('')
p('== [Q0/Q1] 方向扫：逐方向的 $I(u)$（增益由割线解出，使 $J=D_{\rm tgt}$）==')

for t, r in ((80.00, 1), (56.66, 1), (34.41, 2)):
    if r == 1:
        cands = [Uth[:, [i]] for i in range(n)] + \
                [np.linalg.qr(rng.normal(0, 1, (n, 1)))[0] for _ in range(600)]
    else:
        cands = [Uth[:, :r]] + [np.linalg.qr(rng.normal(0, 1, (n, r)))[0] for _ in range(300)]
    ndb0, nv0, fb0 = P['NDB'][0], P['NVALID'][0], P['NFALLBACK'][0]
    vals = []
    for G in cands:
        P['GUESS'][0] = 0.0
        I, D, lg = P['read_pt'](G, t, None)
        if np.isfinite(I):
            vals.append((I, G))
    if not vals:
        p('  $D=%.2f$ $r=%d$：有效 0 条 $\\Rightarrow$ 不判' % (t, r))
        continue
    vals.sort(key=lambda z: z[0])
    Imin = vals[0][0]
    Is = np.array([q[0] for q in vals])
    p('  --- $D=%.2f$，$r=%d$：扫 %d 条，有效 %d 条，$\min I=%.6f$（$-I_{\rm unc}=%+.2e$ bit）---'
      % (t, r, len(cands), len(vals), Imin, Imin - SDP[t]))
    for e in EPS:
        band = [q for q in vals if q[0] <= Imin + e]
        mx = max((gap(a[1], b[1]) for i, a in enumerate(band) for b in band[i + 1:]), default=0.0)
        p('      $\epsilon=%g$ bit：带内 %d/%d 条，带内最大两两角距 %.2f$^\\circ$' % (e, len(band), len(vals), mx))
    p('')
    p('== [Q2] 带内（$\epsilon=10^{-6}$）两两角距 $>5^\circ$ 的前 6 条方向 $\\to$ 各自独立细化 fine=1500 ==')
    band = [q for q in vals if q[0] <= Imin + EPS[0]]
    sel = []
    for q in band:
        if all(gap(q[1], s[1]) > 5.0 for s in sel):
            sel.append(q)
        if len(sel) >= 6:
            break
    pol = []
    for I, G in sel:
        P['GUESS'][0] = 0.0
        I2, D2, lg2, G2 = P['polish'](G, t, None, 1500)
        if np.isfinite(I2):
            pol.append((I2, G2, gap(G, G2)))
        else:
            p('      细化失败（残差门外）$-$ 原始 $I=%.6f$' % I)
    p('      细化成功 %d/%d 条' % (len(pol), len(sel)))
    for i, (I2, G2, mv) in enumerate(pol):
        p('      #%d $I=%.6f$（$-I_{\rm unc}=%+.2e$）；与 $\Theta$ 前 %d 平面主角度 %s；细化位移 %.2f$^\\circ$'
          % (i + 1, I2, I2 - SDP[t], r, np.array2string(pa_vs_th(G2, r), precision=1), mv))
    if len(pol) >= 2:
        dmax, pmax, qmax = 0.0, None, None
        for i, a in enumerate(pol):
            for b in pol[i + 1:]:
                gg = gap(a[1], b[1])
                if gg > dmax:
                    dmax, pmax, qmax = gg, a, b
        dv = abs(pmax[0] - qmax[0])
        p('      最远两两：角距 %.2f$^\\circ$，值差 %.2e bit' % (dmax, dv))
        if dmax >= 5.0 and dv <= 1e-6:
            p('      [Q2] 判决：**秩-$r$ 支撑在 $10^{-6}$ bit 内不唯一（数值级）**')
        elif dmax < 1.0:
            p('      [Q2] 判决：**未见不唯一**（细化读数两两 $<1^\\circ$）')
        else:
            p('      [Q2] 判决：**未决**（角距 %.2f$^\\circ$、值差 %.2e bit $-$ 未同时满足 $\ge5^\\circ$ 与 $\le10^{-6}$）'
              % (dmax, dv))
    else:
        p('      [Q2] 判决：细化成功不足 2 条 $\\Rightarrow$ 不判')
    ndb1, nv1, fb1 = P['NDB'][0], P['NVALID'][0], P['NFALLBACK'][0]
    p('      [Q0] 本格 DARE %d 次，有效读数 %d，回退 %d' % (nd1 - ndb0, nv1 - nv0, fb1 - fb0))
    p('')

p('== 对照：$\Theta$ 首轴/前 %d 平面自身的 $I$ 与 $\min I$ 的差（同一口径，一次求值）==' % 2)
for t, r in ((80.00, 1), (56.66, 1), (34.41, 2)):
    P['GUESS'][0] = 0.0
    I, D, lg = P['read_pt'](Uth[:, :r], t, None)
    if np.isfinite(I):
        p('  $D=%.2f$ $r=%d$：$\Theta$ 前 %d 列 $I=%.6f$，$\min I$ 未在此行重算（见上）' % (t, r, r, I))
    else:
        p('  $D=%.2f$ $r=%d$：$\Theta$ 前 %d 列在残差门外，无读数' % (t, r, r))
p('')
p('用时 %.0f s | DARE %d | 有效读数 %d | 回退 %d | 异常 %s'
  % (time.time() - T0, P['NDB'][0], P['NVALID'][0], P['NFALLBACK'][0], P['EXC'] or '无'))
txt = '\n'.join(out) + '\n'
open('p0/e167_out.txt', 'w', encoding='utf-8').write(txt)
open('p0/e167_run.txt', 'w', encoding='utf-8').write(txt)
