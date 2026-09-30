# -*- coding: utf-8 -*-
r"""e168：**分支感知**的等值带扫（先回答"我的单方向价格曲线是不是多值"，再判支撑唯一性）。

触发（e167 的部分输出，本车道自己的机器问题）：
 `p0/e167_levelset_support.py` 在 $D{=}80,r{=}1$ 从 $\Theta$ 首轴出发：单发割线给 $I=1.653054$，
 同一方向 Powell 细化后给 $I=1.603105$（$=$ SDP 值），而**方向位移 $0.00^\circ$**。
 $\Rightarrow$ 同一方向、同一 $D$ 有两个不同 $I$ $-$ $-$ 只可能是 $D(\lg)$ 非单调（多个增益夹住同一代价）。
 若是这样，则 e162b 的"$\Theta$ 首轴比自由端贵 $4.99e{-}02$ bit"、e165 的"离开角 $62.5^\circ$"、以及 C `c89b` 的同款数，
 全部依赖"落在哪个分支" $-$ $-$ **这条不查清，前面所有 $\Delta$/离开角都不许当结论**（#43）。

机器：对固定方向 $u$，$C=10^{0.5\,\lg}u^{\mathsf T}$，`e162.make_cur` 给 $(I(\lg),D(\lg))$。
逐方向在 $\lg\in[-2,12]$ 打网格（101 点），找 $D-t$ 的**全部**变号，每个变号用二分夹到 $|D-t|\le10^{-9}$。
于是"一个方向 $\to$ 若干条读数"，分支显式，不依赖 GUESS 携带（那是 e167 的病根）。

判据（跑前写死；#46 每条判决自带读数条数）：
 [B0] 打印每格：方向数、有 $\ge1$ 个交点的方向数、交点总数、**有 $\ge2$ 个交点的方向数**、非有限点数。
 [B1] 分支判决：若存在方向有 $\ge2$ 个交点且这些交点的 $I$ 差 $>10^{-3}$ bit $\Rightarrow$ 判"**$D(\lg)$ 非单调：同一方向同一代价有多个增益**
      $-\ -$ 单发割线/GUESS 携带的读数必须标分支"；否则判"单调（每方向至多一交点，$\lg$ 网格分辨率内）"。
 [B2] 单调性侧检：打印 $\lg$ 网格上 $D$ 的**上升段计数**（$D$ 应随增益降；上升段 $>0$ 即非单调的直接证据）。
 [B3] 支撑唯一性（分支清洗后）：全体交点取 $\min I$，带内 $\epsilon\in\{10^{-6},10^{-4},10^{-2}\}$ bit 的方向数与最大两两角距。
      判决：带内存在两条**方向角距 $\ge5^\circ$** 且值差 $\le10^{-6}$ bit $\Rightarrow$ "秩-$r$ 近最优支撑不唯一（数值级）"；
      全部两两 $<1^\circ$ $\Rightarrow$ "未见不唯一"；其余 $\Rightarrow$ "未决"（打具体数）。
 [B4] 口径：$I$ 随增益升、$\min I\ge I_{\rm unc}$（我的可行集 $\subseteq$ 凸集）；这里给的是**近最优**集，不证精确最优。
格：$(D,r)\in\{(80,1),(56.66,1),(34.41,2)\}$；方向：$\Theta$ 前 $r$ 列 $+$ 240 随机（$r\ge2$ 用 120 个随机正交标架）。
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
rng = np.random.default_rng(90210)
GRID = np.linspace(-2.0, 12.0, 101)
EPS = [1e-6, 1e-4, 1e-2]


def gap(Pa, Pb):
    sv = np.linalg.svd(Pa.T @ Pb, compute_uv=False)
    return float(np.degrees(np.arccos(np.clip(sv.min(), -1.0, 1.0))))


def bisect_hit(cur, t, lo, hi):
    """在 $[\ell o,\,hi]$ 内二分 $D(\lg)-t$ 到 $10^{-9}$；返回 $(\lg, I, D)$，无解返回 None。"""
    flo = cur(lo)[1] - t
    fhi = cur(hi)[1] - t
    if not (np.isfinite(flo) and np.isfinite(fhi)) or flo * fhi > 0:
        return None
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        fm = cur(mid)[1] - t
        if not np.isfinite(fm):
            return None
        if flo * fm <= 0:
            hi, fhi = mid, fm
        else:
            lo, flo = mid, fm
        if hi - lo < 1e-12:
            break
    lg = 0.5 * (lo + hi)
    I, D = cur(lg)
    return (lg, I, D) if abs(D - t) <= 1e-9 else None


p(r'== 锚点 $J_C=%.6f$；$I_{\rm unc}$（e163 凸 SDP）：%s ==' % (JC, SDP))
p(r'   $\Theta$ 特征值（降序）：%s ；$\lg$ 网格 $[%g,%g]$ 共 %d 点' % (np.array2string(wth[::-1], precision=4),
                                                                       GRID[0], GRID[-1], len(GRID)))
p('')
p(r'== [B0/B1/B2] 逐方向的全交点扫描 ==')

ALL = {}
for t, r in ((80.00, 1), (56.66, 1), (34.41, 2)):
    nd = 240 if r == 1 else 120
    cands = [Uth[:, :r]] + [np.linalg.qr(rng.normal(0, 1, (n, r)))[0] for _ in range(nd)]
    pts = []
    nbad = 0
    nup = 0
    ndir_hit = 0
    ndir_multi = 0
    maxdI = 0.0
    for G in cands:
        Zn = G / np.linalg.norm(G, axis=0)
        cur = P['make_cur'](Zn, r)
        Ds = []
        for lg in GRID:
            try:
                I, D = cur(lg)
            except Exception:
                I, D = np.nan, np.nan
            Ds.append(D)
        Ds = np.array(Ds)
        nbad += int(np.sum(~np.isfinite(Ds)))
        nup += int(np.sum(np.diff(Ds[np.isfinite(Ds)]) > 0))
        crossings = []
        for i in range(len(GRID) - 1):
            a, b = Ds[i], Ds[i + 1]
            if np.isfinite(a) and np.isfinite(b) and (a - t) * (b - t) <= 0 and a != b:
                got = bisect_hit(cur, t, GRID[i], GRID[i + 1])
                if got is not None:
                    crossings.append(got)
        if crossings:
            ndir_hit += 1
            Is = np.array([q[1] for q in crossings])
            if len(crossings) >= 2:
                ndir_multi += 1
                maxdI = max(maxdI, float(Is.max() - Is.min()))
        for lg, I, D in crossings:
            pts.append((I, Zn, lg))
    p(r'  --- $D=%.2f$，$r=%d$：方向 %d 条 $\Rightarrow$ 有交点 %d 条，交点共 %d 个，'
      r'**有 $\ge2$ 个交点的方向 %d 条**（这些方向上的 $I$ 跨度最大 %.4f bit）---'
      % (t, r, len(cands), ndir_hit, len(pts), ndir_multi, maxdI))
    p(r'      [B2] $\lg$ 网格上 $D$ 的上升段计数 %d（$D$ 应随增益单调降；$>0$ 即非单调的直接证据）；非有限点 %d'
      % (nup, nbad))
    if ndir_multi >= 1 and maxdI > 1e-3:
        p(r'      [B1] 判决：**$D(\lg)$ 非单调 $-\ -$ 同一方向、同一代价有多个增益（$I$ 跨度最大 %.4f bit），'
          r'单发割线/携带 GUESS 的读数必须标分支**' % maxdI)
    elif ndir_multi == 0:
        p(r'      [B1] 判决：**单调（本网格分辨率内每方向至多一交点）**')
    else:
        p(r'      [B1] 判决：有多交点方向但 $I$ 跨度 $\le10^{-3}$ bit（最大 %.4f）$\Rightarrow$ 分支效应可忽略' % maxdI)
    ALL[(t, r)] = pts
    p('')

p(r'== [B3] 分支清洗后的秩-$r$ 等值带（全体交点，含同一方向的多支）==')
for (t, r), pts in ALL.items():
    if not pts:
        p(r'  $D=%.2f$ $r=%d$：0 条交点 $\\Rightarrow$ 不判' % (t, r))
        continue
    pts.sort(key=lambda z: z[0])
    Imin = pts[0][0]
    p(r'  --- $D=%.2f$，$r=%d$：交点 %d 个，$\min I=%.6f$（$-I_{\rm unc}=%+.2e$ bit）---'
      % (t, r, len(pts), Imin, Imin - SDP[t]))
    for e in EPS:
        band = [q for q in pts if q[0] <= Imin + e]
        mx = max((gap(a[1], b[1]) for i, a in enumerate(band) for b in band[i + 1:]), default=0.0)
        nbr = len(set(round(q[2], 6) for q in band))
        p(r'      $\epsilon=%g$ bit：带内 %d/%d 个交点（跨 %d 个增益值），带内最大两两方向角距 %.2f$^\circ$'
          % (e, len(band), len(pts), nbr, mx))
    band = [q for q in pts if q[0] <= Imin + EPS[0]]
    sel = []
    for q in band:
        if all(gap(q[1], s[1]) > 5.0 for s in sel):
            sel.append(q)
        if len(sel) >= 8:
            break
    if len(sel) >= 2:
        dmax = max(gap(a[1], b[1]) for i, a in enumerate(sel) for b in sel[i + 1:])
        dv = max(a[0] for a in sel) - min(a[0] for a in sel)
        p(r'      [B3] 带内 $\ge5^\circ$ 抽样 %d 条：最大角距 %.2f$^\circ$，值跨度 %.2e bit' % (len(sel), dmax, dv))
        if dmax >= 5.0 and dv <= 1e-6:
            p(r'      [B3] 判决：**秩-$r$ 近最优支撑在 $10^{-6}$ bit 内不唯一（数值级）**')
        elif dmax < 1.0:
            p(r'      [B3] 判决：**未见不唯一**')
        else:
            p(r'      [B3] 判决：**未决**（角距 %.2f$^\circ$ 但值跨度 %.2e bit $>$ $10^{-6}$）' % (dmax, dv))
    else:
        p(r'      [B3] 判决：带内 $\ge5^\circ$ 抽样不足 2 条 $\Rightarrow$ 未见不唯一（但功率有限，报数不判）')
    p('')

p(r'== [B4] 口径：以上是**近最优**集（$\min I\ge I_{\rm unc}$，下界由包含关系白送）；不证精确最优 ==')
p('')
p('用时 %.0f s | DARE %d | 有效读数 %d | 回退 %d | 异常 %s'
  % (time.time() - T0, P['NDB'][0], P['NVALID'][0], P['NFALLBACK'][0], P['EXC'] or '无'))
txt = '\n'.join(out) + '\n'
open('p0/e168_out.txt', 'w', encoding='utf-8').write(txt)
open('p0/e168_run.txt', 'w', encoding='utf-8').write(txt)
