# -*- coding: utf-8 -*-
r"""e172：并列集的**归一化角距** $-$ $-$ "支撑是否唯一"的正面判定（C 板 73-A 的命门）。

为什么 e165/e171 都没判成：
 * e165 的角度列用未归一化向量取 $\sigma=\lvert u\cdot v\rvert$，$\sigma>1$ 被 clip $\\Rightarrow$ 角朝 $0$ 塌、$\sigma<1$ 时又被读大 $-\ -$ 整列不可信；
   而且它只打"并列**条数**"，没打并列集**内部**的角散布 $\\Rightarrow$ 条数多 $\ne$ 支撑不唯一（可能全是同一平面的不同标架）。
 * e171 改成宽窗二分 $+$ 归一化，但随机方向扫在 $D{=}80,r{=}1$ 只有 $30/801$ 条方向可定价、$\min I$ 比 Powell 高 $10^{-2}$ bit
   $\\Rightarrow$ 方向采样太稀，$\epsilon=10^{-6}$ 带内永远只剩 1 条 $\\Rightarrow$ **功率不足，不判**（这轮照预注册的 [S3] 老实写了不判）。
本轮直接做**优化器并列集**：25 起点 Powell（与 e165 同机器、同容差），把每条读数的支撑**列归一化**后两两算最大主角度，
再看"值并列（$\le10^{-6}$ bit）的那几条，方向是不是同一个"。这才是 C 73-A 要的那件事。

判据（跑前写死；#46 判决自带读数条数）：
 [U0] 逐格打印原始有效 $k/25$ 与细化后读数总数；残差门 $10^{-7}$（e162 的 RES）。
 [U1] 并列集 $T:=\\{i:\ I_i\le \min I+10^{-6}\\}$；打印 $|T|$、$T$ 内**最大两两角距**（度，列归一化后的最大主角度），
      并打印 $T$ 内每条与最优那条的角距。
 [U2] 判决：$|T|\ge2$ 且 $T$ 内最大两两角距 $\ge5^\circ$ $\\Rightarrow$ 判"**近最优秩-$r$ 支撑不唯一（数值级）**"；
      $|T|\ge2$ 且全部两两 $<1^\circ$ $\\Rightarrow$ 判"**未见不唯一**（起点不同、收敛到同一支撑）"；
      $|T|=1$ $\\Rightarrow$ 判"**不判**"（单条并列，无信息）。
 [U3] 附带口径：$\Delta_r$ 仍只是**上界**（我的可行集 $\subseteq$ 凸集 $\\Rightarrow$ 下界 $\ge0$）；"不唯一"到 $10^{-6}$ bit 为止，
      不升格成"存在两个精确最优点" $-\ -$ 那要闭式或严格凸性。
格：$(D,r)\in\\{(34.41,2),(80.00,1),(56.66,1)\\}$（$34.41,2$ 是 C 路由②最关心的格；两秩-1 格用来对账 e162b/e165 的角）。
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
rng = np.random.default_rng(4242)
TOL_TIE = 1e-6


def unit(G):
    G = np.asarray(G, float).reshape(n, -1)
    return G / np.linalg.norm(G, axis=0)


def gap(Pa, Pb):
    sv = np.linalg.svd(unit(Pa).T @ unit(Pb), compute_uv=False)
    return float(np.degrees(np.arccos(np.clip(sv.min(), -1.0, 1.0))))


def vs_th(G, r):
    sv = np.linalg.svd(Uth[:, :r].T @ unit(G), compute_uv=False)
    return np.degrees(np.arccos(np.clip(sv, -1.0, 1.0)))


def cell(t, r, starts, coarse=700, fine=3000, ntop=6):
    raw = []
    for x0 in starts:
        P['GUESS'][0] = 0.0
        I, D, lg, G = P['polish'](x0, t, None, coarse)
        if np.isfinite(I):
            raw.append((I, abs(D - t), G))
    raw.sort(key=lambda z: z[0])
    ref = []
    for I, D, G in raw[:ntop]:
        P['GUESS'][0] = 0.0
        I2, D2, lg2, G2 = P['polish'](G, t, None, fine)
        if np.isfinite(I2):
            ref.append((I2, abs(D2 - t), G2))
    return raw, sorted(raw + ref, key=lambda z: z[0])


p(r'== 锚点 $J_C=%.6f$；$I_{\rm unc}$（e163 凸 SDP，与 C `c89b` 六位相符）：%s ==' % (JC, SDP))
p(r'   并列容差 $10^{-6}$ bit；角距一律**列归一化**后取最大主角度 $-\ -$ e165 的角度列作废的根因已修')
p('')
p(r'== [U0/U1/U2] 优化器并列集的角散布 ==')

for t, r in ((34.41, 2), (80.00, 1), (56.66, 1)):
    starts = [Uth[:, :r]] + [np.linalg.qr(rng.normal(0, 1, (n, r)))[0] for _ in range(24)]
    raw, allv = cell(t, r, starts)
    p(r'  --- $D=%.2f$，$r=%d$：原始有效 $k=%d/%d$，含细化共 %d 条读数 ---' % (t, r, len(raw), len(starts), len(allv)))
    if not allv:
        p(r'      [U2] 判决：0 条 $\\Rightarrow$ 不判')
        p('')
        continue
    Is = np.array([q[0] for q in allv])
    Imin = Is[0]
    p(r'      $\min I=%.6f$（$\Delta_r=%+.2e$ bit，[U3] 上界口径）' % (Imin, Imin - SDP[t]))
    tie = [q for q in allv if q[0] <= Imin + TOL_TIE]
    p(r'      [U1] 并列集 $|T|=%d$（$\le\min I+10^{-6}$ bit）；最优那条与 $\Theta$ 前 %d 平面主角度 %s'
      % (len(tie), r, np.array2string(vs_th(tie[0][2], r), precision=2)))
    if len(tie) >= 2:
        gmax = 0.0
        lines = []
        for i, a in enumerate(tie):
            for b in tie[i + 1:]:
                g = gap(a[2], b[2])
                gmax = max(gmax, g)
                lines.append((g, b[0] - a[0]))
        near = ', '.join('%.2f$^\\circ$/%+.1e' % (g, dv) for g, dv in sorted(lines, reverse=True)[:8])
        p(r'      [U1] $T$ 内两两角距（最大 8 对，格式 角距/值差）：%s' % near)
        p(r'      [U1] $T$ 内每条与最优的角距：%s'
          % np.array2string(np.array([gap(q[2], tie[0][2]) for q in tie]), precision=2))
        if gmax >= 5.0:
            p(r'      [U2] 判决：**近最优秩-$r$ 支撑不唯一（数值级）** $-\ -$ 存在 $|T|=%d$ 条读数互差 $\le10^{-6}$ bit '
              r'而支撑角距达 %.2f$^\circ$' % (len(tie), gmax))
        elif gmax < 1.0:
            p(r'      [U2] 判决：**未见不唯一** $-\ -$ $|T|=%d$ 条读数全部收敛到同一支撑（两两 $<1^\circ$，最大 %.2f$^\circ$）'
              % (len(tie), gmax))
        else:
            p(r'      [U2] 判决：**未决**（$|T|=%d$，最大角距 %.2f$^\circ$ 落在 $1^\circ$–$5^\circ$ 之间）' % (len(tie), gmax))
    else:
        p(r'      [U2] 判决：**不判**（$|T|=1$，单条并列无信息）')
    # 非并列的次优簇：给"面有多平"的第二层证据
    if len(Is) >= 3:
        rest = [(q[0] - Imin, gap(q[2], tie[0][2])) for q in allv[1:6]]
        p(r'      次优 5 条（值差 bit / 与最优支撑角距）：%s'
          % ', '.join('%.1e/%.1f$^\\circ$' % z for z in rest))
    p('')

p(r'== [U3] 口径重申：以上只到"$10^{-6}$ bit 内并列"，不升格为"存在两个精确最优点"；')
p(r'   后者的否证需要唯一性定理或闭式，本轮没有。C 的 73-A 归约仍然成立：严格情形 $\\Rightarrow$ 支撑唯一性问题。==')
p('')
p('用时 %.0f s | DARE %d | 有效读数 %d | 回退 %d | 异常 %s'
  % (time.time() - T0, P['NDB'][0], P['NVALID'][0], P['NFALLBACK'][0], P['EXC'] or '无'))
txt = '\n'.join(out) + '\n'
open('p0/e172_out.txt', 'w', encoding='utf-8').write(txt)
open('p0/e172_run.txt', 'w', encoding='utf-8').write(txt)
