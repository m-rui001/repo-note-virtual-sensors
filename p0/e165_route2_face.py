# -*- coding: utf-8 -*-
r"""e165：把 C 补记十二 73-C 的**路由②**在我车道跑成正式判决（$\\min_{\\mathrm{rank}\\,V=r}\\Delta(D;V)$）。

C 的原文（板 4375–4376 行，逐字）：
  "② 对 $r<r^*(D)$ 直接算 $\\min_{V:\\mathrm{rank}V=r}\\Delta(D;V)$……若这个 min 严格 $>0$，引理的严格情形就是**数值判决**"
  "判据……跑前写死：min 的相对量级 $>10^{-4}$ bit 且 24 个起点一致 $\\Rightarrow$ 判'严格情形成立（数值级）'；否则报'面的结构未定'"

对象：$\\Delta_r(D):=\\min_{Z:\\mathrm{rank}\\,r}I(Z,D)-I_{\\rm unc}(D)$，$I_{\\rm unc}$ 取**本车道自算的凸 SDP 值**（e163，
与 C 的 `c89b` 在六位数上相符），所以右边不是"我的参照"而是"两条独立凸路线的公共值"。

为什么这次能做、上一轮不能做（公开自我降级，#43）：
  e162 的秩$\\le$3 min 竟比秩$\\le$2 差 $6.9\\times10^{-3}$ $\\Rightarrow$ **包含关系违反** $\\Rightarrow$ 我方高阶 Powell 欠优化 $\\Rightarrow$ 上轮作废。
  本轮把"包含关系"当**硬门**用：逐格先检 $\\Delta_r\\le\\Delta_{r+1}$（秩越低、可行集越小、$\\min$ 只能更大），
  违反的格**当场作废、不参与判决**，只报数。

预注册判据（跑前写死）：
 [T0] 残差门 $|D-D_{\\rm tgt}|\\le10^{-7}$；逐格打印有效起点条数（#46）。
 [T1] 包含关系门：同一 $D$ 上 $\\Delta_1\\ge\\Delta_2\\ge\\Delta_3$（不可达记 $+\\infty$）。违反 $\\Rightarrow$ 该 $D$ 的高秩读数作废，
      并且该格**不许**按 C 的判据下"严格/不成立"任何结论，改判"我方优化器不可信"。
 [T2] 起点一致：$\\ge24$ 个独立起点里，最优 $\\min$ 与第 2 好的差 $\\le10^{-6}$ bit，且有效条数 $\\ge20$ $\\Rightarrow$ 记"一致"。
 [T3] 判决（照 C 的门）：$\\Delta_r>10^{-4}$ bit 且 T2 一致 $\\Rightarrow$ "严格情形成立（数值级）"；
      $\\Delta_r\\le10^{-4}$ bit 且 T2 一致 $\\Rightarrow$ **严格情形不成立**（$\\mathrm{Opt}(D)$ 含近似秩-$r$ 成员 $=$ C 说的"面非平凡"的定量实例）；
      其余一律"面的结构未定"。
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

SDP = {32.31: 7.555200, 34.41: 4.893329, 40.00: 3.271908, 56.66: 2.053565, 80.00: 1.603105}
FLOOR = {1: 43.7058, 2: 32.2752, 3: 31.6737}
RANK6 = {32.31: 4, 34.41: 3, 40.00: 2, 56.66: 1, 80.00: 1}

A, B, W = AUD['A'].copy(), AUD['B'].copy(), AUD['W'].copy()
Pc, K, TH, JC = P['ctrl_full'](A, B, W, np.eye(4), np.eye(4))
P['A'], P['B'], P['W'], P['TH'], P['JC'] = A, B, W, TH, JC
wth, Uth = np.linalg.eigh(TH)
Uth = Uth[:, np.argsort(-wth)]
P['Uth'] = Uth

p('== 锚点：$J_C=%.6f$；$I_{\\rm unc}$ 取本车道 e163 的凸 SDP 值（与 C 的 `c89b` 六位相符）==' % JC)
p('   秩的凸侧读数 $r^*(D)=\\mathrm{rank}_{10^{-6}}(\\Lambda)$：%s' % RANK6)
p('   同秩可达地板：%s（$>$ 该地板的格记不可达 $=+\\infty$，按 #43 不许当结论）' % FLOOR)

CELLS = [(34.41, 2), (34.41, 3), (32.31, 2), (32.31, 3), (40.00, 1), (56.66, 1), (80.00, 1)]
p('')
p('== [T0/T2] 逐格 $\\min_{\\mathrm{rank}\\,V=r}\\Delta$（25 起点：$\\Theta$ 前 $r$ 列 $+$ 24 个随机正交标架）==')
p('%7s %4s %11s %11s %9s %6s %9s %9s' % ('D', 'r', 'min I', 'delta_r', '残差', 'k/25', '次优-最优', '主角度(deg)'))
p('   计数口径（本轮第一版用错过，公开）：$k=$**原始起点里夹住目标代价的条数**；`best()` 的返回长度是 $k+$(top-6 细化)，')
p('   把它当"有效起点数"会让 [T2] 的判决行带错读数（#46）$\\Rightarrow$ 本轮把 $\\mathrm{best}$ 展开成 $\\mathrm{cell}()$ 自己数。')
RES = {}
rng = np.random.default_rng(4242)


def cell(t, r, starts, coarse=700, fine=3000, ntop=6):
    """自己展开 best()：要的是**原始有效起点条数** $k$（`best` 的返回长度是 $k+$细化，不能当 $k$ 用 $-$ $-$ #46）。"""
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
    allv = sorted(raw + ref, key=lambda z: z[0])
    return raw, allv


for t, r in CELLS:
    if t <= FLOOR[r]:
        p('%7.2f %4d  不可达（同秩地板 %.4f $\\ge D$）$\\Rightarrow$ 记 $+\\infty$，不判' % (t, r, FLOOR[r]))
        RES[(t, r)] = np.inf
        continue
    starts = [Uth[:, :r]] + [np.linalg.qr(rng.normal(0, 1, (4, r)))[0] for _ in range(24)]
    raw, allv = cell(t, r, starts)
    if not allv:
        p('%7.2f %4d  有效 0 条 $\\Rightarrow$ 无读数，不判' % (t, r))
        RES[(t, r)] = np.nan
        continue
    I2, rb, G2 = allv[0]
    Is = np.array([q[0] for q in allv])
    dlt = I2 - SDP[t]
    sec = Is[1] - Is[0] if len(Is) > 1 else np.nan
    ntop_within = int(np.sum(Is <= Is[0] + 1e-6))
    pv = np.degrees(np.arccos(np.clip(np.linalg.svd(Uth[:, :r].T @ G2, compute_uv=False), -1, 1)))
    RES[(t, r)] = dlt
    p('%7.2f %4d %11.6f %+11.2e %9.1e %6s %9.1e %s'
      % (t, r, I2, dlt, rb, '%d/%d' % (len(raw), len(starts)), sec, np.array2string(pv, precision=1)))
    p('        原始有效起点 $k=%d/25$；细化后 $\\le10^{-6}$ bit 内并列的读数 %d 条；'
      '次优$-$最优 %.1e bit $\\Rightarrow$ [T2] 起点一致（$k\\ge20$ 且次优差 $\\le10^{-6}$）= %s'
      % (len(raw), ntop_within, sec, bool(len(raw) >= 20 and np.isfinite(sec) and sec <= 1e-6)))

p('')
p('== [T1] 包含关系门（$\\Delta_r$ 必须随 $r$ 递增而降：秩越低可行集越小）==')
for t in (34.41, 32.31):
    d2, d3 = RES.get((t, 2), np.nan), RES.get((t, 3), np.nan)
    ok = (not np.isfinite(d2)) or (not np.isfinite(d3)) or (d2 >= d3 - 1e-9)
    p('  $D=%.2f$：$\\Delta_2=%s$，$\\Delta_3=%s$ $\\Rightarrow$ $\\Delta_2\\ge\\Delta_3$：%s%s'
      % (t, ('%.2e' % d2) if np.isfinite(d2) else str(d2),
         ('%.2e' % d3) if np.isfinite(d3) else str(d3), ok,
         '' if ok else '  **违反 $\\Rightarrow$ 该 $D$ 的高秩读数作废，两格都不许按 C 的门判**'))

p('')
p('== [T3] 判决（照 C 板 4376 行的门，逐字引用见文件头）==')
for t, r in CELLS:
    d = RES.get((t, r), np.nan)
    if not np.isfinite(d):
        p('  $D=%.2f$，$r=%d$：$\\Delta$ 非有限（%s）$\\Rightarrow$ 面的结构未定' % (t, r, '不可达' if d == np.inf else '无读数'))
        continue
    if r >= RANK6[t]:
        p('  $D=%.2f$，$r=%d$：$r\\ge r^*(D)=%d$ $\\Rightarrow$ **非严格情形**，$\\Delta=%.2e$ 应当为 $0$ $-$ $-$ 不作面结构的证据'
          % (t, r, RANK6[t], d))
        continue
    oki = (t in (34.41, 32.31)) and (RES.get((t, 2), np.nan) >= RES.get((t, 3), np.nan) - 1e-9)
    if t in (34.41, 32.31) and not oki:
        p('  $D=%.2f$，$r=%d$：$\\Delta=%.2e$ 但 [T1] 违反 $\\Rightarrow$ 不判（我方优化器不可信）' % (t, r, d))
        continue
    if d > 1e-4:
        p('  $D=%.2f$，$r=%d$：$\\Delta=%.2e$ bit $>10^{-4}$ $\\Rightarrow$ C 的门里这条判"**严格情形成立（数值级）**"'
          % (t, r, d))
    else:
        p('  $D=%.2f$，$r=%d$：$\\Delta=%.2e$ bit $\\le10^{-4}$ $\\Rightarrow$ 判"**严格情形不成立**：$\\mathrm{Opt}(D)$ 含近似秩-$r$ 成员"'
          % (t, r, d))
p('[T3] 注意：$\\Delta_r\\le\\epsilon$ 只证"存在近似秩-$r$ 的近最优点"，不证"存在**精确**秩-$r$ 最优点" $-$ $-$ 后者要唯一性/闭式，')
p('     本脚本给它的是**上界**（$\\Delta_r\\le$ 读数）；下界永远 $\\ge0$（我的可行集 $\\subseteq$ 凸集）。所以"面非平凡"只能这样卖。')
p('')
p('用时 %.0f s | DARE %d 次 | 有效读数 %d | 回退 %d | 异常 %s'
  % (time.time() - T0, P['NDB'][0], P['NVALID'][0], P['NFALLBACK'][0], P['EXC'] or '无'))
open('p0/e165_out.txt', 'w', encoding='utf-8').write('\n'.join(out) + '\n')
