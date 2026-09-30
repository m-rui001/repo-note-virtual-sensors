# -*- coding: utf-8 -*-
r"""e174：秩-2 格 $(D{=}34.41)$ 的**盆地逃逸**测支撑唯一性 $-\ -$ 随机起点在这格只有 $1/25$ 有效，功率不够。

背景（全部有落盘出处）：
 * e172 用"仅列归一化"的 $\arccos\sigma_{\min}$ 当子空间距离 $\\Rightarrow$ 量的其实是标架内部的斜程度 $-\ -$
   铁证是同一条读数与**自己**比得 $63.32^\circ$（`p0/e172_out.txt` 第 9 行）。它的 $r=2$ 的 $65.21^\circ$ 作废。
 * e173 换正交标架后的 Grassmann 主角度，并实测代价对 $\mathrm{O}(r)$ 右乘不变（相对差 $2.9\times10^{-15}$，$R=I_r$ 各向同性）
   $\\Rightarrow$ $r\ge2$ 的"标架不同而值相同"是恒等式。但 e173 的 $r=2$ 格 $|T|=1$（$k=1/25$）$\\Rightarrow$ 不判。
 * e173 的 $r=3$ 格有功率（$k=24/25$、30 条）：并列 2 条、支撑只差 $0.16^\circ$ $\\Rightarrow$ 无支撑不唯一的迹象。

本轮把 $r=2$ 格做**有功率**的检验：先定出最优支撑 $S^*$，再从 $S^*$ 沿 Grassmann 方向**受控位移** $\theta_0$
（$5^\circ$ 到 $80^\circ$）与**内部形状**扰动（对角加权）各若干条重跑 Powell，看它们收敛回哪里、值差多少。
判"不唯一"要的是"从别处也能落到同样好的值"，逃逸重跑比随机撒点直接得多。

口径更正（跑前写死，并打印浮点地板为证）：e173 的自检门 $10^{-9}$ **度** 在 float64 下不可达
（$\arccos$ 在 $1$ 附近导数发散：$\sim\!10^{-16}$ 的奇异值误差就给出 $\sim10^{-6}$ 度的角）。
本轮 [W0] 用 $10^{-3}$ 度，并把自检实测值打出来——这是**数值标定**更正，不是看结果挑门。

预注册判据（#46：每条判决自带有效读数条数）：
 [W0] 度量自检：$\text{gap}(G,G)<10^{-3}$ 度、$\text{gap}(\mathrm{qr}(GQ),G)<10^{-3}$ 度（打印实际量级）。
      失败则本轮全部只报数、不判。
 [W1] $S^*:=$ 全部读数里 $\min I$ 那条的支撑；打印它与 $\Theta$ 前 2 平面的主角度、$\Delta_2$、有效条数。
 [W2] 逃逸集 $E:=\{$ 重跑后 $I\le I(S^*)+10^{-6}$ 且其支撑与 $S^*$ 主角度 $<0.5^\circ$ 之外重 $\\}$ 的补：
      存在 $I\le I(S^*)+10^{-6}$ 且主角度 $\ge5^\circ$ 的读数 $\\Rightarrow$ 判"**近最优秩-2 支撑不唯一（数值级）**"；
      否则全部 $<1^\circ$ $\\Rightarrow$ 判"**单盆地：$10^{-6}$ bit 口径下未见支撑不唯一**"；
      $1^\circ$–$5^\circ$ $\\Rightarrow$ "**未决**"。每条判决打印背后的有效读数条数。
 [W3] 附带报"回落率"：按初始位移 $\theta_0$ 分桶，统计收敛回 $S^*$（$<0.5^\circ$）的比例 $-\ -$ 盆地的形状。
 [W4] 口径：$\Delta_2$ 只是上界；不升格为"精确唯一"，唯一性定理仍欠（C 的 73-A 归约成立）。
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

SDP = {34.41: 4.893329}
A, B, W = AUD['A'].copy(), AUD['B'].copy(), AUD['W'].copy()
Pc, K, TH, JC = P['ctrl_full'](A, B, W, np.eye(4), np.eye(4))
P['A'], P['B'], P['W'], P['TH'], P['JC'] = A, B, W, TH, JC
wth, U0 = np.linalg.eigh(TH)
Uth = U0[:, np.argsort(-wth)]
P['Uth'] = Uth
n, r, tgt = 4, 2, 34.41
rng = np.random.default_rng(20260930)
TOL_TIE = 1e-6
GATE_W0 = 1e-3      # 度；float64 的 arccos 地板远高于 1e-9，见 [W0] 打印


def basis(G):
    G = np.asarray(G, float).reshape(n, -1)
    Q, Rq = np.linalg.qr(G / np.linalg.norm(G, axis=0))
    d = np.sign(np.diag(Rq))
    d[d == 0] = 1.0
    return Q * d


def angles(Ga, Gb):
    sv = np.linalg.svd(basis(Ga).T @ basis(Gb), compute_uv=False)
    return np.degrees(np.arccos(np.clip(sv, -1.0, 1.0)))


def maxgap(Ga, Gb):
    return float(angles(Ga, Gb).max())


def shape(G):
    """内部形状 $=$ Gram 的归一化特征值 $-\ -$ 支撑之外的物理自由度（$\mathrm{O}(r)$ 轨道上不变）。"""
    Gr = np.asarray(G, float).T @ np.asarray(G, float)
    ev = np.linalg.eigvalsh(Gr / np.linalg.norm(Gr))
    return ev[::-1]


p(r'== 锚点 $J_C=%.6f$；格 $D=%.2f$，$r=%d$；$I_{\rm unc}=%.6f$ ==' % (JC, tgt, r, SDP[tgt]))

# ---------- [W0] 度量自检（含浮点地板） ----------
p('')
p(r'== [W0] 度量自检：自身角距与 $\mathrm{O}(r)$ 轨道角距的浮点地板 ==')
tests = [Uth[:, :2], np.linalg.qr(rng.normal(0, 1, (n, 2)))[0]]
selfs, orbs = [], []
for G in tests:
    selfs.append(maxgap(G, G))
    Q = np.linalg.qr(rng.normal(0, 1, (2, 2)))[0]
    orbs.append(maxgap(G @ Q, G))
# 再用真实 Powell 输出（标架很斜）测一次
P['GUESS'][0] = 0.0
I0, D0, lg0, G0 = P['polish'](Uth[:, :2], tgt, None, 900)
if np.isfinite(I0):
    selfs.append(maxgap(G0, G0))
    Q = np.linalg.qr(rng.normal(0, 1, (2, 2)))[0]
    orbs.append(maxgap(G0 @ Q, G0))
    p(r'   原始 Powell 输出 $G_*$ 的 Gram 特征值（归一）：%s $-\-$ 非 $[1,1]$ 即标架不正交，'
      r'e172 的 $63.32^\circ$ 就是这么来的' % np.array2string(shape(G0), precision=3))
p(r'   新口径实测：自身角距最大 %.2e 度 $-$ 轨道右乘角距最大 %.2e 度 $\\Rightarrow$ 门取 $10^{-3}$ 度'
  % (max(selfs), max(orbs)))
W0OK = (max(selfs) < GATE_W0 and max(orbs) < GATE_W0)
p(r'   [W0] 判：自检 %s（%d 条，浮点地板 %.1e 度 $<$ 门 $10^{-3}$ 度）'
  % ('**通过**' if W0OK else '**失败 $\\Rightarrow$ 只报数**', len(selfs), max(selfs)))

# ---------- [W1] 定最优支撑 ----------
p('')
p(r'== [W1] 最优支撑 $S^*$（$\Theta$ 首起点 $+\ 20$ 条随机正交标架，coarse=700 $\\to$ fine=3000）==')
starts = [Uth[:, :2]] + [np.linalg.qr(rng.normal(0, 1, (n, r)))[0] for _ in range(20)]
raw = []
for x0 in starts:
    P['GUESS'][0] = 0.0
    I, D, lg, G = P['polish'](x0, tgt, None, 700)
    if np.isfinite(I):
        raw.append((I, abs(D - tgt), G))
raw.sort(key=lambda z: z[0])
ref = []
for I, D, G in raw[:8]:
    P['GUESS'][0] = 0.0
    I2, D2, lg2, G2 = P['polish'](G, tgt, None, 3000)
    if np.isfinite(I2):
        ref.append((I2, abs(D2 - tgt), G2))
allv = sorted(raw + ref, key=lambda z: z[0])
p(r'   原始有效 $k=%d/%d$，含细化共 %d 条' % (len(raw), len(starts), len(allv)))
if not allv or not W0OK:
    p(r'   [W2] 判决：**不判**（%s）' % ('无有效读数' if not allv else '[W0] 自检失败'))
    sys.exit(0)
Ibest, rb, Gbest = allv[0]
p(r'   $\min I=%.6f$（$\Delta_2=%+.2e$ bit，[W4] 上界口径）；最差残差 $|D-t|=%.2e$'
  % (Ibest, Ibest - SDP[tgt], max(q[1] for q in allv)))
p(r'   $S^*$ 与 $\Theta$ 前 2 平面主角度 %s $^\circ$；$S^*$ 的 Gram 特征值 %s'
  % (np.array2string(angles(Gbest, Uth[:, :2]), precision=2),
     np.array2string(shape(Gbest), precision=3)))

# ---------- [W2/W3] 盆地逃逸 ----------
p('')
p(r'== [W2/W3] 从 $S^*$ 受控位移逃逸：Grassmann 旋转 $\theta_0\in\{5,10,20,40,60,80\}^\circ$（各 4 条）'
  r'$+$ 内部形状扰动 $\mathrm{diag}(a,1/a)$，$a\in\{2,5,20,60\}$（各 3 条）==')
Qb = basis(Gbest)
esc = []
thetas = [5.0, 10.0, 20.0, 40.0, 60.0, 80.0]
for th in thetas:
    for _ in range(4):
        N = rng.normal(0, 1, (n, r))
        Np = N - Qb @ (Qb.T @ N)
        Np = Np / np.linalg.norm(Np, axis=0)
        G0 = np.cos(np.radians(th)) * Qb + np.sin(np.radians(th)) * Np
        th0 = maxgap(G0, Gbest)
        P['GUESS'][0] = 0.0
        I, D, lg, G = P['polish'](G0, tgt, None, 3000)
        if np.isfinite(I):
            esc.append(('rot%.0f' % th, th0, I, abs(D - tgt), maxgap(G, Gbest), G))
for a in [2.0, 5.0, 20.0, 60.0]:
    for _ in range(3):
        G0 = Gbest @ np.diag([a, 1.0 / a])
        th0 = maxgap(G0, Gbest)
        P['GUESS'][0] = 0.0
        I, D, lg, G = P['polish'](G0, tgt, None, 3000)
        if np.isfinite(I):
            esc.append(('shape%g' % a, th0, I, abs(D - tgt), maxgap(G, Gbest), G))
p(r'   逃逸重跑有效 %d 条（共发起 %d 条）' % (len(esc), len(thetas) * 4 + 4 * 3))

vals = np.array([q[2] for q in esc])
resd = np.array([q[3] for q in esc])
mgs = np.array([q[4] for q in esc])
good = (resd <= 1e-7)
tieE = np.where(good & (vals <= Ibest + TOL_TIE))[0]
p(r'   残差 $|D-t|$ 最大 %.2e（**代价单位**，非 bit；`read_pt` 已按 $\mathrm{RES}=10^{-7}$ 卡死）$\\Rightarrow$ '
  r'有效 %d 条 $-$ 其中落入 $\min I+10^{-6}$ **bit** 并列带的 %d 条'
  % (resd.max() if len(resd) else float('nan'), int(good.sum()), len(tieE)))
if len(tieE):
    lines = ', '.join('%.2f$^\\circ$/%+.1e/%s' % (mgs[i], vals[i] - Ibest, esc[i][0]) for i in tieE[:10])
    p(r'   并列条（与 $S^*$ 主角度 / 值差 / 来源）：%s' % lines)
band4 = np.where(good & (vals <= Ibest + 1e-4))[0]
p(r'   放宽到 $10^{-4}$ bit：并列 %d 条，最大主角度 %.2f$^\circ$' % (len(band4), mgs[band4].max() if len(band4) else float('nan')))

far = tieE[mgs[tieE] >= 5.0] if len(tieE) else np.array([], dtype=int)
mid = tieE[(mgs[tieE] >= 1.0) & (mgs[tieE] < 5.0)] if len(tieE) else np.array([], dtype=int)
near = tieE[mgs[tieE] < 1.0] if len(tieE) else np.array([], dtype=int)
if len(far):
    p(r'   [W2] 判决：**近最优秩-2 支撑不唯一（数值级）** $-\ -$ %d 条读数（含逃逸重跑）落在 $\min I+10^{-6}$ bit 内'
      r'而支撑与 $S^*$ 差 $\ge5^\circ$（最大 %.2f$^\circ$）' % (len(far), mgs[far].max()))
elif len(tieE) and len(mid) == 0:
    p(r'   [W2] 判决：**单盆地（$10^{-6}$ bit 口径下未见支撑不唯一）** $-\ -$ 有效读数 %d 条，'
      r'并列带 %d 条与 $S^*$ 的支撑角距全部 $<1^\circ$（最大 %.2f$^\circ$）'
      % (int(good.sum()), len(near), mgs[near].max() if len(near) else 0.0))
elif len(mid):
    p(r'   [W2] 判决：**未决** $-\ -$ %d 条并列读数落在 $1^\circ$–$5^\circ$（最大 %.2f$^\circ$）'
      % (len(mid), mgs[mid].max()))
else:
    p(r'   [W2] 判决：**不判**（并列带 0 条）')

# 回落率：按初始位移分桶
p('')
p(r'== [W3] 回落率（收敛回 $S^*$ 的 $0.5^\circ$ 邻域）按初始位移 $\theta_0$ 分桶 ==')
for th in thetas:
    sel = [q for q in esc if q[0] == 'rot%.0f' % th]
    if sel:
        back = sum(1 for q in sel if q[4] < 0.5)
        p(r'   $\theta_0\approx%.0f^\circ$（实测初位移 %.1f$^\circ$）：%d 条有效，回落 %d 条 $-$ '
          r'终值与 $S^*$ 最大角距 %.2f$^\circ$，值差中位 %+.1e bit'
          % (th, np.mean([q[1] for q in sel]), len(sel), back,
             max(q[4] for q in sel), np.median([q[2] - Ibest for q in sel])))
sh = [q for q in esc if q[0].startswith('shape')]
if sh:
    p(r'   形状扰动（支撑同、内部加权不同）：%d 条有效，支撑角距最大 %.2f$^\circ$，值差最大 %+.1e bit'
      % (len(sh), max(q[4] for q in sh), max(q[2] - Ibest for q in sh)))

p('')
p(r'== [W4] 口径：$\Delta_2$ 只到上界；"单盆地"是**数值级**陈述（$10^{-6}$ bit、Powell、40 条逃逸），')
p(r'   不是唯一性定理 $-\ -$ 后者需要严格凸性或闭式，本轮没有。C 的 73-A 归约照旧成立。==')
p('')
p('用时 %.0f s | DARE %d | 有效读数 %d | 回退 %d | 异常 %s'
  % (time.time() - T0, P['NDB'][0], P['NVALID'][0], P['NFALLBACK'][0], P['EXC'] or '无'))
txt = '\n'.join(out) + '\n'
open('p0/e174_out.txt', 'w', encoding='utf-8').write(txt)
open('p0/e174_run.txt', 'w', encoding='utf-8').write(txt)
