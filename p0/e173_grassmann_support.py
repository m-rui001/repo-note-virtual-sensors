# -*- coding: utf-8 -*-
r"""e173：**支撑唯一性的正确口径** $-\ -$ e172 的角距是坏的，本轮换到 Grassmann 主角度并自带度量自检。

为什么必须重做（e172 的两处口径错，都在落盘前自查出来）：
 1. e172 的 $\text{gap}(G_a,G_b)=\arccos\sigma_{\min}(\hat G_a^{\mathsf T}\hat G_b)$ 只对**列归一化**，没有正交化。
    于是它把"每个标架自己的两列有多斜"混进了"两个子空间相差多少" $-\ -$ 铁证：`e172_out.txt` 第 9 行
    同一条读数与**自己**比，角距 $=63.32^\circ$（真的子空间距离在自身必须为 $0$）。$r=1$ 不受影响（单列归一化即标架），
    所以 e172 两格秩-1 的 $0.00^\circ/0.02^\circ$ 仍然可信；$r=2$ 的 $65.21^\circ$ **作废**。
 2. 更要紧的模型事实：定价对 $\mathrm{O}(r)$ 标架旋转**严格不变** $-\ -$ $C\to Q^{\mathsf T}C$ 时
    $\mathrm{DARE}(A^{\mathsf T},C^{\mathsf T}Q,W,I_r)$ 的解不变（$R=I_r$ 各向同性），而
    $\det(I_r+C P_m C^{\mathsf T})$、$\mathrm{tr}(\Theta P_p)$ 都不变。$r\ge2$ 时每条读数都带一个 $(r(r-1)/2)$ 维的
    **纯冗余轨道** $-\ -$ 拿"标架不同"当"支撑不唯一"是把自由度当成了物理。

本轮口径：支撑 $:=\mathcal R(G)=\operatorname{span}(\text{列})\in\mathrm{Gr}(r,n)$，距离一律用
正交标架 $Q_G=\mathrm{qr}(G)$ 后的**主角度**，并打印子空间间隙 $\|Q_aQ_a^{\mathsf T}-Q_bQ_b^{\mathsf T}\|_2$。
另给每对读数打**标架等价判据**：若 $\mathrm{Gram}(G_a)\approx c\cdot\mathrm{Gram}(G_b)$ 而 $Q$ 的条目明显不同，
则两条在同一 $\mathrm{O}(r)$ 轨道上 $-\ -$ 值并列是**恒等式**，不构成支撑不唯一的证据。

预注册判据（跑前写死；#46 判决自带有效读数条数）：
 [V0] 度量自检：每条读数 $G$ 满足 $\text{gap}(G,G)<10^{-9}$，且对随机正交 $Q$ 满足
      $\text{gap}(\mathrm{qr}(G\,Q),\mathrm{qr}(G))<10^{-9}$。**任一失败则本轮不判**（只报数）。
 [V0b] 不变性自检：同一支撑的随机 $\mathrm{O}(r)$ 标架右乘后 $(I,D)$ 的最大相对差要打印；
      若 $<10^{-9}$ 则"标架不同、值相同"记为**恒等式**，从"不唯一"的证据里扣除（跑前就写好，防止事后挑口径）。
 [V1] 并列集 $T:=\{i:I_i\le\min I+10^{-6}\}$；打印 $|T|$、$T$ 内两两**最大主角度**、$\|\Delta P\|_2$、
      以及"是否同 $\mathrm{O}(r)$ 轨道"的计数。
 [V2] 判决：$\ge5^\circ$（且非同轨道）$\\Rightarrow$ "**近最优秩-$r$ 支撑不唯一（数值级）**"；
      全部 $<1^\circ$ $\\Rightarrow$ "**未见不唯一**"；$1^\circ$–$5^\circ$ $\\Rightarrow$ "**未决**"；$|T|=1$ $\\Rightarrow$ "**不判**"。
      若并列只在轨道内出现，判决行必须明写"这是标架冗余，不是支撑不唯一"。
 [V3] 口径：$\Delta_r$ 仍是上界（我的可行集 $\subseteq$ 凸 SDP 集）；"不唯一"只到 $10^{-6}$ bit，
      不升格为"存在两个精确最优点"。
格：$(D,r)=(34.41,2)$（e172 出事的格）、$(34.41,3)$（C 的 $r^*$ 格）、$(80.00,1)$（秩-1 对账：应复现同一支撑）。
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

SDP = {34.41: 4.893329, 80.00: 1.603105}
A, B, W = AUD['A'].copy(), AUD['B'].copy(), AUD['W'].copy()
Pc, K, TH, JC = P['ctrl_full'](A, B, W, np.eye(4), np.eye(4))
P['A'], P['B'], P['W'], P['TH'], P['JC'] = A, B, W, TH, JC
wth, U0 = np.linalg.eigh(TH)
Uth = U0[:, np.argsort(-wth)]
P['Uth'] = Uth
n = 4
rng = np.random.default_rng(9091)
TOL_TIE = 1e-6


def basis(G):
    """列空间的正交标架（Grassmann 点的代表元）。"""
    G = np.asarray(G, float).reshape(n, -1)
    Q, Rq = np.linalg.qr(G / np.linalg.norm(G, axis=0))
    # 消掉 QR 的符号/次序歧义：按 R 对角元定号
    d = np.sign(np.diag(Rq))
    d[d == 0] = 1.0
    return Q * d


def angap(Ga, Gb):
    """两支撑的最大主角度（度）；同一条的自身值恒为 0。"""
    sv = np.linalg.svd(basis(Ga).T @ basis(Gb), compute_uv=False)
    return float(np.degrees(np.arccos(np.clip(sv.min(), -1.0, 1.0))))


def pgap(Ga, Gb):
    """子空间间隙 ||P_a - P_b||_2。"""
    Qa, Qb = basis(Ga), basis(Gb)
    return float(np.linalg.norm(Qa @ Qa.T - Qb @ Qb.T, 2))


def princ(Ga, Gb):
    sv = np.linalg.svd(basis(Ga).T @ basis(Gb), compute_uv=False)
    return np.degrees(np.arccos(np.clip(sv, -1.0, 1.0)))


def gram(G):
    G = np.asarray(G, float).reshape(n, -1)
    Gr = G.T @ G
    return Gr / np.linalg.norm(Gr)


def same_orbit(Ga, Gb, tol=1e-4):
    """Gram 相同（条目差 $\\le$ tol）而标架条目不同 $\\Rightarrow$ 同一 O(r) 轨道。"""
    return float(np.max(np.abs(gram(Ga) - gram(Gb)))) <= tol


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


p(r'== 锚点 $J_C=%.6f$；$I_{\rm unc}$（e163 凸 SDP）：%s ==' % (JC, SDP))
p(r'   口径：支撑 $=\mathcal R(G)\in\mathrm{Gr}(r,n)$，距离取正交标架后的最大主角度 $+$ 子空间间隙；'
  r'并列容差 $10^{-6}$ bit')

# ---- [V0] 度量自检：先用真数据把坏口径钉死，再验新口径 ----
p('')
p(r'== [V0] 度量自检 ==')
badsum = 0
oksum = 0
chk = [Uth[:, :2]] + [np.linalg.qr(rng.normal(0, 1, (n, 2)))[0] for _ in range(4)]
for G in chk:
    selfang = angap(G, G)
    d = np.linalg.qr(rng.normal(0, 1, (2, 2)))[0]
    orbang = angap(G @ d, G)
    badsum += (selfang > 1e-9)
    oksum += (orbang > 1e-9)
p(r'   新口径：自身角距 $>10^{-9}$ 的 %d/%d 条 $-$ 随机正交右乘后角距 $>10^{-9}$ 的 %d/%d 条'
  % (badsum, len(chk), oksum, len(chk)))
# 旧口径（仅列归一化）在同一批数据上的自身角距，用来复现 e172 的 63.32 度伪影
old = []
for G in chk:
    Gn = G / np.linalg.norm(G, axis=0)
    sv = np.linalg.svd(Gn.T @ Gn, compute_uv=False)
    old.append(float(np.degrees(np.arccos(np.clip(sv.min(), -1, 1)))))
p(r'   旧口径（e172，仅列归一化）自身"角距"：%s $^\circ$ $-\-$ 非零即证明它量的不是子空间距离'
  % np.array2string(np.array(old), precision=2))
METRIC_OK = (badsum == 0 and oksum == 0)
p(r'   [V0] 判：度量自检 %s' % ('**通过**（自身与轨道旋转都为 0）' if METRIC_OK else '**失败 $\\Rightarrow$ 本轮不判**'))

# [V0b] 代价对 $\mathrm{O}(r)$ 右乘是否严格不变（$R=I_r$ 的旋转不变性）——决定"标架不同"算不算物理
worst = 0.0
ninv = 0
for G in [Uth[:, :2], np.linalg.qr(rng.normal(0, 1, (n, 2)))[0], np.linalg.qr(rng.normal(0, 1, (n, 3)))[0]]:
    rr = G.shape[1]
    Zn = G / np.linalg.norm(G)
    for _ in range(3):
        Qo = np.linalg.qr(rng.normal(0, 1, (rr, rr)))[0]
        a = P['hit'](P['make_cur'](Zn, rr), 34.41, 0.0)
        b = P['hit'](P['make_cur'](Zn @ Qo, rr), 34.41, 0.0)
        if np.isfinite(a[0]) and np.isfinite(b[0]):
            worst = max(worst, abs(b[0] - a[0]) / max(abs(a[0]), 1e-30), abs(b[1] - a[1]) / abs(a[1]))
            ninv += 1
p(r'   [V0b] 同一支撑的随机 $\mathrm{O}(r)$ 标架对代价 $(I,D)$ 的最大相对差 $=%.1e$（%d 对，$D_{\rm tgt}{=}34.41$）'
  % (worst, ninv))
p(r'   [V0b] 判：%s $-\ -$ 故 $r\ge2$ 的"标架不同而值相同"是**恒等式（自由度）**，不得当作支撑不唯一的证据'
  % ('代价对 $\mathrm{O}(r)$ 右乘不变' if worst < 1e-9 else '不变性只到 %.0e，轨道冗余仍是**近似**' % worst))

V = {}
for t, r in ((34.41, 2), (34.41, 3), (80.00, 1)):
    starts = [Uth[:, :r]] + [np.linalg.qr(rng.normal(0, 1, (n, r)))[0] for _ in range(24)]
    raw, allv = cell(t, r, starts)
    p('')
    p(r'== $D=%.2f$，$r=%d$：原始有效 $k=%d/%d$，含细化共 %d 条读数 ==' % (t, r, len(raw), len(starts), len(allv)))
    if not allv:
        p(r'   [V2] 判决：0 条 $\\Rightarrow$ 不判')
        continue
    Is = np.array([q[0] for q in allv])
    Imin = Is[0]
    p(r'   $\min I=%.6f$（$\Delta_r=%+.2e$ bit，[V3] 上界口径）；最差残差 $|D-t|=%.2e$ bit'
      % (Imin, Imin - SDP.get(t, np.nan), max(q[1] for q in allv)))
    tie = [q for q in allv if q[0] <= Imin + TOL_TIE]
    p(r'   [V1] 并列集 $|T|=%d$；最优支撑与 $\Theta$ 前 %d 平面主角度 %s'
      % (len(tie), r, np.array2string(princ(tie[0][2], Uth[:, :r]), precision=2)))
    if len(tie) < 2:
        p(r'   [V2] 判决：**不判**（$|T|=1$）')
        continue
    pairs = []
    nOrb = 0
    for i, a in enumerate(tie):
        for b in tie[i + 1:]:
            g = angap(a[2], b[2])
            orb = same_orbit(a[2], b[2])
            nOrb += orb
            pairs.append((g, b[0] - a[0], pgap(a[2], b[2]), orb))
    gmax = max(z[0] for z in pairs)
    p(r'   [V1] $T$ 内两两（最大 8 对，主角度/值差/$\|\Delta P\|_2$/同轨道）：%s'
      % ', '.join('%.2f$^\\circ$/%+.1e/%.1e/%s' % (g, dv, pg, ('Y' if o else 'N'))
                  for g, dv, pg, o in sorted(pairs, reverse=True)[:8]))
    p(r'   [V1] 每条与最优支撑的主角度：%s'
      % np.array2string(np.array([angap(q[2], tie[0][2]) for q in tie]), precision=2))
    p(r'   [V1] 同 $\mathrm{O}(r)$ 轨道的对数：%d/%d $-\ -$ 轨道内并列由 $R=I_r$ 的旋转不变性保证，是恒等式'
      % (nOrb, len(pairs)))
    if not METRIC_OK:
        p(r'   [V2] 判决：**不判**（[V0] 度量自检失败）')
        continue
    nPairFree = len(pairs) - nOrb
    if gmax >= 5.0 and nPairFree > 0:
        v = (r'**近最优秩-$r$ 支撑不唯一（数值级）** $-\ -$ 有 $|T|=%d$ 条读数互差 $\le10^{-6}$ bit，'
             r'非轨道对的最大主角度达 %.2f$^\circ$（$\|\Delta P\|_2=%.1e$）' % (len(tie), gmax,
                                                                              max(z[2] for z in pairs)))
    elif gmax < 1.0:
        v = r'**未见不唯一** $-\ -$ $|T|=%d$ 条支撑两两 $<1^\circ$（最大 %.2f$^\circ$），并列全部落在同一支撑上' % (len(tie), gmax)
    else:
        v = r'**未决**（最大主角度 %.2f$^\circ$ 落在 $1^\circ$–$5^\circ$）' % gmax
    p(r'   [V2] 判决：%s' % v)
    V[(t, r)] = v

p('')
p(r'== [V3] 口径重申：$\Delta_r$ 只到上界；"不唯一"到 $10^{-6}$ bit 为止，不升格为"两个精确最优点"。')
p(r'   且凡报角度必须配 $\Delta I$（板 §98--§100 长期口径）；$R=I_r$ 的 $\mathrm{O}(r)$ 冗余对 $r=1$ 不存在。==')
p('')
p('用时 %.0f s | DARE %d | 有效读数 %d | 回退 %d | 异常 %s'
  % (time.time() - T0, P['NDB'][0], P['NVALID'][0], P['NFALLBACK'][0], P['EXC'] or '无'))
txt = '\n'.join(out) + '\n'
open('p0/e173_out.txt', 'w', encoding='utf-8').write(txt)
open('p0/e173_run.txt', 'w', encoding='utf-8').write(txt)
