# -*- coding: utf-8 -*-
r"""e162：秩-1 **射线族上的六位数对账** $+$ C 的 73-C 路线②（$\min_{\mathrm{rank}\le r}\Delta$）首次执行。

板上原文（读全才写）：
* C 补记十一 `c89b`（板 4264 起）：锚点 $D{=}56.66$：$I(S^*)=2.053565$、$I_{\rm cone}=2.148417$、$\Delta I=+9.49\!\times\!10^{-2}$；
  $D{=}80.00$：$1.603105/1.653054/+5.00\!\times\!10^{-2}$；$D{=}34.41$：$4.893329/4.895405/+2.08\!\times\!10^{-3}$。
  C 写明"$k{=}1$ 时锥是一根**射线**，锥内最优形状 $B=[1]$ 恒等 $\\Rightarrow$ 这两格是**精确值**，不含优化器不确定性"。
* 我方 e152（`p0/e152_out.txt` 4–5 行）：同两格自由最优 $2.063303/1.608043$，$\Theta$ 轴 $2.1484/1.6531$，离开角 $5.13°/4.17°$。
$\\Rightarrow$ $\Theta$ 分支对到 $10^{-4}$，但我的"自由最优"比 C 的 $S^*$ **高** $9.7/4.9\,\mathrm{mbit}$。
  秩-1 时双方可行集都是同一族射线（形状被模长吃掉、只剩方向；尺度由 $J=D$ 定）$\\Rightarrow$
  这 $9.7\,\mathrm{mbit}$ **不含模型差别**，只可能是优化质量。谁对要用更狠的搜索判。

冒烟测试已确认（写在跑前）：秩-1 $\Theta$ 首轴在我方定价下 $D{=}56.66\\to I=2.148417$、$D{=}80\\to1.653054$，
与 C 的 $I_{\\rm cone}=2.148417/1.653054$ **六位一致** $\\Rightarrow$ 双方代价轴/常数没有系统性差别，剩下的只是最优端。

维度约定：$Z\\in\\mathbb R^{n\\times r}$，$C=\\sqrt{s}Z^{\\mathsf T}\\in\\mathbb R^{r\\times n}$，$P_m=\\mathrm{DARE}(A^{\\mathsf T},C^{\\mathsf T},W,I_r)$，
$I=\\tfrac12\\log_2\\det(I_r+CP_mC^{\\mathsf T})$，$D=J_C+\\mathrm{tr}(\\Theta P_p)$。秩-$r$ 全体设计 $=$ 全体 $n\\times r$ 矩阵模去整体尺度
（尺度由 $J=D$ 解出）$\\Rightarrow$ 我的形状搜索就是路线②的 $\\min_V$。

预注册判据（跑前写死）：
 [K0] 有效性门：判决只许用 $|D-D_{\rm target}|\le10^{-7}$ 的点；每行打印其背后的有效读数条数与最差残差（#46），为 0 不判。
      DARE 抛出的异常**计数并打印键名**，绝不静默成"不可达"（板账 #43）。
 [K1] $|\min_{\rm mine} I-C|\le10^{-5}$ $\\Rightarrow$ 射线族六位数锁死 $\\Rightarrow$ C 的 `c89b` 两格 $\Delta I$ 独立复现成立、
      其 72-B 的否证成立；**我方 e152 的 $I$ 与离开角同时作废（欠优化）**。
 [K2] $\min_{\rm mine} I<C-10^{-5}$ $\\Rightarrow$ **挑刺成功**：C 的 $S^*$ 不是该 $D$ 的无约束最优，"$k{=}1$ 精确值"不成立。
 [K3] $C+10^{-5}<\min_{\rm mine} I\le e152+10^{-5}$ $\\Rightarrow$ C 对、e152 欠优化 $\\Rightarrow$ 自我更正并报新角。
 [K4] $\min_{\rm mine} I>e152+10^{-5}$ $\\Rightarrow$ 新搜索没覆盖旧搜索 $\\Rightarrow$ 实现问题，物理全部押后。
 [M1] 路线②（$D{=}34.41$，$r^*{=}3$）：$\Delta_r=\min_{\mathrm{rank}\le r}I-I_{\rm unc}$，$I_{\rm unc}$ 自算不借 C。
      判据照 C 的 73-C 原文（$\min$ 量级 $>10^{-4}$ bit 且起点一致，这里取 $\\ge20/24$ 落在 $\min+10^{-4}$ 内），否则只报数。
 [M2] 同格加两笔对账：我的 $I_{\rm unc}$ vs C 的 $4.893329$；锥内（支撑锁 $\Theta$ 前 3 平面、形状自由）vs C 的 $4.895405$。
      $k\ge2$ 锥内形状不唯一 $\\Rightarrow$ 只有 $k{=}1$ 两格能免优化器对账 $-$ 这条要写进板。
"""
import sys
import time
import numpy as np
from scipy.linalg import solve_discrete_are
from scipy.optimize import minimize
sys.stdout.reconfigure(encoding='utf-8')

out = []
T0 = time.time()


def p(*a):
    s = ' '.join(str(x) for x in a)
    out.append(s)
    print(s)
    sys.stdout.flush()


_src = open('p0/exp_c_audit.py', encoding='utf-8').read().split("print(r'== E53")[0]
_ns = {'__name__': 'p'}
exec(compile(_src, 'p0/exp_c_audit.py[preamble]', 'exec'), _ns)

LOG_LO, LOG_HI, RES = -7.0, 7.0, 1e-7
NDB = [0]                    # DARE 调用总数
NVALID = [0]
NFALLBACK = [0]
EXC = {}


def sym(X):
    return 0.5 * (X + X.T)


def ctrl_full(Ax, Bx, Wx, Qt, Rx):
    Pc = sym(solve_discrete_are(Ax, Bx, sym(Qt), Rx))
    K = np.linalg.solve(Rx + Bx.T @ Pc @ Bx, Bx.T @ Pc @ Ax)
    return Pc, K, sym(K.T @ (Rx + Bx.T @ Pc @ Bx) @ K), float(np.trace(Wx @ Pc))


A, B, W = _ns['A'].copy(), _ns['B'].copy(), _ns['W'].copy()
Pc, K, TH, JC = ctrl_full(A, B, W, np.eye(4), np.eye(4))
n = 4
wth, Uth = np.linalg.eigh(TH)
Uth = Uth[:, np.argsort(-wth)]


def make_cur(Zn, r):
    Ir = np.eye(r)

    def cur(lg):
        NDB[0] += 1
        C = (10.0 ** (0.5 * lg)) * Zn.T
        Pm = sym(solve_discrete_are(A.T, C.T, W, Ir))
        Sm = C @ Pm @ C.T + Ir
        Lk = Pm @ C.T @ np.linalg.inv(Sm)
        Pp = sym(Pm - Lk @ Sm @ Lk.T)
        I = 0.5 * np.log(np.linalg.det(Ir + C @ Pm @ C.T)) / np.log(2.0)
        return I, JC + float(np.trace(TH @ Pp))
    return cur


def hit(cur, target, guess):
    """解 $\\lg s$ 使 $D=$ target：割线为主，失败回退全区间的二分（回退计数打印）。"""
    try:
        f0 = cur(guess)[1] - target
        f1 = cur(guess + 0.05)[1] - target
        lg1 = guess + 0.05
        for _ in range(12):
            if abs(f0) <= RES:
                return cur(guess) + (guess,)
            if abs(f1) <= RES:
                return cur(lg1) + (lg1,)
            if f1 == f0:
                break
            step = -f1 * (lg1 - guess) / (f1 - f0)
            step = float(np.clip(step, -0.5, 0.5))
            guess, f0 = lg1, f1
            lg1 = lg1 + step
            f1 = cur(lg1)[1] - target
        if abs(f1) <= RES:
            return cur(lg1) + (lg1,)
    except Exception as e:
        key = type(e).__name__ + ':' + str(e)[:50]
        EXC['guard:' + key] = EXC.get('guard:' + key, 0) + 1
        return np.nan, np.nan, np.nan
    # 回退：全区间二分（说明初值不可用）
    NFALLBACK[0] += 1
    try:
        lo, hi = LOG_LO, LOG_HI
        Dlo = cur(lo)[1]
        Dhi = cur(hi)[1]
    except Exception as e:
        key = type(e).__name__ + ':' + str(e)[:50]
        EXC['fb:' + key] = EXC.get('fb:' + key, 0) + 1
        return np.nan, np.nan, np.nan
    if not (np.isfinite(Dlo) and np.isfinite(Dhi)) or not (Dhi < target < Dlo):
        return np.nan, np.nan, np.nan
    for _ in range(50):
        mid = 0.5 * (lo + hi)
        if cur(mid)[1] > target:
            lo = mid
        else:
            hi = mid
        if abs(cur(0.5 * (lo + hi))[1] - target) <= RES:
            break
    lg = 0.5 * (lo + hi)
    return cur(lg) + (lg,)


GUESS = [0.0]


def read_pt(Zr, target, V=None):
    if V is not None:
        Zr = V @ Zr
    fro = float(np.linalg.norm(Zr))
    if not np.isfinite(fro) or fro <= 0:
        return np.nan, np.nan, np.nan
    Zn = Zr / fro
    r = Zn.shape[1]
    I, D, lg = hit(make_cur(Zn, r), target, GUESS[0])
    if not (np.isfinite(I) and np.isfinite(D)):
        return np.nan, np.nan, np.nan
    if abs(D - target) > RES:
        return np.nan, np.nan, np.nan
    GUESS[0] = lg
    NVALID[0] += 1
    return I, D, lg


def polish(x0, target, V=None, maxfev=900):
    rows = V.shape[1] if V is not None else n

    def f(x):
        I, D, _ = read_pt(x.reshape(rows, -1), target, V)
        return 1e3 if np.isnan(I) else I

    res = minimize(f, np.asarray(x0, float).ravel(), method='Powell',
                   options={'maxfev': maxfev, 'xtol': 1e-9, 'ftol': 1e-12})
    G = res.x.reshape(rows, -1)
    I, D, lg = read_pt(G, target, V)
    return I, D, lg, G


def best(target, starts, V=None, coarse=250, fine=900, ntop=3):
    res = []
    for x0 in starts:
        GUESS[0] = 0.0
        I, D, lg, G = polish(x0, target, V, coarse)
        if np.isfinite(I):
            res.append((I, abs(D - target), G))
    res.sort(key=lambda z: z[0])
    for I, D, G in list(res[:ntop]):
        I2, D2, lg2, G2 = polish(G, target, V, fine)
        if np.isfinite(I2):
            res.append((I2, abs(D2 - target), G2))
    res.sort(key=lambda z: z[0])
    return res


def frames(r, k, rng, d=None):
    if d is None:
        return [np.linalg.qr(rng.normal(0, 1, (n, r)))[0] for _ in range(k)]
    return [rng.normal(0, 1, (d, r)) for _ in range(k)]


CRAB = {56.66: 2.053565, 80.00: 1.603105, 34.41: 4.893329}
CONE_C = {56.66: 2.148417, 80.00: 1.653054, 34.41: 4.895405}
E152 = {56.66: (2.063303, 5.13), 80.00: (1.608043, 4.17)}
rng = np.random.default_rng(20260930)

p('== 锚点 $J_C=%.6f$ spec$\\Theta$=%s ==  [K0] 门 $|D-D_{\\rm tgt}|\\le10^{-7}$；异常计数见尾部 ==' % (JC, np.array2string(wth[::-1], precision=4)))
p('')
p('== [A] 秩-1 射线（$\\Theta$ 四轴 + 12 随机方向，粗 250 + 细 900）==')
VERD = []
for t in (56.66, 80.00):
    starts = [Uth[:, [i]] for i in range(n)] + frames(1, 12, rng)
    res = best(t, starts)
    if not res:
        p('  D=%.2f 有效 0 条 $\\Rightarrow$ 不判（#46）' % t)
        continue
    Ibest, rb, Gbest = res[0]
    u = Gbest[:, 0] / np.linalg.norm(Gbest[:, 0])
    ang = np.degrees(np.arccos(min(1.0, abs(float(u @ Uth[:, 0])))))
    agree = sum(1 for x in res[:16] if x[0] <= Ibest + 1e-4)
    GUESS[0] = 0.0
    Ia, Da, _ = read_pt(Uth[:, [0]], t)
    p('  D=%.2f 有效 %d 条 | $\\min I$=%.6f（最差残差 %.1e；起点一致 %d/%d）' % (t, len(res), Ibest, max(x[1] for x in res), agree, min(16, len(res))))
    p('        vs C $I(S^*)$=%.6f $\\Rightarrow$ 差 %+.2e bit ；vs e152 %.6f $\\Rightarrow$ 差 %+.2e bit' % (CRAB[t], Ibest - CRAB[t], E152[t][0], Ibest - E152[t][0]))
    p('        离开角 %.2f°（e152 %.2f°）| $\\Theta$ 首轴 $I$=%.6f vs C $I_{\\rm cone}$=%.6f $\\Rightarrow$ 差 %+.2e bit' % (ang, E152[t][1], Ia, CONE_C[t], Ia - CONE_C[t]))
    VERD.append((t, Ibest, CRAB[t], E152[t][0], ang, Ia - CONE_C[t]))

p('')
p('== [K] 判决（有效读数 %d 条 / DARE %d 次 / 回退 %d 次）==' % (NVALID[0], NDB[0], NFALLBACK[0]))
for t, Im, Cc, Ee, ang, dth in VERD:
    if abs(Im - Cc) <= 1e-5:
        tag = '[K1] 六位数锁死 $\\Rightarrow$ C 该格成立，**我方 e152 的 $I$/角度作废（欠优化）**'
    elif Im < Cc - 1e-5:
        tag = '[K2] 挑刺成功 $\\Rightarrow$ C 的 $S^*$ 不是无约束最优，"$k{=}1$ 精确值"须重跑'
    elif Im <= Ee + 1e-5:
        tag = '[K3] C 对、e152 欠优化（差 %.2e bit）$\\Rightarrow$ 自我更正' % (Ee - Im)
    else:
        tag = '[K4] 新搜索未覆盖旧搜索 $\\Rightarrow$ 实现问题，物理押后'
    p('  D=%.2f：%s' % (t, tag))
    p('          新离开角 %.2f°（原 %.2f°）；$\\Theta$ 分支六位对账差 %+.2e bit' % (ang, E152[t][1], dth))

p('')
p('== [M] 路线②（$D{=}34.41$，$r^*{=}3$）：$\\Delta_r=\\min_{\\mathrm{rank}\\le r}I-I_{\\rm unc}$，参照自算 ==')
t = 34.41
res3 = best(t, [Uth[:, :3]] + frames(3, 10, rng))
if not res3:
    p('  参照 0 条 $\\Rightarrow$ 整段不判（#46）')
else:
    Iunc = res3[0][0]
    agree = sum(1 for x in res3[:11] if x[0] <= Iunc + 1e-4)
    p('  $I_{\\rm unc}$（$r{=}3$ 自算，%d 条有效，起点一致 %d/11）=%.6f vs C $I(S^*)$=%.6f $\\Rightarrow$ 差 %+.2e bit'
      % (len(res3), agree, Iunc, CRAB[t], Iunc - CRAB[t]))
    stc = [np.eye(3)] + frames(1, 6, rng, d=3) + frames(2, 6, rng, d=3) + frames(3, 4, rng, d=3)
    resc = best(t, stc, V=Uth[:, :3], fine=1200)
    if resc:
        Ic = resc[0][0]
        p('  [M2] 锥内 $I_{\\rm cone}$=%.6f vs C %.6f $\\Rightarrow$ 差 %+.2e bit；我方 $\\Delta I=%+.2e$ vs C $+2.08\\times10^{-3}$'
          % (Ic, CONE_C[t], Ic - CONE_C[t], Ic - Iunc))
    else:
        p('  [M2] 锥内 0 条有效 $\\Rightarrow$ 不判')
    for r in (1, 2):
        st = ([Uth[:, :r]] if r > 1 else [Uth[:, [i]] for i in range(n)]) + frames(r, 20, rng)
        res = best(t, st)
        if not res:
            p('     r=%d 有效 0 条 $\\Rightarrow$ 不判' % r)
            continue
        mn = res[0][0]
        ag = sum(1 for x in res[:21] if x[0] <= mn + 1e-4)
        d = mn - Iunc
        p('     r=%d（$r^*{=}3$，$r<r^*$）$\\min_{\\mathrm{rank}\\le r} I=%.6f $\\Rightarrow$ $\\Delta_r=%+.2e$ bit；起点一致 %d/21' % (r, mn, d, ag))
        p('        [M1] 按 C 的 73-C 判据 $\\Rightarrow$ %s' % ('**严格情形成立（数值级）**' if (d > 1e-4 and ag >= 20) else '面的结构未定（只报数）'))
p('  D=32.31：$r^*{=}4\\Rightarrow$ 参照需**满秩**形状搜索（16 参数），本轮预算外 $\\Rightarrow$ 只报"未算"，不猜。')
p('  D=56.66/80.00：$r^*{=}1\\Rightarrow$ 无 $r<r^*$ 的格 $\\Rightarrow$ 路线②在这两格为空判据。')
p('')
p('用时 %.0f s | DARE %d 次 | 有效读数 %d | 回退 %d | 异常键 %s' % (time.time() - T0, NDB[0], NVALID[0], NFALLBACK[0], EXC if EXC else '无'))
open('p0/e162_out.txt', 'w', encoding='utf-8').write('\n'.join(out) + '\n')
