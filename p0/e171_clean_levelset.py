# -*- coding: utf-8 -*-
r"""e171：分支无关、归一化正确的秩-$r$ **等值带**扫描（把 e165/e162b/e167 三角度读数一次对齐）。

为什么要重做（三处我自己的读数互相矛盾，全部有明确根因）：
 * e165 $D{=}80,r{=}1$ 报离开角 $62.5^\circ$，e162b 同格报 $6.02^\circ$，e167 细化后又报 $[0.]$；
 * e167 的"位移 $0.00^\circ$"是**我的 bug**：角度用的两个向量没归一化，$|\cos|\gg1$ 被 clip 成 1 $\\Rightarrow$ 恒读 $0^\circ$；
   e165 的 `pa_vs_th` 同样未归一化 $\\Rightarrow$ **那边的角整体偏小（朝 $0$ 塌）**，$62.5^\circ$ 可信、$[0.]$ 存疑；
 * e162/e165 的增益用割线 $+$ 携带 GUESS 求 $J=D_{\rm tgt}$，e169 已证 $D(\lg)$ 对固定方向**单调降**（$\Theta$ 首轴 $0/48$ 上升段）
   $\\Rightarrow$ 本轮**改用宽窗二分**（$[-14,10]$），每个方向的 $I(u)$ 成为**确定函数**，不再有"落在哪一支"的问题。

判据（跑前写死；#46 每条判决自带读数条数）：
 [S0] 一条方向"有效" $:=$ 粗网格 $\lg\in[-6,10]$（33 点）上存在相邻两点夹住 $t$（$D$ 降故穿点唯一），且二分后 $|D-t|\le10^{-9}$；
      逐格打印扫向数/有效数/无穿点数/求解失败数。**窗为什么不能开到 $-14$**：小增益端 `solve_discrete_are` 对随机方向抛
      `Failed to find a finite solution`（第一版死在 $\lg=-14$），故改自适应找括号 $-\ -$ 这是**模型性质**（秩-1 在弱测量下无有限 stabilizing 根）
      而非代码问题，本轮把它单独计数。
 [S1] 单调性侧检：20 条随机方向各打 13 点子网格，$D$ 上升段计数 $\\Rightarrow$ 若非 0，本脚本的二分前提失效，只报数不判。
 [S2] 等值带 $\epsilon\in\{10^{-6},10^{-4},10^{-2}\}$ bit：带内方向数 $+$ 带内**最大两两角距**（列已归一化）。
 [S3] 判决：$\epsilon=10^{-6}$ 带内存在两条角距 $\ge5^\circ$ 的方向 $\\Rightarrow$ 判"**近最优秩-$r$ 支撑不唯一（数值级）**"；
      带内全部两两 $<1^\circ$ $\\Rightarrow$ 判"**未见不唯一**"；带内只有 1 条或 0 条 $\\Rightarrow$ 判"**本扫功率不足，不判**"。
 [S4] 同时打印最优方向与 $\Theta$ 前 $r$ 平面的主角度（归一化后）$-\ -$ 这就是与 e162b $6.02^\circ$、e165 $62.5^\circ$ 对账的那个数。
 [S5] 顺带用"增益饱和端 $\min D$"再测一次同秩可达地板 $-\ -$ 与 e162b [P]（$s=10^7$、120 标架）应是同一个数（不同机制，互为复核）。
网格：$r=1$ 用 $\Theta$ 四轴 $+\,800$ 随机单位射线；$r=2$ 用 $\Theta$ 前 2 列 $+\,400$ 随机正交 2-标架；$r=3$ 用 250 个。
格：$(D,r)\in\{(80,1),(56.66,1),(34.41,2),(34.41,3)\}$。
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
Pc, K0, TH, JC = P['ctrl_full'](A, B, W, np.eye(4), np.eye(4))
P['A'], P['B'], P['W'], P['TH'], P['JC'] = A, B, W, TH, JC
wth, U0 = np.linalg.eigh(TH)
Uth = U0[:, np.argsort(-wth)]
P['Uth'] = Uth
n = 4
rng = np.random.default_rng(20260930)
LO, HI = -14.0, 10.0
EPS = [1e-6, 1e-4, 1e-2]


def sym(X):
    return 0.5 * (X + X.T)


def cur_of(Zn, r):
    Ir = np.eye(r)

    def cur(lg):
        P['NDB'][0] += 1
        C = (10.0 ** (0.5 * lg)) * Zn.T
        Pm = sym(P['solve_discrete_are'](A.T, C.T, W, Ir))
        Sm = C @ Pm @ C.T + Ir
        Lk = Pm @ C.T @ np.linalg.inv(Sm)
        Pp = sym(Pm - Lk @ Sm @ Lk.T)
        I = 0.5 * np.log(np.linalg.det(Ir + C @ Pm @ C.T)) / np.log(2.0)
        return I, JC + float(np.trace(TH @ Pp))
    return cur


def unit(G):
    return G / np.linalg.norm(G, axis=0)


def gap(Pa, Pb):
    sv = np.linalg.svd(unit(Pa).T @ unit(Pb), compute_uv=False)
    return float(np.degrees(np.arccos(np.clip(sv.min(), -1.0, 1.0))))


def vs_th(G, r):
    sv = np.linalg.svd(Uth[:, :r].T @ unit(G), compute_uv=False)
    return np.degrees(np.arccos(np.clip(sv, -1.0, 1.0)))


COARSE = np.linspace(-6.0, 10.0, 33)


def D_of(cur, lg):
    try:
        I, D = cur(lg)
    except Exception:
        return np.nan
    return D if np.isfinite(D) else np.nan


def solve_gain(cur, t):
    """粗网格找**下降穿越**再二分（小增益端 DARE 会 `Failed to find a finite solution`，
    所以不能用固定窗端 $[-14,10]$ $-$ $-$ e171 第一版死在 $\lg=-14$，本轮改自适应找括号）。
    返回 $(\lg,I,D,\min D_{\rm窗})$；窗内无穿点返回 None。"""
    vals = np.array([D_of(cur, lg) for lg in COARSE])
    dmin = float(np.nanmin(vals)) if np.isfinite(vals).any() else np.inf
    for i in range(len(COARSE) - 1):
        if np.isfinite(vals[i]) and np.isfinite(vals[i + 1]) and vals[i] > t > vals[i + 1]:
            lo, hi = COARSE[i], COARSE[i + 1]
            for _ in range(50):
                mid = 0.5 * (lo + hi)
                Dm = D_of(cur, mid)
                if not np.isfinite(Dm):
                    return None
                if Dm > t:
                    lo = mid
                else:
                    hi = mid
                if hi - lo < 1e-12:
                    break
            lg = 0.5 * (lo + hi)
            I, D = cur(lg)
            return (lg, I, D, dmin) if abs(D - t) <= 1e-9 else None
    return None


p(r'== 锚点 $J_C=%.6f$；$\rho(A)=%.4f$；$\lg$ 粗窗 $[%g,%g]$（33 点）；$I_{\rm unc}$（e163 凸 SDP）：%s =='
  % (JC, max(np.abs(np.linalg.eigvals(A))), COARSE[0], COARSE[-1], SDP))
p(r'   $\Theta$ 特征值（降序）：%s' % np.array2string(wth[::-1], precision=4))
p('')
p(r'== [S1] 单调性侧检：20 条随机秩-1 方向，13 点子网格上 $D(\lg)$ 的上升段计数 ==')
sub = np.linspace(-6.0, 10.0, 13)
rup = 0
rcnt = 0
rbad = 0
for _ in range(20):
    u = rng.normal(0, 1, (n, 1))
    Ds = np.array([D_of(cur_of(unit(u), 1), lg) for lg in sub])
    fin = np.isfinite(Ds)
    rbad += int(np.sum(~fin))
    rup += int(np.sum(np.diff(Ds[fin]) > 0))
    rcnt += max(0, int(np.sum(fin)) - 1)
p(r'   上升段合计 %d（可用步数 %d，非有限 %d）$\\Rightarrow$ 单调降假设 %s'
  % (rup, rcnt, rbad, '成立' if rup == 0 else '**失效**'))
p('')
p(r'== [S0/S2/S3/S4] 逐格等值带 ==')

for t, r in ((80.00, 1), (56.66, 1), (34.41, 2), (34.41, 3)):
    nd = 800 if r == 1 else (400 if r == 2 else 250)
    cands = [Uth[:, :r]] + [np.linalg.qr(rng.normal(0, 1, (n, r)))[0] for _ in range(nd)]
    pts = []
    flrs = []
    nun = 0
    nfail = 0
    for G in cands:
        Zn = unit(G)
        try:
            got = solve_gain(cur_of(Zn, r), t)
        except Exception:
            nfail += 1
            continue
        if got is None:
            nun += 1
            continue
        lg, I, D, dmin = got
        pts.append((I, Zn))
        flrs.append(dmin)
    p(r'  --- $D=%.2f$，$r=%d$：扫 %d 条 $-$ $-$ 有效 %d 条，窗内无穿点 %d 条，求解失败 %d 条；'
      r'同秩可达地板（增益饱和端 $\min D$）$=%.4f$ ---'
      % (t, r, len(cands), len(pts), nun, nfail, min(flrs) if flrs else np.nan))
    if not pts:
        p(r'      [S3] 判决：0 条有效读数 $\\Rightarrow$ 不判')
        p('')
        continue
    pts.sort(key=lambda z: z[0])
    Imin = pts[0][0]
    p(r'      $\min I=%.6f$（$-I_{\rm unc}=%+.2e$ bit）；[S4] 最优方向与 $\Theta$ 前 %d 平面主角度 %s'
      % (Imin, Imin - SDP[t], r, np.array2string(vs_th(pts[0][1], r), precision=2)))
    for e in EPS:
        band = [q for q in pts if q[0] <= Imin + e]
        mx = max((gap(a[1], b[1]) for i, a in enumerate(band) for b in band[i + 1:]), default=0.0)
        p(r'      [S2] $\epsilon=%g$ bit：带内 %d/%d 条，带内最大两两角距 %.2f$^\circ$'
          % (e, len(band), len(pts), mx))
    band = [q for q in pts if q[0] <= Imin + EPS[0]]
    if len(band) == 1:
        p(r'      [S3] 判决：**本扫功率不足，不判**（$\epsilon=10^{-6}$ 带内只有 1 条）')
    else:
        sel = []
        for q in band:
            if all(gap(q[1], s[1]) > 5.0 for s in sel):
                sel.append(q)
            if len(sel) >= 8:
                break
        dmax = max(gap(a[1], b[1]) for i, a in enumerate(sel) for b in sel[i + 1:]) if len(sel) >= 2 else 0.0
        dv = max(a[0] for a in band) - min(a[0] for a in band)
        p(r'      [S3] 带内 $\ge5^\circ$ 抽样 %d 条：最大角距 %.2f$^\circ$，带内值跨度 %.2e bit'
          % (len(sel), dmax, dv))
        if dmax >= 5.0:
            p(r'      [S3] 判决：**近最优秩-$r$ 支撑在 $10^{-6}$ bit 内不唯一（数值级）**')
        elif dmax < 1.0:
            p(r'      [S3] 判决：**未见不唯一**')
        else:
            p(r'      [S3] 判决：**未决**（带内最大角距 %.2f$^\circ$ $<5^\circ$）' % dmax)
    p('')

p(r'== [S5] 与 e162b/e165 的对账（同一格的最优角，三种求法）==')
p(r'   e162b [B]：$D{=}80,r{=}1$ 离开角 $6.02^\circ$（割线 $+$ 携带 GUESS，方向已归一）')
p(r'   e165 [T2]：同格 $62.5^\circ$（Powell 多起点，`pa_vs_th` **未归一化** $\\Rightarrow$ 角偏小方向塌缩，此处却读大 $\\Rightarrow$ 需本轮定标）')
p(r'   e167：细化位移 $0.00^\circ$ 是**未归一化 clip 伪影**，已作废')
p(r'   本轮 [S4] 的数在上面的表里，且是**宽窗二分 $+$ 归一化**的唯一口径 $-\ -$ 以后引用它。')
p('')
p('用时 %.0f s | DARE %d | 异常 %s' % (time.time() - T0, P['NDB'][0], P['EXC'] or '无'))
txt = '\n'.join(out) + '\n'
open('p0/e171_out.txt', 'w', encoding='utf-8').write(txt)
open('p0/e171_run.txt', 'w', encoding='utf-8').write(txt)
