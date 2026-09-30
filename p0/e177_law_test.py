# -*- coding: utf-8 -*-
r"""e177：**在 $2/3$ 支的阈值邻域把 $\Delta_2$ 与 $\Lambda$ 的整个谱同时加密** $-\ -$ 审 C 的"定量律未建立"能不能被推进。

来源（逐字，板 4548 之后 78-C）：
  C："与我的最好 $\Delta$ 上界**同向**……$\Rightarrow$ **$\lambda_2/\lambda_1$ 是这条现象的驱动量**"
  C："**但不许写成律**：比值 $\Delta_{\min}/(\lambda_2/\lambda_1)$ 从 $0.32$ 一路涨到 $26$（不是常数，也不是单一幂律：两端局部斜率 $2.8$ 与 $7.5$）"
  C 可担保的那句："严格格的损失随 $S^\star$ 的次主特征值一起消失（11 点单调），定量律未建立"

我 e176 已经把 C 那一列在我机器上复算到 **11/11 格**、且**对象判明**（C 的 $S^\star$ = 我车道的 KKT 阵 $\Lambda$），
并给出 $2/3$ 支的**外推检验 (a)**：归零的是 $\lambda_3/\lambda_1$（降 545 倍），而 $\lambda_2/\lambda_1$ 只降 $1.2$ 倍。
本轮问的是 C 明确留着没答的那一半：**"律"到底存不存在**。$2/3$ 支是我的地盘（C 的 $\Delta$ 在 $r{=}1$ 那支有 53 倍种子分歧，我那支没有），
而且我有**两台独立种子**的读数（e175 seed 777、本脚本另一条流）可以先把**我自己的噪声底量出来**再谈斜率。

预注册判据（跑前写死，本文件先落盘再运行）：
 [Z0] **噪声底标定（本轮的关键，先把尺做出来再量东西）**：与 e175 的 5 个公共格（$33.00/33.50/34.00/34.20/34.41$）逐格比 $\Delta_2$，
      $N:=\max_i\{\max(q_i,1/q_i)\}$，$q_i=$ 本脚本/e175。$N$ 就是"我方 $\Delta$ 的种子噪声倍数"。
      判法：$N\le2$ $\\Rightarrow$ 本车道 $\Delta$ 足以分辨幂指数 $\pm0.3$；$2<N\le5$ $\\Rightarrow$ 只分辨"数量级"；$N>5$ $\\Rightarrow$ **我的 $\Delta$ 曲线也降级为散点**（与 C 78-A 同判，公开自我降级）。
 [Z1] 单调门：$\Delta_2$ 与 $\lambda_3/\lambda_1$ 沿 $D\uparrow$ 逐段不升（容差 $10\%$）$\\Rightarrow$ 支持 e176 的 (a)。
 [Z2] **模型门（穷举三个候选，判决行必须三行齐全）**：在 $\lambda_3>0$ 且 $\Delta_2>10^{-8}$ 的格上
      M1 $\Delta_2/\lambda_3$（绝对谱，线性）；M2 $\Delta_2/(\lambda_3/\lambda_1)$（比值，线性）；M3 $\Delta_2/(\lambda_3/\lambda_1)^2$（比值，平方）。
      每个模型算极差 $\mathrm{spr}=\max/\min$。以 $N$ 为分辨率：**$\mathrm{spr}\le 2N$ 记"该模型在噪声内是常数"**。
        (a) 恰一个模型过 $\\Rightarrow$ **律成立（数值级）**，写清是哪一个；
        (b) 多个过 $\\Rightarrow$ 判"噪声内不可分辨"，并印它们两两之比随 $D$ 的变化（哪个更平）；
        (c) 全不过 $\\Rightarrow$ 印逐段局部斜率 $p_i=\Delta\ln\Delta_2/\Delta\ln(\lambda_3/\lambda_1)$，
             若 $\max p_i-\min p_i\le 2\ln N$ $\\Rightarrow$ "单一幂律、指数 $p=$ 中位数 $\pm\ln N$"；否则**只报散点** $-\ -$ **C 的"律未建立"在我这支也成立**。
 [Z3] 条数门（#46）：每条判决行印有效格数；不可达/无解/非有限计入"退出"，不填分子分母。
口径：$\Delta_2$ 是**上界**（两台机器都已认领）；$\lambda$ 列无搜索、无二分。
      四元组 $(D,r,r^\star,\delta_{\rm sw})$ 逐格印；$\delta_{\rm sw}=|D-34.41|$。
"""
import sys
import re
import time
import numpy as np
import cvxpy as cp
sys.stdout.reconfigure(encoding='utf-8')

out = []
T0 = time.time()


def p(*a):
    s = ' '.join(str(x) for x in a)
    out.append(s)
    print(s)
    sys.stdout.flush()


def sym(X):
    return 0.5 * (X + X.T)


# ---- 复用 e175 的前缀（定价阵、cell/basis/angles、老 SDP 表），把它自己的打印吞掉 ----
P2 = {'__name__': 'p2'}
_src = open('p0/e175_same_cell_delta.py', encoding='utf-8').read().split('RES = {}')[0]
_src = _src.replace("sys.stdout.reconfigure(encoding='utf-8')", "pass")
_stdout = sys.stdout
_buf = __import__('io').StringIO()
_buf.reconfigure = lambda **kw: None
sys.stdout = _buf
exec(compile(_src, 'p0/e175_same_cell_delta.py[prefix]', 'exec'), P2)
sys.stdout = _stdout
A, B, W, TH, JC, n, ln2 = P2['A'], P2['B'], P2['W'], P2['TH'], P2['JC'], P2['n'], P2['ln2']
SDP0, RANK6 = P2['SDP'], P2['RANK6']
cell, angles, Uth, FLOOR = P2['cell'], P2['angles'], P2['Uth'], P2['FLOOR']


def sdp(D):
    Pv = cp.Variable((n, n), symmetric=True)
    Pi = cp.Variable((n, n), symmetric=True)
    PT = A @ Pv @ A.T + W
    lmi = cp.bmat([[Pv - Pi, Pv @ A.T], [A @ Pv, PT]]) >> 0
    prob = cp.Problem(cp.Minimize(0.5 * (np.log(np.linalg.det(W)) - cp.log_det(Pi)) / ln2),
                      [Pv >> 0, Pi >> 0, cp.trace(TH @ Pv) + JC <= D, Pv << PT, lmi])
    prob.solve(solver=cp.CLARABEL, tol_gap_abs=1e-10, tol_gap_rel=1e-10, tol_feas=1e-10)
    if prob.value is None or Pv.value is None:
        return None
    v = sym(Pv.value)
    Lam = sym(np.linalg.inv(v) - np.linalg.inv(A @ v @ A.T + W))
    return float(prob.value), np.linalg.eigvalsh(Lam)[::-1]


DSW = [33.00, 33.25, 33.50, 33.75, 34.00, 34.10, 34.20, 34.30, 34.41]
E175 = open('p0/e175_out.txt', encoding='utf-8').read()
MYD2 = {}
for m in re.finditer(r'^\s*(\d+\.\d+)\s+(\d)\s+(\d)\s+([\d.]+)\s+([+-][\d.eE+-]+)\s+', E175, re.M):
    if int(m.group(3)) == 2:
        MYD2[float(m.group(1))] = float(m.group(5))
assert len(MYD2) == 5, ('e175 的秩-2 行数应为 5，实际 %d' % len(MYD2))
p('== e175 的 $\\Delta_2$（解析、不转写）：%s ==' % {k: '%.3e' % v for k, v in MYD2.items()})

SDP = {}
LAM = {}
RES = {}
COV = {}
p('')
p('== 逐格：凸侧 $I_{\\rm unc}/\\Lambda$ 谱（无搜索）$+$ 非凸侧 $\\Delta_2$（25 起点，seed 20260930，与 e175 的 777 不同流）==')
p('%7s %6s %7s %10s %10s %10s %10s %11s %11s %7s %s'
  % ('D', '$\\delta_{sw}$', '$r^\\star$', '$\\lambda_1$', '$\\lambda_2$', '$\\lambda_3$',
     '$\\lambda_3/\\lambda_1$', '$\\Delta_2$', '本/e175', 'k/25', '残差'))
rng = np.random.default_rng(20260930)
for t in DSW:
    r_ = sdp(t)
    if r_ is None:
        p('%7.2f  SDP 无解 $\\Rightarrow$ 退出' % t)
        continue
    SDP[t], ev = r_
    LAM[t] = ev
    r6 = int(np.sum(ev / ev[0] > 1e-6))
    r = r6 - 1
    if t <= FLOOR.get(r, 0.0) or r < 1:
        p('%7.2f %6.2f %7d  不可达/无严格格 $\\Rightarrow$ 记 $+\\infty$' % (t, abs(t - 34.41), r6))
        RES[t] = np.inf
        continue
    starts = [Uth[:, :r]] + [np.linalg.qr(rng.normal(0, 1, (n, r)))[0] for _ in range(24)]
    raw, allv = cell(t, r, starts)
    if not allv:
        p('%7.2f %6.2f %7d  有效 0 条 $\\Rightarrow$ 无读数' % (t, abs(t - 34.41), r6))
        RES[t] = np.nan
        continue
    I2, rb, G2 = allv[0]
    RES[t] = I2 - SDP[t]
    COV[t] = len(raw)
    q = (RES[t] / MYD2[t]) if t in MYD2 else np.nan
    p('%7.2f %6.2f %7d %10.4f %10.4f %10.4e %11.3e %11.3e %7s %5s %8.1e'
      % (t, abs(t - 34.41), r6, ev[0], ev[1], max(ev[2], 0.0), max(ev[2], 0.0) / ev[0],
         RES[t], ('%.2f' % q) if np.isfinite(q) else '——', '%d/25' % len(raw), rb))

# ---- [Z0] 噪声底 ----
p('')
p('== [Z0] 种子噪声底：公共格 $q=$ 本脚本/e175 ==')
qs = [RES[t] / MYD2[t] for t in sorted(MYD2) if np.isfinite(RES.get(t, np.nan)) and RES.get(t, 0) > 1e-12]
for t in sorted(MYD2):
    d = RES.get(t, np.nan)
    p('  $D=%.2f$：e177 %s / e175 %.3e $= $ %s' % (t, ('%.3e' % d) if np.isfinite(d) else str(d), MYD2[t],
                                               ('%.3f' % (d / MYD2[t])) if np.isfinite(d) and d > 0 else '——'))
N = max([max(x, 1.0 / x) for x in qs]) if qs else np.inf
p('  有效公共格 %d 个 $\\Rightarrow$ $N=%.2f$（我方 $\\Delta$ 的种子噪声倍数）$-$ $-$ #46 条数' % (len(qs), N))
p('  [Z0] %s' % ('**$N\\le2$：可分辨幂指数 $\\pm0.3$**' if N <= 2 else
                ('**$2<N\\le5$：只分辨数量级**' if N <= 5 else
                 '**$N>5$：我的 $\\Delta_2$ 曲线降级为散点（与 C 78-A 同判）**')))

# ---- [Z1] 单调 ----
p('')
p('== [Z1] 单调门（沿 $D\\uparrow$ 逐段不升，容差 $10\\%$）==')
good = [t for t in DSW if np.isfinite(RES.get(t, np.nan)) and RES[t] > 1e-12]
badd = [(good[i], good[i + 1]) for i in range(len(good) - 1) if RES[good[i + 1]] > RES[good[i]] + 1e-6]
badl = []
for i in range(len(good) - 1):
    l3a = max(LAM[good[i]][2], 0.0) / LAM[good[i]][0]
    l3b = max(LAM[good[i + 1]][2], 0.0) / LAM[good[i + 1]][0]
    if l3b > l3a * 1.10:
        badl.append((good[i], good[i + 1]))
p('  $\\Delta_2$ 违反段 %s（有效 %d 格）；$\\lambda_3/\\lambda_1$ 违反段 %s' % (badd if badd else '无', len(good), badl if badl else '无'))
Z1 = not badd and not badl
p('  [Z1] %s $\\Rightarrow$ e176 的分支 (a)（第 $r^*$ 个特征值随阈值归零）%s'
  % ('**通过**' if Z1 else '**不通过**', '在 %d 个点上继续成立' % len(good) if Z1 else '**收回**'))

# ---- [Z2] 模型门：三个候选，判决行必须三行齐全 ----
p('')
p('== [Z2] 模型门（分辨率由 [Z0] 的 $N$ 给定：过线判据 $\\mathrm{spr}\\le 2N$）==')
p('   M1 $\\Delta_2/\\lambda_3$（绝对谱线性）  M2 $\\Delta_2/(\\lambda_3/\\lambda_1)$（比值线性）  M3 $\\Delta_2/(\\lambda_3/\\lambda_1)^2$（比值平方）')
USE = [t for t in good if max(LAM[t][2], 0.0) > 0 and RES[t] > 1e-8]
MODELS = {}
for nm, f in (('M1', lambda t: RES[t] / max(LAM[t][2], 1e-300)),
              ('M2', lambda t: RES[t] / (max(LAM[t][2], 0.0) / LAM[t][0])),
              ('M3', lambda t: RES[t] / (max(LAM[t][2], 0.0) / LAM[t][0]) ** 2)):
    vals = {t: f(t) for t in USE}
    MODELS[nm] = vals
    spr = max(vals.values()) / min(vals.values())
    p('  %s：参与格 %d $-$ $-$ 逐格 %s $\\Rightarrow$ spr $=%.2f$，$2N=%.2f$ $\\Rightarrow$ %s'
      % (nm, len(vals), ' '.join('%.2e' % vals[t] for t in sorted(vals)), spr, 2 * N,
         '**过线（噪声内常数）**' if spr <= 2 * N else '不过线'))
passed = [nm for nm in MODELS if max(MODELS[nm].values()) / min(MODELS[nm].values()) <= 2 * N]
if len(passed) == 1:
    p('  [Z2] 分支 **(a)** $\\Rightarrow$ **律成立（数值级）：$\\Delta_2\\propto$ %s 的自变量**，有效格 %d，$N=%.2f$'
      % (passed[0], len(USE), N))
elif len(passed) > 1:
    p('  [Z2] 分支 **(b)** $\\Rightarrow$ 噪声内不可分辨（%s 同时过线），有效格 %d $-$ $-$ 不升格为律'
      % ('/'.join(passed), len(USE)))
else:
    p('  [Z2] 分支 **(c)** $\\Rightarrow$ 三个模型都不过线（$N=%.2f$）$-$ $-$ 改查单一幂律的指数' % N)
    ps = []
    for i in range(len(USE) - 1):
        ta, tb = USE[i], USE[i + 1]
        la = max(LAM[ta][2], 0.0) / LAM[ta][0]
        lb = max(LAM[tb][2], 0.0) / LAM[tb][0]
        if la > 0 and lb > 0:
            ps.append((np.log(RES[tb]) - np.log(RES[ta])) / (np.log(lb) - np.log(la)))
    pv = list(ps)
    p('   逐段局部斜率 $p_i=\\Delta\\ln\\Delta_2/\\Delta\\ln(\\lambda_3/\\lambda_1)$：%s（段数 %d）'
      % (' '.join('%.2f' % v for v in pv), len(pv)))
    if pv and (max(pv) - min(pv)) <= 2 * np.log(N):
        p('  [Z2] (c-1) $\\Rightarrow$ **单一幂律，指数 $p=%.2f\\pm%.2f$**（极差 %.2f $\\le 2\\ln N=%.2f$），有效格 %d'
          % (float(np.median(pv)), np.log(N), max(pv) - min(pv), 2 * np.log(N), len(USE)))
    else:
        p('  [Z2] (c-2) $\\Rightarrow$ **只报散点：C 的"定量律未建立"在我这支同样成立**（斜率极差 %.2f $>2\\ln N=%.2f$），有效格 %d'
          % ((max(pv) - min(pv)) if pv else np.nan, 2 * np.log(N), len(USE)))

p('')
p('== [Z3] 条数与退出 ==')
p('  网格 %d 格；$\\Delta_2$ 有效 %d 格；不可达/无读数 %d 格；进 [Z2] 的格 %d'
  % (len(DSW), len(good), len(DSW) - len(good), len(USE)))
p('  逐格四元组 $(D,r,r^*,\\delta_{sw})$：%s' % ' '.join('(%.2f,%d,%d,%.2f)' % (t, RANK6.get(t, 2) - 1, RANK6.get(t, 2), abs(t - 34.41)) for t in good))
p('')
p('用时 %.0f s | DARE %d 次 | 有效读数 %d | 回退 %d | 异常 无'
  % (time.time() - T0, P2['NDB'][0] if 'NDB' in P2 else -1, sum(1 for t in DSW if np.isfinite(RES.get(t, np.nan))),
     P2['NFALLBACK'][0] if 'NFALLBACK' in P2 else -1))
open('p0/e177_out.txt', 'w', encoding='utf-8').write('\n'.join(out) + '\n')

