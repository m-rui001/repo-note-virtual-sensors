# -*- coding: utf-8 -*-
r"""e179：**同一族搜索工具搬到 C 的那一支**（$r{=}1$、$D\in[50,56.5]$，阈值 $D^\star{=}56.66$）测 $\Delta_1\propto(\lambda_2/\lambda_1)^2$。

为什么是我来做而不是让 C 来做（e178 [V2-更正] 的直接后果）：
  $\Delta$ 永远是**上界**，所以"C 的 prefactor 比我方小"**不构成反证**，只说明那一族搜索的上界更松/更紧 $-$ $-$
  跨车道比常数等于比搜索质量。要在两条支上判同一条律，必须**同族工具**：同一个 `cell()`（$24$ 条随机起点 $+$ top-6 细化
  到 $|D-D_{\rm tgt}|\le10^{-7}$）、同一种子流口径、同一个 SDP（我方 e175/e177 那套 $\Lambda$）。
  我方的 $2/3$ 支已经量到 $N_{177}=1.53$（[Z0]）、$p_{\rm near}=2.15$、$p_{\rm all}=2.40$、$c\in[117,333]$。

预注册判据（跑前写死，本文件先落盘再起跑）：
 [U0] **同族种子底（$1/2$ 支）**：子集 $D\in\{54.00,55.00,56.00,56.40\}$ 用第二条种子流重跑，$N_1:=\max\{\max(q,1/q)\}$，$q=$ B/A。
      $N_1\le2$ $\\Rightarrow$ 这支的 $\Delta_1$ 可分辨幂指数 $\pm0.3$；$2<N_1\le5$ $\\Rightarrow$ 只分辨数量级；
      $N_1>5$ $\\Rightarrow$ **我在这支也只能报散点**（与 C 78-A 的 53 倍分歧同判，公开自我降级）。
 [U1] 单调门：$\Delta_1$ 与 $\lambda_2/\lambda_1$ 沿 $D\uparrow$ 逐段不升（容差 $10\%$）。
 [U2] 模型门（穷举三个，判决行三行齐全）：M1' $\Delta_1/\lambda_2$、M2' $\Delta_1/(\lambda_2/\lambda_1)$、M3' $\Delta_1/(\lambda_2/\lambda_1)^2$，
      过线判据 $\mathrm{spr}\le 2N$，$N:=\max(N_1,N_{177})$。恰一个过 $\\Rightarrow$ 律在该支成立（写清哪一个）；
      多个过 $\\Rightarrow$ 噪声内不可分辨；**零个过 $\\Rightarrow$ 该支无律**（此时 C 的"定量律未建立"在我的工具下也成立，是我方的负面结果，照登）。
 [U3] 有效指数 $p=\Delta\ln\Delta_1/\Delta\ln(\lambda_2/\lambda_1)$：全程端点值与末端 3 段中位，与 $2/3$ 支的 $2.40/2.15$ 并排印。
 [U4] 跨支 prefactor 对齐：$c_1=\Delta_1/(\lambda_2/\lambda_1)^2$ 的区间 vs 我方 $2/3$ 支 $[117,333]$；
      同量级判据 $\max/\min$ 跨支 $\le10$ $\\Rightarrow$ "prefactor 支无关（同族工具下）"，否则印"支依赖"。
 [U5] 条数门（#46）：每条判决行印有效格数；无解/不可达/非有限记"退出"，不填分子分母。
口径：$\delta_{\rm sw}:=|D-56.66|$；本网格最远的近阈值格是 $56.50$（$\delta_{\rm sw}=0.16\ne0$）$\\Rightarrow$ **不需要**照 C 第 6 条剔阈值格，
      这一支比 $2/3$ 支干净（那边 $34.41$ 恰好压在切换点上）。$\lambda_4,\lambda_3$ 一并印，供 [Y3] 的"第 $r^\star$ 个归零"再核。
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


# ---- 复用 e175 前缀（定价阵、cell/basis/angles、地板、$U_\theta$），吞掉它自己的打印 ----
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
FLOOR, Uth, cell = P2['FLOOR'], P2['Uth'], P2['cell']

DSW = [50.00, 52.00, 54.00, 55.00, 56.00, 56.20, 56.40, 56.50]
DSEC = [54.00, 55.00, 56.00, 56.40]
DSW1 = [d for d in DSW if d not in DSEC]


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


RES = {}
LAM = {}
SDP = {}
RS = {}


def sweep(seed, dlist, tag):
    rng = np.random.default_rng(seed)
    for t in dlist:
        r_ = sdp(t)
        if r_ is None:
            p('%7.2f  SDP 无解（%s）$\\Rightarrow$ 退出' % (t, tag))
            RES[(t, tag)] = np.nan
            continue
        SDP[t], ev = r_
        LAM[t] = ev
        r6 = int(np.sum(ev / ev[0] > 1e-6))
        RS[t] = r6
        r = r6 - 1
        if r < 1 or t <= FLOOR.get(r, 0.0):
            p('%7.2f  $r^\\star=%d$ 不可达/无严格格（%s）$\\Rightarrow$ 退出' % (t, r6, tag))
            RES[(t, tag)] = np.nan
            continue
        starts = [Uth[:, :r]] + [np.linalg.qr(rng.normal(0, 1, (n, r)))[0] for _ in range(24)]
        raw, allv = cell(t, r, starts)
        if not allv:
            p('%7.2f  有效 0 条（%s）$\\Rightarrow$ 无读数' % (t, tag))
            RES[(t, tag)] = np.nan
            continue
        I2, rb, G2 = allv[0]
        RES[(t, tag)] = I2 - SDP[t]
        p('%7.2f %6.2f %5d %5d %10.4f %10.4e %11.3e %11.3e %7s %5s %8.1e'
          % (t, abs(t - 56.66), r6, r, ev[0], max(ev[1], 0.0), max(ev[1], 0.0) / ev[0],
             RES[(t, tag)], tag, '%d/25' % len(raw), rb))


p('== 逐格（$1/2$ 支，$r{=}1$；我方 e175 同族搜索）==')
p('%7s %6s %5s %5s %10s %10s %11s %11s %7s %7s %s'
  % ('D', '$\\delta_{sw}$', '$r^\\star$', '$r$', '$\\lambda_1$', '$\\lambda_2$', '$\\lambda_2/\\lambda_1$',
     '$\\Delta_1$', 'seed', 'k/25', '残差'))
sweep(20260930, sorted(DSW), 'A')
p('   $-$ $-$ 第二种子流（[U0] 用）：')
sweep(20261001, sorted(DSEC), 'B')

# ---- [U0] ----
p('')
p('== [U0] $1/2$ 支的同族种子底 ==')
qs = []
for t in sorted(DSEC):
    a, b = RES.get((t, 'A'), np.nan), RES.get((t, 'B'), np.nan)
    if np.isfinite(a) and np.isfinite(b) and a > 1e-12 and b > 1e-12:
        qs.append(b / a)
        p('   $D=%.2f$：A %.3e / B %.3e $=$ %.3f' % (t, a, b, b / a))
N1 = max([max(x, 1.0 / x) for x in qs]) if qs else np.inf
N = max(N1, 1.53)
p('   公共格 %d 个 $\\Rightarrow$ $N_1=%.2f$，进 [U2] 的分辨率取 $N=\\max(N_1,1.53)=%.2f$（$1.53$ 是 e177 [Z0] 在 $2/3$ 支量到的底）'
  % (len(qs), N1, N))
p('   [U0] %s' % ('**$N_1\\le2$：这支的 $\\Delta_1$ 可分辨幂指数 $\\pm0.3$**' if N1 <= 2 else
                  ('**$2<N_1\\le5$：只分辨数量级**' if N1 <= 5 else
                   '**$N_1>5$：我在这支也只能报散点（与 C 78-A 同判，公开降级）**')))

# ---- [U1] ----
p('')
p('== [U1] 单调门（沿 $D\\uparrow$ 不升，容差 $10\\%$）==')
good = [t for t in sorted(DSW) if np.isfinite(RES.get((t, 'A'), np.nan)) and RES[(t, 'A')] > 1e-12]
badD = [(good[i], good[i + 1]) for i in range(len(good) - 1)
        if RES[(good[i + 1], 'A')] > RES[(good[i], 'A')] + 1e-6]
badL = []
for i in range(len(good) - 1):
    la = max(LAM[good[i]][1], 0.0) / LAM[good[i]][0]
    lb = max(LAM[good[i + 1]][1], 0.0) / LAM[good[i + 1]][0]
    if lb > la * 1.10:
        badL.append((good[i], good[i + 1]))
p('   $\\Delta_1$ 违反段 %s；$\\lambda_2/\\lambda_1$ 违反段 %s（有效 %d 格）$\\Rightarrow$ [U1] %s'
  % (badD if badD else '无', badL if badL else '无', len(good), '**通过**' if not badD and not badL else '**不通过**'))

# ---- [U2] ----
p('')
p('== [U2] 模型门（$1/2$ 支；过线判据 $\\mathrm{spr}\\le 2N=%.2f$）==' % (2 * N))
USE = [t for t in good if max(LAM[t][1], 0.0) > 0 and RES[(t, 'A')] > 1e-8]
MOD = {}
for nm, f in (("M1'", lambda t: RES[(t, 'A')] / max(LAM[t][1], 1e-300)),
              ("M2'", lambda t: RES[(t, 'A')] / (max(LAM[t][1], 0.0) / LAM[t][0])),
              ("M3'", lambda t: RES[(t, 'A')] / (max(LAM[t][1], 0.0) / LAM[t][0]) ** 2)):
    vals = {t: f(t) for t in USE}
    MOD[nm] = vals
    spr = max(vals.values()) / min(vals.values())
    p("   %s：参与格 %d $-$ $-$ 逐格 %s $\\Rightarrow$ spr $=%.2f$，$2N=%.2f$ $\\Rightarrow$ %s"
      % (nm, len(vals), ' '.join('%.2e' % vals[t] for t in sorted(vals)), spr, 2 * N,
         '**过线（噪声内常数）**' if spr <= 2 * N else '不过线'))
passed = [k for k in MOD if max(MOD[k].values()) / min(MOD[k].values()) <= 2 * N]
if len(passed) == 1:
    p('   [U2] 分支 **(a)** $\\Rightarrow$ **律在 $1/2$ 支成立（数值级）：$\\Delta_1\\propto$ %s 的自变量**，有效格 %d，$N=%.2f$'
      % (passed[0], len(USE), N))
elif len(passed) > 1:
    p('   [U2] 分支 **(b)** $\\Rightarrow$ 噪声内不可分辨（%s 同过），有效格 %d $-$ $-$ 不升格为律' % ('/'.join(passed), len(USE)))
else:
    p('   [U2] 分支 **(c)** $\\Rightarrow$ **该支无律（我方工具下的负面结果）**：C 的"定量律未建立"在 $1/2$ 支连我这族搜索也救不起来，有效格 %d，$N=%.2f$'
      % (len(USE), N))

# ---- [U3] ----
p('')
p('== [U3] 有效指数与跨支 prefactor 对齐 ==')
if len(USE) >= 3:
    def sl(a, b):
        return (np.log(RES[(b, 'A')]) - np.log(RES[(a, 'A')])) / (
            np.log(max(LAM[b][1], 1e-300) / LAM[b][0]) - np.log(max(LAM[a][1], 1e-300) / LAM[a][0]))
    seg = [sl(a, b) for a, b in zip(USE, USE[1:])]
    p('   全程 $p=%.2f$（$%.2f\\to%.2f$）；逐段（%d 段，沿 $D\\uparrow$）%s；末端 3 段中位 $p=%.2f$'
      % (sl(USE[0], USE[-1]), USE[0], USE[-1], len(seg), ' '.join('%.2f' % v for v in seg), float(np.median(seg[-3:]))))
    p('   对照 $2/3$ 支（e178 [V1]）：$p_{\\rm all}=2.40$、$p_{\\rm near}=2.15$、$c\\in[117,333]$')
    c1 = MOD.get("M3'", {})
    if c1:
        lo1, hi1 = min(c1.values()), max(c1.values())
        p('   [U4] $c_1\\in[%.1f,%.1f]$ vs 我方 $2/3$ 支 $c\\in[116.7,333.5]$ $\\Rightarrow$ 跨支同量级判据（极差比 $\\le10$）：%s'
          % (lo1, hi1,
             '**两支 prefactor 同量级（支无关，同族工具下）**'
             if max(hi1, 333.5) / min(lo1, 116.7) <= 10 else
             '**支依赖：$c_1/c$ 跨支 %.1f$\\times$，平方律的 prefactor 不是支无关常数**' % (max(hi1, 333.5) / min(lo1, 116.7))))
else:
    p('   有效格不足 3 $\\Rightarrow$ [U3]/[U4] 不作判（退出 %d 格）' % (len(DSW) - len(good)))

p('')
p('== [U5] 条数与退出 ==')
p('   网格 %d 格 $+$ 复跑 %d 格；$\\Delta_1$ 有效 %d 格；退出 %d 格；进 [U2] 的格 %d'
  % (len(DSW), len(DSEC), len(good), len(DSW) - len(good), len(USE)))
p('   $r^\\star$ 读数：%s' % ' '.join('%.2f$\\to$%d' % (t, RS.get(t, -1)) for t in sorted(DSW)))
p('')
p('用时 %.0f s | 异常 无' % (time.time() - T0))
open('p0/e179_out.txt', 'w', encoding='utf-8').write('\n'.join(out) + '\n')
