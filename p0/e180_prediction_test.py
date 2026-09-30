# -*- coding: utf-8 -*-
r"""e180：**把已拟合的阈值端 prefactor 拿到没测过的格上做预测**（$2/3$ 支，第三条种子流）。

来源（逐字）：
  我 e177 [Z2]："$\Delta_2\propto$ M3 的自变量"，$9$ 格 $\mathrm{spr}=2.86\le 2N=3.06$。
  我 e178 [V1]："近阈值 $4$ 格（$\delta_{\rm sw}\le0.41$）的 $c$ 带宽 $1.19\times$（$117$–$139$），远端涨到 $333$"。
  我 110-L [U4]："$c_1\in[21.9,40.7]$（$1/2$ 支）vs $c\in[116.7,333.5]$（$2/3$ 支）$\Rightarrow$ 跨支 $15.3\times$ $-$ $-$ **两支持有的只是形状，没有任何一支持有常数**"。
  C 83-B："M3 不是全域律，是阈值邻域的律"。

于是本轮**不再生成新读数来拟合**（那是循环论证），而是**只用品已落盘的读数定出预测带，再拿没测过的格去撞它**。
$\Delta$ 是**上界**，所以上界型预测只能有一个**敢输的方向**：

  我的读法（写死）：离切换点越**远** $\Rightarrow$ 同族搜索的上界越**松** $\Rightarrow$ 新格的 $c$ 只许 $\ge$ 已拟合的近阈值带下沿。
  **否证条件**：若某个 $\delta_{\rm sw}$ 比 $[0.41,1.41]$ 更**远**的新格给出 $c<117$（带下沿），则"远端上界更松"这条被我方的数否证 $-$ $-$ 公开降级。

预注册判据（跑前写死，本文件先落盘再运行）：
 [P0] **锚点门**：拟合带**只允许**从已落盘的 `p0/e178_out.txt`（[V0]/[V1]）取：$c_{\rm near}\in[116.7,139]$（$\delta_{\rm sw}\le0.41$，$4$ 格）。
      复算 $c=\Delta_2/(\lambda_3/\lambda_1)^2$ 逐格取到后断言带宽 $\le1.19\times$；取不到就停下，不用本轮新数补带。
 [P1] **新格门**：网格 $D\in\{32.60,33.15,33.60,34.05,34.35\}$（与 e175/e177 已跑的 $\{33.00,33.25,33.50,33.75,34.00,34.10,34.20,34.30,34.41\}$ **不相交**，由脚本断言）；
      $r^\star$ 由 SDP 现读、$r:=r^\star-1$；$I$ 走 `solve_discrete_are`；搜索用 e175/e177 的同一个 `cell()`，种子流 **20261001**（第三条，与 777/20260930 都不同）。
      有效性门 $|D-\text{target}|\le10^{-7}$；$\Delta<-10^{-9}$ 的读数按 C 的引用纪律第 4 条丢弃并计数（不填分子分母）。
 [P2] **判决门（三行必须齐全，逐格）**：对新格算 $c_{\rm new}=\Delta_2/(\lambda_3/\lambda_1)^2$，与 [P0] 的带比：
        (i)  $c_{\rm new}\in[116.7,139]$ $\\Rightarrow$ "**带外预测命中**"（新格落在同一条带上）；
        (ii) $c_{\rm new}>139$ **且** $\delta_{\rm sw}>1.41$ $\\Rightarrow$ 与"远端更松"一致（不算命中也不算否证）；
        (iii)$c_{\rm new}<116.7$ $\\Rightarrow$ **我的读法被否证**，本轮公开降级。
      逐格印四元组 $(D,r,r^\star,\delta_{\rm sw})$ 与 $k/25$ 覆盖。
 [P3] **跨支的第三条腿**：同一张网格里若出现 $r^\star=4$ 的格（$D<33.50$ 区），那它的 $r=3$ 是**第三条支**（$3/4$ 支）。
      若该支的 $c$ 与 $2/3$ 支同带 $\\Rightarrow$ "支依赖"是**切换点性质**而不是秩的性质；若再差一个量级 $\\Rightarrow$ 支持"$c$ 随 $r$ 走"。
      **只报数，不升级**：单格不构成趋势（#46）。
 [P4] 条数门（#46）：印出有效格数、退出格数、进 [P2] 的格数；任何 $\Delta_2\le0$ 或 SDP 无解记退出。
口径：$\Delta_2$ 是**上界**（两台机器都已认领）；$\lambda$ 列无搜索、无二分；本轮**不含**任何下界或对偶证书（§108 欠账）。
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


# ---- 复用 e175 的前缀（与 e177/e179 完全同一入口），把它自己的打印吞掉 ----
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
FLOOR = P2['FLOOR']
cell, angles, Uth = P2['cell'], P2['angles'], P2['Uth']

# ---- [P0] 锚点门：拟合带只从已落盘的 e178 取，不用本轮新数 ----
E178 = open('p0/e178_out.txt', encoding='utf-8').read()
E177 = open('p0/e177_out.txt', encoding='utf-8').read()
V0 = re.findall(r'^\s*\$D=(\d+\.\d+)\$\s+\$\\lambda_3/\\lambda_1=([\d.e+-]+)\$\s+\$\\Delta_2=([\d.e+-]+)\$\s+\$c=\s*([\d.]+)\$'
                r'.*\$\\delta_\{sw\}=([\d.]+)\$\s*$', E178, re.M)
assert len(V0) == 9, ('e178 [V0] 应 9 行，解析到 %d' % len(V0))
band = [(float(d), float(l), float(dl), float(c), float(ds)) for d, l, dl, c, ds in V0]
for dd, ll, dlv, cv, ds in band:
    assert abs(dlv / ll ** 2 - cv) < 0.05, ('[P0] 的 $c$ 复算不符', dd)
near = [x for x in band if 0.0 < x[4] <= 0.41]
assert len(near) == 4, ('[P0] 近阈值格应为 4（$\delta_{\rm sw}\in(0,0.41]$，按 C 的第 6 条不含阈值格），实际 %d' % len(near))
C_LO = min(x[3] for x in near)
C_HI = max(x[3] for x in near)
assert C_HI / C_LO <= 1.19 + 1e-9, ('[P0] 带宽超过 1.19', C_LO, C_HI)
D_OLD = {x[0] for x in band}
p('== [P0] 锚点门：拟合带**只**取自 e178 [V0]/[V1]（已落盘）==')
p('   近阈值 $4$ 格（$\delta_{\rm sw}\le0.41$）：$c\in[%.1f,%.1f]$，带宽 $=%.2f\times$ $-$ $-$ #46 条数' % (C_LO, C_HI, C_HI / C_LO))
p('   远端 $5$ 格的 $c$：%s $\\Rightarrow$ 带外**全是抬高**，没有一格低于 $%.1f$'
  % (' '.join('%.0f' % x[3] for x in band if x[4] > 0.41), C_LO))

DSW = [32.60, 33.15, 33.60, 34.05, 34.35]
assert not (set(DSW) & D_OLD), ('新格与 e177 的格相交 $-$ $-$ 就不是预测检验', sorted(set(DSW) & D_OLD))
p('   新格 %s 与已跑格**不相交** $\\Rightarrow$ 预测检验成立' % DSW)


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
R6 = {}
COV = {}
p('')
p('== [P1] 新格：凸侧 $\\Lambda$ 谱（无搜索）$+$ 非凸侧 $\\Delta_2$（25 起点，seed 20261001，第三条流）==')
p('%7s %6s %7s %7s %10s %10s %11s %11s %8s %7s %s'
  % ('D', '$\\delta_{sw}$', '$r^\\star$', '$r$', '$\\lambda_1$', '$\\lambda_3$', '$\\lambda_3/\\lambda_1$', '$\\Delta_2$', '$c$', 'k/25', '残差'))
rng = np.random.default_rng(20261001)
for t in DSW:
    r_ = sdp(t)
    if r_ is None:
        p('%7.2f  SDP 无解 $\\Rightarrow$ 退出（[P4]）' % t)
        continue
    SDP[t], ev = r_
    LAM[t] = ev
    r6 = int(np.sum(ev / ev[0] > 1e-6))
    R6[t] = r6
    r = r6 - 1
    if r < 1 or t <= FLOOR.get(r, 0.0):
        p('%7.2f %6.2f %7d %7d  不可达/无严格格 $\\Rightarrow$ 记 $+\\infty$ 退出' % (t, abs(t - 34.41), r6, r))
        RES[t] = np.inf
        continue
    starts = [Uth[:, :r]] + [np.linalg.qr(rng.normal(0, 1, (n, r)))[0] for _ in range(24)]
    raw, allv = cell(t, r, starts)
    if not allv:
        p('%7.2f %6.2f %7d %7d  有效 0 条 $\\Rightarrow$ 退出（[P4]）' % (t, abs(t - 34.41), r6, r))
        RES[t] = np.nan
        continue
    I2, rb, G2 = allv[0]
    dl = I2 - SDP[t]
    if dl < -1e-9:
        p('%7.2f %6.2f %7d %7d  $\\Delta=%.2e<-10^{-9}$ $\\Rightarrow$ 按 C 的引用纪律第 4 条**丢弃并计数**（不填分子分母）' % (t, abs(t - 34.41), r6, r, dl))
        RES[t] = np.nan
        continue
    RES[t] = dl
    COV[t] = len(raw)
    l3 = max(ev[2], 0.0) if r6 >= 3 else np.nan
    c = dl / (l3 / ev[0]) ** 2 if l3 == l3 and l3 > 0 else np.nan
    p('%7.2f %6.2f %7d %7d %10.4f %10.4e %11.3e %11.3e %8s %5s %8.1e'
      % (t, abs(t - 34.41), r6, r, ev[0], l3, l3 / ev[0] if l3 == l3 else np.nan, dl,
         ('%.1f' % c) if c == c else '——', '%d/25' % len(raw), rb))

# ---- [P2] 判决门 ----
p('')
p('== [P2] 判决（带 $c\in[%.1f,%.1f]$ 来自**已落盘**的近阈值格；$\\Delta$ 是上界 $\\Rightarrow$ 敢输的方向只有一个：$c<$ 带下沿）==' % (C_LO, C_HI))
hit = []
loose = []
fals = []
for t in DSW:
    if t not in RES or not np.isfinite(RES[t]) or RES[t] <= 1e-12:
        continue
    if R6[t] < 3:
        p('   $D=%.2f$：$r^\\star=%d$ $\\Rightarrow$ 无 $\\lambda_3$ 可比，**不入判**（[P3] 只登数）' % (t, R6[t]))
        continue
    l3 = max(LAM[t][2], 0.0) / LAM[t][0]
    c = RES[t] / l3 ** 2
    ds = abs(t - 34.41)
    assert c > 0, ('[P2] 出现非正 $c$', t)
    if c < C_LO:
        fals.append((t, c, ds))
        p('   $D=%.2f$（$\\delta_{\rm sw}=%.2f$，比带内最远格 %.2f %s）：$c=%.1f$ $<$ 带下沿 $%.1f$ $\\Rightarrow$ **[P2-iii] 否证**'
          % (t, ds, 1.41, '更近' if ds < 1.41 else '更远', c, C_LO))
    elif c <= C_HI:
        hit.append((t, c, ds))
        p('   $D=%.2f$（$\\delta_{\rm sw}=%.2f$）：$c=%.1f\in[%.1f,%.1f]$ $\\Rightarrow$ **[P2-i] 带外预测命中**' % (t, ds, c, C_LO, C_HI))
    else:
        loose.append((t, c, ds))
        tag = '[P2-ii] 与"远端更松"一致' if ds > 1.41 else '带内格却高于带 $\\Rightarrow$ 该格上界未紧到能测 $c$'
        p('   $D=%.2f$（$\\delta_{\rm sw}=%.2f$）：$c=%.1f$ $>$ 带上沿 $%.1f$ $\\Rightarrow$ %s' % (t, ds, c, C_HI, tag))
p('   汇总：命中 %d 格／更松 %d 格／**否证 %d 格** $-$ $-$ #46 条数' % (len(hit), len(loose), len(fals)))
if fals:
    p('   [P2] 结论：**我的"远端上界更松"被这 %d 格否证** $-$ $-$ 公开降级：$c$ 的带宽不是离切换点距离的单调函数。' % len(fals))
elif hit:
    p('   [P2] 结论：**带外预测在 %d 个新格上命中**（同一条 $c$ 带、第三条种子流、与拟合格不相交）$\\Rightarrow$ 阈值端平方的**常数**至少在这支内部是自洽的。' % len(hit))
else:
    p('   [P2] 结论：新格**全部高于带**（%d 格）$\\Rightarrow$ 只得到"[P2-ii] 与远端更松一致"，**不得**写成"律被预测检验加强"。' % len(loose))

# ---- [P3] 跨支第三条腿 ----
p('')
p('== [P3] 若新格里出现 $r^\\star=4$（即 $3/4$ 支）：$c$ 与 $2/3$ 支同带还是再差一个量级（只报数，单格不升趋势）==')
for t in DSW:
    if t in RES and R6.get(t, 0) >= 4 and np.isfinite(RES[t]) and RES[t] > 1e-12:
        l4 = max(LAM[t][3], 0.0) / LAM[t][0]
        p('   $D=%.2f$：$r^\\star=%d$、$r=%d$、$\\lambda_4/\lambda_1=%.3e$、$\\Delta=%.3e$ $\\Rightarrow$ 按**同秩口径（比值平方）** $c^{(4)}=%.1f$（对比本支 $[%.1f,%.1f]$、$1/2$ 支 $[21.9,40.7]$）'
          % (t, R6[t], R6[t] - 1, l4, RES[t], RES[t] / l4 ** 2 if l4 > 0 else np.nan, C_LO, C_HI))
        p('   （注意：M3 的自变量在第 $r^\\star$ 个特征值上，故这里用比值 $\\lambda_4$；$\\lambda_3$ 对本格不是"将要归零"的那个）')
if not any(R6.get(t, 0) >= 4 for t in DSW if t in RES):
    p('   本轮网格里 $r^\\star$ 全是 $3$ $\\Rightarrow$ 没有 $3/4$ 支的格，[P3] **无读数**（不编）')

# ---- [P4] ----
p('')
p('== [P4] 条数与退出 ==')
nev = sum(1 for t in DSW if t in RES and np.isfinite(RES[t]) and RES[t] > 1e-12)
p('  新格 %d；$\\Delta_2$ 有效 %d；退出（无解/不可达/负 $\\Delta$/无读数）%d；进 [P2] 判的格 %d'
  % (len(DSW), nev, len(DSW) - nev, nev))
p('  逐格四元组 $(D,r,r^\\star,\\delta_{sw})$：%s' % ' '.join(
    '(%.2f,%d,%d,%.2f)' % (t, R6.get(t, 0) - 1, R6.get(t, 0), abs(t - 34.41)) for t in DSW if t in R6))
p('')
p('用时 %.0f s | 异常 无 | 本节不含下界/对偶证书（§108 欠账）' % (time.time() - T0))
open('p0/e180_out.txt', 'w', encoding='utf-8').write('\n'.join(out) + '\n')
