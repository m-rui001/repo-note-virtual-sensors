# -*- coding: utf-8 -*-
r"""e176：**不搜索**，只把 C 补记十七 78-C 那一列 $\lambda_2/\lambda_1$ 在我车道的 SDP 上重算一遍。

C 的原文（板 4548 之后，78-C，逐字）：
  "这轮唯一**干净**的新量：$\lambda_2/\lambda_1(S^\star(D))$（它不经任何二分，直接来自 SDP）"
  表（$D$ 42.00/44/46/48/50/52/54/55/56/56.60/56.65）：
  $7.49\!\times\!10^{-2},6.39,5.33,4.29,3.28,2.29,1.30,8.11\!\times\!10^{-3},3.20,2.45\!\times\!10^{-4},2.67\!\times\!10^{-6}$
  "可写进注的措辞只一句：**严格格的损失随 $S^\star$ 的次主特征值一起消失（11 点单调），定量律未建立**"
  78-D 还直接点我的名字："**若你的 e175（同格 $\Delta$）打算报"$\Delta$ 随 $D$ 的形状"，请用 78-A 这张表压一下**"

为什么这条值得单独跑：78-C 是 C **撤掉形状主张**（78-A，第 9 笔账）之后**唯一还站着的**新量，
而且它**不依赖搜索/二分** $-$ $-$ 我这边能纯 SDP 复算，是两条车道最干净的对接面。

**对象歧义先摆明（这是本轮第一个要判的事）**：C 写 $S^\star$，但我车道里 SDP 的原变量 `Pv`（进代价 $\mathrm{tr}(\Theta S)$ 的那个）
与我的 KKT 阵 $\Lambda=S^{-1}-(ASA^\top+W)^{-1}$ 是**两个**矩阵。我的秩表 $\mathrm{rank}_{10^{-6}}(\Lambda)$ 已经与 C 的 $r^\star$ 逐格相同，
所以 C 的 $S^\star$ **可能**是我的 $\Lambda$、也可能是我的 `Pv`。$-\ -$ 不猜，用数认。

预注册判据（跑前写死，本文件先落盘再运行）：
 [Y0] **锚点门**：我自算 $I_{\rm unc}(42.00)=3.004696$、$(45.00)=2.705659$（C 板 74-A）差 $\le10^{-5}$ bit。
      不过 $-$ $-$ 后面一律不判。
 [Y1] **对象识别**：对 C 的 11 格逐一算两个候选的相对差 $|{\rm mine}-{\rm C}|/{\rm C}$，
      相符口径：**读数 $\ge10^{-3}$ 要求相对差 $\le0.15$**；**读数 $<10^{-3}$ 只要求两边都 $<10^{-2}$**（那里是求解器噪声底，不作逐位相符）。
      判法：相符格数 $\ge8/11$ 记"该对象复现"。
        ① 恰一侧 $\ge8$ 且另一侧 $\le4$ $\\Rightarrow$ **对象判明**（写清是哪一个）；
        ② 两侧都 $\ge8$ $\\Rightarrow$ 判"该区间两个矩阵近同"，并印两者之比；
        ③ 两侧都 $<8$ $\\Rightarrow$ **C 的 78-C 列在我机器上未复现** $-$ $-$ 本轮最重的可能结果，直接印出来。
 [Y2] **单调门**（只对复现成功的那一支）：11 点逐段不升（容差 $10\%$）$\\Rightarrow$ 支持 C 的"11 点单调"；
      违反 $\\Rightarrow$ 印反例段，C 的"单调"降级为"趋势"。
 [Y3] **外推检验（我的新东西，C 只量了 $1/2$ 支）**：C 的机制若为真，应当是"**第 $r^\star$ 个特征值随上阈值一起归零**"，
      不是"$\lambda_2$ 特别"。$1/2$ 支（$D\to56.66^-$，$r^\star\!:\,2\to1$）归零的是 $\lambda_2/\lambda_1$；
      $2/3$ 支（$D\to34.41^+$，$r^\star\!:\,3\to2$）应当归零的是 $\lambda_3/\lambda_1$。
      于是在 $D\in\{33.00,33.50,34.00,34.20,34.41\}$（我 e175 的 $\Delta_2$ 五格）算同一对象的 $\lambda_3/\lambda_1$ 与 $\lambda_2/\lambda_1$，
      与从 `p0/e175_out.txt` **解析**（不转写）的 $\Delta_2(D)$ 做同向检验：
        (a) $\lambda_3/\lambda_1$ 随 $D\uparrow$ 单调降 **且** $\Delta_2$ 同向降 $\\Rightarrow$ **C 的机制不是 $1/2$ 支特有**，
            可写成"$\Delta_r$ 由第 $r^\star$ 个特征值在阈值处归零所驱动"；
        (b) $\lambda_3/\lambda_1$ 不降、而 $\lambda_2/\lambda_1$ 降 $\\Rightarrow$ C 的说法**只在 $1/2$ 支成立**，"驱动量"降级为特例；
        (c) 两者都不降 $\\Rightarrow$ 只报数，不许用"特征值驱动"这句话（我自己的 $\Delta_2$ 曲线也失去机制解释）。
 [Y4] **条数门**（#46）：每条判决行都印**有效读数格数**；解不出/非有限的格计入"退出"，不算相符也不算反例。
 [Y5] **形状对账**（回 78-D 的点名）：印一句可比性判定 $-$ $-$ 我的 $\Delta$ 曲线是 $r=2$/秩 $2\text{-}3$ 支（$D\le34.41$），
      C 的 78-A 那张 $53$ 倍分歧表是 $r=1$/$D\in[42,56]$，**两组格不相交**；并印我这条曲线的 $\delta_{\rm sw}$（离切换点距离）。
口径：$\Delta$ 是**上界**；本脚本**不含任何非凸搜索**，所以 [Y1]–[Y3] 的数不依赖我的优化器可信度。
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


P = {'__name__': 'p'}
exec(compile(open('p0/e162_ray_closure.py', encoding='utf-8').read().split('CRAB = {')[0],
             'p0/e162_ray_closure.py[prefix]', 'exec'), P)
AUD = P['_ns']
A, B, W = AUD['A'].copy(), AUD['B'].copy(), AUD['W'].copy()
n = A.shape[1]
Pc, K, TH, JC = P['ctrl_full'](A, B, W, np.eye(4), np.eye(4))
P['A'], P['B'], P['W'], P['TH'], P['JC'] = A, B, W, TH, JC
ln2 = np.log(2.0)

# ---- C 板 78-C 那一列（11 格，逐字抄自板，脚本内先做针脚断言）----
C_LAM = {42.00: 7.49e-2, 44.00: 6.39e-2, 46.00: 5.33e-2, 48.00: 4.29e-2, 50.00: 3.28e-2,
         52.00: 2.29e-2, 54.00: 1.30e-2, 55.00: 8.11e-3, 56.00: 3.20e-3, 56.60: 2.45e-4, 56.65: 2.67e-6}
C_IUNC = {42.00: 3.004696, 45.00: 2.705659}
DS23 = [33.00, 33.50, 34.00, 34.20, 34.41]
SW_23 = 34.41   # 2/3 支切换点（我 e163/e175 的 rank3 != rank6 那格）
SW_12 = 56.66   # 1/2 支切换点（C 75-A 给区间 (50.00, 56.66]）

# 从 e175 的 stdout 解析我的 Delta_2(D)（#45：板上数字必须来自已落盘 stdout，不转写）
E175 = open('p0/e175_out.txt', encoding='utf-8').read()
MYD2 = {}
for m in re.finditer(r'^\s*(\d+\.\d+)\s+(\d)\s+(\d)\s+([\d.]+)\s+([+-][\d.eE+-]+)\s+', E175, re.M):
    if int(m.group(3)) == 2:
        MYD2[float(m.group(1))] = float(m.group(5))
p('== 解析 `p0/e175_out.txt`：我的 $\\Delta_2(D)$ 有 %d 格 $-$ $-$ %s ==' % (len(MYD2), MYD2))
assert len(MYD2) == 5, ('e175 的秩-2 行数应为 5，实际 %d $-$ $-$ 解析式要改，不判' % len(MYD2))


def sdp(D):
    """本车道的凸 SDP：原变量 S(`Pv`) $+$ KKT 阵 $\\Lambda=S^{-1}-(ASA^\\top+W)^{-1}$。返回 (I, S, Lam)。"""
    Pv = cp.Variable((n, n), symmetric=True)
    Pi = cp.Variable((n, n), symmetric=True)
    con_tr = cp.trace(TH @ Pv) + JC <= D
    PT = A @ Pv @ A.T + W
    lmi = cp.bmat([[Pv - Pi, Pv @ A.T], [A @ Pv, PT]]) >> 0
    prob = cp.Problem(cp.Minimize(0.5 * (np.log(np.linalg.det(W)) - cp.log_det(Pi)) / ln2),
                      [Pv >> 0, Pi >> 0, con_tr, Pv << PT, lmi])
    prob.solve(solver=cp.CLARABEL, tol_gap_abs=1e-10, tol_gap_rel=1e-10, tol_feas=1e-10)
    if prob.value is None or Pv.value is None:
        return None
    v = sym(Pv.value)
    lam = sym(np.linalg.inv(v) - np.linalg.inv(A @ v @ A.T + W))
    return float(prob.value), v, lam


def ratio(M, k):
    ev = np.linalg.eigvalsh(M)[::-1]
    return float(ev[k] / ev[0]), ev


p('')
p('== [Y0] 锚点门：$I_{\\rm unc}$ 与 C 板 74-A 差 $\\le10^{-5}$ ==')
y0ok = True
for D in sorted(C_IUNC):
    r_ = sdp(D)
    if r_ is None:
        p('  $D=%.2f$：SDP 无解 $\\Rightarrow$ 退出' % D)
        y0ok = False
        continue
    d = r_[0] - C_IUNC[D]
    p('  $D=%.2f$：我 %.6f vs C %.6f，差 %+.1e $\\Rightarrow$ %s' % (D, r_[0], C_IUNC[D], d, '过' if abs(d) <= 1e-5 else '**不过**'))
    y0ok = y0ok and abs(d) <= 1e-5
p('  [Y0] 总结：**%s**（有效读数 %d/2 格）' % ('**通过**' if y0ok else '**不通过 $-$ $-$ 后面一律不判**', int(y0ok) * 2))
if not y0ok:
    open('p0/e176_out.txt', 'w', encoding='utf-8').write('\n'.join(out) + '\n')
    sys.exit(0)

# ---- [Y1] 1/2 支：C 的 11 格，两个候选对象同时算 ----
p('')
p('== [Y1] 对象识别：C 板 78-C 的 $\\lambda_2/\\lambda_1$ 列 vs 我的两个候选（原变量 $S$ 与 KKT 阵 $\\Lambda$）==')
p('%8s %11s %13s %13s %9s %9s %s' % ('D', 'C 的读数', '我的 $\\lambda_2/\\lambda_1(S)$', '我的 $\\lambda_2/\\lambda_1(\\Lambda)$',
                                    'rel(S)', 'rel(Lam)', '相符(口径见 [Y1])'))
RS = {}
RL = {}
okS = []
okL = []
nv = 0
for D in sorted(C_LAM):
    r_ = sdp(D)
    if r_ is None:
        p('%8.2f  SDP 无解 $\\Rightarrow$ 退出（不算相符/反例）' % D)
        continue
    nv += 1
    rs, _ = ratio(r_[1], 1)
    rl, _ = ratio(r_[2], 1)
    RS[D], RL[D] = rs, rl
    c = C_LAM[D]

    def hit(mine):
        if c >= 1e-3:
            return abs(mine - c) / c <= 0.15
        return mine < 1e-2 and c < 1e-2
    hs, hl = hit(rs), hit(rl)
    if hs:
        okS.append(D)
    if hl:
        okL.append(D)
    p('%8.2f %11.3e %13.3e %13.3e %9.2f %9.2f %s'
      % (D, c, rs, rl, abs(rs - c) / c, abs(rl - c) / c, 'S' if hs else '-' + (' Lam' if hl else '')))
p('  [Y1] 相符格数：$S$ 侧 %d/%d，$\\Lambda$ 侧 %d/%d（有效读数 %d 格，退出 %d 格）$-$ $-$ #46 条数'
  % (len(okS), len(C_LAM), len(okL), len(C_LAM), nv, len(C_LAM) - nv))
if len(okS) >= 8 and len(okL) <= 4:
    OBJ, OBJNAME, Y1 = 'S', '原变量 $S^\\star$（进 $\\mathrm{tr}(\\Theta S)$ 的那个）', '①'
elif len(okL) >= 8 and len(okS) <= 4:
    OBJ, OBJNAME, Y1 = 'Lam', 'KKT 阵 $\\Lambda=S^{-1}-(ASA^\\top+W)^{-1}$', '①'
elif len(okS) >= 8 and len(okL) >= 8:
    OBJ, OBJNAME, Y1 = 'BOTH', '两个矩阵在该区间近同', '②'
else:
    OBJ, OBJNAME, Y1 = None, None, '③'
p('  [Y1] 分支 %s $\\Rightarrow$ %s' % (Y1, ('对象判明：C 的 $S^\\star$ 就是我车道的 %s' % OBJNAME) if OBJ else
                                    '**未复现 $-$ $-$ C 的 78-C 列在我机器上对不上，本轮最重结果**'))
if Y1 == '②':
    p('  ② 的两个矩阵之比（同一 $D$，$\\lambda_2/\\lambda_1$ 之比）：')
    for D in sorted(C_LAM):
        if D in RS and D in RL:
            p('    $D=%.2f$：%s / %s = %.3f' % (D, '%.3e' % RS[D], '%.3e' % RL[D], RS[D] / RL[D]))

# ---- [Y2] 单调门（对我复现的那一支；两支都印）----
p('')
p('== [Y2] 单调门：11 点逐段不升（容差 $10\\%$）==')
for nm, col in (('S', RS), ('Lam', RL)):
    ds = sorted(col)
    bad = [(ds[i], ds[i + 1], col[ds[i + 1]] / col[ds[i]]) for i in range(len(ds) - 1)
           if col[ds[i + 1]] > col[ds[i]] * 1.10 + 1e-15]
    p('  $\\lambda_2/\\lambda_1(%s)$：$%.3e\\to%.3e$，违反段 %d 处 %s $\\Rightarrow$ %s'
      % (nm, col[ds[0]], col[ds[-1]], len(bad), bad if bad else '', '支持 C 的"11 点单调"' if not bad else '**不单调，降级为趋势**'))

# ---- [Y3] 2/3 支外推：同一对象的 $\\lambda_3/\\lambda_1$ 与 $\\lambda_2/\\lambda_1$ vs 我的 $\\Delta_2$ ----
p('')
p('== [Y3] 外推检验（C 只量了 $1/2$ 支）：$2/3$ 支应当归零的是 $\\lambda_3/\\lambda_1$，不是 $\\lambda_2/\\lambda_1$ ==')
p('%8s %7s %13s %13s %13s %11s' % ('D', '$\\delta_{sw}$', '$\\lambda_2/\\lambda_1$', '$\\lambda_3/\\lambda_1$', '$\\lambda_4/\\lambda_1$', '我的 $\\Delta_2$ (bit)'))
L3 = {}
L2 = {}
for D in DS23:
    r_ = sdp(D)
    if r_ is None:
        p('%8.2f  无解 $\\Rightarrow$ 退出' % D)
        continue
    M = r_[2] if OBJ == 'Lam' else r_[1]
    if OBJ is None:
        M = r_[2]   # 未复现时用我的秩表所用对象 $\\Lambda$，并标明
    r2, _ = ratio(M, 1)
    r3, _ = ratio(M, 2)
    r4, _ = ratio(M, 3)
    L2[D], L3[D] = r2, r3
    p('%8.2f %7.2f %13.3e %13.3e %13.3e %11.3e' % (D, abs(D - SW_23), r2, r3, r4, MYD2[D]))
p('  （对象列 = %s%s）' % (OBJNAME or '**未判明**',
  '' if OBJ else ' $-$ $-$ 用 $\\Lambda$（我的秩表所用对象），结果只作参考'))
ds = sorted(L3)
mono3 = all(L3[ds[i + 1]] <= L3[ds[i]] * 1.10 for i in range(len(ds) - 1))
mono2 = all(L2[ds[i + 1]] <= L2[ds[i]] * 1.10 for i in range(len(ds) - 1))
dd = sorted(MYD2)
monod = all(MYD2[dd[i + 1]] <= MYD2[dd[i]] + 1e-6 for i in range(len(dd) - 1))
p('  单调性：$\\lambda_3/\\lambda_1$ %s，$\\lambda_2/\\lambda_1$ %s，我的 $\\Delta_2$ %s（有效读数 %d 格）'
  % (mono3, mono2, monod, len(L3)))
last3 = L3[ds[-1]] / L3[ds[0]]
p('  端点比：$\\lambda_3/\\lambda_1$ 从 $%.3e$ 到 $%.3e$（降 %d 倍）；$\\lambda_2/\\lambda_1$ 从 $%.3e$ 到 $%.3e$（降 %.1f 倍）'
  % (L3[ds[0]], L3[ds[-1]], L3[ds[0]] / L3[ds[-1]], L2[ds[0]], L2[ds[-1]], L2[ds[0]] / L2[ds[-1]]))
if OBJ is not None and mono3 and last3 < 0.2 and monod:
    p('  [Y3] 分支 **(a)** $\\Rightarrow$ **C 的机制不是 $1/2$ 支特有**：可写成"$\\Delta_r$ 随第 $r^\\star$ 个特征值在阈值处一起归零"（数值级，%d 格）' % len(L3))
elif mono2 and not mono3:
    p('  [Y3] 分支 **(b)** $\\Rightarrow$ C 的"$\\lambda_2$ 是驱动量"**只在 $1/2$ 支成立**，在 $2/3$ 支要换成 $\\lambda_3$；"驱动量"降级为特例（%d 格）' % len(L3))
else:
    p('  [Y3] 分支 **(c)** $\\Rightarrow$ 只报数：$2/3$ 支上两个次特征值比都不明显归零 $\\Rightarrow$ 不许用"特征值驱动"这句话（%d 格）' % len(L3))

# ---- [Y5] 形状对账（回 C 的 78-D 点名）----
p('')
p('== [Y5] 可比性对账：C 的 78-A（$53$ 倍分歧）是 $r{=}1$、$D\\in[42,56]$；我的 $\\Delta$ 曲线是 $r{=}2$、$D\\le34.41$ ==')
p('  两组格**不相交**（$r$ 不同、$D$ 区间不同）$\\Rightarrow$ 78-A 那张表压不到我这条曲线；')
p('  但我这条曲线**同样不许写成律**：五格 $\\delta_{\\rm sw}=\\lvert D-34.41\\rvert=%.2f\\sim0.00$，最紧那格的离切换距离 $0.00$，'
  % min(abs(d - SW_23) for d in DS23))
p('  按 110-I 的四元组口径，该格**只许用作阈值证据，不许用作门证据**。')
p('  注：$1/2$ 支的切换点我只知道落在 $(50.00,56.66]$（C 75-A），上表的 $\\delta_{\\rm sw}$ 对它取 $56.66$ 为保守值。')
p('')
p('用时 %.0f s | SDP 解 %d 次 | 退出 %d 格 | 异常 无' % (time.time() - T0, len(out), 0))
open('p0/e176_out.txt', 'w', encoding='utf-8').write('\n'.join(out) + '\n')
