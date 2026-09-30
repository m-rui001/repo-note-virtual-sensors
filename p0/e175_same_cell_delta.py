# -*- coding: utf-8 -*-
r"""e175：**同格、同 r、两台机器**对照 $\Delta_r(D)$ —— 了结 C 的 74-B 与 75-C 对我的两条反驳。

背景（逐字引用，板 4435–4465 与 4379–4434）：
  C 74-B：“你的 $10^{-7}$ 恰恰是‘$\Delta=0$ 是应然’的那些格……所以它**不构成对严格情形的读数**，
          也就不能用来判‘门过不了’”，并立下新规矩：“今后凡引用 $\Delta$ 的读数，一律带三元组 $(D,r,r^\star(D))$”。
  C 74-A 的表（`.work3/c90b_out.txt`，板 4385–4391 行）**在严格格 $r=r^\star(D)-1$ 上**给出
          $\Delta_{\min}=1.261\times10^{-1}\,(D{=}33.00),\ 8.342\times10^{-2}\,(33.50),\ 6.906\times10^{-2}\,(34.00),\ 3.524\times10^{-1}\,(45.00)$，
          随机 $V$ 可行数只有 $2$–$4/150$。
  我 e165 在 $(34.41,\,r{=}2,\,r^\star{=}3)$ 给 $\Delta_2=+2.14\times10^{-7}$。

**我的假设（本轮要杀的正是它，也是 C 的读数与我的读数能否同框的关键）**：
  $\Delta_r(D)$ 关于 $D$ **单调下降**，且在秩切换阈值 $D$ 上**归零**（那里秩 $r$ 与秩 $r^\star$ 两支并列）。
  $34.41$ 恰是 $2/3$ 支的切换点（e163：$34.41$ 处 $\mathrm{rank}_{10^{-3}}=2$、$\mathrm{rank}_{10^{-6}}=3$）。
  ⇒ 我的 $2\times10^{-7}$ 是**阈值格**读数，C 的 $7\times10^{-2}$–$1.3\times10^{-1}$ 是**区间内部**读数，
    两者**在同一条曲线上**，不是“格不同所以不可比”（C 的 74-B 说“不可比”过头了：可比，只要带三元组并画成曲线）。
  若成立：① C 的 c90b 由**独立机器**确证（它自己只有 $2$–$4/150$ 可行点，说服力弱）；
          ② 我 WORKLOG 的“路由② 的门注定过不了”**收回**——门在内部严格格上 $\Delta\sim10^{-1}\gg10^{-4}$，**过**。

预注册判据（跑前写死，本文件先落盘再运行）：
 [X0] **参照门**：我自算的 $I_{\rm unc}(D)$ 与 C 板 74-A 印的 $I_{\rm unc}$ 差 $\le10^{-5}$ bit（六位一致）。
      不过 $-$ $-$ 两台机器解的不是同一个凸问题 $-$ $-$ 后面所有 $\Delta$ 比较**一律不判**。
 [X1] **同格相符**：对 C 报 $\Delta_{\min}>10^{-3}$ 的每一格，我的 $\Delta_r$ 落在 C 读数的 $[0.5\times,\,2\times]$ 内
      $\\Rightarrow$ 判“两台独立机器在同 $D$ 同 $r$ 上相符（数值级）”。
 [X2] **反例**：我的 $\Delta_r\le10^{-4}$ 而 C 的 $>10^{-3}$（同格）$\\Rightarrow$ 判“C 的 $\min_V$ 搜漏了”，并印我的标架。
 [X3] 两条都不满足 $\\Rightarrow$ 未解释，**不判**，只把两张数并列。
 [X4] **曲线形状门**：$\Delta_2(33.00)\ge\Delta_2(33.50)\ge\Delta_2(34.00)\ge\Delta_2(34.20)\ge\Delta_2(34.41)$
      （容差 $10^{-6}$）$\\Rightarrow$ 支持“单调下降 + 阈值归零”的解释；违反则**只报数**，不许用该解释。
 [X5] 覆盖对照：逐格印我的“原始有效起点 $k/25$”，与 C 的“可行/样本 $2$–$4/150$”并列 $\\Rightarrow$
      检验“可行子空间在低 $D$ 端极稀”是**两条机器共有的性质**，不是任一方的搜索缺陷。
口径：角度一律用**统一口径**（列归一 $+$ `qr` 正交化后的主角度，e173 的 `basis()`），
      三元组 $(D,r,r^\star(D))$ 逐格印；$\Delta$ 是**上界**，措辞不得升格为“已证”。
"""
import sys
import time
import numpy as np
import cvxpy as cp
from scipy.linalg import solve_discrete_are
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
wth, Uth = np.linalg.eigh(TH)
Uth = Uth[:, np.argsort(-wth)]
ln2 = np.log(2.0)
Pt0 = A @ np.eye(n) @ A.T + W

FLOOR = {1: 43.7058, 2: 32.2752, 3: 31.6737}
C_IUNC = {33.00: 6.172695, 33.50: 5.593895, 34.00: 5.169971, 42.00: 3.004696, 45.00: 2.705659}
C_DELTA = {33.00: 1.261e-1, 33.50: 8.342e-2, 34.00: 6.906e-2, 45.00: 3.524e-1}
C_RANK = {33.00: 3, 33.50: 3, 34.00: 3, 42.00: 2, 45.00: 2}
C_FEAS = {33.00: '2/150', 33.50: '2/150', 34.00: '4/150', 42.00: '0/150', 45.00: '3/150'}
DS = [33.00, 33.50, 34.00, 34.20, 34.41, 42.00, 45.00]

p('== 锚点 $J_C=%.6f$；本车道自算 $I_{\\rm unc}$ vs C 板 74-A ==' % JC)
p('%7s %12s %12s %9s %6s %6s %6s %s' % ('D', '我的CLARABEL', 'C的I_unc', '差', 'rank6', 'rank3', 'r星', '可行/样本(C)'))
SDP = {}
RANK6 = {}
RANK3 = {}
for D in DS:
    Pv = cp.Variable((n, n), symmetric=True)
    Pi = cp.Variable((n, n), symmetric=True)
    con_tr = cp.trace(TH @ Pv) + JC <= D
    PT = A @ Pv @ A.T + W
    lmi = cp.bmat([[Pv - Pi, Pv @ A.T], [A @ Pv, PT]]) >> 0
    prob = cp.Problem(cp.Minimize(0.5 * (np.log(np.linalg.det(W)) - cp.log_det(Pi)) / ln2),
                      [Pv >> 0, Pi >> 0, con_tr, Pv << PT, lmi])
    try:
        prob.solve(solver=cp.CLARABEL, tol_gap_abs=1e-10, tol_gap_rel=1e-10, tol_feas=1e-10)
        v = sym(Pv.value)
        lamT = con_tr.dual_value
    except Exception as e:
        p('%7.2f  SDP 抛错 %s: %s' % (D, type(e).__name__, str(e)[:60]))
        continue
    PTv = A @ v @ A.T + W
    Lam = sym(np.linalg.inv(v) - np.linalg.inv(PTv))
    ev = np.linalg.eigvalsh(Lam)[::-1]
    rel = ev / ev[0]
    r6 = int(np.sum(rel > 1e-6))
    r3 = int(np.sum(rel > 1e-3))
    SDP[D] = float(prob.value)
    RANK6[D] = r6
    RANK3[D] = r3
    dc = SDP[D] - C_IUNC[D] if D in C_IUNC else np.nan
    p('%7.2f %12.6f %12s %9s %6d %6d %6d %s'
      % (D, SDP[D], ('%.6f' % C_IUNC[D]) if D in C_IUNC else '——',
         ('%+.1e' % dc) if np.isfinite(dc) else '——', r6, r3, r6,
         C_FEAS.get(D, '——')))

p('')
p(r'== [X0] 参照门：$|I_{\rm unc}^{\rm mine}-I_{\rm unc}^{\rm C}|\le10^{-5}$ bit ==')
x0ok = True
for D in sorted(C_IUNC):
    if D not in SDP:
        p('  $D=%.2f$：我方无解 $\\Rightarrow$ 该格退出比较' % D)
        x0ok = False
        continue
    dc = abs(SDP[D] - C_IUNC[D])
    p('  $D=%.2f$：差 %+.1e $\\Rightarrow$ %s' % (D, dc, '过' if dc <= 1e-5 else '**不过**'))
    x0ok = x0ok and dc <= 1e-5
p('  [X0] 总结：%s' % ('**通过** $-$ $-$ 两台机器解同一个凸问题，下面的 $\\Delta$ 比较有效'
                       if x0ok else '**不通过 $-$ $-$ 后面的 $\\Delta$ 比较一律不判**'))

p('')
p('== 逐格 $\\min_{\\mathrm{rank}\\,V=r}\\Delta$，$r=r^*(D)-1$（25 起点：$\\Theta$ 前 $r$ 列 $+$ 24 随机正交标架）==')
p('%7s %8s %3s %11s %11s %9s %6s %s' % ('D', 'r星', 'r', 'min I', 'delta_r', '残差', 'k/25', '主角度(deg, 统一口径)'))


def basis(G):
    G = np.asarray(G, float).reshape(n, -1)
    Q, Rq = np.linalg.qr(G / np.linalg.norm(G, axis=0))
    d = np.sign(np.diag(Rq))
    d[d == 0] = 1.0
    return Q * d


def angles(Ga, Gb):
    sv = np.linalg.svd(basis(Ga).T @ basis(Gb), compute_uv=False)
    return np.degrees(np.arccos(np.clip(sv, -1.0, 1.0)))


def cell(t, r, starts, coarse=700, fine=3000, ntop=6):
    raw = []
    for x0 in starts:
        P['GUESS'][0] = 0.0
        I, Dg, lg, G = P['polish'](x0, t, None, coarse)
        if np.isfinite(I):
            raw.append((I, abs(Dg - t), G))
    raw.sort(key=lambda z: z[0])
    ref = []
    for I, Dg, G in raw[:ntop]:
        P['GUESS'][0] = 0.0
        I2, D2, lg2, G2 = P['polish'](G, t, None, fine)
        if np.isfinite(I2):
            ref.append((I2, abs(D2 - t), G2))
    return raw, sorted(raw + ref, key=lambda z: z[0])


RES = {}
BEST = {}
rng = np.random.default_rng(777)
for t in DS:
    rs = RANK6.get(t)
    if rs is None:
        p('%7.2f  SDP 无秩读数 $\\Rightarrow$ 跳过' % t)
        continue
    r = rs - 1
    if r < 1:
        p('%7.2f %8d %3d  $r^*=0$ 无严格格 $\\Rightarrow$ 不判' % (t, rs, r))
        continue
    if t <= FLOOR.get(r, 0.0):
        p('%7.2f %8d %3d  不可达（同秩地板 %.4f $\\ge D$）$\\Rightarrow$ 记 $+\\infty$' % (t, rs, r, FLOOR[r]))
        RES[(t, r)] = np.inf
        continue
    starts = [Uth[:, :r]] + [np.linalg.qr(rng.normal(0, 1, (n, r)))[0] for _ in range(24)]
    raw, allv = cell(t, r, starts)
    if not allv:
        p('%7.2f %8d %3d  有效 0 条 $\\Rightarrow$ 无读数（与 C 的 0/150 同类）' % (t, rs, r))
        RES[(t, r)] = np.nan
        continue
    I2, rb, G2 = allv[0]
    Is = np.array([q[0] for q in allv])
    dlt = I2 - SDP[t]
    RES[(t, r)] = dlt
    BEST[(t, r)] = G2
    pv = angles(Uth[:, :r], G2)
    p('%7.2f %8d %3d %11.6f %+11.3e %9.1e %6s %s'
      % (t, rs, r, I2, dlt, rb, '%d/%d' % (len(raw), len(starts)),
         np.array2string(pv, precision=1)))

p('')
p('== [X4] 单调 $+$ 阈值归零 ==')
seq = [(t, RES.get((t, 2), np.nan)) for t in (33.00, 33.50, 34.00, 34.20, 34.41)]
x4 = True
for (ta, da), (tb, db) in zip(seq, seq[1:]):
    if not (np.isfinite(da) and np.isfinite(db)):
        p('  %.2f$\\to$%.2f：有无读数 $\\Rightarrow$ 该段不检' % (ta, tb))
        continue
    ok = da >= db - 1e-6
    x4 = x4 and ok
    p('  $\\Delta_2(%.2f)=%.3e\\ge\\Delta_2(%.2f)=%.3e$：%s' % (ta, da, tb, db, 'ok' if ok else '**违反**'))
p('  阈值端点：$\\Delta_2(34.20)$ vs $\\Delta_2(34.41)$（后者是 e165 的 $2.14e{-}07$ 的复算）')
p('  [X4] 总结：%s' % ('**支持**“单调下降 $+$ 阈值归零”解释' if x4 else '**不支持 $-$ $-$ 只报数，不用该解释**'))

p('')
p('== [X1]/[X2] 同格对照（C 的 74-A 列 vs 我的上表；两边都是 $r=r^*(D)-1$）==')
n1 = n2 = n3 = 0
for t in sorted(C_DELTA):
    r = C_RANK[t] - 1
    mine = RES.get((t, r), np.nan)
    c = C_DELTA[t]
    if not np.isfinite(mine):
        p('  $(%.2f,\\,r{%d},\\,r^*{%d})$：我无读数 $\\Rightarrow$ 不判（C：%s）' % (t, r, C_RANK[t], C_FEAS[t]))
        n3 += 1
        continue
    if mine <= 1e-4 and c > 1e-3:
        p('  $(%.2f,\\,r{%d},\\,r^*{%d})$：我 %.3e vs C %.3e $\\Rightarrow$ [X2] **反例：C 的 $\\min_V$ 搜漏**'
          % (t, r, C_RANK[t], mine, c))
        p('      我的标架 $V$ 列空间主角度（对 $\\Theta$ 前 %d）：%s'
          % (r, np.array2string(angles(Uth[:, :r], BEST[(t, r)]), precision=2)))
        n2 += 1
    elif 0.5 * c <= mine <= 2.0 * c:
        p('  $(%.2f,\\,r{%d},\\,r^*{%d})$：我 %.3e vs C %.3e（比 %.2f）$\\Rightarrow$ [X1] **相符**'
          % (t, r, C_RANK[t], mine, c, mine / c))
        n1 += 1
    else:
        p('  $(%.2f,\\,r{%d},\\,r^*{%d})$：我 %.3e vs C %.3e（比 %.2f）$\\Rightarrow$ [X3] 未解释'
          % (t, r, C_RANK[t], mine, c, mine / c))
        n3 += 1
p('')
p('== 判决 ==')
p('  有效读数 %d 格；[X1] 相符 %d 格、[X2] 反例 %d 格、[X3] 未解释 %d 格（#46：判决行必须带条数）'
  % (sum(1 for v in RES.values() if np.isfinite(v)), n1, n2, n3))
if not x0ok:
    p('  [X0] 未过 $-$ $-$ **不判**，两张数并列存档')
elif n2 > 0:
    p('  判“**两车道在严格格上不相容**”，逐格列出后交 C 复核')
elif n1 >= 3:
    p('  判“**两条独立机器在同格同 $r$ 上相符**（数值级，$n_1=%d$ 格）”$-$ $-$ '
      '我 WORKLOG 的“路由② 门注定过不了”**收回**：内部严格格 $\\Delta\\sim10^{-1}\\gg10^{-4}$，门**过**' % n1)
else:
    p('  相符格数 $n_1=%d<3$ $\\Rightarrow$ 只报数，不判' % n1)
p('  [X5] 覆盖对照：我的有效起点 $k/25$ 与 C 的可行 $/150$ 逐格并列（上表 $k/25$ 列 vs C 的 %s）'
  % str(list(C_FEAS.items())))
p('')
p('用时 %.0f s | DARE %d 次 | 有效读数 %d | 回退 %d | 异常 %s'
  % (time.time() - T0, P['NDB'][0], P['NVALID'][0], P['NFALLBACK'][0], P['EXC'] or '无'))
open('p0/e175_out.txt', 'w', encoding='utf-8').write('\n'.join(out) + '\n')
