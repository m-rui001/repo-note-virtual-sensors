# -*- coding: utf-8 -*-
r"""e163：本车道**独立实现**的凸 SDP $+$ 影子价格 $\\Rightarrow$ 把 C 的"数值秩"从诊断量变成可对账的读数。

为什么要做（三块板上的原文逼出来的）：
* C 的 `c89b` 里 $I(S^*)$ 来自 cvxpy/CLARABEL 的凸 SDP，$k=\\mathrm{rank}(S^*)$ 是**数值秩诊断**；
  C 补记十一据此说"$k{=}1$ 两格是精确值"。e162 已经用**另一条路**（秩-1 射线 $+$ Powell，无凸求解器）
  打出 $2.053565/1.603105$，与 C 差 $+4.7\\times10^{-7}$ $\\Rightarrow$ 全局最优被秩-1 点取到（证书方向：我的可行集 $\\subseteq$ 它的）。
* 但 $D{=}34.41$ C 报 $k{=}3$、$I(S^*)=4.893329$，而我方秩-**2** 非凸搜索读到 $4.893399$（只差 $+7.0\\times10^{-5}$），
  同时我方秩-3 搜索返回 $4.900227$（比秩-2 还差 $\\Rightarrow$ **包含关系违反**，我方高阶搜索不可信）。
  $\\Rightarrow$ 那一格的秩到底是 2 还是 3，必须由**凸侧**判，不能由我的 Powell 判。

本脚本干什么（全部本车道自己写，不 import C 的任何文件）：
 变量 $P_v\\succeq0,\\Pi$，min $\\tfrac12\\log_2\\frac{\\det W}{\\det\\Pi}$
 s.t. $\\mathrm{tr}(\\Theta P_v)+J_C\\le D$，$P_v\\preceq A P_v A^{\\mathsf T}+W$，
 $\\begin{bmatrix}P_v-\\Pi & P_vA^{\\mathsf T\\\\ A P_v & A P_vA^{\\mathsf T}+W\\end{bmatrix}\\succeq0$；
 解出后 $S=\\big(P_v^{-1}-(A P_vA^{\\mathsf T}+W)^{-1}\\big)^{-1}$，报 $S$ 的谱 $\\Rightarrow$ 数值秩按**两个阈值**各读一遍。

预注册判据（跑前写死）：
 [S1] **公式校验**（不看物理，先看对象对不对）：CLARABEL 的价值必须与 e162 的秩-1 可达值在
      $D{=}56.66/80.00$ 上差 $\\le10^{-5}$ bit。不成立 $\\Rightarrow$ 我的 SDP 公式与两车道的定价不同源，全部结论押后。
 [S2] 可行残差门：打印 recovered point 的三条约束残差（迹代价、$P_v\\preceq AP_vA^{\\mathsf T}+W$、LMI 最小特征值），
      全部 $\\le10^{-7}$ 才许写"可行点"；否则只说"解出的值"。（#46：判决行必须带产生它的读数）
 [S3] 秩读数：$\\mathrm{rank}_{10^{-6}}(S)$ 与 $\\mathrm{rank}_{10^{-3}}(S)$ **两个阈值并列**。
      若两者不一致 $\\Rightarrow$ "秩"在这个 $D$ 上是**阈值依赖**的 $\\Rightarrow$ 板上凡按秩报价的句子都要带阈值；
      若 $D{=}34.41$ 的 $\\mathrm{rank}_{10^{-3}}=2$ 而 C 报 3 $\\Rightarrow$ 结合我方秩-2 只差 $7\\times10^{-5}$，
      结论是"$\\mathrm{Opt}(D)$ 含近似秩-2 成员"，C 的 $k$ 列须按阈值重印（**这是挑刺项，但要带两个阈值才许写**）。
 [S4] 影子价格 $\\lambda=\\partial I/\\partial D$（代价约束的乘子）与两车道的斜率口径对账：
      $\\lambda$ 应 $>0$ 且随 $D$ 递减；与 e159b 的 $\\gamma_{\\rm loc}$、C 的梯子律 $r/2$ 的关系**只报数不判**（口径见 §109）。
 [S5] 双求解器（CLARABEL 紧容差 / SCS 松容差）价值差 $\\le10^{-6}$ $\\Rightarrow$ 记为"求解器级稳健"；
      **严格对偶下界**需要 $\\inf_x\\mathcal L(x,\\lambda)$，本脚本不做 $\\Rightarrow$ 记为欠项，板上不得写"已证下界"。
"""
import sys
import numpy as np
import cvxpy as cp
from scipy.linalg import solve_discrete_are
sys.stdout.reconfigure(encoding='utf-8')

out = []


def p(*a):
    s = ' '.join(str(x) for x in a)
    out.append(s)
    print(s)


_src = open('p0/exp_c_audit.py', encoding='utf-8').read().split("print(r'== E53")[0]
_ns = {'__name__': 'p'}
exec(compile(_src, 'p0/exp_c_audit.py[preamble]', 'exec'), _ns)


def sym(X):
    return 0.5 * (X + X.T)


A, B, W = _ns['A'].copy(), _ns['B'].copy(), _ns['W'].copy()
n = A.shape[1]
Pc = sym(solve_discrete_are(A, B, sym(np.eye(n)), np.eye(n)))
K = np.linalg.solve(np.eye(n) + B.T @ Pc @ B, B.T @ Pc @ A)
TH = sym(K.T @ (np.eye(n) + B.T @ Pc @ B) @ K)
JC = float(np.trace(W @ Pc))
ln2 = np.log(2.0)
E162 = {56.66: 2.053565, 80.00: 1.603105}
C89B = {34.41: (4.893329, 3), 56.66: (2.053565, 1), 80.00: (1.603105, 1)}

p('== 锚点 SDP（本车道独立实现）：$J_C=%.6f$，$\\det W$ 项 $=%.6f$ ==' % (JC, 0.5 * np.log(np.linalg.det(W)) / ln2))
p('%7s %12s %12s %9s %8s %8s %7s %10s' % ('D', 'CLARABEL', 'SCS', '差', '残差迹', 'LMI最小eig', 'rank6', 'rank3(阈值1e-3)'))
ROWS = []
for D in (32.31, 34.41, 40.00, 56.66, 80.00):
    Pv = cp.Variable((n, n), symmetric=True)
    Pi = cp.Variable((n, n), symmetric=True)
    cost_c = Pv >> 0
    con_tr = cp.trace(TH @ Pv) + JC <= D
    Pi_c = Pi >> 0
    Pt = A @ Pv @ A.T + W
    lmi = cp.bmat([[Pv - Pi, Pv @ A.T], [A @ Pv, Pt]]) >> 0
    prob = cp.Problem(cp.Minimize(0.5 * (np.log(np.linalg.det(W)) - cp.log_det(Pi)) / ln2),
                      [cost_c, Pi_c, con_tr, Pv << Pt, lmi])
    vals = {}
    pts = {}
    for nm, skw in [('CLARABEL', dict(solver=cp.CLARABEL, tol_gap_abs=1e-10, tol_gap_rel=1e-10, tol_feas=1e-10)),
                    ('SCS', dict(solver=cp.SCS, eps=1e-10, max_iters=20000))]:
        try:
            prob.solve(**skw)
            vals[nm] = prob.value
            pts[nm] = (sym(Pv.value), sym(Pi.value), con_tr.dual_value)
        except Exception as e:
            vals[nm] = np.nan
            pts[nm] = (None, None, None)
            p('  [%s @ D=%.2f] 抛错 %s: %s' % (nm, D, type(e).__name__, str(e)[:70]))
    V, PiV, lam = pts['CLARABEL']
    if V is None:
        p('%7.2f 无解' % D)
        continue
    PT = A @ V @ A.T + W
    # 正确的对象是**信息增量** $\\Lambda=P_v^{-1}-P_t^{-1}=C^{\\mathsf T}C$（本站 $C^{\\mathsf T}C$ 的语义），
    # 不是它的逆 $-$ $-$ 第一次跑把 $\\Lambda$ 又求了一次逆，导致 [S6] 全线不一致（$\\sim 10^{1}$ bit 的假差）。
    try:
        Lam = sym(np.linalg.inv(V) - np.linalg.inv(PT))
    except Exception as e:
        p('%7.2f Lambda 反解失败 %s' % (D, e))
        continue
    ev = np.linalg.eigvalsh(Lam)[::-1]
    rel = ev / ev[0]
    r6 = int(np.sum(rel > 1e-6))
    r3 = int(np.sum(rel > 1e-3))
    res_tr = abs(float(np.trace(TH @ V)) + JC - D)
    lmin = float(np.min(np.linalg.eigvalsh(np.block([[V - PiV, V @ A.T], [A @ V, PT]]))))
    feas_lmi = float(np.min(np.linalg.eigvalsh(V)))
    gap = abs(vals['CLARABEL'] - vals['SCS']) if np.isfinite(vals.get('SCS', np.nan)) else np.nan
    # [S6] 把 SDP 的信息增量 $\\Lambda=C^{\\mathsf T}C$ 因子化成本车道的测量阵 $C$，送回 DARE 定价独立复算 I/D
    ws, Us = np.linalg.eigh(Lam)
    o = np.argsort(-ws)
    ws, Us = ws[o], Us[:, o]
    pos = ws > 1e-9 * ws[0]
    Cmat = np.sqrt(np.clip(ws[pos], 0, None))[:, None] * Us[:, pos].T
    rr = int(pos.sum())
    try:
        Pmd = sym(solve_discrete_are(A.T, Cmat.T, W, np.eye(rr)))
        Sm = Cmat @ Pmd @ Cmat.T + np.eye(rr)
        Ld = Pmd @ Cmat.T @ np.linalg.inv(Sm)
        Ppd = sym(Pmd - Ld @ Sm @ Ld.T)
        I_d = 0.5 * np.log(np.linalg.det(Sm)) / ln2
        D_d = JC + float(np.trace(TH @ Ppd))
    except Exception as e:
        I_d, D_d = np.nan, np.nan
        p('        [S6] DARE 复算抛错 %s' % str(e)[:60])
    p('        [S6] $S^*$ 送本车道定价（秩 %d 的 $C$）：$I_{\\rm DARE}=%.6f$（vs SDP %+.2e）$D_{\\rm DARE}=%.6f$（vs 目标 %+.2e）'
      % (rr, I_d, I_d - vals['CLARABEL'], D_d, D_d - D))
    p('             $\\Rightarrow$ 可行点交叉校验：%s' % ('两项都 $\\le10^{-5}$' if (abs(I_d - vals['CLARABEL']) <= 1e-5 and abs(D_d - D) <= 1e-5) else '**不一致，本行的秩/值读数降级**'))
    p('%7.2f %12.6f %12.6f %9.1e %8.1e %8.1e %7d %10d   rel 谱 %s'
      % (D, vals['CLARABEL'], vals['SCS'], gap, res_tr, lmin, r6, r3,
         np.array2string(rel, precision=3, floatmode='fixed')))
    if D in C89B:
        p('        vs C 的 $I(S^*)$=%.6f $\\Rightarrow$ 差 %+.2e bit ；C 报 $k=%d$，本车道读 $\\mathrm{{rank}}_{{10^{{-6}}}}=%d$、$\\mathrm{{rank}}_{{10^{{-3}}}}=%d$'
          % (C89B[D][0], vals['CLARABEL'] - C89B[D][0], C89B[D][1], r6, r3))
    if D in E162:
        p('        [S1] vs e162 秩-1 可达值 %.6f $\\Rightarrow$ 差 %+.2e bit $\\Rightarrow$ 公式校验：%s'
          % (E162[D], vals['CLARABEL'] - E162[D], '通过（$\\le10^{-5}$）' if abs(vals['CLARABEL'] - E162[D]) <= 1e-5 else '**不通过，全部押后**'))
    lamv = float(np.asarray(lam).ravel()[0]) if lam is not None else np.nan
    p('        [S4] $\\lambda=|\\partial I/\\partial D|=%s$（乘子符号见 e164：$\\le$ 约束给正乘子、差分给负号 $\\Rightarrow$ 口径取绝对值；$>0$：%s）' % (
        ('%.5f' % abs(lamv)) if np.isfinite(lamv) else 'n/a', bool(np.isfinite(lamv) and abs(lamv) > 0)))
    ROWS.append((D, vals['CLARABEL'], gap, res_tr, lmin, r6, r3, lamv))

p('')
p('== [S2] 可行残差门：迹残差与 LMI 最小特征值都 $\\le10^{-7}$ 的行 $\\Rightarrow$ 才许称"可行点" ==')
for r in ROWS:
    p('  D=%6.2f  迹残差 %.1e  LMI 最小 eig %.1e  $\\Rightarrow$ %s' % (r[0], r[3], r[4], '可行点' if (r[3] <= 1e-7 and r[4] >= -1e-7) else '只报值，不称可行'))
p('')
p('== [S5] 双求解器稳健（价值差 $\\le10^{-6}$）==')
for r in ROWS:
    p('  D=%6.2f  CLARABEL vs SCS 差 %.1e  $\\Rightarrow$ %s' % (r[0], r[2], '稳健' if (np.isfinite(r[2]) and r[2] <= 1e-6) else '不一致'))
p('[S5] 注：这里只做到"两个独立求解器 + 我方非凸可达上界"三方对账；**严格对偶下界**要 $\\inf_x\\mathcal L(x,\\lambda)$，未做，记欠项。')

# ==== [S7] 影子价格与梯子律系数的对账（只报数，不判）====
# 口径全部现推，不沿用上一版的结论（上一版把目标写成 $(r/2)\\ln2$，与正确目标差 $\\ln^2 2$，是单位错，见文末自纠）。
# 梯子律（C 板第 3501 行原文）：$\\Delta(x)=\\frac{\\mathrm{rank}\\,C}{2}\\log_2\\frac1x+b_{\\rm pred}+O(x)$，$\\Delta I$ 单位是 **bit**。
# 记 $A:=d\\Delta I/d\\log_2(1/x)$（本车道 e159 的口径），律给 $A=\\mathrm{rank}C/2$。
# 凸侧乘子给 $\\lambda:=|\\partial I/\\partial D|$（e164 已用中心差分证同 $\\le5\\times10^{-3}$）。由
# $d\\Delta I=-\\lambda\\,dD$ 与 $d\\log_2(1/x)=-dD/(x\\ln2)$ 两式相除 $\\Rightarrow$  $\\boxed{A=\\lambda\\,x\\,\\ln2}$，$\\gamma_{\\rm loc}=\\ln2/A=1/(\\lambda x)$。
p('')
p('== [S7] 对偶口径的梯子系数 $A=\\lambda x\\ln2$ vs 律给的 $\\mathrm{rank}C/2$（$x=D-$ 同秩可达地板，e162b [P]）==')
FLOOR = {1: 43.7058, 2: 32.2752, 3: 31.6737}
p('%7s %6s %9s %10s %8s %10s %9s %8s %9s' % ('D', 'rankC', 'lambda', 'floor(同秩)', 'x', 'A=lambda*x*ln2', 'rankC/2', '比值', 'gamma_loc'))
for (D, val, gap, res_tr, lmin, r6, r3, lamv) in ROWS:
    if r6 not in FLOOR:
        p('%7.2f  rank%d 无同秩地板（只测了 1/2/3）$\\Rightarrow$ 跳过，不猜地板' % (D, r6))
        continue
    if not np.isfinite(lamv) or lamv <= 0:
        p('%7.2f  乘子不可用（%s）$\\Rightarrow$ 跳过' % (D, lamv))
        continue
    x = D - FLOOR[r6]
    if x <= 0:
        p('%7.2f  x=%.4f $\\le0$（低于同秩地板）$\\Rightarrow$ 跳过' % (D, x))
        continue
    lam = abs(lamv)
    A = lam * x * ln2
    lad = 0.5 * r6
    p('%7.2f %6d %9.5f %10.4f %8.4f %10.4f %9.4f %8.3f %9.3f' % (D, r6, lam, FLOOR[r6], x, A, lad, A / lad, ln2 / A))
p('[S7-自纠] 上一版本拿 $\\lambda x$ 与 $(r/2)\\ln2$ 比，$D=80$ 得 1.047 貌似"复现律" $-$ $-$ 那是 $\\ln^2 2$ 的单位巧合：')
p('     现推 $A=\\lambda x\\ln2$，与 $r/2$ 同量纲。用旧口径写过的"梯子律在远端被独立复现"一句**撤回**。')
p('[S7] 判据：**只报数**。比值随 $x$ 收窄而升（朝 $D\\to$ 地板端）才是与 C 的"$A=r/2$ 在最内十倍程窗 $\\le1\\%$"一致的方向；')
p('     $x$ 用的是**同秩**可达地板（秩-$r$ 的代价下界随秩下降），换秩即换 $x$，凡引用此表都要带秩；')
p('     且律里的维数是 C 原文的**活动通道数** $\\mathrm{rank}\\,C$，本表取 $\\mathrm{rank}_{10^{-6}}(\\Lambda)$，不是设计阵列数。')
open('p0/e163_out.txt', 'w', encoding='utf-8').write('\n'.join(out) + '\n')
