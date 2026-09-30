# -*- coding: utf-8 -*-
r"""e164：把 e163 的 $\\lambda$ 从"求解器给的乘子"降级为"可独立复核的差分读数"。

为什么要做：e163 [S7] 用的是 cvxpy 的 `dual_value`，那是一个**规范化之后**的量。
我要在板上写 $\\gamma_{\\rm loc}=1/\\lambda$ 与 $\\lambda=\\partial I/\\partial D$，
就得先证明"乘子 == 价值对右端项的导数"，而不是靠求解器贴的标签。
$D=80$ 那格读出 $\\lambda=0.01000$ $-$ $-$ 过分整齐的数最可疑，先差分再说话。

坑（第一次跑就踩到，如实记在这里）：$D=32.31-h$ 会掉到 SDP 自身的代价地板之下 $\\Rightarrow$ CLARABEL 抛
`SolverError`。所以步长必须**自适应**，并且把地板本身量出来。

做法（同一套 SDP 公式，右端项当参数；全部本车道自写）：
 [F0] 二分 SDP 的可达地板 $D_{\\min}$（价值 $\\to0$ 的那一端）。它必须 $\\le$ e162b 的秩-3 地板 31.6737，
      因为我的秩-$r$ 可行集是 SDP 可行集的子集。
 [F1] 中心差分 $\\lambda_{\\rm fd}=[I(D+h)-I(D-h)]/(2h)$，$h$ 取"不越地板"的最大值与 $2.0$ 的较小者，
      再用 $h/4$ 复核；两 $h$ 相对差 $\\le5\\%$ 才算差分收敛。
 [F2] 乘子与差分（细 $h$）相对差 $\\le2\\%$ $\\Rightarrow$ 允许写 $\\lambda=\\partial I/\\partial D$；
      否则板上只许写"乘子读数"，$\\gamma_{\\rm loc}=1/\\lambda$ 一并撤回。
 [F3] $\\lambda$ 随 $D$ 递减？$I(D)$ 的二阶差分符号？（$\\partial I/\\partial D<0$ 且绝对值递减 $\\Rightarrow$ 凸）
"""
import sys
import numpy as np
import cvxpy as cp
from scipy.linalg import solve_discrete_are
sys.stdout.reconfigure(encoding='utf-8')

out = []
EXC = {}


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
SOLVE = dict(solver=cp.CLARABEL, tol_gap_abs=1e-10, tol_gap_rel=1e-10, tol_feas=1e-10)


def sdp(D, want_dual=False):
    try:
        Pv = cp.Variable((n, n), symmetric=True)
        Pi = cp.Variable((n, n), symmetric=True)
        Pt = A @ Pv @ A.T + W
        con_tr = cp.trace(TH @ Pv) + JC <= D
        prob = cp.Problem(cp.Minimize(0.5 * (np.log(np.linalg.det(W)) - cp.log_det(Pi)) / ln2),
                          [con_tr, Pv >> 0, Pi >> 0, Pv << Pt,
                           cp.bmat([[Pv - Pi, Pv @ A.T], [A @ Pv, Pt]]) >> 0])
        prob.solve(**SOLVE)
        v = prob.value
        if v is None or not np.isfinite(v):
            EXC['value-nonfinite@%.4f' % D] = EXC.get('value-nonfinite@%.4f' % D, 0) + 1
            return (np.nan, np.nan) if want_dual else np.nan
        if want_dual:
            d = con_tr.dual_value
            d = float(np.asarray(d).ravel()[0]) if d is not None else np.nan
            return float(v), d
        return float(v)
    except Exception as e:
        k = type(e).__name__ + ':' + str(e)[:50]
        EXC[k] = EXC.get(k, 0) + 1
        return (np.nan, np.nan) if want_dual else np.nan


p('== 锚点（$J_C=%.6f$）：$I(D)$ 是"给定代价预算下的最小率"，故 $\\partial I/\\partial D<0$ ==' % JC)

p('')
p('== [F0] 低代价端：秩-1 在本车道不可达的 $D$，SDP 仍给有限值 $-$ $-$ 这本身就是"秩必须升高"的凸侧证据 ==')
p('  （第一次跑我把地板方向写反了：$I(D)$ 随 $D$ **递减**，$I\\to0$ 在 $D\\to$ 大的一端，')
p('    而"代价最小端"对应 $I\\to+\\infty$，所以二分"价值 $\\to0$"是无意义的操作，改成直接对照读数。）')
for D, r1f in [(32.31, 43.7058), (34.41, 43.7058)]:
    v = sdp(D)
    p('  $D=%.2f$：本车道**秩-1** 可达地板 %.4f $>D$ $\\Rightarrow$ 秩-1 族在该格为空；SDP（秩不限）价值 %.6f 有限'
      % (D, r1f, v))
p('  $\\Rightarrow$ 两车道口径自洽的证据：低代价格的 $S^*$ 必然用掉更高秩（e163 读 rank$_{10^{-6}}$=4/3），')
p('    而高代价格（$D=56.66,80$）秩-1 就已取到 SDP 值 $-$ $-$ 见 e163 [S1]。')
p('  [F0-注] 这里**不**宣称量出了"SDP 的代价地板"：秩不限的 $D_{\\min}$ 需要 $\\inf$ 秩-$n$ 设计的极限，本脚本未做，记欠项。')

DS = [32.31, 34.41, 40.00, 56.66, 80.00]
FLOOR_RANK = {1: 43.7058, 2: 32.2752, 3: 31.6737}
p('')
p('== [F1/F2] 乘子 vs 中心差分（$h$ 自适应：下界不得越出本车道同秩地板，否则 $I$ 无定义）==')
p('   **符号口径先立**：$I(D)$ 递减 $\\Rightarrow$ 中心差分为负，而 cvxpy 对 $\\le$ 约束给**正**乘子。')
p('   故板上口径必须是 $\\lambda:=|\\partial I/\\partial D|=-\\text{(带号乘子)}$；下表"乘子"列直接取绝对值，')
p('   并单列核对符号（若某格乘子与差分同号 $\\Rightarrow$ 该格口径作废）。')
p('%7s %10s %8s %10s %8s %10s %9s %8s' % ('D', '|乘子|', 'h(粗)', 'fd(粗)', 'h(细)', 'fd(细)', '相对差', '符号相反'))
ROWS = []
for D in DS:
    val, mult = sdp(D, want_dual=True)
    # 可用步长：上界无限，下界不得越出本格的"同秩可达地板"（否则 $I$ 无定义、价值读数无意义）
    fl = np.nan
    for r, f in FLOOR_RANK.items():
        if D > f:
            fl = f if np.isnan(fl) else max(fl, f)
    hbig = 2.0 if np.isnan(fl) else max(1e-3, min(2.0, 0.8 * (D - fl)))
    hs = [hbig, hbig / 4.0]
    fds = []
    for h in hs:
        a = sdp(D + h)
        b = sdp(D - h)
        fds.append((a - b) / (2 * h))
    if not np.isfinite(mult) or not np.isfinite(fds[1]):
        p('%7.2f 乘子 %s 差分 %s $\\Rightarrow$ 不可判（异常计数见文末）' % (D, mult, fds[1]))
        continue
    am = abs(mult)
    rel = abs(am - abs(fds[1])) / abs(fds[1])
    sign_ok = bool(mult * fds[1] < 0)
    p('%7.2f %10.5f %8.4f %10.5f %8.4f %10.5f %9.2e %8s' % (D, am, hs[0], fds[0], hs[1], fds[1], rel, sign_ok))
    ROWS.append((D, am, hs[0], fds[0], hs[1], fds[1], rel, abs(fds[0] - fds[1]) / abs(fds[1]), sign_ok))

p('')
p('== [F2] 结算（判决行自带产生它的读数条数，#46）==')
n1 = sum(1 for r in ROWS if r[7] <= 0.05)
n2 = sum(1 for r in ROWS if r[6] <= 0.02)
ns = sum(1 for r in ROWS if r[8])
p('  有效格 %d/%d；两 $h$ 之间的相对差 $\\le$ 5%% 的格 %d 个；|乘子| 与 |差分| 相对差 $\\le$ 2%% 的格 %d 个；符号相反 %d 个'
  % (len(ROWS), len(DS), n1, n2, ns))
p('  [F2] 判定：%s' % ('全部通过 $\\Rightarrow$ 板上可写 $\\lambda=\\lvert\\partial I/\\partial D\\rvert$、$\\gamma_{\\rm loc}=1/\\lambda$'
                        if (len(ROWS) == len(DS) and n2 == len(DS) and ns == len(DS)) else
                        '未全通过（有效 %d，相对差通过 %d，符号通过 %d）$\\Rightarrow$ 只许写"乘子读数"，'
                        '$\\gamma_{\\rm loc}=1/\\lambda$ 的推论降级' % (len(ROWS), n2, ns)))
p('')
p('== [F3] 单调性与曲率（全用差分口径，不用乘子标签）==')
lam = [abs(r[5]) for r in ROWS]
p('  $\\lambda_{\\rm fd}=|\\partial I/\\partial D|$（细 $h$）：%s' % np.array2string(np.array(lam), precision=5))
p('  随 $D$ 严格递减？ %s（$\\Rightarrow$ 预算越宽、边际信息越便宜，与"率是代价的减函数"自洽）'
  % bool(all(lam[i] > lam[i + 1] for i in range(len(lam) - 1))))
vals = [sdp(D) for D in DS]
p('  $I(D)$：%s' % np.array2string(np.array(vals), precision=6))
sec = [(vals[i + 1] - 2 * vals[i] + vals[i - 1]) / (DS[i + 1] - DS[i - 1]) ** 2 for i in range(1, len(DS) - 1)]
p('  二阶差分（按网格宽度归一）%s $\\Rightarrow$ 全部非负（凸且降）：%s'
  % (np.array2string(np.array(sec), precision=5, floatmode='fixed'), bool(all(x >= -1e-12 for x in sec))))
i80 = DS.index(80.0)
p('[F3] 单查 $D=80$ 的整齐读数：|乘子| %.5f vs |差分| %.5f，相对差 %.2e $\\Rightarrow$ %s'
  % (ROWS[i80][1], ROWS[i80][5], ROWS[i80][6], '不是标签伪影' if ROWS[i80][6] <= 0.02 else '**是伪影，撤回**'))
p('')
p('异常计数（绝不静默）：%s' % (dict(EXC) if EXC else '无'))
open('p0/e164_out.txt', 'w', encoding='utf-8').write('\n'.join(out) + '\n')
