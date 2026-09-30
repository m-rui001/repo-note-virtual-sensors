r"""E96e = 把式 (18) 的**无约束精确前沿**当证书用：受限族的最小代价不可能低于它。

理由只有一条，但很硬：$\mathrm{DI}_{\text{free}}(D)$ 是对**全体 Borel 因果随机核**（1510.04214 的
$\Gamma$，Theorem 1 证明线性高斯三阶段结构在其中最优，且 (18) 的 max-det SDP 精确求出它）的最小信息。
把传感族缩小成"平面 / 秩 $r$ / 任务受限"只能让最小代价**变大**，所以
$$\boxed{\;D_{\mathcal V}(I)\ \ge\ D^*_{\text{free}}(I)\quad\forall\,\mathcal V\;}$$
其中右端可由 SDP 精确算到任意位数。这一条同时替换我车道原来自己造的代价侧下界 $41.485$
（`p0/e95_kh_inverse_bound.py`），并且**不需要任何新定理**——它是已发表结论的直接推论。

== 判据（跑前写死）==
 [Z1] 反解单调性：$D\uparrow$ 时 $\mathrm{DI}_{\text{free}}(D)\downarrow$，网格上必须严格单调，
      否则 SDP 读数不可用（$\pm10^{-3}$ 内的平台视为数值噪声，不算违反）。
 [Z2] 每个报告点必须同时给出 SDP 状态（optimal/nan）与求解器残余；任何 `nan` 行不进板子。
 [Z3] 与 E96c 的 $D^*_{\text{free}}(3)=42.0402$ 对账，差 $\le10^{-3}$。
 [Z4] 关键判决量：蓝线代价 $43.711$ 处的无约束所需率 $\mathrm{DI}_{\text{free}}(43.711)$。
      若 $>3$ ⇒ 蓝线**在任何**机制族里都不可行，Remark 2 直接结案（最强的可能结果）；
      若 $\le3$ ⇒ 蓝线不被排除，Remark 2 只能表述为"蓝线低于已知受限族构造值、
      但高于无约束精确前沿"。**这条判据不许事后改。**
"""
import sys
sys.stdout.reconfigure(encoding='utf-8')
import numpy as np
import cvxpy as cp
from scipy.optimize import brentq

np.set_printoptions(precision=5, suppress=True, linewidth=170)
_src = open('p0/exp_c_audit.py', encoding='utf-8').read().split("print(r'== E53")[0]
_ns = {'__name__': 'p'}
exec(compile(_src, 'p0/exp_c_audit.py[preamble]', 'exec'), _ns)
A, W, TH, JC, n, sym = _ns['A'], _ns['W'], _ns['TH'], _ns['JC'], _ns['n'], _ns['sym']
ln2 = np.log(2.0)
LOGDET_W = float(np.linalg.slogdet(W)[1])


def Ptil(P):
    return sym(A @ P @ A.T + W)


def free_sdp(D, solver=cp.CLARABEL):
    P = cp.Variable((n, n), symmetric=True)
    Pi = cp.Variable((n, n), symmetric=True)
    cons = [P >> 0, Pi >> 0, P << Ptil(P), cp.trace(TH @ P) + JC <= D,
            cp.bmat([[P - Pi, P @ A.T], [A @ P, Ptil(P)]]) >> 0]
    prob = cp.Problem(cp.Minimize(-0.5 * cp.log_det(Pi) + 0.5 * LOGDET_W), cons)
    prob.solve(solver=solver)
    if P.value is None:
        return np.nan, None, prob.status
    return prob.value / ln2, sym(P.value), prob.status


def Dstar_of_rate(I0, lo=None, hi=200.0):
    """沿 D 二分反解 DI_free(D)=I0；先扫出变号区间。"""
    if lo is None:
        lo = JC + 0.05
    g = lambda D: free_sdp(D)[0] - I0
    a, b = lo, hi
    ga, gb = g(a), g(b)
    if not (np.isfinite(ga) and np.isfinite(gb)) or ga * gb > 0:
        # 扫网格找变号
        grid = np.concatenate([JC + np.geomspace(0.05, 5.0, 12), np.linspace(8.0, 200.0, 30)])
        vals = []
        for D in grid:
            v = free_sdp(D)[0]
            vals.append(v - I0 if np.isfinite(v) else np.nan)
        vals = np.array(vals)
        ok = np.where(np.isfinite(vals))[0]
        for i in range(len(ok) - 1):
            j, k = ok[i], ok[i + 1]
            if vals[j] * vals[k] < 0:
                a, b = grid[j], grid[k]
                break
        else:
            return np.nan, 'no bracket'
    try:
        D = brentq(g, a, b, xtol=1e-5, rtol=1e-12, maxiter=200)
    except Exception as e:
        return np.nan, f'brentq:{type(e).__name__}'
    v, P, st = free_sdp(D)
    return D, f'{st} 核验 DI={v:.5f} 区间[{a:.3f},{b:.2f}]'


print('== (1) 无约束精确 converse：给定代价 D 所需的最小率（bit/sample）==')
print('   注：D 必须 > 地板 JC=%.4f，否则 DI=+inf（原文 §V 的垂直渐近线）' % JC)
rows = []
for D in (42.0402, 43.0, 43.711, 43.80, 45.4537, 50.6023):
    v, P, st = free_sdp(D)
    sv = None
    if P is not None:
        S = sym(np.linalg.inv(P) - np.linalg.inv(Ptil(P)))
        s = np.sort(np.linalg.svd(S, compute_uv=False))[::-1]
        sv = int((s > 1e-3 * s.max()).sum())
    rows.append((D, v, st, sv))
    print(f'   D={D:9.4f}  DI_free={v:8.5f} bit   rank(0.1%)={sv}   status={st}')
prev = None
mono = True
for D, v, st, sv in rows:
    if prev and np.isfinite(v) and np.isfinite(prev[1]) and v > prev[1] + 1e-3:
        mono = False
    prev = (D, v)
print(f' [Z1] 单调性 {"通过" if mono else "违反——SDP 读数不可用"}'
      f'   [Z2] 全部 status=optimal: {all(r[2]=="optimal" for r in rows)}')

print('\n== (2) 无约束精确 converse：给定率 I 所需的最小代价（对任何受限族都是下界）==')
tab = []
for I0 in (1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0):
    D, info = Dstar_of_rate(I0)
    tab.append((I0, D, info))
    print(f'   I={I0:4.1f} bit   D*_free={D:9.4f}   {info}')
d3 = [D for I, D, _ in tab if I == 3.0][0]
print(f' [Z3] 与 E96c 的 42.0402 差 {d3-42.0402:+.5f}'
      f'   [Z2] 有效行 {sum(1 for _, D, _ in tab if np.isfinite(D))}/{len(tab)}')

print('\n== (3) 受限族构造值 vs 无约束精确下界（Remark 2 的诚实口径）==')
cands = [('无约束（本车道构造 E71/E74）', 42.0427, 3.0),
         ('tail-3 平面（Powell 上界）', 45.4537, 3.0),
         ('tail-2 平面（Powell 上界）', 50.6023, 3.0),
         ('原文蓝线数字化（Remark 2 争议点）', 43.711, 3.0)]
for name, D, I in cands:
    dv = free_sdp(D)[0]
    print(f'   {name:32s} D={D:8.4f} @I={I}  |  无约束在同一 D 只需 {dv:.4f} bit '
          f'=> 该族在该 D 的**率溢价** {I-dv:+.4f} bit')
for name, D, I in cands:
    if np.isfinite(D):
        print(f'   {name:32s}：代价下界 D*_free(3 bit)={d3:.4f} ⇒ 高出精确 converse {D-d3:+.4f}')
bl = free_sdp(43.711)[0]
print(f' [Z4] 蓝线判决：DI_free(43.711) = {bl:.5f} bit —— '
      f'{"蓝线在任何机制族都不可行（Remark 2 结案）" if bl > 3.0 else "蓝线不被排除；Remark 2 只能写成区间陈述"}')

print('\n== (4) 与原文发表锚点的对账（工具链验收用例，必须在偏差内）==')
for D, pub in ((33.0, 6.133), (40.0, 3.266), (80.0, 1.602)):
    v = free_sdp(D)[0]
    print(f'   D={D:6.1f}  原文 {pub:.3f}  本脚本 {v:.4f}  差 {v-pub:+.4f}')
