r"""E96f = 墙律 $2^{-2I/r}$ 的指数 $r$ 是否**就是原文 SDP 输出的 rank(SNR)**——在已发表精确前沿上核对。

我车道 §25 的墙律是"经验规律"（$r{=}1/2/3$ 三个族各自拟合斜率 $-2\log_2 2/r$）。
E96e 顺手给了第一个**不受我构造影响**的检验台：式 (18) 的 $\mathrm{DI}_{\text{free}}(D)$ 是
已证明精确的（Theorem 1/2），它的秩也是它自己输出的（$r_D=\mathrm{rank}(\mathrm{SNR}_D)$，0.1% 口径）。
于是可以直接问：
$$\frac{d\,\mathrm{DI}}{d\,\log_2(1/x)},\qquad x=D-\mathrm{tr}(WS)$$
在数值上是否 $\approx r_D/2$。若是，墙律的"指数由秩选"这句话**不再依赖我的构造族**，
而是在文献自己的曲线上成立——这对定理化是有用的，但也**很可能是高斯率失真的常识**
（"斜率计活跃维数"是逆水 filling 的标准读法），所以：

== 判据（跑前写死）==
 [U1] 逐窗拟合：把 $D$ 网格按 $r_D$ 分段，每段内用最小二乘拟合斜率 $s$，报告 $2s$ 与 $r_D$。
      只有 $|2s-r_D|/r_D\le5\%$ 才算支持"指数=秩"。
 [U2] 反例优先：只要出现 $r_D$ 相同而 $2s$ 明显不同（$>10\%$）的相邻窗口，立刻记为否证。
 [U3] **查文献义务**（跑完就必须做，不许跳过）：这条陈述若能在
      Tanaka–Mohajerin Esfahani–Mitter 的注水/秩讨论、或 Sahai–Mitter / 高斯 RD 的斜率讨论里
      找到等价表述，则本车道只保留"在 (18) 前沿上的数值核对 + 口径"，**不得**宣称新定理。

本脚本同时给一张可直接进注记的表：$I\in\{1.5,\dots,4\}$ 上的 $D^*_{\text{free}}$（任何受限族的下界）。
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


def sdp(D):
    P = cp.Variable((n, n), symmetric=True)
    Pi = cp.Variable((n, n), symmetric=True)
    cons = [P >> 0, Pi >> 0, P << Ptil(P), cp.trace(TH @ P) + JC <= D,
            cp.bmat([[P - Pi, P @ A.T], [A @ P, Ptil(P)]]) >> 0]
    prob = cp.Problem(cp.Minimize(-0.5 * cp.log_det(Pi) + 0.5 * LOGDET_W), cons)
    prob.solve(solver=cp.CLARABEL)
    if P.value is None:
        return np.nan, np.nan, 'nan', prob.status
    Pv = sym(P.value)
    S = sym(np.linalg.inv(Pv) - np.linalg.inv(Ptil(Pv)))
    s = np.sort(np.linalg.svd(S, compute_uv=False))[::-1]
    rk = int((s > 1e-3 * s.max()).sum()) if s.size and s.max() > 0 else 0
    return prob.value / ln2, rk, s, prob.status


print('== (1) 精确前沿的 (D, x, DI, rank, 奇异值) 表：x 从 1e-2 到 60 跨四个秩段 ==')
xs = np.unique(np.concatenate([np.geomspace(1e-2, 60.0, 46), JC + np.array([0.05, 0.1, 0.3, 0.7, 1.5, 3.0, 5.0, 8.0]) - JC]))
rows = []
for x in np.sort(xs):
    D = JC + x
    I, rk, s, st = sdp(D)
    if not np.isfinite(I):
        continue
    rows.append((x, D, I, rk, s))
    if len(rows) % 6 == 0 or x == xs.min():
        print(f'   x={x:9.4f}  D={D:9.4f}  DI={I:9.5f}  rank={rk}  sv={s}')
print(f'   共 {len(rows)} 行有效')

print('\n== (2) 逐段斜率：s = dDI/d log2(1/x)，检验 2s vs rank ==')
X = np.array([r[0] for r in rows]); DI = np.array([r[2] for r in rows])
RK = np.array([r[3] for r in rows])
# 中心差分（对 log2(1/x)）
t = np.log2(1.0 / X)
s_local = np.gradient(DI, t)
seg = []
i = 0
while i < len(rows):
    j = i
    while j + 1 < len(rows) and RK[j + 1] == RK[i]:
        j += 1
    if j > i + 1:
        sl, ic = np.polyfit(t[i:j + 1], DI[i:j + 1], 1)
        seg.append((RK[i], X[i], X[j], sl, ic))
    i = j + 1
print('   rank段  x范围                斜率s      2s        2s/rank    [U1]')
for rk, xlo, xhi, sl, ic in seg:
    rel = abs(2 * sl - rk) / rk
    print(f'   {rk}     [{xlo:8.4f},{xhi:8.4f}]   {sl:8.4f}  {2*sl:8.4f}  {rel*100:7.2f}%'
          f'   {"支持" if rel <= 0.05 else (" borderline" if rel <= 0.10 else "不支持")}')

print('\n   [U2] 同秩分段：把每个秩段再切两半，看斜率是否在同秩下漂移 >10%')
devs = []
i = 0
while i < len(rows):
    j = i
    while j + 1 < len(rows) and RK[j + 1] == RK[i]:
        j += 1
    m = j - i + 1
    if m >= 8:
        for half, (a, b) in enumerate(((i, i + m // 2), (i + m // 2, j + 1))):
            sl, ic = np.polyfit(t[a:b], DI[a:b], 1)
            devs.append(abs(2 * sl - RK[i]) / RK[i])
            print(f'     rank={RK[i]} 段{half}  x∈[{X[a]:7.3f},{X[b-1]:7.3f}]  2s={2*sl:7.4f}'
                  f'  与秩比 {(2*sl-RK[i])/RK[i]*100:+6.1f}%')
    i = j + 1
dev = max(devs) if devs else np.nan
print(f'   [U2] 判决：同秩段内 2s 的最大漂移 {dev*100:.1f}%  ⇒ '
      f'{"触发 >10% 否证阈值：\"墙律指数=rank(SNR)\"不成立，此线终止" if dev > 0.10 else "未触发否证"}')

print('\n== (3) 进注记用的表：D*_free(I)（任何受限族在下 I bit 时的代价下界）==')
print('   做法：不用 brentq 直接调 SDP（近地板处求解器会 SolverError），'
      '改为在单调的 (x,DI) 网格上做反插值；每行再用一次正解核验。')
# 扩展网格到 x=1200 以覆盖 I=1.0
xs2 = np.unique(np.concatenate([X, np.geomspace(60.0, 4000.0, 20)]))
DI2 = []
for x in xs2:
    I, rk, s, st = sdp(JC + x)
    DI2.append(I)
DI2 = np.array(DI2)
ok = np.isfinite(DI2)
xs2, DI2 = xs2[ok], DI2[ok]
o = np.argsort(DI2)                  # np.interp 要求 xp 递增：按 I 升序排
xq, Iq = xs2[o], DI2[o]
def invert(I0):
    """反插值起点 + 割线修正（每步都用正解核验；近地板求解器失败就退回上一点）。"""
    D = JC + float(np.interp(I0, Iq, xq))
    best = (np.inf, D, np.nan, 0, 'nan')
    for _ in range(5):
        I, rk, s, st = sdp(D)
        if not np.isfinite(I):
            break
        if abs(I - I0) < best[0]:
            best = (abs(I - I0), D, I, rk, st)
        h = 0.02 * max(D - JC, 1e-3)
        Ip = sdp(D + h)[0]
        Im = sdp(D - h)[0]
        if not (np.isfinite(Ip) and np.isfinite(Im)):
            break
        d = (Ip - Im) / (2 * h)
        if not np.isfinite(d) or abs(d) < 1e-9:
            break
        Dn = D + (I0 - I) / d
        if not np.isfinite(Dn) or Dn <= JC:
            break
        D = Dn
        if abs(I - I0) < 1e-6:
            break
    return best


for I0 in (1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.5):
    if I0 < Iq.min() - 1e-9 or I0 > Iq.max() + 1e-9:
        print(f'   I={I0:5.2f} bit  超出网格 [{Iq.min():.3f},{Iq.max():.3f}]')
        continue
    err, Dg, I, rk, st = invert(I0)
    print(f'   I={I0:5.2f} bit  D*_free≈{Dg:10.4f}  正解核验 DI={I:.6f}（残差 {I-I0:+.2e}）  rank={rk}  {st}')

print('\n== (4) 与 E24 近地板拟合的口径对账 ==')
s_e24 = 6.501 / np.log2(10.0)
print(f'   E24（近地板，x∈[1e-3,1]）：DI = 6.501·log10(1/x)+6.86 ⇒ 每 decade 6.501 bit')
print(f'   换算成对 log2(1/x) 的斜率 s = 6.501/log2(10) = {s_e24:.4f}  ⇒ 2s = {2*s_e24:.4f}')
head = [r for r in rows if r[0] <= 1.0]
print(f'   本表 x<=1 的前 {len(head)} 行秩读数 {[r[3] for r in head]}'
      f'   对照 2s={2*s_e24:.3f}   [U3] 查文献前不得宣称新定理')
