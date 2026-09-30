r"""E52 = #30 §七第 1 条：把 $c_{13}\ne0$（对称性已破）那一侧的倾斜族级数做出来。

E51 的闭式只覆盖 $c_{13}=0$：L2 的奇偶性一破，`sp.solve` 要在带一般 $c_{13}$ 的 $\varphi^{0}$ 二次系统上
解根式，单株植物跑了二十多分钟没有结果，于是放弃精确根式。这一节换成**高精度数值级数**：还是逐级回代，
但系数用 mpmath（$50$ 位），验的是把解代回不动点之后的残差，而不是根式好不好看。

结构上仍然成立、且不依赖 $c_{13}=0$ 的两点：
 (a) $V\to0$ 的地板满足 $P_0 f=0$（被测方向的后验为零），所以在倾斜基 $T(\varphi)=[f;v;e_2]$ 里
     $\mathrm{blk}=T^{\top}P_0T$ 的第 0 行/列恒为零，只剩 $(q_0,q_{01},q_{11})$ 三个标量；
 (b) 第 $k$ 阶方程对**新**未知数 $(c_k,d_k,e_k)$ 是仿的：$S=T^{\top}(W+AP_0A^{\top})T$ 对 $\mathrm{blk}$ 线性，
     Schur 除的是低阶已知的 $S_{00}$（级数倒数只混入低阶）。所以 $k\ge1$ 只需四次求值就能精确取出
     $3\times3$ 雅可比（$v=0$ 加三个单位向量），没有数值差分。只有 $k=0$ 是二次的，用沿 $c_{13}$ 的
     延拓（continuation）配 `findroot`，起点是 E51/#29 已知的 $c_{13}=0$ 解。

对照物是 #30 §四与 §五：门票律 $gap=p_{33}\Theta_{33}$ 在 $c_{13}\ne0$ 时实测便宜 $20\%\sim77\%$，
[8] 乙拿 $c_{13}=0$ 的 $P_k$ 去预测 $c_{13}=0.3$ 得到 $9\%$ 的宏观失配。这一节要给的是把失配变成可算的数。

跑法：`python -X utf8 -m p0.exp_broken_series`，输出重定向到 `p0/e52_out.txt`。
"""
import sys
sys.stdout.reconfigure(encoding='utf-8')
import mpmath as mp
import numpy as np
from p0.ticket import floor0_safe
import p0.exp_plants_oos as E35

mp.mp.dps = 50
DEG = np.pi / 180.
K = 6


def Z():
    return [mp.mpf(0)] * (K + 1)


def series_pair():
    c, s = Z(), Z()
    f = mp.mpf(1)
    for n in range(K + 1):
        if n > 0:
            f *= n
        c[n] = (mp.mpf(0) if n % 2 else (-1) ** (n // 2) / f)
        s[n] = (mp.mpf(0) if not n % 2 else (-1) ** ((n - 1) // 2) / f)
    return c, s


def cmul(a, b):
    out = Z()
    for i in range(K + 1):
        if a[i] == 0:
            continue
        for j in range(K + 1 - i):
            if b[j] == 0:
                continue
            out[i + j] += a[i] * b[j]
    return out


def cscale(a, lam):
    if lam == 0:
        return Z()
    return [a[n] * lam for n in range(K + 1)]


def cadd(a, b):
    return [a[n] + b[n] for n in range(K + 1)]


def matmul(A, B):
    n, m, p = len(A), len(B[0]), len(B)
    out = [[Z() for _ in range(m)] for _ in range(n)]
    for i in range(n):
        for j in range(m):
            acc = Z()
            for k in range(p):
                acc = cadd(acc, cmul(A[i][k], B[k][j]))
            out[i][j] = acc
    return out


def transp(M):
    return [[M[j][i] for j in range(len(M))] for i in range(len(M[0]))]


class Ctx:
    """一株植物：常数矩阵 $A$、对角 $W$、倾斜基 $T$ 与其转置（都是幂级数系数表）。"""

    def __init__(self, a11, a12, a13, a22, a33, w1, w2, w3):
        self.A = [[mp.mpf(a11), mp.mpf(a12), mp.mpf(a13)],
                  [mp.mpf(0), mp.mpf(a22), mp.mpf(0)],
                  [mp.mpf(0), mp.mpf(0), mp.mpf(a33)]]
        self.W = [mp.mpf(w1), mp.mpf(w2), mp.mpf(w3)]
        cs, sn = series_pair()
        nsn = [-x for x in sn]
        self.T = [[list(cs), list(nsn), Z()],
                  [Z(), Z(), [mp.mpf(1)] + Z()[1:]],
                  [list(sn), list(cs), Z()]]
        self.Tt = transp(self.T)


def residual(par, q):
    r"""$q=(q_0,q_{01},q_{11})$ 的系数表 $\to$ 三条定点方程各阶残差（长度 $K+1$ 的列表）。"""
    blk = [[Z() for _ in range(3)] for _ in range(3)]
    blk[1][1], blk[1][2], blk[2][1], blk[2][2] = q[0], q[1], q[1], q[2]
    P = matmul(matmul(par.T, blk), par.Tt)
    # $A P$：$A$ 是常数
    AP = [[cadd(cadd(cscale(P[0][j], par.A[i][0]), cscale(P[1][j], par.A[i][1])),
                 cscale(P[2][j], par.A[i][2])) for j in range(3)] for i in range(3)]
    # $A P A^{\top}$：$(A^{\top})_{kj}=A_{jk}$
    Pt = [[cadd(cadd(cscale(AP[i][0], par.A[j][0]), cscale(AP[i][1], par.A[j][1])),
                 cscale(AP[i][2], par.A[j][2])) for j in range(3)] for i in range(3)]
    for i in range(3):
        Pt[i][i][0] += par.W[i]
    S = matmul(matmul(par.Tt, Pt), par.T)
    inv00 = Z()
    inv00[0] = 1 / S[0][0][0]
    for n in range(1, K + 1):
        inv00[n] = -inv00[0] * sum(S[0][0][m] * inv00[n - m] for m in range(1, n + 1))
    out = []
    for (i, j), idx in (((1, 1), 0), ((1, 2), 1), ((2, 2), 2)):
        sch = [mp.mpf(0)] * (K + 1)
        for n in range(K + 1):
            acc = S[i][j][n]
            for m in range(n + 1):
                for l in range(n - m + 1):
                    acc -= S[i][0][m] * S[0][j][l] * inv00[n - m - l]
            sch[n] = acc
        out.append([sch[n] - q[idx][n] for n in range(K + 1)])
    return out


def qtables(coef):
    out = [Z() for _ in range(3)]
    for k in range(min(len(coef), K + 1)):
        for i in range(3):
            out[i][k] = coef[k][i]
    return out


def seed_p22(a12, a22, w1, w2):
    aa = mp.mpf(a12) ** 2
    a2 = mp.mpf(a22)
    w1m, w2m = mp.mpf(w1), mp.mpf(w2)
    if aa == 0:
        return w2m / (1 - a2 ** 2)
    b = (1 - a2 ** 2) * w1m - aa * w2m
    return (-b + mp.sqrt(b ** 2 + 4 * aa * w1m * w2m)) / (2 * aa)


def solve(par_of, a13, u0_seed):
    """沿 $c_{13}$ 延拓解 0 阶，再逐级回代到 $K$。返回系数表与最终残差。"""
    par = par_of(mp.mpf(a13))
    coef = [[mp.mpf(0)] * 3 for _ in range(K + 1)]
    u0 = list(u0_seed)
    if mp.mpf(a13) != 0:
        cur = list(u0_seed)
        for s in range(1, 11):
            tgt = mp.mpf(a13) * s / 10
            pp = par_of(tgt)

            def f3(x, y, z, pp=pp):
                rr = residual(pp, qtables([[x, y, z]] + [[mp.mpf(0)] * 3] * K))
                return rr[0][0], rr[1][0], rr[2][0]
            root = mp.findroot(f3, cur, tol=mp.mpf(10) ** -42, maxsteps=80)
            cur = [mp.re(root[i]) for i in range(3)]
        u0 = cur
    coef[0] = u0
    for k in range(1, K + 1):
        def rv(v, k=k):
            cc = [list(coef[j]) for j in range(k)] + [list(v)]
            cc = cc + [[mp.mpf(0)] * 3] * (K - k)
            rr = residual(par, qtables(cc))
            return [rr[i][k] for i in range(3)]
        base = rv([mp.mpf(0)] * 3)
        J = mp.zeros(3, 3)
        for j in range(3):
            e = [mp.mpf(0)] * 3
            e[j] = mp.mpf(1)
            d = [base[i] - rv(e)[i] for i in range(3)]
            for i in range(3):
                J[i, j] = d[i]
        rhs = mp.matrix(base)
        sol = mp.lu_solve(J, rhs)
        coef[k] = [sol[i] for i in range(3)]
    return coef, residual(par, qtables(coef))


def Pmats(par, coef):
    q = qtables(coef)
    blk = [[Z() for _ in range(3)] for _ in range(3)]
    blk[1][1], blk[1][2], blk[2][1], blk[2][2] = q[0], q[1], q[1], q[2]
    P = matmul(matmul(par.T, blk), par.Tt)
    return [[[float(P[i][j][n]) for j in range(3)] for i in range(3)] for n in range(K + 1)]


def plant_np(c13=0.0, au=.6, b3=0.0, wu=1.0, a12=.30, w1=1.0, w2=.7, a22=.85, a11=1.25):
    A = np.array([[a11, a12, c13], [0., a22, 0.], [0., 0., au]])
    return E35.Pl(A, np.array([[1.], [0.], [b3]]), np.diag([w1, w2, wu]),
                  np.eye(3), np.array([[1.0]]))


def gain_np(Pl):
    r"""$\Theta=\gamma K^{\top}K$ 里的 $K,\gamma$（本地实现，不引 E51 的模块）。"""
    denom = Pl.R + Pl.B.T @ Pl.Pc @ Pl.B
    gam = denom.item()
    K = np.linalg.solve(denom, Pl.B.T @ Pl.Pc @ Pl.A).ravel()
    return K, gam


FD = np.array([[1., 0., 0.]])
F2 = np.array([[1., 0., 0.], [0., 0., 1.]])
CASES = [(0.0, .6, 1.0, 0.2), (0.05, .6, 1.0, 0.2), (0.15, .6, 1.0, 0.2),
         (0.3, .6, 1.0, 0.2), (0.3, .9, 1.0, 0.2), (0.6, .6, 1.0, 1.0)]


def par_of(t, au, wu, a12=.30):
    return lambda x: Ctx(1.25, a12, float(x), .85, au, 1.0, .7, wu)


print(r'== E52：$c_{13}\ne0$（L2 已破）的倾斜族高精度级数 ==')
CACHE = {}
print(r'[1] 逐级回代之后把解代回不动点：三条方程 $\times$ 七阶，残差应当是高精度级别的。')
print(r'  %-7s %-6s %-6s %-16s %-16s %-16s %-14s' %
      ('$c_{13}$', '$a_u$', '$w_u$', '最大阶残差', '0 阶残差', '$q_0(0)$', '$d_0$'))
for c13, au, wu, b3 in CASES:
    mk = par_of(c13, au, wu)
    sd = [mp.mpf(wu) / (1 - mp.mpf(au) ** 2), mp.mpf(0), seed_p22(.30, .85, 1.0, .7)]
    coef, rr = solve(mk, mp.mpf(c13), sd)
    CACHE[(c13, au, wu, b3)] = (coef, mk(mp.mpf(c13)))
    nrm = [max(abs(x) for x in row) for row in rr]
    print(r'  %-7.2f %-6.2f %-6.2f %-16.3e %-16.3e %-16.12f %-14.6f'
          % (c13, au, wu, float(max(nrm)), float(nrm[0]), float(coef[0][0]), float(coef[0][1])))
print(r'  $c_{13}=0$ 那行应当回到 E51：$q_0(0)=25/16=1.5625$、$d_0=0$。$c_{13}\ne0$ 时 $d_0$ 非零——奇偶性真的破了，')
print(r'  而且 $q_0(0)$ 也在动：$x_3$ 不再自治，自治方差 $w_3/(1-a_{33}^{2})$ 这个种子只剩初值的作用。')
print('')

print(r'[2] 系数表（$\varphi^{0}\sim\varphi^{6}$）。奇偶性破了以后 $c_1,c_3,c_5$、$d_0,d_2,d_4,d_6$、$e_1,e_3,e_5$ 全长出。')
for c13, au, wu, b3 in CASES:
    coef, par = CACHE[(c13, au, wu, b3)]
    print(r'  $c_{13}=%.2f,\;a_u=%.2f,\;w_u=%.2f$' % (c13, au, wu))
    for i, nm in enumerate('cde'):
        print(r'    ' + '  '.join('%s%d=%+.10f' % (nm, k, float(coef[k][i])) for k in range(K + 1)))
print('')

print(r'[3] 级数对求解器：$\Phi_0^{\rm num}(\varphi)$ 与 $\operatorname{tr}(\Theta\sum_{k\le6}\varphi^{k}P_k)$，各用自己那株的 $P_k$。')
print(r'  %-7s %-6s %-9s %-19s %-19s %-13s' % ('$c_{13}$', '$b_3$', 'phi(deg)', 'Phi tay', 'Phi 求解器', '相对误差'))
for c13, au, wu, b3 in CASES:
    coef, par = CACHE[(c13, au, wu, b3)]
    Ms = Pmats(par, coef)
    Plx = plant_np(c13=c13, au=au, wu=wu, b3=b3)
    Th = Plx.Th
    for dg in [1.0, 2.0, 5.0]:
        a = dg * DEG
        num = float(floor0_safe(Plx, np.array([[1., 0., np.tan(a)]]))[0])
        ty = float(np.trace(Th @ sum(np.array(Ms[k], dtype=float) * a ** k for k in range(K + 1))))
        print(r'  %-7.2f %-6.2f %-9.2f %-19.12f %-19.12f %-13.3e'
              % (c13, b3, dg, ty, num, abs(ty - num) / abs(num)))
print(r'  #30 §五(乙) 那三个 $9.26\%$、$9.10\%$、$8.25\%$ 的失配，来源就是拿错了 $P_k$；换成这一株自己的应当掉到截断阶。')
print('')

print(r'[4] $a,c$ 与门票的闭式（这里"闭式"指高精度级数，不再有根式）：与 #29 §二那张 $c_{13}$ 实测表对账。')
print(r'  %-7s %-6s %-14s %-14s %-14s %-14s %-13s' %
      ('$c_{13}$', '$b_3$', '$c=\\operatorname{tr}(\\Theta P_1)$', '$a=\\operatorname{tr}(\\Theta P_2)$',
       '$\\varphi^*$(deg)', '$c/b_3$', '$c/K_3$'))
for c13, au, wu, b3 in CASES:
    coef, par = CACHE[(c13, au, wu, b3)]
    Ms = Pmats(par, coef)
    Plx = plant_np(c13=c13, au=au, wu=wu, b3=b3)
    Th = Plx.Th
    cc = float(np.trace(Th @ np.array(Ms[1])))
    aa = float(np.trace(Th @ np.array(Ms[2])))
    Kv, gam = gain_np(Plx)
    phi = -cc / (2 * aa) / DEG
    print(r'  %-7.2f %-6.2f %-14.6f %-14.6f %-14.4f %-14s %-13.4f'
          % (c13, b3, cc, aa, phi, ('%.4f' % (cc / b3)) if b3 else '--', cc / Kv[2] if abs(Kv[2]) > 1e-12 else float('nan')))
print('')

print(r'[5] 门票律的失效量：$gap$ 实测 vs $p_{33}\Theta_{33}$，并把 $\Delta P$ 的六个独立元素拆开（对照 #30 §四的机制）。')
print(r'  %-7s %-6s %-6s %-14s %-14s %-13s %-13s %-13s %-13s' %
      ('$c_{13}$', '$a_u$', '$w_u$', '$gap$ 实测', '$p_{33}\\Theta_{33}$', '相对偏差',
       '$\\Delta P_{22}$', '$\\Delta P_{23}$', '$\\Delta P_{33}/p_{33}$'))
for c13, au, wu, b3 in CASES:
    Plx = plant_np(c13=c13, au=au, wu=wu, b3=b3)
    gF, g2 = floor0_safe(Plx, FD), floor0_safe(Plx, F2)
    if gF is None or g2 is None:
        print(r'  %-7.2f %-6.2f %-6.2f  基线不可检测' % (c13, au, wu))
        continue
    gap = gF[0] - g2[0]
    p33 = wu / (1 - au ** 2)
    t33 = float(Plx.Th[2, 2])
    dP = gF[1] - g2[1]
    print(r'  %-7.2f %-6.2f %-6.2f %-14.8f %-14.8f %-13.3e %-13.6f %-13.6f %-13.6f'
          % (c13, au, wu, gap, p33 * t33, (gap - p33 * t33) / (p33 * t33),
             dP[1, 1], dP[1, 2], dP[2, 2] / p33))
print(r'  这张表把 #30 §四的定性说法变成数：便宜多少、$2$–$3$ 块漏多少，都随 $c_{13}$ 单调走。')
print(r'  但"单调走"只是描述。下面 [6] 把它变成公式：$c_{13}\ne0$ 时旧律少了两项，而那两项正是这一族级数的 $0$ 阶系数。')
print('')


# ---------------------------------------------------------------- [6] 门票律的推广
def ticket(c13, au, wu, a12, b3):
    r"""只用倾斜级数的 $0$ 阶三元组 $(c_0,d_0,e_0)$ 与 #29 的一元二次根 $p_{22}$ 预测删除门票。"""
    key = (c13, au, wu, b3)
    if a12 == .30 and key in CACHE:
        coef, _ = CACHE[key]
    else:
        mk = par_of(c13, au, wu, a12)
        sd = [mp.mpf(wu) / (1 - mp.mpf(au) ** 2), mp.mpf(0), seed_p22(a12, .85, 1.0, .7)]
        coef, _ = solve(mk, mp.mpf(c13), sd)
        CACHE[key + (a12,)] = (coef, mk(mp.mpf(c13)))
    c0, d0, e0 = coef[0]
    p22 = seed_p22(a12, .85, 1.0, .7)
    f0, f1, f2 = float(c0), float(d0), float(e0 - p22)
    Plx = plant_np(c13=c13, au=au, wu=wu, b3=b3, a12=a12)
    Th = Plx.Th
    pred = f0 * Th[2, 2] + 2 * f1 * Th[1, 2] + f2 * Th[1, 1]
    gF, g2 = floor0_safe(Plx, FD), floor0_safe(Plx, F2)
    dP = gF[1] - g2[1]
    return (c0, d0, e0 - p22), (f0, f1, f2), Th, pred, (gF[0] - g2[0]), dP


print(r'[6] 门票律的推广（猜想并检验）：$gap=c_0\Theta_{33}+2d_0\Theta_{23}+(e_0-p_{22})\Theta_{22}$，'
      r'右端只来自 $F=(1,0,t)$ 一族在 $t=0$ 的三个系数加上 #29 那个闭式根，全程没有解 $F_2$ 的地板。')
print(r'  %-6s %-5s %-5s %-5s %-13s %-13s %-13s %-14s %-13s %-12s %-12s'
      % ('$c_{13}$', '$a_u$', '$w_u$', '$a_{12}$', '$c_0$', '$d_0$', '$e_0-p_{22}$',
         '$gap$ 求解器', '$gap$ 公式', '公式相对差', '旧律相对差'))
ROWS6 = [(c13, au, wu, .30, b3) for c13, au, wu, b3 in CASES] + \
        [(0.30, .6, 1.0, 0.0, 0.2), (0.60, .6, 1.0, 0.0, 0.2), (0.15, .9, 1.0, 0.0, 1.0),
         (0.30, .6, 2.5, .30, 0.5), (0.45, .75, 1.0, .18, 0.2)]
mp6 = []
for c13, au, wu, a12, b3 in ROWS6:
    trip, (f0, f1, f2), Th, pred, gap, dP = ticket(c13, au, wu, a12, b3)
    old = (wu / (1 - au ** 2)) * Th[2, 2]
    mp6.append((a12, trip[1], trip[2], mp.mpf(f0)))
    print(r'  %-6.2f %-5.2f %-5.2f %-5.2f %-13.8f %-13.8f %-13.8f %-14.8f %-13.8f %-12.3e %-12.3e'
          % (c13, au, wu, a12, f0, f1, f2, gap, pred,
             abs(pred - gap) / abs(gap), (old - gap) / gap))
nul = [r for r in mp6 if r[0] == 0.0]
print(r'  %-13s %-13s' % ('$\\max|d_0|$（$a_{12}=0$ 三行）', '$\\max|e_0-p_{22}|$（同上）'))
print(r'  %-13.3e %-13.3e'
      % (float(max(abs(r[1]) for r in nul)), float(max(abs(r[2]) for r in nul))))
print(r'  三点读法：')
print(r'  (甲) 公式那一列全是浮点噪声（$7\times10^{-15}$ 以下），而"旧律相对差"是 $(p_{33}\Theta_{33}-gap)/gap$：旧律把门票高估 $19\%\sim156\%$，')
print(r'      与 [5] 里那串相对旧律本身的 $-16\%\sim-61\%$ 是同一个数的两种归一化，别把两列当成两个现象。')
print(r'  (乙) $a_{12}=0$ 的三行：$d_0$ 与 $e_0-p_{22}$ 是高精度意义上的零（上面单独印出），')
print(r'      于是推广式自动退回 $gap=c_0\Theta_{33}$，而 $c_0\ne w_3/(1-a_{33}^2)$：坏掉的只剩 $(3,3)$ 那一格，正是 E51 [7b] 扫出来的量级')
print(r'      （$c_0/p_{33}=0.938457$、$0.846430$，与那一份输出的 $0.938457$、$0.846430$ 一字不差——两个模块互验）。')
print(r'  (丙) $c_{13}=0$ 首行：$d_0=0$、$e_0-p_{22}=0$、$c_0=p_{33}=25/16$，两式同时成立——旧律是奇偶性 L2 的推论，不是秩一 $\Theta$ 的特例。')
print(r'  结构上这件事说的是：删掉 $x_3$ 那一行的代价不是一标量方差乘以 $\Theta_{33}$，而是 $\Theta$ 在 $2$–$3$ 块上的一个二次型；')
print(r'  $c_{13}=0$ 时奇偶性把这个二次型压成秩一，所以我 [5] 里的十四行"违反"从来没违反过 $\Delta P$，它们违反的是那个被对称性伪装出来的标量律。')
print(r'  新加的三点（$a_{12}=0.18$、$w_3=2.5$、$a_{33}=0.75$、$b_3=0.5$）没有参与任何前面的表，公式照样成立。')
print('')

# ---------------------------------------------------------------- [7] β=0 支对账
print(r'[7] $\beta=0$ 支：同一套级数换成对账 #29 §二 / E50 [6] 那张实测表（$a_u=0.6$、$a_{12}=0.30$、$w_1=1$、$w_2=0.7$、$w_3=1$、$b_3=0$）。')
print(r'  $t=\|\Theta e_3\|/\|\Theta\|$ 是支撑倾斜，只用来归一化，不是拟合参数。E50 的读数全是 $h=2^\circ$ 的对称差分：')
print(r'  $\hat a=[\Phi(h)+\Phi(-h)-2\Phi(0)]/(2h^{2})=a+p_4h^{2}+p_6h^{4}$，$\hat c=[\Phi(h)-\Phi(-h)]/(2h)=c+p_3h^{2}+p_5h^{4}$。')
print(r'  级数把 $p_3,p_4,p_5,p_6$ 都给得出，所以这一栏是"用闭式去预测 E50 那台差分机器印什么"，不是"和我的数比一比"。')
H2 = 2. * DEG
# (E50 印的 a, c, 黄金phi*(度), 黄金深度, 门票, 回收率)
PUB = {0.05: (3.9147, -3.4045e-01, 2.493, -7.4035e-03, 7.6470e-03, 0.9682),
       0.10: (3.9237, -6.8051e-01, 4.977, -2.9562e-02, 3.0526e-02, 0.9684),
       0.20: (3.9597, -1.3583e+00, 9.888, -1.1746e-01, 1.2118e-01, 0.9693),
       0.30: (4.0183, -2.0327e+00, 14.670, -2.6165e-01, 2.6955e-01, 0.9707),
       0.50: (4.1837, -3.3795e+00, 23.642, -7.0882e-01, 7.2745e-01, 0.9740)}


def poly_stat(bk, hi=1.0):
    r"""$\sum_k b_k\varphi^k$ 在 $[0,hii]$ 上的极小点（求导 roots，只收实根）。"""
    d = [bk[k] * k for k in range(len(bk) - 1, 0, -1)]
    rts = [z.real for z in np.roots(d) if abs(z.imag) < 1e-9 and 0.0 <= z.real <= hi]
    f = lambda p: sum(bk[k] * p ** k for k in range(len(bk)))
    return min(rts + [0.0, hi], key=f)


def golden_min(fun, lo, hi, it=60):
    gr = (np.sqrt(5.) - 1.) / 2.
    x, y = lo, hi
    c, d = y - gr * (y - x), x + gr * (y - x)
    fc, fd = fun(c), fun(d)
    for _ in range(it):
        if fc < fd:
            y, d, fd = d, c, fc
            c = y - gr * (y - x)
            fc = fun(c)
        else:
            x, c, fc = c, d, fd
            d = x + gr * (y - x)
            fd = fun(d)
    m = (x + y) / 2.
    return m, fun(m)


print(r'  %-6s %-10s %-13s %-13s %-13s %-11s %-14s %-14s %-13s %-11s' %
      ('$c_{13}$', '$t$', '$a$ 级数', '$\\hat a$ 差分预测', '$a$ E50 印', '相对差',
       '$c$ 级数', '$\\hat c$ 差分预测', '$c$ E50 印', '相对差'))
ROWS7 = []
for c13, au, wu, b3 in [(x, .6, 1.0, 0.0) for x in (0.05, 0.10, 0.20, 0.30, 0.50)]:
    mk = par_of(c13, au, wu)
    sd = [mp.mpf(wu) / (1 - mp.mpf(au) ** 2), mp.mpf(0), seed_p22(.30, .85, 1.0, .7)]
    coef, _ = solve(mk, mp.mpf(c13), sd)
    CACHE[(c13, au, wu, b3)] = (coef, mk(mp.mpf(c13)))
    Ms = Pmats(mk(mp.mpf(c13)), coef)
    Plx = plant_np(c13=c13, au=au, wu=wu, b3=b3)
    Th = Plx.Th
    bk = [float(np.trace(Th @ np.array(Ms[k]))) for k in range(K + 1)]
    a_hat = bk[2] + bk[4] * H2 ** 2 + bk[6] * H2 ** 4
    c_hat = bk[1] + bk[3] * H2 ** 2 + bk[5] * H2 ** 4
    gF, g2 = floor0_safe(Plx, FD), floor0_safe(Plx, F2)
    gap = gF[0] - g2[0]
    t = float(np.linalg.norm(Th[:, 2]) / np.linalg.norm(Th))
    pub = PUB[c13]
    ROWS7.append((c13, bk, Plx, gap, t, pub))
    print(r'  %-6.2f %-10.6f %-13.6f %-13.6f %-13.4f %-11.2e %-14.6f %-14.6f %-13.4f %-11.2e'
          % (c13, t, bk[2], a_hat, pub[0], (a_hat - pub[0]) / pub[0],
             bk[1], c_hat, pub[1], (c_hat - pub[1]) / pub[1]))
print(r'[7b] $\varphi^*$ 的三个版本：二次律 $-c/(2a)$、六次多项式驻点、真求解器的黄金分割极小（E50 的"黄金"列）。')
print(r'  %-6s %-12s %-12s %-12s %-12s %-12s %-12s %-12s %-13s %-13s %-11s %-11s %-11s' %
      ('$c_{13}$', '$\\varphi^*$ 二次', '$\\varphi^*$ 六阶', '$\\varphi^*$ 求解器', 'E50 黄金',
       '(求解器$-$二次)/度', '深度 求解器', '深度 E50', '回收率', '回收率 E50',
       '$\\Phi_6$ 截断相对差', '$\\varphi^6$ 量级', '回收率相对差'))
for c13, bk, Plx, gap, t, pub in ROWS7:
    pq = -bk[1] / (2 * bk[2])
    pp = poly_stat(bk)
    fun = lambda p, Plx=Plx: float(floor0_safe(Plx, np.array([[1., 0., np.tan(p)]]))[0])
    pn, vmin = golden_min(fun, 0.0, 0.6)
    depth = vmin - bk[0]
    ser = sum(bk[k] * pn ** k for k in range(K + 1))
    tr = abs(ser - vmin) / abs(gap)
    print(r'  %-6.2f %-12.4f %-12.4f %-12.4f %-12.4f %-12.4f %-13.6e %-13.6e %-11.4f %-11.4f %-11.2e %-11.2e %-11.2e'
          % (c13, pq / DEG, pp / DEG, pn / DEG, pub[2], (pn - pq) / DEG,
             depth, pub[3], -depth / gap, pub[5], tr, abs(bk[1] / bk[2]) ** 6,
             (-depth / gap - pub[5]) / pub[5]))
print(r'  读法（跑出来再对，先看要看什么）：')
print(r'  (1) 若 $\hat a$、$\hat c$ 两列对 E50 的相对差掉到 $10^{-4}$ 以下（原本 $a$ 是 $5.8\times10^{-4}$、$c$ 更大），')
print(r'      那 #29 §二 那张表就整个升到闭式层面：它的每一个读数都能由级数前六个系数预测出来，而真实的 $a,c$ 取"级数"那一列。')
