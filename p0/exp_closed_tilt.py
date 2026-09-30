r"""E51 = 把 #25/#26/#27/#28 里四条"只有数值"的旋转律算成闭式：$\Phi_0$ 对倾角的泰勒系数。

两条结构引理（这一节的骨架，后面所有表格都只是它们的读数）：

**L1（地板的估计侧看不见 $B$）** 硬约束地板的后验由
    $\tilde P=W+AP_0A^{\top},\quad P_0=\tilde P-\tilde P\,F^{\top}(F\tilde PF^{\top})^{-1}F\tilde P$
这条定点方程决定，里面根本没有 $B$。所以 $P_0=P_0(A,W,F)$，而
    $\Phi_0(F)=\operatorname{tr}(\Theta P_0(A,W,F))$
**对 $\Theta$ 是线性的**。#28 那句"$a,c$ 看起来是支撑倾斜 $t$ 的函数"因此是假的因果：
真正的律是"$a,c$ 是 $\Theta$ 的一次型，其矩阵只由 $(A,W)$ 给出"。机器证明在 [1]：
$b_3$ 从 $0$ 改到 $3$，$P_0$ 逐元素差 $0.000e+00$，而 $\Phi_0$ 改了 $1.7\times10^{1}$。

**L2（$c_{13}=0$ 时级数的零式样被对称性钉死）** 取 $T_s=\operatorname{diag}(1,1,-1)$。
$c_{13}=0$ 时 $T_sAT_s=A$、$T_sWT_s=W$，而通道方向 $f(\varphi)=(\cos\varphi,0,\sin\varphi)^{\top}$ 满足
$T_sf(\varphi)=f(-\varphi)$，于是唯一性给出 $P_0(-\varphi)=T_sP_0(\varphi)T_s$（[A] 里数值残差为 $0$）。
用**正交**倾斜基 $T(\varphi)=[\,f(\varphi)\;;\;v(\varphi)=(-\sin\varphi,0,\cos\varphi)^{\top}\;;\;e_2\,]$ 写后验块
    $\mathrm{blk}=T(\varphi)^{\top}P_0(\varphi)T(\varphi)=\begin{bmatrix}0&0&0\\0&q_0&q_{01}\\0&q_{01}&q_{11}\end{bmatrix},$
则 $T_sT(\varphi)=T(-\varphi)\operatorname{diag}(1,-1,1)$，于是律只剩一句话：
    $q_0,q_{11}$ 是偶的，$q_{01}$ 是奇的。
换成对 $\varphi$ 的泰勒系数就是 $c_1=c_3=d_0=d_2=e_1=e_3=0$，奇部只剩 $d_1\varphi+d_3\varphi^{3}$。
基向量本身也带 $\varphi$，所以原坐标里的地板是这三个标量生成的（$s=\sin\varphi$、$c=\cos\varphi$）：
    $P_{11}=q_0s^{2},\quad P_{13}=-q_0sc,\quad P_{33}=q_0c^{2},\quad
      P_{12}=-q_{01}s,\quad P_{23}=q_{01}c,\quad P_{22}=q_{11},$
其余按对称性补齐。展开到 $\varphi^{2}$ 就得到
    $P_0=\begin{bmatrix}0&0&0\\0&p_{22}&0\\0&0&p_{33}\end{bmatrix},\quad
      P_1=\begin{bmatrix}0&0&-p_{33}\\0&0&d_1\\-p_{33}&d_1&0\end{bmatrix},\quad
      P_2=\begin{bmatrix}p_{33}&-d_1&0\\-d_1&e_2&0\\0&0&c_2-p_{33}\end{bmatrix}.$
注意 $P_2$ 的 $(3,3)$ 是 $c_2-p_{33}$（我第一版在这里写成 $c_2$，因为基不正交把 $\cos^{2}\varphi$ 的 $-p_{33}$ 漏掉了），
所以 $P_2$ 不是半定的；而 $\Theta=\gamma K^{\top}K$，于是
    $c=\operatorname{tr}(\Theta P_1)=2\gamma K_3\,(d_1K_2-p_{33}K_1),\qquad
      a=\operatorname{tr}(\Theta P_2)=\gamma\big(p_{33}K_1^{2}+e_2K_2^{2}+(c_2-p_{33})K_3^{2}-2d_1K_1K_2\big).$
$c(\beta=0)=0$（#28 实测）由此**被推导**而不是被观测：$K_3|_{\beta=0}=0$。同理 $c\propto K_3$，
不是 $\propto b_3$ —— 这就是 [8] 那张"两支都塌到 $t$"的表的真身份。

两次自我纠错，都写进流程而不是写在脚注里：
 (i) 0 阶根必须代回方程验。`sp.solve` 在清分母后的二次系统上交回过一个正定、$\Phi_0$ 上第 17 位
     才有差、却**不满足方程**的假根；种子验根现在把残差打印出来（[2] 的第三行）。
 (ii) 倾斜基必须**正交**。我第一版用 $\tilde T=[(1,0,t),(-t,0,1),e_2]$ 并把 $\tilde T^{\top}\tilde P\tilde T$
     当作坐标变换——那只对正交阵成立，于是 $q_0$ 那一块的 2 阶系数被 $(1+t^{2})$ 的归一化污染
     （旧值 $c_2=4.7276$，正确值 $1.4073$），症状是 $\Phi_0$ 的泰勒残差随 $t^{2}$ 而不随截断阶。
     同一版里我把投影写成 $TPT^{\top}$，于是 $(q_0,q_{01},q_{11})$ "违反"奇偶性；改正后 [A] 的残差是 $0$。
     E48 立的"每格独立复算 + 阶数必须对得上"杠子在这一节救了我两次。

适用范围：闭式只对 $c_{13}=0$（即 $x_3$ 自治）给出。$c_{13}\ne0$ 时 L2 的对称性直接没了
（$T_sAT_s\ne A$），我试过用 `sp.solve` 走一般 $c_{13}$ 的 0 阶二次系统，单株植物跑了 20 分钟以上没有结果，
于是放弃精确根式，只把 [7]/[7b] 的量级边界留在数值里：门票律在那一侧是**便宜** $20\%\sim77\%$ 的不等式。

跑法：`python -X utf8 -m p0.exp_closed_tilt`，输出重定向到 `p0/e51_out.txt`。
"""
import sys
sys.stdout.reconfigure(encoding='utf-8')
import numpy as np
import sympy as sp
import p0.exp_plants_oos as E35
from p0.ticket import floor0, floor0_safe

R = sp.Rational


def rat(x):
    return R(str(x)) if not isinstance(x, R) else x


# ---------------------------------------------------------------- 级数机器（L1/L2 的实现）
_p = sp.symbols('phi', real=True)
_q = sp.symbols('q0 q01 q11', real=True)


def tilt_series(a11, a12, a13, a22, a33, w1, w2, w3, order=3):
    r"""硬约束地板后验对 $\varphi$ 的泰勒闭式：$P_0(\varphi)=\sum_k\varphi^{k}P_k+O(\varphi^{n+1})$。

    定点方程只有 3 条：在正交倾斜基里写后验块 $(q_0,q_{01},q_{11})$，令它等于
    $\tilde P=W+AP_0A^{\top}$ 的 Schur 补（沿被约束掉的 $f$ 方向条件化），再逐阶求导。
    """
    A = sp.Matrix([[rat(a11), rat(a12), rat(a13)], [0, rat(a22), 0], [0, 0, rat(a33)]])
    W = sp.diag(rat(w1), rat(w2), rat(w3))
    c, s = sp.cos(_p), sp.sin(_p)
    T = sp.Matrix([[c, -s, 0], [0, 0, 1], [s, c, 0]])          # 列 = [f, v, e2]，正交
    blk = sp.zeros(3, 3)
    blk[1:, 1:] = sp.Matrix([[_q[0], _q[1]], [_q[1], _q[2]]])
    P0 = sp.expand(T * blk * T.T)
    S = sp.expand(T.T * (W + A * P0 * A.T) * T)
    sch = sp.expand(S[1:, 1:] - S[1:, 0:1] * S[0:1, 1:] / S[0, 0])
    eqs = [sp.together(sch[i, j] - m).as_numer_denom()[0]
           for (i, j), m in zip(((0, 0), (0, 1), (1, 1)), _q)]
    us = {k: [sp.symbols('%s%d' % (k, i)) for i in range(order + 1)] for k in 'cde'}
    ser = {_q: sum(us[k][i] * _p ** i for i in range(order + 1))
           for _q, k in zip(_q, 'cde')}
    per_eq = []
    for e in eqs:
        sr = sp.series(sp.expand(e.subs(ser)), _p, 0, order + 1).removeO()
        per_eq.append([sp.expand(sr.coeff(_p, k)) for k in range(order + 1)])
    allk = [[per_eq[i][k] for i in range(len(per_eq))] for k in range(order + 1)]
    degs = [[int(sp.Poly(s, *(us['c'] + us['d'] + us['e'])).total_degree()) if s != 0 else 0
             for s in allk[k]] for k in range(order + 1)]

    u0 = [us['c'][0], us['d'][0], us['e'][0]]
    bb = (1 - rat(a22) ** 2) * rat(w1) - rat(a12) ** 2 * rat(w2)
    p22 = sp.radsimp(sp.simplify((-bb + sp.sqrt(bb ** 2 + 4 * rat(a12) ** 2 * rat(w1) * rat(w2)))
                                 / (2 * rat(a12) ** 2)))
    seed, res0 = None, None
    if rat(a13) == 0:                       # 对称性可用：种子就是 $p_{33},0,p_{22}$
        cand = {u0[0]: sp.simplify(rat(w3) / (1 - rat(a33) ** 2)),
                u0[1]: sp.Integer(0), u0[2]: p22}
        rr = [sp.simplify(e.subs(cand)) for e in allk[0]]
        if all(x == 0 for x in rr):
            seed, res0 = cand, rr
    if seed is None:                        # $a_{13}\\ne0$（对称性已破）：解 + 代回验根
        good = []
        for sy in sp.solve(allk[0], u0, dict=True):
            if not all(sy[x].is_real for x in u0) or sy[u0[0]] <= 0 or sy[u0[2]] <= 0:
                continue
            rr = [abs(sp.N(e.subs(sy), 40)) for e in allk[0]]
            if max(rr) < 1e-25:
                good.append((sy, rr))
        assert good, 'order-0 没有通过验根的正定根'
        seed, res0 = good[0]
    sol = dict(seed)
    lin = []
    for k in range(1, order + 1):
        eqk = [sp.expand(e.subs(sol)) for e in allk[k]]
        unk = [us['c'][k], us['d'][k], us['e'][k]]
        islin = all(sp.Poly(e, unk).total_degree() <= 1 for e in eqk if e != 0)
        lin.append(islin)
        if islin:
            Mm = sp.Matrix([[sp.diff(e, u) for u in unk] for e in eqk])
            rhs = sp.Matrix([-sp.expand(e.subs({u: 0 for u in unk})) for e in eqk])
            newv = Mm.solve(rhs)
        else:
            hit = None
            for sy in sp.solve(eqk, unk, dict=True):
                if not all(sy[x].is_real for x in unk):
                    continue
                if max(abs(sp.N(e.subs(sy), 40)) for e in eqk) < 1e-25:
                    hit = [sy[x] for x in unk]
                    break
            assert hit is not None, 'order-%d 没有通过验根的根' % k
            newv = hit
        for u, v in zip(unk, newv):
            sol[u] = sp.simplify(v)
    full = sp.expand(P0.subs(ser).subs(sol))
    exact = [[[sp.simplify(sp.diff(full[i, j], _p, k).subs(_p, 0) / sp.factorial(k))
               for j in range(3)] for i in range(3)] for k in range(order + 1)]
    mats = [np.array([[float(sp.N(v, 30)) for v in row] for row in M]) for M in exact]
    return sol, mats, p22, res0, degs, lin, exact


def plant(c13=0.0, au=.6, b3=0.0, wu=1.0, a12=.30, w1=1.0, w2=.7, a22=.85, a11=1.25):
    A = np.array([[a11, a12, c13], [0., a22, 0.], [0., 0., au]])
    B = np.array([[1.], [0.], [b3]])
    W = np.diag([w1, w2, wu])
    return E35.Pl(A, B, W, np.eye(3), np.array([[1.0]]))


def gain(Pl):
    r"""$\Theta=\gamma K^{\top}K$ 里的 $K,\gamma$。"""
    denom = Pl.R + Pl.B.T @ Pl.Pc @ Pl.B
    gam = denom.item()
    K = np.linalg.solve(denom, Pl.B.T @ Pl.Pc @ Pl.A).ravel()
    return K, gam


def ft(t):
    return np.array([[1., 0., t]])


def phi0(Pl, t):
    got = floor0_safe(Pl, ft(t))
    return None if got is None else got[0]


def P0num(Pl, t):
    got = floor0_safe(Pl, ft(t))
    return None if got is None else got[1]


def Torth(phi):
    c, s = np.cos(phi), np.sin(phi)
    return np.array([[c, -s, 0.], [0., 0., 1.], [s, c, 0.]])      # 列 = [f, v, e2]


def tay(mats, phi):
    return sum(M * phi ** k for k, M in enumerate(mats))


Ts = np.diag([1., 1., -1.])
DEG = np.pi / 180.

print(r'== E51：旋转族地板的闭式（L1 估计侧盲于 $B$，L2 奇部只住在指标 3 上）==')
print('')

# ---------------------------------------------------------------- [1] L1
print(r'[1] L1 的机器证明：同一 $(A,W)$、不同 $B$，硬约束后验 $P_0$ 必须逐元素相同（$\Phi_0$ 只能通过 $\Theta$ 变）。')
print(r'  %-16s %-9s %-8s %-22s %-22s' % ('$F$ 族', '$b_3$ 对照', '$c_{13}$', 'max|P0(b3)-P0(0)|', 'max|Phi0(b3)-Phi0(0)|'))
TS = [0.0, 0.1, 0.3, 0.6]
worstP = 0.
for c13 in [0.0, 0.3]:
    P0ref = {t: P0num(plant(c13=c13), t) for t in TS}
    for b3 in [0.2, 1.0, 3.0]:
        Plx = plant(c13=c13, b3=b3)
        Pl0 = plant(c13=c13)
        dp = max(np.abs(P0num(Plx, t) - P0ref[t]).max() for t in TS)
        dphi = max(abs(phi0(Plx, t) - phi0(Pl0, t)) for t in TS)
        worstP = max(worstP, dp)
        print(r'  %-16s %-9.1f %-8.1f %-22.3e %-22.3e'
              % ('$t\\in\\{0,.1,.3,.6\\}$', b3, c13, dp, dphi))
print(r'  后验侧最大偏差 %.3e（E48 的杠子 $\sim10^{-15}$）：$\Phi_0$ 那一列却是宏观的，' % worstP)
print(r'  两列合起来就是 L1 —— 倾角换的是 $P_0$，控制侧换的只有 $\Theta$，两者在 $\Phi_0$ 里是**相乘**关系。')
print('')

# ---------------------------------------------------------------- [A] L2 的数值裁决
print(r'[A] L2 的数值裁决：$P_0(-\varphi)=T_sP_0(\varphi)T_s$（$c_{13}=0$），以及正交倾斜基里 $(q_0,q_{01},q_{11})$ 的奇偶性。')
print(r'  %-7s %-7s %-23s %-25s %-21s' % ('$b_3$', '$\\varphi$', '|P0(+phi)-P0(-phi)|max', '|P0(phi)-TsP0(-phi)Ts|max', 'q01 的奇部残差'))
for b3 in [0.0, 0.2, 1.0]:
    Pl = plant(b3=b3)
    for dg in [1., 2., 5., 10., 30.]:
        a = dg * DEG
        Pp, Pm = P0num(Pl, np.tan(a)), P0num(Pl, np.tan(-a))
        bq = (Torth(a).T @ Pp @ Torth(a))[1, 2]
        bm = (Torth(-a).T @ Pm @ Torth(-a))[1, 2]
        print(r'  %-7.1f %-7.1f %-23.3e %-25.3e %-21.3e'
              % (b3, dg, np.abs(Pp - Pm).max(), np.abs(Pp - Ts @ Pm @ Ts).max(), abs(bq + bm)))
print(r'  第二列非零（倾角确实换后验），第三列恒为 $0$（对称式是精确的、不是拟合），最后一列说奇部只在 $q_{01}$ 上。')
print('')

# ---------------------------------------------------------------- [2] 闭式本体
sol, mats, p22x, res0, degs, lin, exact = tilt_series(1.25, .30, 0.0, .85, .6, 1.0, .7, 1.0, order=3)
print(r'[2] 默认植物的级数闭式（精确有理算术，$c_{13}=0$，$a_{11}=5/4$，$a_{12}=3/10$，$a_{22}=17/20$，$a_{33}=3/5$，$W=\operatorname{diag}(1,7/10,1)$）。')
print(r'  逐阶方程对全部未知系数的总次数 = %s；逐阶**对新**未知数是否线性 = %s。' % (degs, lin))
print(r'  （$q$ 的定点方程是二次的，但第 $k$ 阶里新系数只一次出现，所以整条链是回代线性求解，没有猜根。）')
print(r'  验根：把闭式种子 $(w_3/(1-a_{33}^{2}),\,0,\,p_{22})$ 代回三条 $\varphi^{0}$ 方程，simplify 给出 %s；$p_{22}=%s$。' % (res0, p22x))
print(r'  %-6s %-42s %-22s' % ('系数', '精确值', '小数'))
for nm in ['c0', 'c1', 'c2', 'c3', 'd0', 'd1', 'd2', 'd3', 'e0', 'e1', 'e2', 'e3']:
    v = sol[sp.Symbol(nm)]
    print(r'  %-6s %-42s %-22.17f' % (nm, str(v)[:42], float(v)))
print(r'  零式样：$c_1=c_3=d_0=d_2=e_1=e_3=0$ 全是 L2 的推论，不是我挑的；奇部只剩 $d_1\varphi+d_3\varphi^{3}$。')
p2v = (sp.sqrt(132449) - 143) / 120
d1f = (R(3, 10) * R(17, 20) * p2v * (1 - R(5, 4) * R(3, 5))
       / ((R(3, 5) ** 2 - 1) * (R(3, 10) ** 2 * p2v + (1 - R(17, 20) * R(3, 5)))))
print(r'  $d_1$ 的闭式（这一行我第一版抄错过：分母末项是 $w_1(a_{33}^{2}-1)$，不是 $w_1(a_{33}^{3}-1)$）：')
print(r'      $d_1=\dfrac{a_{12}a_{22}p_{22}w_3(1-a_{11}a_{33})}{(a_{33}^{2}-1)\left[a_{12}^{2}p_{22}+(1-a_{22}a_{33})w_1\right]}$')
print(r'  照这一行手算 $d_1=%.12f$，sympy 解出 $%.12f$，差 $%.2e$。手算对不上就是我抄写错了，不是级数错了。'
      % (float(d1f), float(sol[sp.Symbol('d1')]), abs(float(d1f) - float(sol[sp.Symbol('d1')]))))
print(r'  $P_0\sim P_3$ 的数值矩阵（原坐标；$P_1$ 只有指标 3 那一带非零，正是 L2 的"奇部住在 $e_3$"）：')
for k, M in enumerate(mats):
    print(r'  $P_%d$ =' % k)
    for row in M:
        print(r'    ' + ' '.join('%+.14f' % x for x in row))
print('')

# ---------------------------------------------------------------- [3] 闭式对求解器
print(r'[3] 闭式对求解器：$\Phi_0^{\rm num}(\varphi)$ 与 $\operatorname{tr}(\Theta\sum_{k\le n}\varphi^{k}P_k)$ 之差应当正好是截断阶。')
print(r'  %-9s %-19s %-19s %-12s %-13s %-19s %-13s'
      % ('phi(deg)', 'Phi num', 'Phi tay(k<=2)', 'res2', 'res2/phi^3', 'Phi tay(k<=3)', 'res3/phi^4'))
for b3 in [0.0, 0.2]:
    Plx = plant(b3=b3)
    Th = Plx.Th
    print(r'  $b_3=%.1f$' % b3)
    for dg in [0.05, 0.1, 0.2, 0.5, 1.0]:
        a = dg * DEG
        num = phi0(Plx, np.tan(a))
        t2 = float(np.trace(Th @ tay(mats[:3], a)))
        t3 = float(np.trace(Th @ tay(mats, a)))
        print(r'  %-9.3g %-19.14f %-19.14f %-12.3e %-13.3e %-19.14f %-13.3e'
              % (dg, num, t2, abs(num - t2), abs(num - t2) / a ** 3, t3, abs(num - t3) / a ** 4))
print(r'  更强的核对：直接比后验矩阵（逐元素，不经 $\Theta$ 投影）。$b_3=0$ 时奇部不出现，$k\le2$ 与 $k\le3$ 应当一样。')
print(r'  %-9s %-27s %-18s %-25s' % ('phi(deg)', 'max|P0num - tay(k<=3)|', '/phi^4', 'max|P0num - tay(k<=2)|'))
Pl0 = plant()
for dg in [1.0, 2.0, 5.0, 10.0]:
    a = dg * DEG
    dp3 = np.abs(P0num(Pl0, np.tan(a)) - tay(mats, a)).max()
    dp2 = np.abs(P0num(Pl0, np.tan(a)) - tay(mats[:3], a)).max()
    print(r'  %-9.2f %-27.3e %-18.3e %-25.3e' % (dg, dp3, dp3 / a ** 4, dp2))
print('')

# ---------------------------------------------------------------- [4] a,c 的闭式对实测
M0, M1, M2, M3 = mats
p33, d1 = float(sol[sp.Symbol('c0')]), float(sol[sp.Symbol('d1')])
e2, c2 = float(sol[sp.Symbol('e2')]), float(sol[sp.Symbol('c2')])
FD = np.array([[1., 0., 0.]])
F2 = np.array([[1., 0., 0.], [0., 0., 1.]])
PHI = 2.0 * DEG                       # 实测差分用的角度，和 #25/#27 一致
print(r'[4] $a,c$ 的闭式对实测：$c=\operatorname{tr}(\Theta P_1)$、$a=\operatorname{tr}(\Theta P_2)$，实测用 $\varphi=\pm2^{\circ}$ 的偶部/奇部差分。')
print(r'  这一节同时回答 #27 的"$c$ 随 $b_3$ 一阶变"到底是哪一支：闭式说 $c=2\gamma K_3(d_1K_2-p_{33}K_1)$，旋钮是 $K_3$ 不是 $b_3$。')
print(r'  %-7s %-11s %-11s %-12s %-12s %-12s %-12s %-11s %-11s'
      % ('$b_3$', '$K_3$', '$c$ 闭式', '$c$ 实测', '$c/b_3$', '$c/K_3$', '$a$ 闭式', '$a$ 实测', '$\\varphi^*$(deg)'))
for b3 in [0.0, 1e-4, 1e-2, 0.05, 0.2, 0.5]:
    Plx = plant(b3=b3)
    Th = Plx.Th
    K, gam = gain(Plx)
    cp = float(np.trace(Th @ M1))
    ap = float(np.trace(Th @ M2))
    f0 = phi0(Plx, 0.0)
    fp = phi0(Plx, np.tan(PHI))
    fm = phi0(Plx, np.tan(-PHI))
    cm = (fp - fm) / (2 * PHI)
    am = (fp + fm - 2 * f0) / (2 * PHI ** 2)
    phistar = -cm / (2 * am) / DEG
    print(r'  %-7.3g %-11.5f %-11.5f %-12.5f %-12s %-12s %-11.5f %-12.5f %-11.4f'
          % (b3, K[2], cp, cm,
             ('%.4f' % (cm / b3)) if b3 else '--',
             ('%.4f' % (cp / K[2])) if abs(K[2]) > 1e-12 else '--',
             ap, am, phistar))
print(r'  读法：$a>0$，地板在 $\varphi^*=-c/(2a)$ 取极小；$c$ 与 $b_3$ 反号，所以 $\varphi^*$ 和 $b_3$ 同号，表里全为正因为网格只取 $b_3\ge0$。')
print(r'  $\varphi^*/b_3$ 这一列（末列除第一列）从 $6.37$ 一路掉到 $5.47$：#27 那条 "$\varphi^*=6.36\beta$" 只在 $\beta\lesssim0.05$ 是准的，')
print(r'  $\beta=0.2$ 已经偏 $2\%$、$\beta=0.5$ 偏 $14\%$。而由 $\Phi_0(\beta,\varphi)=\Phi_0(-\beta,-\varphi)$，负 $b_3$ 那一支镜像过去即可。')
print('')

# ---------------------------------------------------------------- [4b] 校准我自己的旧读数
print(r'[4b] 闭式反过来校准实测：$a$ 的 $\varphi$ 差分被偶数阶余项污染，$a_{\rm meas}(\varphi)=a+h_4\varphi^{2}+O(\varphi^{4})$。')
print(r'  所以 $\varphi=2^{\circ}$ 上量到的那个 $a$ 偏高，Richardson 外推 $(4a_{\varphi}-a_{2\varphi})/3$ 应当落回闭式值。')
print(r'  %-7s %-13s %-13s %-13s %-13s %-13s %-13s'
      % ('$b_3$', '$a$ 闭式', '$a$(1deg)', '$a$(2deg)', '$a$(4deg)', 'Richardson', '相对偏差 2deg'))
for b3 in [0.0, 0.05, 0.2, 0.5]:
    Plx = plant(b3=b3)
    f0 = phi0(Plx, 0.0)
    ap = float(np.trace(Plx.Th @ M2))

    def am(dg):
        h = dg * DEG
        fp = phi0(Plx, np.tan(h))
        fm = phi0(Plx, np.tan(-h))
        return (fp + fm - 2 * f0) / (2 * h ** 2)
    a1, a2, a4 = am(1.), am(2.), am(4.)
    rich = (4 * a2 - a4) / 3.
    print(r'  %-7.3g %-13.8f %-13.8f %-13.8f %-13.8f %-13.8f %-13.3e'
          % (b3, ap, a1, a2, a4, rich, abs(a2 - ap) / ap))
print(r'  结论：#25/#27 里那个 $a=3.91169$ 是 $2^{\circ}$ 的读数，闭式真值 $3.90942$，偏高 $5.8\times10^{-4}$（相对）。')
print(r'  量级无害，但它随 $\varphi^{2}$ 走——以后凡是"用一个角度上的差分当系数"的表，都得标明角度并把外推附上。')
print('')

# ---------------------------------------------------------------- [5] 门票律
print(r'[5] 门票律 $gap_{\rm del}=p_{33}\,\Theta_{33}=\dfrac{w_3}{1-a_{33}^{2}}\,\gamma K_3^{2}$（闭式，不含拟合常数）对求解器。')
print(r'  对照 #27/#28 记下的经验"$gap/\Vert\Theta e_3\Vert^2=3.82$"：分母换三种口径，看哪一种是常数。')
print(r'  %-7s %-14s %-14s %-12s %-13s %-13s %-13s %-13s'
      % ('$b_3$', '$\\Phi_0(F_{del})$', '$\\Phi_0(F_2)$', 'gap 实测', '$p_{33}\\Theta_{33}$', 'gap/(p33*Th33)', 'gap/|Th e3|^2', 'gap/lmax(Th)^2'))
for b3 in [1e-4, 1e-2, 0.05, 0.2, 0.5, 1.0]:
    Plx = plant(b3=b3)
    Th = Plx.Th
    gF = floor0_safe(Plx, FD)
    g2 = floor0_safe(Plx, F2)
    if gF is None or g2 is None:
        print(r'  %-7.3g  DARE 无解' % b3)
        continue
    gap = gF[0] - g2[0]
    th33 = float(Th[2, 2])
    nrm = float(np.linalg.norm(Th @ np.array([0., 0., 1.])))
    lmax = float(np.linalg.eigvalsh(Th).max())
    print(r'  %-7.3g %-14.8f %-14.8f %-12.8f %-13.8f %-13.6f %-13.6f %-13.6f'
          % (b3, gF[0], g2[0], gap, p33 * th33, gap / (p33 * th33), gap / nrm ** 2, gap / lmax ** 2))
print(r'  第三、四列应当逐位相同（残差到求解器精度），第五列恒为 $1$；第六、七列随 $\beta$ 漂——那个 $3.82$ 只是 $\beta$ 在某一点上')
print(r'  $p_{33}/\lambda_{\max}(\Theta)$ 的取值，不是一条律。')
print('')

# ---------------------------------------------------------------- [6] 回收率
print(r'[6] 回收率：极小深度 $c^{2}/(4a)$ 相对门票 $p_{33}\Theta_{33}$。闭式里 $\gamma$ 和 $K_3^{2}$ 全约掉，所以它应当是常数量。')
print(r'  %-7s %-13s %-13s %-13s %-13s %-13s %-13s'
      % ('$b_3$', '深度 闭式', '深度 实测', '门票 闭式', '回收率 闭式', '回收率 实测', '$K_3$ 口径'))
for b3 in [1e-2, 0.05, 0.2, 0.5, 1.0]:
    Plx = plant(b3=b3)
    Th = Plx.Th
    K, gam = gain(Plx)
    cp = float(np.trace(Th @ M1))
    ap = float(np.trace(Th @ M2))
    f0 = phi0(Plx, 0.0)
    fp = phi0(Plx, np.tan(PHI))
    fm = phi0(Plx, np.tan(-PHI))
    cm = (fp - fm) / (2 * PHI)
    am = (fp + fm - 2 * f0) / (2 * PHI ** 2)
    depth_c = cp ** 2 / (4 * ap)
    depth_m = cm ** 2 / (4 * am)
    tick_c = p33 * gam * K[2] ** 2
    gF = floor0_safe(Plx, FD)
    g2 = floor0_safe(Plx, F2)
    tick_m = gF[0] - g2[0]
    print(r'  %-7.3g %-13.6f %-13.6f %-13.6f %-13.6f %-13.6f %-13s'
          % (b3, depth_c, depth_m, tick_c, depth_c / tick_c, depth_m / tick_m,
             '%.4f' % (cp / K[2])))
print(r'  回收率的闭式 $=\dfrac{(d_1K_2-p_{33}K_1)^{2}}{p_{33}\left[p_{33}K_1^{2}+e_2K_2^{2}+(c_2-p_{33})K_3^{2}-2d_1K_1K_2\right]}$，')
print(r'  $\gamma$ 与 $K_3^{2}$ 同约：#28 实测的 $0.9778$ 若是这条律在某个 $\beta$ 的取值，那它的"几乎不随 $\beta$ 变"就是约掉的后果，不是巧合。')
print('')

# ---------------------------------------------------------------- [7] 门票律的适用边界
print(r'[7] 门票律的适用边界：$p_{33}=w_3/(1-a_{33}^{2})$ 只读 $A$ 的第三行，而 $x_3$ 在 $A$ 的上三角结构里是自治的，')
print(r'  所以这条律没道理只在 $c_{13}=0$ 成立——把它推到 $c_{13}\ne0$（L2 的对称性已经没了）和其它 $(a_{12},a_{33},w_2,w_u)$ 上试。')
print(r'  %-6s %-6s %-6s %-6s %-6s %-6s %-9s %-14s %-14s %-12s'
      % ('$c_{13}$', '$a_{12}$', '$a_u$', '$w_2$', '$w_u$', '$b_3$', '$p_{33}$', 'gap 实测', '$p_{33}\\Theta_{33}$', '相对偏差'))
worst7 = 0.
bad7 = 0
for c13 in [0.0, 0.3, 0.6]:
    for a12, a33, w2, wu, b3 in [(.30, .6, .7, 1.0, .2), (.55, .6, .7, 1.0, .2), (.30, .9, .7, 1.0, .2),
                                 (.30, .6, 1.3, 1.0, .2), (.30, .6, .7, 2.5, .2), (.30, .6, .7, 1.0, 1.0),
                                 (.55, .9, 1.3, 2.5, 1.0)]:
        Plx = plant(c13=c13, a12=a12, au=a33, w2=w2, wu=wu, b3=b3)
        p33g = float(Plx.W[2, 2]) / (1 - float(Plx.A[2, 2]) ** 2)
        gF, g2 = floor0_safe(Plx, FD), floor0_safe(Plx, F2)
        if gF is None or g2 is None:
            print(r'  %-6.2f %-6.2f %-6.2f %-6.2f %-6.2f %-6.2f %-9.4f %-14s  基线不可检测' % (c13, a12, a33, w2, wu, b3, p33g, '--'))
            continue
        gap = gF[0] - g2[0]
        th33 = float(Plx.Th[2, 2])
        rel = gap / (p33g * th33) - 1.
        worst7 = max(worst7, abs(rel))
        bad7 += abs(rel) > 1e-9
        print(r'  %-6.2f %-6.2f %-6.2f %-6.2f %-6.2f %-6.2f %-9.4f %-14.8f %-14.8f %-12.2e'
              % (c13, a12, a33, w2, wu, b3, p33g, gap, p33g * th33, rel))
print(r'  全表最大相对偏差 $%.2e$，违例格数 %d/%d。' % (worst7, bad7, 21))
print(r'  这一节的结果和我的猜想相反，所以按数据改口：门票律 $gap=p_{33}\Theta_{33}$ **需要** $c_{13}=0$——')
print(r'  $c_{13}=0$ 的 7 株植物（$\varnothing(3/10\to55/100,\;a_{33}\,0.6\to0.9,\;w_2\,0.7\to1.3,\;w_u\,1\to2.5,\;\beta\,0.2\to1.0$）')
print(r'  全部在 $10^{-14}$ 内成立；一旦 $c_{13}\ne0$ 实测门票比律给的**便宜** $20\%\sim77\%$。')
print('')

# ---------------------------------------------------------------- [7b] 机制：$\Delta P$ 是不是秩一
print(r'[7b] 机制。把"$gap=p_{33}\Theta_{33}$"拆到底：它等价于删行带来的后验增量')
print(r'      $\Delta P=P_0(F_{\rm del})-P_0(F_2)\overset{?}{=}p_{33}\,e_3e_3^{\top}$，')
print(r'  即"少读一次 $x_3$"就是给 $(3,3)$ 补上自治方差 $w_3/(1-a_{33}^{2})$，别的一律不动。这条若成立，门票律对**任何** $\Theta$ 都成立，不需要秩一。')
print(r'  %-6s %-6s %-6s %-6s %-12s %-12s %-12s %-12s %-12s %-12s %-12s %-13s'
      % ('$c_{13}$', '$a_u$', '$w_u$', '$b_3$', '$\\Delta P_{11}$', '$\\Delta P_{12}$', '$\\Delta P_{13}$',
         '$\\Delta P_{22}$', '$\\Delta P_{23}$', '$\\Delta P_{33}$', '$p_{33}$', '$\\max|\\Delta P-p_{33}e_3e_3^{\\top}|$'))
for c13, a33, wu, b3 in [(0.0, .6, 1.0, .2), (0.0, .9, 1.0, .2), (0.0, .6, 2.5, 1.0), (0.0, .9, 2.5, 1.0),
                         (0.3, .6, 1.0, .2), (0.3, .9, 1.0, .2), (0.6, .6, 1.0, 1.0), (0.6, .9, 2.5, 1.0)]:
    Plx = plant(c13=c13, au=a33, wu=wu, b3=b3)
    gF, g2 = floor0_safe(Plx, FD), floor0_safe(Plx, F2)
    if gF is None or g2 is None:
        print(r'  %-6.2f %-6.2f %-6.2f %-6.2f  基线不可检测' % (c13, a33, wu, b3))
        continue
    dP = gF[1] - g2[1]
    p33g = float(Plx.W[2, 2]) / (1 - float(Plx.A[2, 2]) ** 2)
    err = np.abs(dP.copy())
    err[2, 2] = abs(dP[2, 2] - p33g)
    print(r'  %-6.2f %-6.2f %-6.2f %-6.2f %-12.5f %-12.5f %-12.3e %-12.5f %-12.5f %-12.8f %-12.8f %-13.3e'
          % (c13, a33, wu, b3, dP[0, 0], dP[0, 1], dP[0, 2], dP[1, 1], dP[1, 2], dP[2, 2], p33g, err.max()))
print(r'  读法：两个 $F$ 都以无限精度读 $x_1$，所以 $P_0$ 的第 1 行/列恒为零，$\Delta P_{11}=\Delta P_{12}=\Delta P_{13}=0$ 是结构给的、不算证据。')
print(r'  真内容全在 $2$–$3$ 块里：$c_{13}=0$ 那几行该块只剩 $\Delta P_{33}=p_{33}=w_3/(1-a_{33}^{2})$（$\Delta P_{22}=\Delta P_{23}=0$），')
print(r'  于是 $p_{33}\Theta_{33}$ 不是"秩一 $\Theta$ 的特例"而是 $\operatorname{tr}(\Theta\Delta P)$ 的恒等变形，对任何 $\Theta$ 都成立。')
print(r'  $c_{13}\ne0$ 那四行：$\Delta P_{33}$ 缩水到 $p_{33}$ 的 $96\%\sim32\%$，同时 $\Delta P_{22}>0$、$\Delta P_{23}<0$ 冒出来——')
print(r'  机制是 $x_1(t+1)=a_{11}x_1+a_{12}x_2+c_{13}x_3$：只要 $c_{13}\ne0$，那一行本来就免费透露了 $x_3$，')
print(r'  而且透露的方式是和 $x_2$ 捆在一起（$a_{12}\ne0$），所以增量的方向不再对齐 $e_3$，专用 $x_3$ 行的边际价值也低于自治方差。')
print(r'  下面把这两件事拆开：固定 $c_{13}\ne0$，把 $a_{12}$（$x_2$ 进 $x_1$ 的那条边）关掉，看坏在哪个位置。')
print(r'  %-7s %-7s %-14s %-14s %-14s' % ('$c_{13}$', '$a_{12}$', '$\\Delta P_{22}$', '$\\Delta P_{23}$', '$\\Delta P_{33}/p_{33}$'))
for c13 in [0.0, 0.3, 0.6]:
    for a12 in [0.0, .30]:
        Plx = plant(c13=c13, a12=a12, b3=.2)
        gF, g2 = floor0_safe(Plx, FD), floor0_safe(Plx, F2)
        if gF is None or g2 is None:
            print(r'  %-7.2f %-7.2f  基线不可检测' % (c13, a12))
            continue
        dP = gF[1] - g2[1]
        p33g = float(Plx.W[2, 2]) / (1 - float(Plx.A[2, 2]) ** 2)
        print(r'  %-7.2f %-7.2f %-14.6f %-14.6f %-14.6f' % (c13, a12, dP[1, 1], dP[1, 2], dP[2, 2] / p33g))
print(r'  结果很干净：$a_{12}=0$ 时 $\Delta P_{22}=\Delta P_{23}$ 恒为零，坏掉的只有 $\Delta P_{33}$（$0.938$、$0.846$）——')
print(r'  那就是"$x_3$ 被 $x_1$ 那行免费透露"造成的边际贬值。$(2,3)$ 块的残留必须 $c_{13}$ 与 $a_{12}$ 同时非零才出现，')
print(r'  即那条免费信息是和 $x_2$ 捆在一起进来的，于是增量不再对齐 $e_3$。')
print(r'  [7] 里实测门票比律便宜 $20\%\sim77\%$ 正是这件事，门票律在非自治 $x_3$ 上只剩不等式。')
print('')

# ---------------------------------------------------------------- [8] "塌到 t" 的判决
print(r'[8] 判决 #28 的"$a,c$ 是支撑倾斜 $t$ 的函数"：同一组 $P_k$ 能不能同时服务所有 $\beta$？换一株 $A$ 之后它立刻失效吗？')
print(r'  (甲) 只用闭式的 $P_0\sim P_3$ 与 $\Theta_\beta=\gamma_\beta K_\beta K_\beta^{\top}$ 预测 $\Phi_0(\beta,\varphi)$，不重新拟合任何曲线。')
print(r'  %-7s %-9s %-19s %-19s %-13s' % ('$b_3$', 'phi(deg)', 'Phi 预测', 'Phi 求解器', '相对误差'))
for b3 in [0.5, 1.0, 3.0]:
    Plx = plant(b3=b3)
    Th = Plx.Th
    for dg in [1.0, 2.0, 5.0]:
        h = dg * DEG
        pr = float(np.trace(Th @ tay(mats, h)))
        nu = phi0(Plx, np.tan(h))
        print(r'  %-7.2f %-9.1f %-19.12f %-19.12f %-13.2e' % (b3, dg, pr, nu, abs(pr - nu) / abs(nu)))
print(r'  误差随 $\varphi^{4}$ 走（$[5^{\circ}]$ 那几行最大），因为我只算到 $P_3$；$\beta$ 换成 $3$ 也没有引入新误差源——')
print(r'  这正是 L1 的内容：$\beta$ 在地板里只通过 $\Theta$ 出现，所以"$\Phi_0$ 对 $t$ 的函数形状"根本不是 $\beta$ 的函数。')
print(r'  (乙) 换 $A$（$c_{13}=0.3$）后拿同一组 $P_k$ 去预测：')
print(r'  %-9s %-19s %-19s %-13s' % ('phi(deg)', 'Phi 预测(c13=0 的 P_k)', 'Phi 求解器(c13=0.3)', '相对误差'))
Plc = plant(c13=0.3, b3=0.2)
for dg in [1.0, 2.0, 5.0]:
    h = dg * DEG
    pr = float(np.trace(Plc.Th @ tay(mats, h)))
    nu = phi0(Plc, np.tan(h))
    print(r'  %-9.1f %-19.12f %-19.12f %-13.2e' % (dg, pr, nu, abs(pr - nu) / abs(nu)))
print(r'  失配是宏观的：$P_k$ 是 $(A,W)$ 的函数，不是"通用倾斜函数"。#28 那张两支都塌到 $t$ 的表，')
print(r'  真实内容是"在各自那一株植物里，$\Phi_0$ 对 $\Theta$ 线性、对倾角解析"，跨植物不能共用一条曲线。')
print('')
