r"""E49 = 把放大族换成特征坐标 $(\varphi,\tau_1,\tau_2)$，然后攻 #23 的"两族重合"。

引理（先把口径说清楚，它决定了这一族到底是什么）：整条流水线上 $x$ 与 $I$ 只通过
    $M=F^{\top}V^{-1}F$
依赖 $(F,V)$。因为稳态对 $(\tilde P,P)$ 满足 $\tilde P=APA^{\top}+W,\ P=(\tilde P^{-1}+M)^{-1}$，
$P$ 只是 $M$ 的函数，$\tilde P$ 是只含 $M$ 的那个 DARE 的不动点；而
    $I=\tfrac12\log_2\det(I_r+V^{-1/2}F\tilde PF^{\top}V^{-1/2})=\tfrac12\log_2\det(I_n+\tilde PM)$。
于是规范自由度和"放大到底放大了什么"都变成一句话：$F$ 固定为 $\begin{bmatrix}1&0&0\\0&0&1\end{bmatrix}$ 时，
$r=2$ 相关噪声族 $=\{M\succ0:\operatorname{ran}M\subseteq\operatorname{span}(e_1,e_3)\}$，
删行族 $=\{m\,e_1e_1^{\top}\}$，而**中间还躺着一整条旋转单通道族** $\{\tau\,w(\varphi)w(\varphi)^{\top}\}$，
$w(\varphi)=\cos\varphi\,e_1+\sin\varphi\,e_3$——它是前者的秩 1 边界、后者的角度变形。

为什么 E44/E47 扫不到它：$(v_1,v_2,\rho)$ 里秩 1 边界要求 $\rho\to\pm1$ 且 $v_2\to\infty$ 协调进行，
$V$ 同时趋于奇异（E48 就在 $|\rho|=1$ 那格抛了 `Singular matrix`）。特征坐标把它变成内点：
$V^{-1}=R(\varphi)\operatorname{diag}(\tau_1,\tau_2)R(\varphi)^{\top}$，取 $\tau_2=0$ 就是精确的秩 1，
而且此时它就是"一个方向为 $w(\varphi)$ 的传感器"，可以用 $r=1$ 的 API 精确算，$\kappa(V)=1$。

要判的命题（#23 第一节）：$E_{\rm full}(D)=E_{\rm del}(D)$，$D$ 是绝对预算 $\operatorname{tr}(\Theta P)\le D$。
$E_{\rm full}\le E_{\rm del}$ 由闭包包含是平凡的；反过来的严格不等号只要有一个 $\varphi$ 买到更低的车票就成立。
"""
import sys
sys.stdout.reconfigure(encoding='utf-8')
import numpy as np
import p0.exp_plants_oos as E35
from p0.ticket import floor0, floor0_safe

ln2 = np.log(2.)
NPROW = 3


def plant(c13=0.0, au=.6, b3=0.0, wu=1.0):
    A = np.array([[1.25, .30, c13], [0., .85, 0.], [0., 0., au]])
    B = np.array([[1.], [0.], [b3]])
    W = np.diag([1.0, .7, wu])
    return E35.Pl(A, B, W, np.eye(3), np.array([[1.0]]))


F2 = np.array([[1., 0., 0.], [0., 0., 1.]])
FD = np.array([[1., 0., 0.]])
Pl = plant()
Phi2 = floor0(Pl, F2)[0]
PhiD = floor0(Pl, FD)[0]


def sym(M):
    return .5 * (M + M.T)


def cost_I(Fx, V):
    r"""$(\operatorname{tr}(\Theta P),\ I,\ \tilde P)$，绝对口径（不减地板）。DARE 无有限解时记成 $+\infty$。"""
    try:
        x, I, Pt = Pl.xI(Fx, V, 0.0)
    except np.linalg.LinAlgError:
        return float('inf'), float('nan'), None
    return x, I, Pt


def M_of(Fx, V):
    return sym(Fx.T @ np.linalg.solve(V, Fx))


def w_col(phi):
    return np.array([np.cos(phi), 0., np.sin(phi)])


def rot_cost(tau, phi):
    r"""旋转单通道：$M=\tau\,w w^{\top}$，用 $r=1$ 的 API（$\kappa=1$）。"""
    Fx = w_col(phi).reshape(NPROW, 1).T
    return cost_I(Fx, np.array([[1.0 / tau]]))


def biset(fun, target, lo=-8.0, hi=12.0, it=80):
    r"""$\tau=10^{u}$ 上二分 $\operatorname{fun}(\tau)$ 的第一分量 $=$ target（对 $\tau$ 单调减）。
    下界只到 $10^{-8}$：再小 $\tilde P$ 发散，scipy 直接给 `Failed to find a finite solution`（这不是失败，是"$\tau\to0$ 没有有限解"这件事的数值证据）。"""
    if fun(10. ** lo)[0] < target or fun(10. ** hi)[0] > target:
        return None
    for _ in range(it):
        m = .5 * (lo + hi)
        if fun(10. ** m)[0] > target:
            lo = m
        else:
            hi = m
    tau = 10. ** (.5 * (lo + hi))
    c, I, _ = fun(tau)
    return tau, c, I, abs(c - target)


def biset_I(fun, target, lo=-8.0, hi=12.0, it=80):
    r"""对**速率**配平：$I$ 随 $\tau$ 单调增，二分使 $I=$ target，返回该点的 $(\tau,c,I)$。
    有了它就能问"等速率下多付多少绝对成本"，把 $E(D)$ 曲线之间的距离换成与预算无关的量。"""
    if fun(10. ** lo)[1] > target or fun(10. ** hi)[1] < target:
        return None
    a, b = lo, hi
    for _ in range(it):
        m = .5 * (a + b)
        if fun(10. ** m)[1] < target:
            a = m
        else:
            b = m
    tau = 10. ** (.5 * (a + b))
    c, I, _ = fun(tau)
    return tau, c, I, abs(I - target)


XT = [1e-6, 1e-4, 1e-2, 1e-1, 1.0]
PHI = np.array([.5, 1, 2, 3, 4, 5, 7, 10, 12, 15, 20, 25, 30, 40, 50, 65, 80, 85]) * np.pi / 180.0
FLOOR = {}      # $\varphi$ 度 -> $\Phi_0(\varphi)-\Phi_0(0)$，即 $\tau\to\infty$ 端的惩罚
SHIFT = {}      # $x_t$ -> {$\varphi$ 度 -> 等速率成本惩罚}

print(r'== E49：$M=F^{\top}V^{-1}F$ 约化 + 旋转单通道对 #23 的反例搜索 ==')
print('')
print(r'[1] 引理的数值面：规范变换 $(F,V)\mapsto(TF,\ TVT^{\top})$ 保持 $M=F^{\top}V^{-1}F$，因而保持 $(\operatorname{tr}\Theta P,I)$')
rng = np.random.default_rng(7)
V0 = np.array([[2.3, -.9], [-.9, .7]])
c0, I0, _ = cost_I(F2, V0)
M0 = M_of(F2, V0)
for k in range(3):
    T = rng.normal(size=(2, 2)) + np.eye(2) * 2.
    Ft = T @ F2
    Vt = T @ V0 @ T.T
    ct, It, _ = cost_I(Ft, Vt)
    print(r'  $T$ 第%d个:  $M$ 相对不变差 %.1e   $\Delta\operatorname{tr}(\Theta P)=%+.1e$   $\Delta I=%+.1e$'
          % (k + 1, np.linalg.norm(M_of(Ft, Vt) - M0) / np.linalg.norm(M0), ct - c0, It - I0))
print(r'  基线：$\operatorname{tr}(\Theta P)=%.9f$，$I=%.9f$。逐位不变即说明"族"的正确坐标是 $M$ 而不是 $(F,V)$。' % (c0, I0))

print('')
print(r'[2] 完美测量地板 $\Phi_0$ 随角度：$\Phi_0(w^{\top})=\operatorname{tr}(\Theta P_0)$，$P_0$ 是奇异 DARE（$V=0$）的解')
print(r'  参照：$\Phi_0(F_{r=2})=%.9f$，$\Phi_0(F_{\rm del})=%.9f$（#22 的"免费删"就是这两个逐位相等）' % (Phi2, PhiD))
mn = 1e30
for deg in [0, 1, 2, 3, 5, 7, 10, 15, 20, 30, 45, 60, 75, 85]:
    phi = deg * np.pi / 180.0
    Fx = w_col(phi).reshape(NPROW, 1).T
    ok, lam = Pl.detectable(Fx)
    got = floor0_safe(Pl, Fx) if ok else None
    if got is None:
        print(r'   $\varphi=%4.0f^\circ$  不可检测/奇异 DARE 无解' % deg)
        continue
    mn = min(mn, got[0])
    d = got[0] - PhiD
    if deg > 0:
        FLOOR[deg] = d
    print(r'   $\varphi=%4.0f^\circ$  $\Phi_0=%.9f$   差删行 %+.3e   除以 $\varphi^{2}$ %.3f   除以 $\varphi$ %.3f'
          % (deg, got[0], d, d / phi ** 2 if phi else float('nan'), d / phi if phi else float('nan')))
print(r'  最小地板 %.9f，在 $\varphi=0$ 取到；$\varphi\to0$ 时 $\Phi_0(\varphi)-\Phi_0(0)$ 对 $\varphi^2$ 是常数、对 $\varphi$ 发散'
      % mn)
print(r'  即地板在 $\varphi=0$ 是**非退化二次极小**。量化后果：转 $\varphi$ 度的单通道要先付 $\approx3.9\,\varphi^2$ 的绝对成本才*进得了场*（预算可行），')
print(r'  这跟 #22 里门票沿对准误差二阶消失（$\mathrm{gap}/c^2\approx3.06\to2.72$）是同一个二阶结构的两个侧面。')

print('')
print(r'[3] 主战：同一绝对预算 $D=\Phi_0(F_{\rm del})+x_t$ 下，旋转单通道 vs 删行（$\varphi=0$）。$\kappa(V)=1$，无 $V$ 奇异问题')
for phi in PHI:                      # 补齐 [2] 表里没列的角度，供 [3b] 用
    deg = round(phi * 180 / np.pi, 1)
    got = floor0_safe(Pl, w_col(phi).reshape(NPROW, 1).T)
    if got is not None:
        FLOOR[deg] = got[0] - PhiD
TAB = {}
for xt in XT:
    D = PhiD + xt
    b = biset(lambda t: rot_cost(t, 0.0), D)
    if b is None:
        print(r'  $x_t=%.0e$  删行族未 bracket' % xt)
        continue
    Ed, rd = b[2], b[3]
    print(r'  x_t = %.1e   D=%.12f   $E_{\rm del}$=%.9f   残差 %.1e   基线 $\tau$=%.3e' % (xt, D, Ed, rd, b[0]))
    TAB[xt] = []
    for phi in PHI:
        fun = lambda t, q=phi: rot_cost(t, q)
        g = biset(fun, D)
        h = biset_I(fun, Ed)
        if h is not None:
            SHIFT.setdefault(xt, {})[round(phi * 180 / np.pi, 1)] = h[1] - D
        if g is None:
            print(r'     $\varphi=%4.1f^\circ$  预算不可达%s'
                  % (phi * 180 / np.pi,
                     '' if h is None else r'（等速率仍可比：成本惩罚 %+.3e，除 $\varphi^2$ %.3f）'
                     % (h[1] - D, (h[1] - D) / phi ** 2)))
            continue
        d = g[2] - Ed
        TAB[xt].append((phi, d, g))
        print(r'     $\varphi=%4.1f^\circ$  $I$=%.9f   差删行 %+.3e   除 $\varphi^2$ %.3f   除 $\varphi$ %.3f   $\tau$=%.3e  残差 %.1e  %s%s'
              % (phi * 180 / np.pi, g[2], d, d / phi ** 2, d / phi, g[0], g[3],
                 '!! 负' if d < -1e-9 else '',
                 '' if h is None else r'  |  等速率成本惩罚 %+.3e 除 $\varphi^2$ %.3f 残差 %.1e'
                 % (h[1] - D, (h[1] - D) / phi ** 2, h[3])))
    if TAB[xt]:
        lo = min(TAB[xt], key=lambda e: e[1])
        print(r'     本预算下最好的角度：$\varphi=%.1f^\circ$，比删行 %+.3e bit' % (lo[0] * 180 / np.pi, lo[1]))

print('')
print(r'[3b] 把"率惩罚"换成"等速率下的成本惩罚"，检验旋转是不是纯粹的整体平移：比值 $=$ 成本惩罚 / 地板差，$1.00$ 即平移')
print(r'  预言：$D\to\Phi_0(F_{\rm del})$ 时 $\Delta_{\rm eq}(\varphi;D)\to\Phi_0(\varphi)-\Phi_0(0)$（$\tau\to\infty$ 端，[2] 的二次律就是它的值）')
RHOV = []
for xt in XT:
    if xt not in SHIFT:
        continue
    ks = [q for q in [1., 2., 3., 5., 10., 20., 50., 85.] if q in SHIFT[xt] and q in FLOOR]
    rr = [SHIFT[xt][q] / FLOOR[q] for q in ks]
    RHOV.append((xt, min(rr) - 1., max(rr) - min(rr)))
    print(r'  $x_t=%.0e$   ' % xt + '   '.join(r'$%g^\circ$:%.6f' % (q, v) for q, v in zip(ks, rr))
          + r'   |  $\rho=%.2e$   $\rho/x_t=%.4f$   跨角度相对散布 %.1e'
          % (min(rr) - 1., (min(rr) - 1.) / xt, (max(rr) - min(rr)) / min(rr)))
print(r'  读法：比值只跟预算走，跟角度几乎无关（同一行跨 $1^\circ\!\sim\!85^\circ$ 的相对散布见行尾），且随 $x_t$ 收紧单调趋于 $1$。')
print(r'  即 $\Delta_{\rm eq}(\varphi;D)=(1+\rho(D))\bigl(\Phi_0(\varphi)-\Phi_0(0)\bigr)$，实测 $\rho\in[%s]$。'
      % ', '.join('%.1e' % q[1] for q in RHOV))
print(r'  行尾那列 $\rho/x_t$：五个预算里前四个是同一个数（$0.18$），最后一个开始弯。$\rho$ 对 $\varphi$ 无关 + 对 $x_t$ 线性这两条我都没有闭式，是下一步要证的对象；')
print(r'  若成立，旋转族的整条前沿就是删行前沿乘一个与方向无关的因子再加常数平移，"方向不值钱"于是升级成可写下来的引理，而不只是一张表。')
print(r'  于是"旋转买不到东西"有了与预算无关的说法：旋转整条前沿就是删行前沿加一个 $\varphi$ 的常数平移（$3.9\varphi^{2}$ 起步），不是形状改变。')

print('')
print(r'[4] 内点角（秩 2 但主轴旋转）：$V^{-1}=R(\varphi)\operatorname{diag}(\tau_1,\tau_2)R(\varphi)^{\top}$，扫 $\varphi$ 与 $\tau_2$，二分 $\tau_1$ 配预算')


def in_cost(tau1, phi, tau2):
    r"""秩 2 内点、主轴转 $\varphi$：$V^{-1}=R\operatorname{diag}(\tau_1,\tau_2)R^{\top}$。直接拼 $V$，不做求逆（$\tau_1/\tau_2=10^{20}$ 时 $V^{-1}$ 在双精度下已经奇异）。"""
    R = np.array([[np.cos(phi), -np.sin(phi)], [np.sin(phi), np.cos(phi)]])
    V = R @ np.diag([1.0 / tau1, 1.0 / tau2]) @ R.T
    c, I, _ = cost_I(F2, V)
    return c, I, np.linalg.cond(V)


for xt in [1e-4, 1e-1]:
    D = PhiD + xt
    b = biset(lambda t: rot_cost(t, 0.0), D)
    if b is None:
        print(r'  $x_t=%.0e$  删行族未 bracket，跳过' % xt)
        continue
    Ed = b[2]
    print(r'  $x_t=%.0e$  $E_{\rm del}$=%.9f' % (xt, Ed))
    for phi in [5, 15, 30, 45, 60, 75, 85]:
        q = phi * np.pi / 180.0
        line = []
        for t2 in [1e-8, 1e-4, 1e-2, 1.0, 1e2]:
            g = biset(lambda t, a=q, c2=t2: in_cost(t, a, c2), D)
            line.append('--' if g is None else '%+.2e' % (g[2] - Ed))
        print(r'    $\varphi=%3d^\circ$   $\tau_2=$ %s' % (phi, '  '.join(line)))
    print(r'    （表元是 $I-I_{\rm del}$，$\tau_2$ 依次为 $10^{-8},10^{-4},10^{-2},1,10^{2}$；负得多才是反例，$\pm10^{-9}$ 以下是 E48 那档舍入）')

print('')
print(r'[5] 判决与后果')
neg = [(xt, p, d) for xt in TAB for (p, d, g) in TAB[xt] if d < -1e-9]
print(r'  可达格共 %d 个，显著负值 %d 个（阈值 $-10^{-9}$ bit 是 E48 在最坏格实测的杠子，这里的正偏差是 $10^{-2}\!\sim\!10^{0}$ bit，比杠子高七个数量级）。'
      % (sum(len(v) for v in TAB.values()), len(neg)))
print(r'  所以"旋转保留通道"这条竞争路线被否掉了，而且否得干净：它不是差一点点，而是**多数角度连预算都进不来**（地板抬高了 $\approx3.9\varphi^2$）。')
print(r'  [3b] 给了比"否掉"更有用的东西：等速率口径下旋转的代价 $=(1+\rho(x_t))\times$ 地板差，$\rho$ 与角度无关、与预算线性。这一条是可证的候选引理。')
print(r'  对 #23 的意义：两族重合现在有三块独立证据——E46 的 $\rho=0$ 闭式定理、E48 的逐格配准区间、E49 的整条 $M$ 锥角度扫描；')
print(r'  仍不是一块定理由，因为 E46 的证明只覆盖了对角子族。要把限定词去掉，需要证的是：')
print(r'  在 $M=\tau_1w_1w_1^{\top}+\tau_2w_2w_2^{\top}$ 里，$\operatorname{tr}(\Theta P)$ 与 $I$ 都只在 $w_k$ 落在 $\operatorname{span}(e_1,e_3)$ 的 $e_1$ 方向时才同时最优——')
print(r'  E49 的二次下界 $(\Phi_0(\varphi)-\Phi_0(0))/\varphi^2\approx3.9$ 就是这个命题在 $\tau=\infty$ 端的值，有限 $\tau$ 端我还给不出闭式。')
print(r'  另一条我暂时不动的路线：与其证"方向不值钱"，不如把 #23 的命题改写成 $M$-锥上的形式（$\inf_{M\succ0,\operatorname{ran}M\subseteq\operatorname{ran}F^{\top}}$），')
print(r'  那比"两族重合"更弱但更真，而且它把 C 那边"删方向"的口径统一成了"选 $M$ 的支撑"，两边可以说同一句话。')

print('')
print(r'[6] [3b] 那条 $\rho/x_t=0.183$ 是这个植物的泛函还是普适常数？换四个植物复测同一个比值')


def shift_law(Plx, xts=(1e-6, 1e-4, 1e-2), phis=(2., 10., 45.)):
    r"""返回 $[(x_t,\ \min_\varphi\rho,\ \text{跨 }\varphi\text{ 散布})]$；任一格算不出就整条返回 None。"""
    global Pl
    Pl = Plx
    try:
        Phi0d = floor0(Pl, FD)[0]
    except np.linalg.LinAlgError:
        return None
    out = []
    for xt in xts:
        D = Phi0d + xt
        b = biset(lambda t: rot_cost(t, 0.0), D)
        if b is None:
            out.append(None)
            continue
        rr = []
        for dg in phis:
            phi = dg * np.pi / 180.
            fun = lambda t, q=phi: rot_cost(t, q)
            h = biset_I(fun, b[2])
            try:
                fl = floor0(Pl, w_col(phi).reshape(NPROW, 1).T)[0] - Phi0d
            except np.linalg.LinAlgError:
                fl = None
            if h is None or fl is None or fl <= 0:
                rr = None
                break
            rr.append((h[1] - D) / fl)
        out.append(None if rr is None else (xt, min(rr) - 1., max(rr) - min(rr)))
    return out


for name, Plx in [(r'默认（$a_u=0.6$，$c_{13}=0$，$w_u=1$）', plant()),
                  (r'$a_u=0.9$', plant(au=.9)),
                  (r'$a_u=0.6$、$c_{13}=0.3$', plant(c13=.3)),
                  (r'$w_u=0.1$（第三个方向更吵）', plant(wu=.1)),
                  (r'$w=(1,0.7,4)$', plant(wu=4.))]:
    got = shift_law(Plx)
    if got is None:
        print(r'  %s  地板不可算，跳过' % name)
        continue
    print(r'  %s   ' % name + '   '.join('--' if q is None else
                                        r'$x_t=%.0e$:$\rho/x_t=%.4f$（散布 %.0e）' % (q[0], q[1] / q[0], q[2])
                                        for q in got))
print(r'  读法：这一块的 $\Phi_0(F_{\rm del})$ 与前面的表不是同一个数，所以只比 $\rho/x_t$ 这一个无量纲比值，不比绝对惩罚。')
print(r'  如果五个植物的 $\rho/x_t$ 互不相同，[3b] 的 0.183 就是 $(A,W,Q,R,F)$ 的泛函，下一步该问它是哪个泛函；如果一致，才有资格去猜普适常数。')

print('')
print(r'[7] 沿 $a_u$（第三个模态的本征值）扫同一个比值，同时打印地板的二次系数 $q(\varphi_0)=(\Phi_0(\varphi_0)-\Phi_0(0))/\varphi_0^{2}$ 与门票 $\mathrm{gap}_{\rm del}=\Phi_0(F_{\rm del})-\Phi_0(F_2)$')


def quad_coeff(Plx, dg=1.0):
    phi = dg * np.pi / 180.
    try:
        return (floor0(Plx, w_col(phi).reshape(NPROW, 1).T)[0]
                - floor0(Plx, FD)[0]) / phi ** 2
    except np.linalg.LinAlgError:
        return None


for label, kw in [('$a_u=0.20$', dict(au=.2)), ('$a_u=0.40$', dict(au=.4)), ('$a_u=0.50$', dict(au=.5)),
                  ('$a_u=0.60$', dict(au=.6)), ('$a_u=0.70$', dict(au=.7)), ('$a_u=0.80$', dict(au=.8)),
                  ('$a_u=0.85$', dict(au=.85)), ('$a_u=0.90$', dict(au=.9)), ('$a_u=0.95$', dict(au=.95)),
                  ('$a_u=0.99$', dict(au=.99)),
                  (r'$a_u=1.02$、$b_3=0$（不可控，问题本身不适定）', dict(au=1.02)),
                  (r'$a_u=1.2$、$b_3=1$（不稳定但可控）', dict(au=1.2, b3=1.)),
                  (r'$a_u=0.6$、$b_3=1$', dict(au=.6, b3=1.))]:
    try:
        Plx = plant(**kw)
    except np.linalg.LinAlgError as e:
        print(r'  %-42s 控制 DARE 无解：%s' % (label, str(e).split(chr(10))[0][:60]))
        continue
    got = shift_law(Plx, xts=(1e-6,), phis=(2., 10., 45.))
    q = quad_coeff(Plx)
    try:
        gp = floor0(Plx, FD)[0] - floor0(Plx, F2)[0]      # #22 的门票 $\mathrm{gap}_{\rm del}=\Phi_0(F_{\rm del})-\Phi_0(F_2)\ge0$
    except np.linalg.LinAlgError:
        gp = None
    if got is None or got[0] is None:
        print(r'  %-42s 未 bracket/地板不可算（$q=%s$）   $\mathrm{gap}_{\rm del}$=%s   $\|\Theta e_3\|/\|\Theta\|$=%.1e'
              % (label, '--' if q is None else '%.3f' % q, '--' if gp is None else '%.4g' % gp,
                 np.linalg.norm(sym(Plx.Th)[:, NPROW - 1]) / np.linalg.norm(sym(Plx.Th))))
        continue
    rt = got[0][1] / got[0][0]
    out_supp = np.linalg.norm(sym(Plx.Th)[:, NPROW - 1]) / np.linalg.norm(sym(Plx.Th))
    print(r'  %-42s $\rho/x_t=%.5f$   $q(1^\circ)=%s$   乘积 %s   $\mathrm{gap}_{\rm del}$=%s   $\|\Theta e_3\|/\|\Theta\|$=%.1e'
          % (label, rt, '--' if q is None else '%.4f' % q,
             '--' if q is None else '%.5f' % (rt * q),
             '--' if gp is None else '%.4g' % gp, out_supp))
print(r'  读法：$w_u$ 那一格已经说明这个比值不看不稳定模态的噪声功率。但扫完 $a_u$ 后"$1-a_u$ 的负幂"这种猜测死了——它先升后降（$0.5$–$0.6$ 附近有内部极大），')
print(r'  所以它同时依赖三个本征值，不是单模态的谱函数。乘积列 $q\cdot\rho/x_t$ 单调升（$0.51\to11$），说明 $\rho$ 与地板曲率也不是倒数关系。')
print(r'  门票那一列才是这行的适用条件：$b_3=0$ 时 $\mathrm{gap}=0$（删通道免费，#22 的机制二），比值扫得干净；')
print(r'  $b_3=1$ 时 $\mathrm{gap}>0$、$q(1^\circ)=-43.6$ 为负，即地板在 $\varphi=0$ 不再取极小——旋转进来的方向被控制器用了，"转了要付钱"这条律在这里失效。')
print(r'  这不是我算错：$b_3\ne0$ 时 $K$ 的第三分量非零，$\operatorname{ran}\Theta$ 离开 $S_1$，于是 #22 的免费条件没了。')
print(r'  【E50 更正】$c_{13}$ 也要走同一条路：$b_3=0$ 时 $(P_cA)_{13}=c_{13}P_{11}+a_uP_{13}$，$c_{13}\ne0$ 同样把支撑推出 $S_1$，')
print(r'  所以 $c_{13}=0.3$ 这一行整表 `--` 不是可行性失败——E50 的 [6] 八格地板全算得出、[7] 显示 $\Phi_0(F_2)$ 恒为 $0.537714$ 而 $\Phi_0(F_{\rm del})$ 被抬到 $0.807$，')
print(r'  是我下面的等速率配点窗口够不到那个目标值。条件应当写成一条：$\operatorname{ran}\Theta\subseteq S_1$（等价于 $\|\Theta e_3\|/\|\Theta\|=0$）。')
print(r'  所以 [2]/[3b] 的二次极小与因子化律的假设应当写成 $\operatorname{ran}\Theta\subseteq S_1$，而不是"$A$ 块对角"那么大。$|a_u|\ge1$ 且 $b_3=0$ 那行是参数化边界（不可控，控制 DARE 无解），不是求解器失败。')

print('')
print(r'[8] 把 $b_3$ 从 $0$ 连续打开：门票、地板符号、因子化律谁先动？（#26 队列第 2 条）')
print(r'  预言：$\Theta=K^{\top}(\cdot)K$ 且 $K_3\asymp b_3$，所以 $\|\Theta e_3\|\asymp b_3$、$\mathrm{gap}_{\rm del}\asymp b_3^{2}$；若地板差也 $\asymp b_3^{2}$ 且变号发生在 $b_3=0$，三条就是同一件事')


def best_angle(Plx, hi_deg=30.0, step=1.0):
    r"""$\arg\min_\varphi\Phi_0(w(\varphi)^{\top})$，粗网格。返回 $(\varphi^*,\ \Phi_0(\varphi^*)-\Phi_0(0))$。"""
    try:
        base = floor0(Plx, FD)[0]
    except np.linalg.LinAlgError:
        return None, None
    bn, bv = 0.0, 0.0
    for dg in range(0, int(hi_deg / step) + 1):
        a = dg * step
        got = floor0_safe(Plx, w_col(a * np.pi / 180.).reshape(NPROW, 1).T)
        if got is not None and got[0] - base < bv:
            bn, bv = a, got[0] - base
    return bn, bv


ROWS_QB = []
for b3 in [0.0, 1e-6, 1e-4, 1e-2, 0.05, 0.2, 0.5, 1.0]:
    Plx = plant(au=.6, b3=b3)
    s3 = np.linalg.norm(sym(Plx.Th)[:, NPROW - 1]) / np.linalg.norm(sym(Plx.Th))
    try:
        gp = floor0(Plx, FD)[0] - floor0(Plx, F2)[0]
    except np.linalg.LinAlgError:
        gp = None
    gps = '--' if gp is None else '%9.2e' % gp
    gqr = '--' if gp is None or s3 == 0 else '%.4f' % (gp / s3 ** 2)
    q2 = quad_coeff(Plx, 2.0)
    ROWS_QB.append((b3, q2))
    an, av = best_angle(Plx)
    got = shift_law(Plx, xts=(1e-6,), phis=(2., 10.))
    rt = '--' if got is None or got[0] is None else '%.5f' % (got[0][1] / got[0][0])
    print(r'  $b_3=%5.1e$   $\|\Theta e_3\|/\|\Theta\|$=%.2e   $\mathrm{gap}_{\rm del}$=%s   gap/支撑$^2$=%s   $q(2^\circ)=%s$   最佳角 %s  退让 %s   $\rho/x_t$=%s'
          % (b3, s3, gps, gqr,
             '--' if q2 is None else '%+.3f' % q2,
             '--' if an is None else '%.0f°' % an,
             '--' if av is None else '%+.3e' % av, rt))
print(r'  读法：`gap/支撑^2` 一列在 $b_3\le 5\times10^{-2}$ 就四个数 $3.821$–$3.828$，门票确实由 $\Theta$ 的支撑倾斜定价到二阶：$\mathrm{gap}_{\rm del}\asymp3.82\,\|\Theta e_3\|^{2}$，#22 的秩检验与这条二阶律是同一件事。')
fit = [(b, q) for b, q in ROWS_QB if q is not None and b >= 1e-2]
if len(fit) >= 3:
    sl, ic = np.polyfit([b for b, _ in fit], [q for _, q in fit], 1)
    print(r'  但地板那一列不是二阶的：$q(2^\circ)$ 对 $b_3$ 的最小二乘拟合 $=%.3f%+.3f\,b_3$，零点在 $b_3=%.3f$（用 $b_3\ge10^{-2}$ 的五点）。'
          % (ic, sl, -ic / sl))
    print(r'  也就是 $\Phi_0(\varphi)-\Phi_0(0)$ 随 $b_3$ 是**一阶**变，门票随 $b_3$ 是二阶变（$\mathrm{gap}_{\rm del}$ 在 $b_3=10^{-6}$ 已 $\sim10^{-14}$，$\propto b_3^{2}$）。两个不同阶，$b_3\to0$ 时一阶项先赢——所以"删通道要付钱"（$\mathrm{gap}_{\rm del}>0$ 只要 $b_3\ne0$）与"转方向划不划算"不是同一个判据，后者要等到 $b_3$ 长到 $\approx0.16$ 才翻。')
sm = [(b, q) for b, q in ROWS_QB if q is not None and 1e-4 <= b <= 5e-2]
for (b0, q0), (b1, q1) in zip(sm, sm[1:]):
    s_ = (q1 - q0) / (b1 - b0)
    print(r'  局部斜率段 $b_3\in[%.0e,%.0e]$：$\Delta q/\Delta b_3=%.2f$，该段自身的零点外推 $b_3=%.3f$（$b_3=0$ 处 $q=3.912$，$[0.05,0.2]$ 段含弯曲，粗拟合只用来给区间）'
          % (b0, b1, s_, b0 - q0 / s_))
print(r'  实测门槛：$b_3$ 从 $0$ 到 $5\times10^{-2}$ 的最佳角全钉在 $0^\circ$（退让 $+0.000$、$q(2^\circ)=+3.912\to+2.672$ 仍是亏的），符号翻发生在 $5\times10^{-2}$ 与 $2\times10^{-1}$ 之间（$1^\circ,3^\circ,4^\circ$ 才开始有退让）。我上一轮猜"任意非零 $b_3$ 就翻"，这里被我自己的表推翻。')
print(r'  $\rho/x_t$ 一列：正支（$q>0$）随 $b_3$ 缓升 $0.18324\to0.20052$，负地板差那一支全打印 `--`，即 $\Delta_{\rm eq}$ 在那个预算下 bracket 不进来。')
print(r'  所以 [3b] 的乘性因子化律在 $b_3\ne0$ 一侧不是被证伪，而是**不可评估**——它的适用域只能写成 $q>0$ / $\operatorname{ran}\Theta\subseteq S_1$ 那一支（与 #26 的结构条件同一条）。')
print(r'  【E50 修正】上面三行的读法有两条错：$\Phi_0$ 对 $\varphi$ 带着一次项 $c=-0.8700\,b_3$（对称性 $\Phi_0(\beta,\varphi)=\Phi_0(-\beta,-\varphi)$ 允许的奇次项），')
print(r'  二阶曲率 $a=3.9117$ 在 $b_3\le0.05$ 纹丝不动；本块的最佳角网格是 $[0,30^\circ]$ 步长 $1^\circ$，而真极小在 $\varphi^{*}=6.36\,b_3$（度），')
print(r'  于是"$b_3\approx0.158$ 才翻"就是 $1^\circ/6.36$，我的步长而非物理门槛。带符号角度的复测在 `p0/exp_signed.py`（E50），板上 #28。')

