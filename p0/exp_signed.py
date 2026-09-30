r"""E50 = 把旋转族的角度放开到带符号，用一条两行的对称性判决 #27 的"$q(2^\circ)$ 随 $b_3$ 一阶变"。

对称性（$c_{13}=0$ 时，$T=\operatorname{diag}(1,1,-1)$ 是植物的相似变换）：
    $A\to TAT=A,\;B\to TB=(1,0,-b_3)^{\top},\;W\to W,\;Q\to Q$，
于是 $b_3\to-b_3$ 的植物与原木是同一条系统、只是 $x_3$ 换了符号，而通道方向随之变成
    $w(\varphi)\to Tw(\varphi)=w(-\varphi)$。
结论是恒等式 $\Phi_0(\beta,\varphi)=\Phi_0(-\beta,-\varphi)$，且在 $\beta=0$ 处 $\Phi_0$ 对 $\varphi$ 是偶函数。

它把 #27 那条"一阶"读数拆成两个互斥的解释：
 (甲) $\Phi_0$ 在 $(\beta,0)$ 处真有 $\varphi$ 的奇次项 $c(\beta)\varphi$（$c$ 对 $\beta$ 为奇），
      那么我按 $\varphi^2$ 归一化看到的 $-24.8\,\beta$ 就是 $c(\beta)/\varphi$，
      而"$\varphi=0$ 是地板极小"这条律对任意非零 $\beta$ 都错——真正的极小在带符号的那一侧，
      我 #25/#26/#27 的角度网格只在 $\varphi\ge0$ 上跑，整个漏掉了这一支；
 (乙) $c\equiv0$，那么 $q(+2^\circ)=q(-2^\circ)$，$q$ 对 $\beta$ 必须是偶的，
      我那条线性拟合只是 $[10^{-4},0.05]$ 这段上的弯曲伪装。
两个判决都在这张表里，跑法：`python -X utf8 -m p0.exp_signed`。
"""
import sys
sys.stdout.reconfigure(encoding='utf-8')
import numpy as np
import p0.exp_plants_oos as E35
from p0.ticket import floor0, floor0_safe

NPROW = 3


def plant(c13=0.0, au=.6, b3=0.0, wu=1.0):
    A = np.array([[1.25, .30, c13], [0., .85, 0.], [0., 0., au]])
    B = np.array([[1.], [0.], [b3]])
    W = np.diag([1.0, .7, wu])
    return E35.Pl(A, B, W, np.eye(3), np.array([[1.0]]))


def wrow(deg):
    a = deg * np.pi / 180.
    return np.array([[np.cos(a), 0., np.sin(a)]])


FD = np.array([[1., 0., 0.]])
BETAS = [0.0, 1e-4, 1e-2, 0.05, 0.2, 0.5]
DEG = [1.0, 2.0, 5.0, 10.0, 45.0]

print(r'== E50：带符号角度 + 对称性判决（推翻或坐实 #27 的"$q$ 随 $b_3$ 一阶"）==')
print('')
print(r'[1] 恒等式 $\Phi_0(\beta,\varphi)=\Phi_0(-\beta,-\varphi)$。若这条就不过，我的求解器和我的对称性论证里已经有一个在骗人。')
print(r'  %-8s %-8s %-14s %-14s %-14s %-14s' % ('$b_3$', '$\\varphi$', '$\\Phi(+b,+\\varphi)$', '$\\Phi(-b,-\\varphi)$', '$\\Phi(+b,-\\varphi)$', '$\\Phi(-b,+\\varphi)$'))
worst = 0.
for b3 in BETAS:
    for dg in DEG:
        p1 = floor0_safe(plant(b3=+b3), wrow(dg))
        p2 = floor0_safe(plant(b3=-b3), wrow(-dg))
        p3 = floor0_safe(plant(b3=+b3), wrow(-dg))
        p4 = floor0_safe(plant(b3=-b3), wrow(dg))
        if None in (p1, p2, p3, p4):
            print(r'  %-8.3g %-8.1f  %-14s %-14s %-14s %-14s  DARE 无解'
                  % (b3, dg, '--', '--', '--', '--'))
            continue
        f1, f2, f3, f4 = p1[0], p2[0], p3[0], p4[0]
        r1 = abs(f1 - f2)
        r2 = abs((f3 - f4))
        worst = max(worst, r1 / abs(f1), r2 / abs(f1))
        print(r'  %-8.3g %-8.1f  %-14.12f %-14.12f %-14.12f %-14.12f  残差 %.2e / %.2e'
              % (b3, dg, f1, f2, f3, f4, r1, r2))
print(r'  跨全部格子的最大相对残差 = %.2e（E48 的杠子在 $10^{-15}$ 档，这条若成立就是解析恒等式而非巧合）' % worst)
print('')

print(r'[2] 把偶部与奇部分开：$a(\beta)=[\Phi(+\varphi)+\Phi(-\varphi)-2\Phi(0)]/(2\varphi^{2})$，$c(\beta)=[\Phi(+\varphi)-\Phi(-\varphi)]/(2\varphi)$，$\varphi=2^\circ$。')
print(r'  二次模型给出的极小在 $\varphi^{*}=-c/(2a)$、深度 $c^{2}/(4a)$；门票 $\mathrm{gap}_{\rm del}=\Phi_0(F_{\rm del})-\Phi_0(F_2)$ 同行列出对照。')
print(r'  %-8s %-13s %-11s %-12s %-12s %-9s %-9s %-12s %-12s %-12s %-9s'
      % ('$b_3$', '$\\Phi_0(0)$', '$a(\\beta)$', '$c(\\beta)$', '$c/b_3$', '$\\varphi^*$', '黄金$\\varphi^*$', '黄金深度', '模型 $c^2/4a$', '门票', '深度/门票'))
F2 = np.array([[1., 0., 0.], [0., 0., 1.]])
PH = 2. * np.pi / 180.


def golden_min(fun, lo, hi, it=70):
    r"""一维黄金分割，返回 $(\varphi_{\min},\ \Phi_{\min})$，$\varphi$ 以度计。"""
    gr = (np.sqrt(5.) - 1.) / 2.
    a, b = lo, hi
    c, d = b - gr * (b - a), a + gr * (b - a)
    fc, fd = fun(c), fun(d)
    for _ in range(it):
        if fc < fd:
            b, d, fd = d, c, fc
            c = b - gr * (b - a)
            fc = fun(c)
        else:
            a, c, fc = c, d, fd
            d = a + gr * (b - a)
            fd = fun(d)
    x = (a + b) / 2.
    return x, fun(x)


rows = []
for b3 in [0.0, 1e-6, 1e-4, 1e-2, 0.05, 0.2, 0.5, 1.0]:
    Plx = plant(b3=b3)
    base = floor0(Plx, FD)[0]
    fp = floor0(Plx, wrow(+2.))[0]
    fm = floor0(Plx, wrow(-2.))[0]
    a_ = (fp + fm - 2 * base) / (2 * PH ** 2)
    c_ = (fp - fm) / (2 * PH)
    star = -c_ / (2 * a_)
    mod = c_ ** 2 / (4 * a_)
    gg = floor0_safe(Plx, F2)
    gap = base - gg[0] if gg is not None else None
    gdeg, gval = golden_min(lambda g: floor0(Plx, wrow(g))[0], -15., 15.)
    rows.append((b3, a_, c_, star, mod, gval - base, gap, gdeg))
    print(r'  %-8.3g %-13.9f %-11.5f %-12.4e %-12s %-9.3f %-9.3f %-12.4e %-12.4e %-12s %-9s'
          % (b3, base, a_, c_, '--' if b3 == 0 else '%.4f' % (c_ / b3), star * 180 / np.pi, gdeg,
             gval - base, mod,
             '--' if gap is None else '%.4e' % gap,
             '--' if (gap is None or gap <= 0) else '%.4f' % (-(gval - base) / gap)))
print('')

print(r'[3] 判决')
a0 = rows[0][1]
print(r'  偶部：$a(0)=%.4f$，其余 $\beta$ 的 $a(\beta)$ 见上表。若 $a$ 几乎是常数，那 #27 说的"$q(2^\circ)$ 随 $\beta$ 一阶降"就整个是奇次项 $c/\varphi$ 冒充的——我那个 $-24.8$ 的斜率应当等于 $c/b_3$ 除以 $2^\circ$（弧度）。' % a0)
print(r'  $c/b_3$ 一列若在各 $\beta$ 之间稳定，就坐实 $c(\beta)$ 对 $\beta$ 是奇且一阶；$\beta=0$ 那行 $c$ 应当恰好是 $0$（对称性给的免费对照）。')
cs = [r[2] for r in rows if r[0] > 0]
bs = [r[0] for r in rows if r[0] > 0]
print(r'  我原先的网格是 $[0,30^\circ]$、步长 $1^\circ$。二次模型说极小在 $\varphi^{*}=%.2f\,\beta$（度）处，于是 $\beta<%.3f$ 时最优点整个落在步长以下、被网格读成 $0^\circ$——#27 那条"门槛 $\approx0.16$"若是这么来的，它就是我的网格伪影，不是物理。'
      % (rows[4][3] * 180 / np.pi / rows[4][0], 1.0 / (rows[4][3] * 180 / np.pi / rows[4][0])))
print(r'  最后一列把两件事放在一起：最优旋转实际买到的深度，与 #22 删通道的门票。若它在小 $\beta$ 上贴着一个常数（尤其 $\approx1$），"往有用方向转"与"把那条通道删掉"就是同一笔钱的两面；若它发散或乱跳，我 #27 里 $3.82$ 与 $3.912$ 相差 $2.3\%$ 就纯属巧合，那条猜测可以直接埋了。')
print('')
print(r'[4] 与 #25/#26/#27 的关系：哪句要改')
print(r'  #25 的"$\varphi=0$ 是地板的非退化二次极小"：在 $\beta=0$ 由对称性成立（[1] 那几行残差恰好 $0$，不是"小"），$\beta\ne0$ 时被奇次项取代——极小仍在、但挪到 $\varphi^{*}\propto\beta$，二阶曲率不动。')
print(r'  #26 的适用条件 $\operatorname{ran}\Theta\subseteq S_1$：它管的是"转了要不要付钱"的*符号*，不是曲率。奇次项的存在说明这条律的失效方式是"最优点平移"，比我在 #26 写的"整个反过来"要温和。')
print(r'  #27 的"$b_3\approx0.158$ 门槛"和"门票二阶、地板一阶"：两条都要撤回，正确形式在 [2]/[3]。')

print('')
print(r'[5] 回收率 $-(\text{黄金深度})/\mathrm{gap}_{\rm del}$ 是植物的泛函还是常数？换 $(a_u,w_u)$ 复测')
print(r'  %-9s %-9s %-8s %-11s %-11s %-11s %-11s %-11s' % ('$a_u$', '$w_u$', '$b_3$', '$a$', '$c/b_3$', '黄金深度', '门票', '回收率'))
for au in [0.2, 0.6, 0.9]:
    for wu in [0.1, 1.0, 4.0]:
        for b3 in [0.01, 0.05]:
            Plx = plant(au=au, wu=wu, b3=b3)
            base = floor0(Plx, FD)[0]
            fp = floor0(Plx, wrow(+2.))[0]
            fm = floor0(Plx, wrow(-2.))[0]
            a_ = (fp + fm - 2 * base) / (2 * PH ** 2)
            c_ = (fp - fm) / (2 * PH)
            gg = floor0_safe(Plx, F2)
            if gg is None:
                print(r'  %-9.2f %-9.2f %-8.2f  控制/估计 DARE 无解' % (au, wu, b3))
                continue
            gap = base - gg[0]
            gdeg, gval = golden_min(lambda g: floor0(Plx, wrow(g))[0], -15., 15.)
            depth = gval - base
            print(r'  %-9.2f %-9.2f %-8.2f %-11.4f %-11.4f %-11.4e %-11.4e %-11.4f'
                  % (au, wu, b3, a_, c_ / b3, depth, gap, -depth / gap if gap > 0 else float('nan')))
print(r'  读法：回收率若在九个植物上都贴着同一个数，那"一条倾斜的单通道几乎等于两条正交通道"就有普适味道，下一步该问它等于什么；')
print(r'  若彼此差得开，$0.9778$ 只是默认植物的泛函，和 #25 的 $0.1832$ 同级别——那也不坏，它给了"删方向"与"转方向"一个可算的换汇比率。')
print('')
print(r'[6] 对称性破缺：$c_{13}\ne0$ 时 $T$ 不再使 $A$ 不变，预测是 $\beta=0$ 处也长出一次项 $c$、极小离开 $0^\circ$。')
print(r'  同时查 #25 那一行 $c_{13}=0.3$ 到底为什么整表 `--`：是 DARE 无解、还是等预算点配不出来。')
print(r'  %-8s %-8s %-11s %-12s %-12s %-9s %-12s %-12s %-14s %-11s' % ('$c_{13}$', '$b_3$', '$a$', '$c$', r'$\varphi^*$', r'黄金$\varphi^*$', '黄金深度', '门票', r'$\Phi_0(F_2)$ 状态', r'$\|\Theta e_3\|/\|\Theta\|$'))
for c13 in [0.0, 0.05, 0.1, 0.2, 0.3, 0.5]:
    for b3 in [0.0, 0.01]:
        Plx = plant(c13=c13, b3=b3)
        base = floor0(Plx, FD)[0]
        fp = floor0(Plx, wrow(+2.))[0]
        fm = floor0(Plx, wrow(-2.))[0]
        a_ = (fp + fm - 2 * base) / (2 * PH ** 2)
        c_ = (fp - fm) / (2 * PH)
        gg = floor0_safe(Plx, F2)
        st = 'ok' if gg is not None else 'DARE 无解'
        gap = base - gg[0] if gg is not None else None
        gdeg, gval = golden_min(lambda g: floor0(Plx, wrow(g))[0], -25., 25.)
        s3 = np.linalg.norm(Plx.Th[:, 2]) / np.linalg.norm(Plx.Th)
        print(r'  %-8.2f %-8.2f %-11.4f %-12.4e %-12.4e %-9.3f %-12.4e %-12s %-14s %-11.3e'
              % (c13, b3, a_, c_, -c_ / (2 * a_) * 180 / np.pi, gdeg, gval - base,
                 '--' if gap is None else '%.4e' % gap, st, s3))
print(r'  读法：$c$ 若在 $\beta=0$ 就随 $c_{13}$ 线性非零，那"$\varphi=0$ 是极小"的完整条件是 $\beta=0$ 且 $c_{13}=0$ 两条，#25 的 [2]/[3] 少写了一条。')
print(r'  门票那一列若给得出数，$c_{13}=0.3$ 在 #25 的 `--` 就不是可行性问题而是我的等速率配点失败——两者要分清。')
print('')
print(r'[7] $c_{13}=0.3$ 的 bracket 诊断：把 $x_t=10^{-6}$ 那一格的地板与目标速率都打出来')
for c13 in [0.0, 0.2, 0.3]:
    Plx = plant(c13=c13)
    for F, lab in [(FD, r'$F_{\rm del}$'), (F2, r'$F_2$'), (wrow(2.), r'$w(2^\circ)$'), (wrow(-2.), r'$w(-2^\circ)$')]:
        got = floor0_safe(Plx, F)
        print(r'  $c_{13}=%.1f$  %-12s $\Phi_0$=%s' % (c13, lab, '--' if got is None else '%.6f' % got[0]))
print('')
print(r'[8] 把 [7] 里那个"$\Phi_0(F_2)$ 对 $c_{13}$ 一动不动"升成可验算的闭式。')
print(r'  断言（$b_3=0$，$a_{12}=A_{12}$、$a_{22}=A_{22}$、$w_k=W_{kk}$）：$F_2$ 的两行是 $e_1,e_3$，硬约束后验只剩 $e_2$ 一个分量，$P_0=p\,e_2e_2^{\top}$。')
print(r'  先验 $\tilde P=A(P_0)A^{\top}+W$ 只用到 $A$ 的第二列 $(a_{12},a_{22},0)^{\top}$，$c_{13}$ 在第三列、永远乘不到，故')
print(r'  $\tilde P_{22}=a_{22}^{2}p+w_2$、$\tilde P_{12}=a_{12}a_{22}p$、$\tilde P_{11}=a_{12}^{2}p+w_1$，而 $x_1$ 已知精确意味着后验要把 $\tilde P_{12}^{2}/\tilde P_{11}$ 扣掉：')
print(r'  $p=a_{22}^{2}p+w_2-\dfrac{a_{12}^{2}a_{22}^{2}p^{2}}{a_{12}^{2}p+w_1}$，整理成一元二次')
print(r'  $$a_{12}^{2}p^{2}+\big[(1-a_{22}^{2})w_1-a_{12}^{2}w_2\big]p-w_1w_2=0 ,\qquad \Phi_0(F_2)=\Theta_{22}\,p .$$')
print(r'  极限对照：$a_{12}\to0$ 时 $p\to w_2/(1-a_{22}^{2})$（解耦的那株）。$\Theta_{22}$ 与 $c_{13}$ 无关是因为 $b_3=0$ 时控制 DARE 的 $(e_1,e_2)$ 块自闭（$A$ 上三角、$B=e_1$）。')
print(r'  %-8s %-8s %-13s %-13s %-13s %-11s %-11s' % (r'$c_{13}$', r'$b_3$', r'$\Phi_0(F_2)$', r'闭式 $\Theta_{22}p$', '残差', r'$p$', r'$p_{a_{12}=0}$'))
for c13, b3 in [(0.0, 0.0), (0.1, 0.0), (0.3, 0.0), (0.5, 0.0), (0.0, 0.01), (0.0, 0.05), (0.3, 0.01), (0.3, 0.05)]:
    Plx = plant(c13=c13, b3=b3)
    got = floor0_safe(Plx, F2)
    a12, a22 = Plx.A[0, 1], Plx.A[1, 1]
    w1, w2 = Plx.W[0, 0], Plx.W[1, 1]
    bb = (1. - a22 ** 2) * w1 - a12 ** 2 * w2
    p = (-bb + np.sqrt(bb ** 2 + 4. * a12 ** 2 * w1 * w2)) / (2. * a12 ** 2)
    if got is None:
        print(r'  %-8.2f %-8.2f  F2 地板 DARE 无解' % (c13, b3))
        continue
    closed = Plx.Th[1, 1] * p
    print(r'  %-8.2f %-8.2f %-13.9f %-13.9f %-13.3e %-11.6f %-11.5f'
          % (c13, b3, got[0], closed, abs(got[0] - closed), p, w2 / (1. - a22 ** 2)))
print(r'  读法：$\beta=0$ 那四行的残差若在 $10^{-15}$ 档，"$\Phi_0(F_2)$ 与 $c_{13}$ 无关"就不是数值巧合，而是上面那三行代数——')
print(r'  删除门票 $\mathrm{gap}(c_{13})=\Phi_0(F_{\rm del})-\Phi_0(F_2)$ 的全部 $c_{13}$ 依赖都住在被删掉的那一行上（$x_3$ 的经 $a_{13}$ 漏进 $x_1$），实测 $\approx3.03\,c_{13}^{2}$，')
print(r'  与 $\beta$ 那一支的 $0.0495\,\beta^{2}$ 同阶不同系数。')
print(r'  最后一件是我没预料到的：$\beta\ne0$ 那四行闭式照样贴到 $10^{-16}$，我原先预测它会偏离。原因是地板的估计侧根本不碰 $B$——')
print(r'  $p$ 只吃 $(a_{12},a_{22},w_1,w_2)$，$\beta$ 与 $c_{13}$ 的影响全部集中在 $\Theta_{22}$。所以闭式比我写的前提大，"$\Phi_0(F_2)$ 与 $c_{13}$ 无关"应当收窄成"$b_3=0$ 时 $\Theta_{22}$ 与 $c_{13}$ 无关"这一条。')
print('')
print(r'[9] 把两支（$b_3$ 与 $c_{13}$）归一到同一个坐标：支撑倾斜 $t=\|\Theta e_3\|/\|\Theta\|$。')
print(r'  %-14s %-8s %-9s %-11s %-11s %-11s %-11s %-9s' % ('分支', r'$b_3$', '$c_{13}$', r'$t$', r'$a$', r'$c/t$', r'门票$/t^2$', '回收率'))
for b3, c13, lab in [(1e-4, 0.0, '$b_3$ 支'), (0.01, 0.0, '$b_3$ 支'), (0.05, 0.0, '$b_3$ 支'),
                     (0.0, 0.05, '$c_{13}$ 支'), (0.0, 0.10, '$c_{13}$ 支'), (0.0, 0.20, '$c_{13}$ 支'), (0.0, 0.30, '$c_{13}$ 支')]:
    Plx = plant(c13=c13, b3=b3)
    base = floor0(Plx, FD)[0]
    fp = floor0(Plx, wrow(+2.))[0]
    fm = floor0(Plx, wrow(-2.))[0]
    a_ = (fp + fm - 2 * base) / (2 * PH ** 2)
    c_ = (fp - fm) / (2 * PH)
    t_ = np.linalg.norm(Plx.Th[:, 2]) / np.linalg.norm(Plx.Th)
    gg = floor0_safe(Plx, F2)
    gap = base - gg[0]
    gdeg, gval = golden_min(lambda g: floor0(Plx, wrow(g))[0], -25., 25.)
    print(r'  %-14s %-8.3g %-9.2f %-9.4e %-11.4f %-11.4f %-11.4f %-9.4f'
          % (lab, b3, c13, t_, a_, c_ / t_, gap / t_ ** 2, -(gval - base) / gap))
print(r'  读法：$a$ 若在两支之间只差千分之几，那"$\varphi$ 方向的曲率"就确实是 $\operatorname{ran}\Theta$ 的几何量而不是旋钮的函数；')
print(r'  $c/t$ 与门票$/t^{2}$ 若也在几个百分点内对齐，#26 的条件 $\operatorname{ran}\Theta\subseteq S_1$ 就是唯一干活的结构量，两支的参数差异只是通往同一个 $t$ 的两条路。')
