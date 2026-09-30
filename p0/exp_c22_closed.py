"""E54：给 C 的 c22 一个解析判决 —— 截距移 $b$ 的闭式，和"0.183 泛函"的身份证。

C 在 12:16 的 `.work3/c22_n1.py` 做了两件事：
  [1] 逐角度拟合 ratio-1 $=(C\\cdot x_t)$，得 $C_\\varphi$ 中位 0.18375（我的 E49[3b] 律的形状成立）；
  [2] 把它换成截距语言：$b(\\varphi)-b(0)=\\tfrac12\\log_2(1+C\\delta)$，$\\delta=\\Phi_0(F_\\varphi)-\\Phi_0(F_0)$。
[2] 是新的说法。但它没有来源——$C$ 仍是拟合出来的数。

我这边有现成的闭式（E28 的推导 + E36 的 Fréchet/Stein 解 $C$），对 $r=1$ Family：
$$I_{\\rm TRV}(x)=\\tfrac12\\log_2\\tfrac1x+\\tfrac12\\log_2\\bigl(c(\\varphi)\\,m_0(\\varphi)\\bigr)+o(1),\\qquad
  m_0=F\\tilde P_0F^\\top,$$
所以**不需要任何假设**就有
$$b(\\varphi)-b(0)=\\tfrac12\\log_2\\frac{c(\\varphi)m_0(\\varphi)}{c(0)m_0(0)}.$$
与 C 的律对照，$C$ 的"泛函身份"立刻读出来：
$$\\boxed{\\;C_0=\\lim_{\\varphi\\to0}\\frac{c(\\varphi)m_0(\\varphi)/(c(0)m_0(0))-1}{\\Phi_0(F_\\varphi)-\\Phi_0(F_0)}\\;}$$
分子分母都是我已有的闭式量（$c$ 由 E36 的 Stein 解、$\Phi_0$ 由奇异 DARE），**零自由参数**。
本文件判三件事：(a) 这个解析 $b$ 移不移动得动 C 的实测；(b) $C_0$ 是否真为常数（若是，C 的律是推论；若否，
它只在 $\\delta$ 小的窗口成立）；(c) $C_0$ 解析值与 0.183 差多少。

口径两条，都是我自己踩过的坑，这里当场自检：$\\Theta$ 来自控制 DARE（$Q=I$，不是 $W$）；
$\\Phi$ 用**后验** $P$，不是预报 $\\tilde P$。
"""
import sys
import numpy as np
sys.stdout.reconfigure(encoding='utf-8')
import p0.exp_plants_oos as E35
import p0.exp_ladder as X37
import p0.exp_C_exact as X36
from p0.ticket import floor0_safe

sym = X36.sym
DEG = np.pi / 180.

A3 = np.array([[1.25, .30, 0.], [0., .85, 0.], [0., 0., .6]])
B3 = np.array([[1.], [0.], [0.]])
W3 = np.diag([1.0, .7, 1.0])
PL = E35.Pl(A3, B3, W3, np.eye(3), np.array([[1.0]]))
JC = PL.UNC


def wdir(p):
    return np.array([np.cos(p), 0., np.sin(p)])


def F_of(p):
    return wdir(p)[None, :]


def pieces(p):
    """$(\\Phi_0,\\ m_0,\\ c)$：地板、地板处的先验投影、单通道代价敏感度。"""
    F = F_of(p)
    phi, Pk, res = floor0_safe(PL, F)
    Pt0 = sym(PL.A @ Pk @ PL.A.T + PL.W)
    m0 = float(F @ Pt0 @ F.T)
    C = X37.C_of(PL, F)
    return phi, m0, float(C[0, 0]), res


ROWS = [0, 2, 5, 10, 20, 45, 85]
MEAS_DB = {2: 0.00063, 5: 0.00394, 10: 0.01582, 20: 0.06440, 45: 0.36312}
C_C22 = 0.18375

print('== E54：截距移的闭式 vs C 的 c22 [2] ==')
print('  植物 $\\beta=0,c_{13}=0$（= C 的 c22 那株）  $\\operatorname{tr}(WP_c)=%.4f$  $\\operatorname{rank}\\Theta=%d'
      % (JC, np.linalg.matrix_rank(PL.Th, tol=1e-9)))

d = {}
for q in ROWS:
    p = q * DEG
    phi, m0, c, res = pieces(p)
    d[q] = (phi, m0, c)
    print('  %3d°  $\\Phi_0$=%11.6f  $m_0=F\\tilde P_0F^\\top$=%10.6f  $c$=%11.6f  DARE残差 %.1e'
          % (q, phi, m0, c, res))

phi0, m00, c0 = d[0]
print('\n  [A] 解析截距移 $\\tfrac12\\log_2\\frac{c\\,m_0}{c_0m_{00}}$  vs  C 实测 $b_\\varphi-b_0$')
print('      同时给 C 的律 $\\tfrac12\\log_2(1+C_{22}\\delta)$（$C_{22}=%.5f$，纯拟合）' % C_C22)
print('  %-5s %-11s %-13s %-13s %-13s %-11s %-11s' %
      ('phi', '$\\delta$', '解析 db', 'C 实测 db', '拟合律 db', '解析/实测', '(cm/cm0c0-1)/δ'))
for q in ROWS[1:]:
    phi, m0, c = d[q]
    delta = phi - phi0
    db_an = .5 * np.log2(c * m0 / (c0 * m00))
    db_fit = .5 * np.log2(1 + C_C22 * delta)
    ratio = (c * m0 / (c0 * m00) - 1.) / delta
    meas = MEAS_DB.get(q)
    print('  %-5d %-11.5f %-13.6f %-13s %-13.6f %-11s %-11.6f'
          % (q, delta, db_an, '%.6f' % meas if meas else '   —   ',
             db_fit, '%.4f' % (db_an / meas) if meas else '  —  ', ratio))

# 奇偶性先定：这株植物 $c_{13}=0=b_3$，对角符号变换 $x_3\to-x_3$ 是系统的对称性，
# 故 $\Phi_0$ 与 $c m_0$ 都应是 $\varphi$ 的**偶**函数 —— 一阶导为 $0$，$C_0$ 是二阶项之比。
print('\n  [B0] 偶性自检（$\varphi=\pm h$ 同 Plant，只换方向符号）')
for hq in (1e-2, 2e-2):
    pa, _, _, _ = pieces(hq)
    mb, _, _, _ = pieces(-hq)
    _, m0a, ca, _ = pieces(hq)
    _, m0b, cb, _ = pieces(-hq)
    print('    $h=%.0e$：$\\Phi_0(\\pm h)$ 差 $%.2e$   $\\log(cm_0)$ 差 $%.2e$'
          % (hq, pa - mb, np.log(ca * m0a) - np.log(cb * m0b)))

print('\n  [B] $C_0$ 的解析极限：$\\bigl[\\log(cm_0)-\\log(c_0m_{00})\\bigr]/\\bigl[\\Phi_0-\\Phi_{00}\\bigr]$ 的 $\\varphi\\to0$')
print('      （分子分母同为 $O(\\varphi^2)$，两点 Richardson 在 $h^2$ 上剥掉 $O(\\varphi^4)$）')
phi85, m085, c85 = d[85]


def pair(h):
    phi_p, m0_p, c_p, _ = pieces(h)
    phi_m, m0_m, c_m, _ = pieces(-h)
    return (.5 * (phi_p + phi_m) - phi0, .5 * (np.log(c_p * m0_p) + np.log(c_m * m0_m)) - np.log(c0 * m00))


hs = [1e-2, 2e-2, 4e-2]
vals = []
for h in hs:
    dd, dg = pair(h)
    vals.append(dg / dd)
    print('    $h=%.2e$ rad：$\\delta=%.3e$  $\\Delta\\log(cm_0)=%.3e$  商 $=%.6f$' % (h, dd, dg, dg / dd))
c0a, c0b = vals[0], vals[1]
c0_lim = c0a + (c0a - c0b) / ((hs[1] / hs[0]) ** 2 - 1.)
print('    Richardson（$h^2$ 步）$\\to$ $C_0^{\\rm closed}=%.6f$' % c0_lim)
print('    C 的拟合中位 $C_{22}=%.5f$   → 相对差 $%.2e$' % (C_C22, c0_lim / C_C22 - 1.))
print('    E49 报的 $0.183$            → 相对差 $%.2e$' % (c0_lim / 0.183 - 1.))
print('    $\\varphi=85°$ 的商 $=%.5f$（漂移 $%.1f\\%%$）—— 这就是 c22 表里"跨角 2e-2"的来源：'
      % ((c * m0 / (c0 * m00) - 1.) / (phi - phi0), 100 * ((c * m0 / (c0 * m00) - 1.) / (phi - phi0) / c0_lim - 1.)))
print('    律 $\\tfrac12\\log_2(1+C\\delta)$ 是 $C_0$ 常数化的一阶截断，$O(\\delta^2)$ 项由上面的漂移直接读出。')

print('\n  [C] 判读（写死，跑出来再核）')
print('    (1) 若 [A] 的"解析/实测"列在 $10^{-3}$ 内为 1：C 的 c22[2] 是我的闭式的推论，不是新律；')
print('        0.183 的泛函身份 = $\\lim[\\Delta\\log(cf\\tilde P_0f^\\top)]/\\Delta\\Phi_0$，开放问题 6 关闭。')
print('    (2) 若 $C_0$ 的商随 $\\varphi$ 漂：C 的律只在 $\\delta$ 小端成立，漂的一阶就是 $\\delta^2$ 项，')
print('        而 c22 自己的表已经显示漂移（跨角 2e-2、跨预算 $x_t$ 6e-2），这条必须并报。')
print('    (3) 若 [A] 整列不对：错在我的 $c$（Richardson 外推）或 C 的 $b$ 拟合（c22 打印 $\\varphi=0$ 的斜率 0.3612，')
print('        它自己的"应 ~0"检查是红的）。两边都要报，不许只报好看的。')

print('\n  [D] 闭式反过来用：不扫描，直接给 $C_0$ 对植物参数的导数（#25/#26 的"对 $w_u$ 盲、被 $a_u$ 拉动"）')
print('      每个植物只解两次奇异 DARE + 两次 Stein，$h=0.02$ rad 对称差分。')


def C0_of(au, wu, h=2e-2):
    Plx = E35.Pl(np.array([[1.25, .30, 0.], [0., .85, 0.], [0., 0., au]]),
                 B3, np.diag([1.0, .7, wu]), np.eye(3), np.array([[1.0]]))
    def brk(pm):
        F = pm
        got = floor0_safe(Plx, F)
        if got is None:
            return None
        phi, Pk, _ = got
        Pt = sym(Plx.A @ Pk @ Plx.A.T + Plx.W)
        return phi, np.log(float(F @ Pt @ F.T * X37.C_of(Plx, F)[0, 0]))
    a0 = brk(F_of(0.))
    ap = brk(F_of(h))
    am = brk(F_of(-h))
    if a0 is None:
        return None, None, '不可行'
    phip, logp = ap
    phim, logm = am
    dd = .5 * (phip + phim) - a0[0]
    dg = .5 * (logp + logm) - a0[1]
    return dg / dd, dd, 'ok'


print('  %-6s %-6s %-13s %-11s %s' % ('$a_u$', '$w_u$', '$C_0$', '$\\delta(h)$', '备注'))
for wu in (0.5, 1.0, 2.0):
    for au in (.3, .45, .6, .75, .85):
        v, dd, tag = C0_of(au, wu)
        print('  %-6.2f %-6.2f %-13s %-11.3e %s'
              % (au, wu, '%.6f' % v if v is not None else '   —   ',
                 dd if dd is not None else 0., tag))
print('  读法：$w_u$ 那一列若三行完全一样，"$0.183$ 对 $w_u$ 盲"就从数值观察升为闭式事实；')
print('      若不同，则 #26 的"盲"是差分窗口太窄。两种都算把线程关掉，不再有第三种读法。')

print('\n  [E] 位似恒等式对 $C_0$ 的作用（条目 19 的 $\\Phi_0(\\varepsilon W)=\\varepsilon\\Phi_0$、$C(\\varepsilon W)=C$ 的直接推论）：')
print('      $\\log\\det C$ 只加常数、$\\delta$ 乘 $\\varepsilon$，所以预言是 $C_0(\\varepsilon W)=C_0(W)/\\varepsilon$，')
print('      即 $\\varepsilon C_0$ 应为常数。这是闭式给出的**可推翻预测**，不是拟合出来的。')


def C0_scaled(eps, au=.6, h=2e-2):
    Plx = E35.Pl(A3, B3, eps * W3, np.eye(3), np.array([[1.0]]))

    def brk(F):
        phi, Pk, _ = floor0_safe(Plx, F)
        Pt = sym(Plx.A @ Pk @ Plx.A.T + Plx.W)
        return phi, np.log(float(F @ Pt @ F.T * X37.C_of(Plx, F)[0, 0]))
    a0, p0 = brk(F_of(0.))
    dp = .5 * (brk(F_of(h))[0] + brk(F_of(-h))[0]) - a0
    dg = .5 * (brk(F_of(h))[1] + brk(F_of(-h))[1]) - p0
    return dg / dp


print('  %-8s %-12s %-12s' % ('$\\varepsilon$', '$C_0(\\varepsilon W)$', '$\\varepsilon C_0$'))
for eps in (.25, .5, 1.0, 2.0, 4.0):
    v = C0_scaled(eps)
    print('  %-8.2f %-12.6f %-12.6f' % (eps, v, eps * v))

print('\n  [F] [D] 的 $w_u$ 残差是 $O(h^2)$ 污染还是真依赖？把 $h$ 收三档并对 $h^2$ 外推')
print('  %-6s %-11s %-11s %-11s %-11s' % ('$w_u$', 'h=4e-2', 'h=2e-2', 'h=1e-2', '$h\\to0$（Richardson）'))
for wu in (0.5, 1.0, 2.0):
    vv = [C0_of(.6, wu, h=h)[0] for h in (4e-2, 2e-2, 1e-2)]
    lim = vv[2] + (vv[2] - vv[1]) / 3.
    print('  %-6.2f %-11.6f %-11.6f %-11.6f %-11.6f' % (wu, vv[0], vv[1], vv[2], lim))
print('  判据：外推值若三行合并到 $10^{-6}$，"对 $w_u$ 盲"是闭式事实；若仍差 $10^{-4}$，')
print('      它就是真依赖，#26 要改写成"盲到 $2\\times10^{-4}$ 相对量"。')
