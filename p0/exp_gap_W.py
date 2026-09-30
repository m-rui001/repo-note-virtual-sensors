r"""E38：梯子的门槛怎样依赖过程噪声 $W$。

动机来自 #18 读 Tzoumas et al.（arXiv:1802.08376）的 Theorem 2：它的 NP-难实例是 $T=1$、
$A=B=C=Q=R=I$，没有动力学、没有率变量，难度全在信息矩阵求和的组合结构里。反过来我的梯子住在
DARE 里，那就必须问清楚：台阶是 $W$ 的规模造出来的，还是 $(A,F)$ 的几何造出来的？

$V=0$ 的预报协方差解 $\tilde P=A\tilde PA^\top+W-A\tilde PA^\top F^\top(F\tilde PF^\top)^{-1}F\tilde PA^\top$，
右端创新项是"二次 $\div$ 一次"，所以 $(\varepsilon W,V)\mapsto\varepsilon\tilde P(W,V/\varepsilon)$。三条推论，
每条都可测：

(1) $\Phi_0(\varepsilon W)=\varepsilon\Phi_0(W)$ **精确**（地板一次齐次）；$\operatorname{gap}_A$ 同样一次齐次。
(2) $C_{ij}=\operatorname{tr}(\Theta\,\partial P/\partial V_{ij})|_{V=0}$ 是**零次**齐次：
    $\partial_V[\varepsilon P(W,V/\varepsilon)]=\partial_{\tilde V}P(W,\tilde V)$。所以 $C$ 与 $\varepsilon$ 无关。
(3) 于是 $\mathrm{gap}_A/\lambda_i(C)\propto\varepsilon$：**横轴单位不动、门槛随 $W$ 缩小而塌向原点**，
    $\varepsilon\to0$ 时整副梯子合并成一条斜率 $\operatorname{rank}C/2$ 的垂直线。
    这就是"$W\to0$ 没有梯子"的定量形式，也正是 Tzoumas 那个 $T=1$ 实例所处的角落。

注意 (2) 与朴素猜想相反：如果连 $V$ 一起缩放才是位似，单独缩放 $W$ 只压门槛、不压通道尺度。
$C$ 的测量必须用**固定**的 $t$ 网格（$C$ 不随 $\varepsilon$ 移动，网格跟着动就等于把深端搬到浅端，
E36 的两点外推会在 $\varepsilon>1$ 处给出负的 $\lambda(C)$，这我已经踩过一次）。

脚本另加一条反向核对：恒等式 $\Phi(\varepsilon W,V)=\varepsilon\Phi(W,V/\varepsilon)$ 直接逐点测。
日志 `p0/e38_out.txt`。
"""
import sys
import numpy as np
sys.stdout.reconfigure(encoding='utf-8')
from scipy.linalg import solve_discrete_are
import p0.exp_plants_oos as E35
import p0.exp_ladder as X37
import p0.exp_C_exact as X36

sym = X36.sym


def rebuild(Pl, Wnew):
    return E35.Pl(Pl.A, Pl.B, Wnew, Pl.Q, Pl.R)


def Pt0_exact(Pl, F, cap=200000):
    r = F.shape[0]
    Pt0 = sym(solve_discrete_are(Pl.A.T, F.T, Pl.W, np.zeros((r, r))))
    Pt0, _Pk, _res, _k = X37.polish_Pt0(Pl, F, Pt0, cap=cap)
    return Pt0


def ladder(Pl, F, label, ts_c=None):
    fx = X37.floor_exact(Pl, F, None)
    if fx['闭式'] is None:
        print('  %s 地板不可解：%s' % (label, fx['闭式err']))
        return None
    phi0 = fx['闭式']
    old = X37.TS_C
    if ts_c is not None:
        X37.TS_C = ts_c
    C = X37.C_of(Pl, F)
    X37.TS_C = old
    fams, e, U, rho = X37.fams_of(Pl, F, C, phi0)
    return dict(phi0=phi0, C=C, elam=e, fams=fams, rho=rho)


print('== E38：门槛对 W 的依赖 ==')
base = None
for nm, pl, F in E35.plants():
    if nm.find('锚点') >= 0:
        base = (nm.replace('$', ''), pl, F)
        break
nm0, Pl0, F0 = base
n, r = Pl0.nn, F0.shape[0]
print('算例：%s，$n=%d$，$r=%d$；$C$ 用固定网格 $t=%s$'
      % (nm0, n, r, ' '.join('%.0e' % q for q in X37.TS_C)))

print('\n[1] 规模 $W\to\varepsilon W$。预言 $\Phi_0,\mathrm{gap}\propto\varepsilon$，$\lambda(C)$ 不变，'
      '故 $\\mathrm{gap}/\lambda\propto\varepsilon$')
print('  %-7s %13s %11s %13s   %-20s %s'
      % ('$\\varepsilon$', '$\\Phi_0$', '$\\Phi_0/\\varepsilon$', '$\\lambda(C)$', '$\\mathrm{gap}_A$',
         '$\\mathrm{gap}_A/\\lambda_i$'))
ref = None
for eps in (1e-3, 1e-2, 1.0, 1e2):
    L = ladder(rebuild(Pl0, eps * Pl0.W), F0, r'$\varepsilon=%g$' % eps)
    if L is None:
        continue
    if ref is None:
        ref = L
    gs = ' '.join('%.3e' % f['gap'] for f in L['fams'] if len(f['A']) == 1)
    rs = ' '.join('%s:%.3f' % (f['name'], f['gap'] / f['w'][0]) for f in L['fams'] if len(f['A']) == 1)
    print('  %-7.0e %13.6f %11.4f %13s   %-20s %s'
          % (eps, L['phi0'], L['phi0'] / eps, ' '.join('%.5g' % q for q in L['elam']), gs, rs))

Pa = Pt0_exact(Pl0, F0)
print('\n  一次齐次的逐点核对（$\\tilde P_0$、$\\Phi$ 恒等式、$C$ 零次齐次）：')
for eps in (1e-3, 1e2):
    Pe = rebuild(Pl0, eps * Pl0.W)
    Pb = Pt0_exact(Pe, F0)
    print('    $\\varepsilon=%-6.0e$ $\\|\\tilde P_0(\\varepsilon W)/\\varepsilon-\\tilde P_0(W)\\|/\\|\\tilde P_0\\|=%.2e$'
          % (eps, np.linalg.norm(Pb / eps - Pa) / np.linalg.norm(Pa)))
    dev = []
    for t in (1e-3, 1e-5):
        V = t * np.eye(r)
        x1, I1, _p1, _q1, _r1 = X37.pv(Pl0, F0, V / eps, 0.0)
        x2, I2, _p2, _q2, _r2 = X37.pv(Pe, F0, V, 0.0)
        dev.append(abs(x2 - eps * x1) / (eps * x1))
    print('    $\\varepsilon=%-6.0e$ 预算恒等式 $x(\\varepsilon W,V)$ vs $\\varepsilon x(W,V/\\varepsilon)$：相对偏差 %s'
          % (eps, ' '.join('%.2e' % q for q in dev)))
Ca = X37.C_of(Pl0, F0)
for eps in (1e-2, 1e2):
    Cb = X37.C_of(rebuild(Pl0, eps * Pl0.W), F0)
    print('    $\\varepsilon=%-6.0e$ $C$ 零次齐次：$\\|C(\\varepsilon W)-C(W)\\|/\\|C(W)\\|=%.2e$'
          % (eps, np.linalg.norm(Cb - Ca) / np.linalg.norm(Ca)))

print('\n[2] 方向 $W=R(\\theta)\\operatorname{diag}(1,10^{-3},10^{-6},10^{-9})R(\\theta)^\\top$'
      '（谱固定，只在 $(e_0,e_1)$ 平面转）')
print('  齐次性管不了方向：看 $\\mathrm{gap}_A/\\lambda_i$ 能否被转到合并')
sp = np.diag([1.0, 1e-3, 1e-6, 1e-9])


def rot(i, j, th):
    R = np.eye(n)
    c, s = np.cos(th), np.sin(th)
    R[i, i] = c; R[i, j] = -s; R[j, i] = s; R[j, j] = c
    return R


print('  %-6s %12s %10s %-22s %s' % ('$\\theta$', '$\\Phi_0$', '$\\kappa(C)$', '$\\lambda(C)$',
                                     '$\\mathrm{gap}_A/\\lambda_i$'))
for deg in (0, 15, 30, 45, 60, 75, 90):
    Wth = rot(0, 1, np.deg2rad(deg)) @ sp @ rot(0, 1, np.deg2rad(deg)).T
    L = ladder(rebuild(Pl0, Wth), F0, r'$\theta=%d$' % deg)
    if L is None:
        continue
    rs = ' '.join('%s:%.3f' % (f['name'], f['gap'] / f['w'][0]) for f in L['fams'] if len(f['A']) == 1)
    kap = L['elam'][0] / L['elam'][-1] if L['elam'][-1] > 0 else np.inf
    print('  %-6s %12.6f %10.2e %-22s %s'
          % ('%d°' % deg, L['phi0'], kap, ' '.join('%.4g' % q for q in L['elam']), rs))

print('\n[3] 反向：把 $\\varepsilon$ 拉到极大，看门槛是否超过可测区间（梯子被推走 $=$ 观测不到垂直段）')
for eps in (1e3, 1e4):
    L = ladder(rebuild(Pl0, eps * Pl0.W), F0, r'$\varepsilon=%g$' % eps)
    if L is None:
        continue
    rs = ' '.join('%s:%.3f' % (f['name'], f['gap'] / f['w'][0]) for f in L['fams'] if len(f['A']) == 1)
    print('  $\\varepsilon=%-6.0e$ $\\Phi_0/\\varepsilon=%.4f$  $\\lambda(C)$=%s   %s'
          % (eps, L['phi0'] / eps, ' '.join('%.5g' % q for q in L['elam']), rs))
