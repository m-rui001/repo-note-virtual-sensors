r"""E51b：把 $c_{13}=0$ 的级数往高阶推（$P_4,P_5$），拿到 $h_4$ 的闭式，看 $\varphi=5^{\circ},10^{\circ}$ 的预测能进到哪。

动机来自 #30 §七第 2 条：[3] 里 $\Phi_0$ 的截断余项除以 $\varphi^{4}$ 是常数（$1.862$、$1.903$），
说明 $k\le3$ 的误差完全由 4 阶项主导；若把级数推到 $k\le5$，同一格子的残差应当改由 $\varphi^{6}$ 主导。
这一节只做两件事：(i) 精确有理算术跑 order=5；(ii) 用新的 $P_4,P_5$ 重算 [3]/[8] 的那两张表。
"""
import sys
sys.stdout.reconfigure(encoding='utf-8')
import numpy as np
import sympy as sp
from p0.exp_closed_tilt import tilt_series, plant, phi0, P0num, tay, ft, DEG

ORDER = 5
print(r'== E51b：$c_{13}=0$ 级数推到 $k\le5$（$P_4,P_5$ 与 $h_4$）==')
print(r'[1] 跑 order=%d 的级数机器（逐步线性回代，无猜根）。' % ORDER)
sol, mats, p22, res0, degs, lin, exact = tilt_series(1.25, .30, 0.0, .85, .6, 1.0, .7, 1.0, order=ORDER)
print(r'  验根残差（simplify 后应当全零）：%s' % str(res0))
print(r'  逐阶对新未知数是否线性：%s' % str(lin))
print(r'  %-6s %-30s %-22s' % ('系数', '精确值', '小数'))
for k in range(ORDER + 1):
    for nm in 'cde':
        v = sol[sp.Symbol('%s%d' % (nm, k))]
        print(r'  %-6s %-30s %-22.17f' % ('%s%d' % (nm, k), str(v)[:30], float(sp.N(v, 30))))
print('')

# ---------------------------------------------------------------- [2] 截断阶是否随阶走
print(r'[2] 截断阶：$\Phi_0$ 的残差除以下一个偶阶幂应当是常数。$k\le3$ 时除 $\varphi^{4}$，$k\le5$ 时除 $\varphi^{6}$。')
print(r'  %-9s %-14s %-14s %-14s %-14s' % ('phi(deg)', 'res(k<=3)', '/phi^4', 'res(k<=5)', '/phi^6'))
for dg in [0.5, 1.0, 2.0, 5.0, 10.0]:
    a = dg * DEG
    Plx = plant()
    num = phi0(Plx, np.tan(a))
    r3 = abs(num - float(np.trace(Plx.Th @ tay(mats[:4], a))))
    r5 = abs(num - float(np.trace(Plx.Th @ tay(mats, a))))
    print(r'  %-9.3g %-14.3e %-14.3e %-14.3e %-14.3e' % (dg, r3, r3 / a ** 4, r5, r5 / a ** 6))
print('')

# ---------------------------------------------------------------- [5] h4 的闭式：旧读数的偏差现在能算
print(r'[2b] $h_4$ 有闭式了：$\Phi_0^{\rm diff}(\varphi)=a+h_4\varphi^{2}+O(\varphi^{4})$ 里的 $h_4=\operatorname{tr}(\Theta P_4)$。')
print(r'  这一行把 #30 §三里"$2^\circ$ 差分把 $a$ 抬高 $5.8\times10^{-4}$"那件事变成可手算的数，不再依赖 Richardson。')
print(r'  %-7s %-19s %-14s %-14s %-14s' % ('$b_3$', '$h_4$ 闭式', '$h_4$ 拟合', '$a$ 闭式', '$a$(2deg) 差分'))
M4 = np.array(mats[4])
for b3 in [0.0, 0.2, 0.5]:
    Plx = plant(b3=b3)
    Th = Plx.Th
    h4 = float(np.trace(Th @ M4))
    f0 = phi0(Plx, 0.0)
    fp = phi0(Plx, np.tan(2 * DEG))
    fm = phi0(Plx, np.tan(-2 * DEG))
    a_meas = (fp + fm - 2 * f0) / (2 * (2 * DEG) ** 2)
    a_closed = float(np.trace(Th @ np.array(mats[2])))
    num = phi0(Plx, np.tan(2 * DEG))
    fit = abs(num - float(np.trace(Th @ tay(mats[:4], 2 * DEG)))) / (2 * DEG) ** 4
    print(r'  %-7.3g %-19.12f %-14.4f %-14.5f %-14.5f' % (b3, h4, fit, a_closed, a_meas))
print(r'  拟合的那列是 [2] 里 $\varphi=2^{\circ}$ 的 $\operatorname{res}/\varphi^{4}$；它与闭式 $h_4$ 应当在求解器精度内一致。')
print(r'  并且 $a$ 差分与闭式之差应当正好是 $h_4\varphi^{2}$：')
for b3 in [0.0, 0.5]:
    Plx = plant(b3=b3)
    Th = Plx.Th
    h4 = float(np.trace(Th @ M4))
    f0 = phi0(Plx, 0.0)
    fp = phi0(Plx, np.tan(2 * DEG))
    fm = phi0(Plx, np.tan(-2 * DEG))
    a_meas = (fp + fm - 2 * f0) / (2 * (2 * DEG) ** 2)
    a_closed = float(np.trace(Th @ np.array(mats[2])))
    print(r'  %-7.3g $a$ 差分 $-$ $a$ 闭式 $=%+.8e$，$h_4\varphi^{2}=%+.8e$，差 $=%+.2e$'
          % (b3, a_meas - a_closed, h4 * (2 * DEG) ** 2, (a_meas - a_closed) - h4 * (2 * DEG) ** 2))
print('')

# ---------------------------------------------------------------- [3] 逐元素（不经 Theta 投影）
print(r'[3] 绕开 $\Theta$：逐元素比后验矩阵，$\max\|P_0^{\rm num}-\sum_{k\le n}\varphi^{k}P_k\|$。')
print(r'  %-9s %-16s %-14s %-16s %-14s' % ('phi(deg)', 'n=3', '/phi^4', 'n=5', '/phi^6'))
Pl0 = plant()
for dg in [1.0, 2.0, 5.0, 10.0, 20.0]:
    a = dg * DEG
    num = P0num(Pl0, np.tan(a))
    d3 = np.abs(num - tay(mats[:4], a)).max()
    d5 = np.abs(num - tay(mats, a)).max()
    print(r'  %-9.2f %-16.3e %-14.3e %-16.3e %-14.3e' % (dg, d3, d3 / a ** 4, d5, d5 / a ** 6))
print('')

# ---------------------------------------------------------------- [4] 收敛半径的样子
print(r'[4] 大角度：级数还能用多远？（同一组 $P_k$，$\varphi$ 一路推到 $40^{\circ}$）')
print(r'  %-9s %-19s %-19s %-13s %-13s' % ('phi(deg)', 'Phi tay(k<=5)', 'Phi 求解器', '绝对差', '相对差'))
for dg in [5.0, 10.0, 15.0, 20.0, 30.0, 40.0]:
    a = dg * DEG
    Plx = plant(b3=0.5)
    num = phi0(Plx, np.tan(a))
    t5 = float(np.trace(Plx.Th @ tay(mats, a)))
    print(r'  %-9.2f %-19.12f %-19.12f %-13.3e %-13.3e' % (dg, t5, num, abs(t5 - num), abs(t5 - num) / abs(num)))
print('')
