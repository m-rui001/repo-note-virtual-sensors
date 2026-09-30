r"""E80 = 审 C 的 tier (i)：静态 Jensen 间隙闭式  ½log2( E[g'^2] / exp E[log g'^2] )。

C 的 note 里 Theorem 的**陈述**是：对任务 g（statement 里写 "g measurable"），
    R_lin - R_task = 1/2 log2 ( E[g'(x)^2] / exp E[log g'(x)^2] ) >= 0，等号 <=> g 仿射，
证明只有一行："h(g(x)) = h(x) + E[log g'(x)]"。

两件事要分开查：
  (A) 公式在他们**选的那族** g = x + eps x^3 上算得对不对（g' = 1+3eps x^2 > 0，是同胚，
      换元公式合法）。他们的数字：0.043561/0.127998/0.312347/0.474891/0.819021。
      => 这一条应该**承认**：如果我的独立复算到 6 位一致，公式在他们族内没错。
  (B) 陈述的**假设**够不够。换元公式要求 g 是 C^1 同胚（单射）。取 g(x)=x^2（非单射，
      且不是仿射）：
        - 公式给：AM/E[GM] = E[4x^2]/exp(E log 4x^2) => 间隙 = 1/2 log2(4/e^{-gamma}) = +1.416 bit > 0；
        - 真相：x -> x^2 丢掉符号位，h(x^2) = h(|x|) = h(x) - ln2，即任务所需率比传状态**少 1 bit**。
      => 若两者矛盾，则 (i) 的"gap >= 0，等号 <=> 仿射"在其陈述的假设类里**为假**，
         需要补"同胚"这条假设。这是可一行修复的**假设漏洞**，不是结果崩塌。
  (C) 顺手把"h(g(x)) = h(x) + E log|g'|"本身数值验证一次（同胚 vs 非同射各一例）。

判据（先写下）：
 [A] 我独立算的五个 eps 值与 C 的表差 < 1e-5 => 公式族内正确（承认）。
 [B] 对 g=x^2：报 公式值(+?)、h(x^2)-h(x)(?)。若公式 > 0 而 h(x^2) < h(x) => 陈述需补假设。
 [C] 换元式在同胚例子上 < 1e-4，在 x^2 上**不成立**（应差 ~ln2），证明我读对了它的前提。
"""
import sys
sys.stdout.reconfigure(encoding='utf-8')
import numpy as np
from scipy.optimize import brentq
ln2 = np.log(2.0)
print('== E80：C 的 tier (i) 闭式：独立复算 + 假设检查 ==')

# 用高密度 Gauss-Legendre 在 (-L,L) 上算 N(0,1) 测度下的积分
L = 12.0
u, w = np.polynomial.legendre.leggauss(4000)
xs = 0.5 * L * (u + 1.0) - 0.5 * L
ws = 0.5 * L * w * np.exp(-0.5 * xs ** 2) / np.sqrt(2 * np.pi)


def make_w(n):
    u, w = np.polynomial.legendre.leggauss(n)
    x = 0.5 * L * (u + 1.0) - 0.5 * L
    return x, 0.5 * L * w * np.exp(-0.5 * x ** 2) / np.sqrt(2 * np.pi)


print('\n[预检] 求积精度：光滑被积函数 log(1+3*0.3*x^2) 与对数奇异 log|2x| 的 n=1000/4000 漂移')
for n in (1000, 4000):
    xa, wa = make_w(n)
    s1 = float(np.sum(wa * np.log(1.0 + 0.9 * xa ** 2)))
    s2 = float(np.sum(wa * np.log(np.abs(2 * xa) + 1e-300)))
    print('    n=%4d | 光滑 %+.9f | 奇异 %+.9f（精确值 (ln4-gamma-ln2)/... = %+.9f）'
          % (n, s1, s2, np.log(2) + (-np.euler_gamma - np.log(2)) / 2))
xs, ws = make_w(4000)


def jensen_gap(gp2):
    am = float(np.sum(ws * gp2))
    gm = float(np.exp(np.sum(ws * np.log(np.maximum(gp2, 1e-300)))))
    return 0.5 * np.log2(am / gm), am, gm


print('\n[A] g = x + eps x^3, g\' = 1 + 3 eps x^2（同胚）：闭式 vs C 表')
print('    注意被积量是 g\'**2（第一次跑我把 g\' 直接传进去了，得到 0.1037@eps=0.3，是脚本错）')
C_TAB = {0.05: 0.043561, 0.1: 0.127998, 0.2: 0.312347, 0.3: 0.474891, 0.6: 0.819021}
worst = 0.0
for eps, ref in C_TAB.items():
    gp = (1.0 + 3 * eps * xs ** 2) ** 2
    g, am, gm = jensen_gap(gp)
    # 解析对照：E[g'^2] = 1 + 6eps + 9eps^2 * E[x^4] = 1 + 6eps + 27 eps^2
    am_ex = 1 + 6 * eps + 27 * eps ** 2
    worst = max(worst, abs(g - ref))
    print('    eps=%.2f  我算 %.6f   C %.6f   差 %+.2e   [AM 数值 %.6f vs 解析 %.6f]'
          % (eps, g, ref, g - ref, am, am_ex))
print('    最大绝对差 = %.2e  => [A] %s' % (worst, '一致，公式在其族内正确' if worst < 1e-5 else '不一致'))


def h_cont(vals, dens_w):
    """由 pushforward 样本的蒙特卡洛密度算微分熵（核方法会偏，这里用解析式对照即可）。"""
    return None


print('\n[B] g = x^2（非单射，非仿射）：闭式给什么，真相是什么')
g2 = (2.0 * xs) ** 2                      # g' = 2x, g'^2 = 4x^2
gap, am, gm = jensen_gap(np.maximum(g2, 1e-300))
print('    公式值 = 1/2 log2(AM/GM) = 1/2 log2(%.4f/%.6f) = %+.4f bit（C 的定理预言：严格正）'
      % (am, gm, gap))
print('    闭式核对：GM 解析值 = 2 e^{-gamma} = %.6f，数值 GM = %.6f（奇异积分求积误差 ~1e-3）'
      % (2 * np.exp(-np.euler_gamma), gm))
gap_exact = 0.5 * (1.0 + np.euler_gamma / ln2)
print('    精确间隙 = 1/2 log2(4 / (2 e^{-gamma})) = 1/2 (1 + gamma/ln2) = %+.6f bit（数值 %.6f）'
      % (gap_exact, gap))
# 真相：x^2 = |x|^2 且 sign 独立 => h(x^2) = h(|x|) = h(x) - ln2 => 任务所需率比传状态少 1 bit
hx = 0.5 * np.log(2 * np.pi * np.e)                 # h(N(0,1)) nats
# h(|x|) 直接积分：密度 2*phi(t), t>0
t = np.linspace(1e-9, 8.0, 400001)
pabs = 2 * np.exp(-0.5 * t ** 2) / np.sqrt(2 * np.pi)
h_abs = -np.trapezoid(pabs * np.log(np.maximum(pabs, 1e-300)), t)
print('    h(x) = %.6f nats, h(|x|) = %.6f nats  => h(x^2) - h(x) = %+.6f nats = %+.4f bit'
      % (hx, h_abs, h_abs - hx, (h_abs - hx) / ln2))
print('    符号判定：公式 %+.4f bit vs 真实 %+.4f bit => %s'
      % (gap, (h_abs - hx) / ln2,
         '**矛盾：(i) 的 "gap>=0，等号<=>仿射" 在其陈述的假设类里为假，需补 g 为同胚**'
         if gap > 0 and h_abs < hx else '不矛盾'))

print('\n[C] 换元式 h(g(x)) = h(x) + E[log|g\'|] 的前提（同胚 vs 非单射）')
for eps in (0.3,):
    gp = 1.0 + 3 * eps * xs ** 2
    Elog = float(np.sum(ws * np.log(gp)))
    gfun = lambda x, e=eps: x + e * x ** 3
    hy_est = hx + Elog
    print('    eps=%.1f 同胚：E log|g\'| = %+.6f，换元式预测 h(g(x)) = %.6f nats' % (eps, Elog, hy_est))
    # 独立构造 y 的密度（单调 => 反演），先用两个不依赖 Jacobian 的量校验这套构造
    yg = np.linspace(gfun(xs).min(), gfun(xs).max(), 400001)
    xin = np.interp(yg, gfun(xs), xs)
    pin = np.exp(-0.5 * xin ** 2) / np.sqrt(2 * np.pi) / (1.0 + 3 * eps * xin ** 2)
    m2_y = np.trapezoid(yg ** 2 * pin, yg)
    m2_x = float(np.sum(ws * gfun(xs) ** 2))
    print('           构造校验：∫p_y - 1 = %+.2e | E[y^2] 由 y 密度 %.6f vs 由 x 空间 %.6f（差 %+.2e）'
          % (np.trapezoid(pin, yg) - 1.0, m2_y, m2_x, m2_y - m2_x))
    hy_num = -np.trapezoid(pin * np.log(np.maximum(pin, 1e-300)), yg)
    print('           数值 h(g(x)) = %.6f nats  => 与换元式差 %+.2e（应 <1e-4）' % (hy_num, hy_num - hy_est))
# 非单射：g=x^2。正确的 pushforward 密度有两个原像；换元式只数了一个
Elog_abs = float(np.sum(ws * np.log(np.maximum(np.abs(2 * xs), 1e-300))))
print('    g=x^2 非单射：换元式预测 h = h(x)+E log|2x| = %.6f nats' % (hx + Elog_abs))
print('           真值 h(x^2) = %.6f nats  => 换元式高估 %+.6f nats（= %+.4f bit）'
      % (h_abs, (hx + Elog_abs) - h_abs, ((hx + Elog_abs) - h_abs) / ln2))
print('           => 换元式在非单射下**不成立**（漏了原像求和），C 那一行证明的前提必须写进定理。')

print('\n================ 判定 ================')
print('[A] 公式在 g=x+eps x^3 上复算到 %.1e => %s'
      % (worst,
         '承认：闭式正确，"verified to six digits" 站得住' if worst < 1e-5 else
         '不一致：需查明谁错（先怀疑我自己的求积/入参）'))
print('[B] g=x^2：公式给精确值 %+.6f bit（1/2(1+gamma/ln2)），真相 %+.6f bit => 符号相反'
      % (gap_exact, (h_abs - hx) / ln2))
print('    => 定理陈述缺假设：需 "g 为 C^1 同胚"（换元式 h(g(x))=h(x)+E log|g\'| 只在单射时成立）')
print('    注：他们的所有数值实验都用 g=x+eps x^3（同胚），所以这是**陈述漏洞**，不是数字错误；')
print('        但 note 里 "equality iff g is affine" 作为一般命题为假（正确说法：iff g\' 几乎必然为常数）。')
print('        非单射任务（平方、取模、量化、饱和）在实践中很常见，补上假设后 (i) 仍然完整成立。')
