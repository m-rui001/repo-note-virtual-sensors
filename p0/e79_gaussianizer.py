r"""E79 = 打 C 的 tier (iv)："匹配（传感器输出）功率下，没有非线性观测映射能提高关于状态的信息率"。

这条对**高斯源**是真的，一句话可证；对**非高斯源**是假的，而且假的量可以闭式算出来：
    I(x; phi(x)+v) = h(phi(x)+v) - h(v) <= 1/2 log2(2 pi e (P+V)) - 1/2 log2(2 pi e V) = 1/2 log2(1+P/V),
P := Var(phi(x)) 是被"匹配功率"钉住的常数；等号成立 <=> phi(x)+v 高斯；v 高斯且与 phi(x) 独立，
由 Cramer 定理 <=> phi(x) 高斯。取 phi = 把 x 的分布运输到高斯的单调映射（Gaussianizer），
它严格单调故可逆（不销毁任何状态信息），于是上界**可达**。
=> 率最大化的传感器 = "把自己的输出高斯化"的那个映射，达到 AWGN 容量 1/2 log2(1+P/V)；
   identity 只达到容量减去源的非高斯性（negentropy）。identity 最优 <=> x 高斯。
C 在 (iv) 里扫的是**单调三次族** phi=c(t+q t^3)：三次增长的像太重尾，造不出 Gaussianizer
（重尾源到高斯需要次多项式增长），所以"族内 identity 最优"是**族的伪影**，不是定理。
更要紧：C 自己的 tier (ii) 用**非线性植物**，状态因此非高斯 —— (iv) 的前提在 (ii) 里不成立，
所以 (ii) 报的 0.14~0.17 bit 里可能混着这个高斯化基线，而不（只）是任务形状。

判据（跑之前先写下）：
 [0] 自检：高斯源 ⇒ Gaussianizer 应退化为线性，identity 率 == 容量 == Gaussianizer 率。
     对不上 = 我的数值实现有问题，后面全部作废。
 [1] 非高斯源（Laplace；两个不同混合度）：报 identity 率、单调三次族最优率（扫 q）、
     Gaussianizer 率、容量 1/2log2(1+P/V)。
     判据：Gaussianizer ≈ 容量（数值相等）且明显 > 三次族最优 => (iv) 作为一般陈述**错**。
 [2] 报"高斯化 − identity"的增益幅度（bit），与 C 在 (ii) 报的 0.14~0.17 bit 对表：
     同量级 => (ii) 的归因必须重做（先减掉高斯化基线，再谈任务形状）。
     另外核对：C 说他的源"kurtosis 最高 4.9"——若 kurt=3 的 Laplace 就已有大增益，
     则 (iv) 的失效阈值比 C 的表述更低。
"""
import sys
sys.stdout.reconfigure(encoding='utf-8')
import numpy as np
from scipy.stats import norm

ln2 = np.log(2.0)
print('== E79：匹配功率下的率最大化传感器 = Gaussianizer（审 C 的 tier (iv)）==')
print('   论断：max_{phi: Var(phi)=P} I(x;phi(x)+v) = 1/2 log2(1+P/V)，')
print('         最大值点 = Gaussianizer phi=sqrt(P) * Phi^-1(F_x(x))（严格单调=>可逆=>不损信息）；')
print('         identity 达到最大 <=> x 高斯。')

NX, NZ, CHUNK = 3001, 4001, 100


def rate(phi, V, lo, hi):
    """I(x; phi(x)+v) = h(z) - h(v)，z 密度 = x 密度与 N(0,V) 的卷积（网格实现）。"""
    xs = np.linspace(lo, hi, NX)
    dx = xs[1] - xs[0]
    px = PDF(xs) * dx
    ph = np.asarray(phi(xs), dtype=float)
    zs = np.linspace(min(ph.min(), lo) - 7 * np.sqrt(V),
                     max(ph.max(), hi) + 7 * np.sqrt(V), NZ)
    dz = zs[1] - zs[0]
    pz = np.zeros(NZ)
    s = np.sqrt(2 * np.pi * V)
    for i0 in range(0, NX, CHUNK):
        sl = slice(i0, min(i0 + CHUNK, NX))
        pz += px[sl] @ np.exp(-0.5 * (zs[None, :] - ph[sl, None]) ** 2 / V) / s
    hz = -np.trapezoid(pz * np.log(np.maximum(pz, 1e-300)), zs)
    return (hz - 0.5 * np.log(2 * np.pi * np.e * V)) / ln2, float(pz.sum() * dz)


def set_src(name):
    """返回 (lo, hi)。同时设置全局 PDF / Fx。源一律归一到 Var=1。"""
    global PDF, FX
    if name == 'gauss':
        PDF, FX = norm.pdf, norm.cdf
        return (-9.0, 9.0)
    if name == 'laplace':
        b = 1.0 / np.sqrt(2.0)                          # Var = 2b^2 = 1
        PDF = lambda t: np.exp(-np.abs(t) / b) / (2 * b)
        FX = lambda t: np.where(t < 0, 0.5 * np.exp(t / b), 1 - 0.5 * np.exp(-t / b))
        return (-12.0, 12.0)
    if name == 'mix18':                                  # 轻混合
        w, s2 = 0.93, 3.0
    else:                                                # mix33：重混合
        w, s2 = 0.9, 9.0
    PDF = lambda t: w * norm.pdf(t) + (1 - w) * norm.pdf(t / np.sqrt(s2)) / np.sqrt(s2)
    FX = lambda t: w * norm.cdf(t) + (1 - w) * norm.cdf(t / np.sqrt(s2))
    return (-20.0, 20.0)


SUM = {}
for name in ('gauss', 'laplace', 'mix18', 'mix33'):
    lo, hi = set_src(name)
    # 归一到 Var = 1：若 f0 的方差是 s2_0，则 g(t)=s0*f0(s0*t), s0=sqrt(s2_0)，才使 Var(g)=1
    ts = np.linspace(lo, hi, 200001)
    f0, F0 = PDF, FX
    p0 = f0(ts)
    Z0 = np.trapezoid(p0, ts)
    var0 = np.trapezoid(ts ** 2 * p0, ts) / Z0
    s0 = np.sqrt(var0)
    PDF = lambda t, f=f0, s0=s0: f(t * s0) * s0
    FX = lambda t, F=F0, s0=s0: F(t * s0)
    p = PDF(ts)
    Z = np.trapezoid(p, ts)                         # 截断质量（应 ~1）
    m2 = np.trapezoid(ts ** 2 * p, ts) / Z
    m4 = np.trapezoid(ts ** 4 * p, ts) / Z
    kurt = m4 / m2 ** 2
    assert abs(m2 - 1.0) < 5e-3, '源归一失败：Var=%.6f' % m2
    print('\n--- 源 %-8s：归一后数值 Var=%.6f（断言 <5e-3），kurtosis=%.3f，截断质量=%.6f'
          % (name, m2, kurt, Z))
    P = m2
    for V in (0.5, 2.0):
        cap = 0.5 * np.log2(1.0 + P / V)
        I_id, n_id = rate(lambda t: t, V, lo, hi)       # identity（功率已匹配为 P）
        assert I_id <= cap + 2e-3, 'identity 率 %.4f 超过容量上界 %.4f => 数值/归一有错' % (I_id, cap)
        # 单调三次族 phi = c(t + q t^3)，c 由功率匹配定，扫 q（范围对齐 C 的 [-0.6,0.6]）
        bestq, bestI = 0.0, -9e9
        for q in np.linspace(-0.6, 0.6, 25):
            u = ts + q * ts ** 3
            c = np.sqrt(P / np.trapezoid(u ** 2 * PDF(ts), ts))
            Iq, _ = rate(lambda t, c=c, q=q: c * (t + q * t ** 3), V, lo, hi)
            if Iq > bestI:
                bestI, bestq = Iq, q
        gz = lambda t: np.sqrt(P) * norm.ppf(np.clip(FX(t), 1e-11, 1 - 1e-11))
        vg = np.trapezoid(gz(ts) ** 2 * PDF(ts), ts)
        I_g, n_g = rate(gz, V, lo, hi)
        # 我自己的网格收敛检查（这正是我审 C 的那一条，先自查）
        NZ0, NX0, CH0 = NZ, NX, CHUNK
        NZ, NX, CHUNK = 2 * NZ0, 2 * NX0, max(1, CH0 // 2)
        I_id2, _ = rate(lambda t: t, V, lo, hi)
        I_g2, _ = rate(gz, V, lo, hi)
        NZ, NX, CHUNK = NZ0, NX0, CH0
        print('    V=%.1f | 容量=%.4f | Gaussianizer(Var=%.4f)=%.4f | 三次族最优(q*=%+.2f)=%.4f '
              '| identity=%.4f | 质量检查 %.6f/%.6f'
              % (V, cap, vg, I_g, bestq, bestI, I_id, n_id, n_g))
        print('        [自查] 网格加倍：identity %+0.1e | Gaussianizer %+0.1e bit'
              % (I_id2 - I_id, I_g2 - I_g))
        print('        [2]   高斯化−identity = %+.4f | 高斯化−三次族最优 = %+.4f | 容量−identity = %+.4f bit'
              % (I_g - I_id, I_g - bestI, cap - I_id))
        SUM[(name, V)] = (cap, I_g, bestI, I_id)

print('\n================ 判定 ================')
cap, Ig, Ib, Ii = SUM[('gauss', 0.5)]
print('[0] 高斯源自检：容量 %.4f vs Gaussianizer %.4f vs identity %.4f（差 %.1e / %.1e）'
      % (cap, Ig, Ii, abs(cap - Ig), abs(cap - Ii)))
print('    => 若三者相等（<1e-3）则数值实现正确，且与"identity 最优 <=> 源高斯"一致。')
for name in ('laplace', 'mix18', 'mix33'):
    for V in (0.5, 2.0):
        cap, Ig, Ib, Ii = SUM[(name, V)]
        print('[1] %-8s V=%.1f：Gaussianizer−容量 = %+.2e（应 ~0，上界可达）；'
              'Gaussianizer−三次族最优 = %+.4f bit => (iv) %s'
              % (name, V, Ig - cap, Ig - Ib,
                 '**被否**' if Ig - Ib > 5e-3 else '未被否'))
print('\n[2] 与 C 的 tier (ii)（0.144~0.173 bit）对表：')
for name in ('laplace', 'mix18', 'mix33'):
    cap, Ig, Ib, Ii = SUM[(name, 2.0)]
    print('    %-8s：高斯化基线 = %+.4f bit（V=2）%s'
          % (name, Ig - Ii, '<= 与 (ii) 同量级，(ii) 的归因需重做'
             if abs(Ig - Ii) > 0.05 else '小于 (ii)，(ii) 不能全部由此解释'))
print('\n判据回顾：Gaussianizer == 容量 > 三次族最优 => (iv) 只在高斯源成立，作为一般陈述为假；')
print('          且单调三次族**够不到**最优，所以 C 的"族内扫描"不构成对 (iv) 的证明。')
