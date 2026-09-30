r"""E81 = 把 C 的口号 "nonlinear sensing buys task-shaped information, not more information"
升级成定理，并顺带把 tier (iv) 的"没有非线性映射能提高率"钉死为**两重的族伪影**。

设定（任务版容量）：任务 g，传感器 z = phi(x) + v，v~N(0,V)，功率约束 Var(phi(x)) <= P。
    sup_phi  I(g(x); phi(x)+v) = 1/2 log2(1 + P/V)，
最优元 = **任务 Gaussianizer**  phi*(x) = sqrt(P) * Phi^{-1}( F_{g(x)}(g(x)) )。
证明三步：(i) I(g;z) <= h(z) - h(v) <= 1/2 log(2pi e(P+V)) - h(v) = 容量；
 (ii) phi* 使 phi*(x) 精确为 N(0,P)（概率积分变换）=> 第一个 <= 等号需要 I(x;z|g)=0，
      即 phi 在 g 的每个原像集（fiber）上为常数 —— phi* 只通过 g(x) 依赖 x，正是如此；
 (iii) 于是两个 <= 同时取等。
推论：
 (a) g 为同胚时 fiber 是单点，I(x;z|g)=0 自动成立，phi* 是 x 的单调映射；
     对高斯 x 它就是线性 => **identity 最优**（= E79 的 [0]，也 = C 的 (iv) 在高斯源那半）。
 (b) g 为 k-to-1 时（平方、取模、饱和、量化），identity 的 fiber 熵 I(x;z|g) > 0，
     最多可达 log2 k bit；一个**非线性且故意不单调**的传感器能把它整个拿回来。
     => 在 C 的 tier (iv) 说"没有非线性映射能提高率"的地方，只要允许"关于任务"的率，
        就有一个 ~1 bit 量级的**正**增益，而且闭式给得出。
 (c) C 扫的族 phi = c(x + q x^3) 是**奇函数**：它既造不出 Gaussianizer（E79），
     也永远无法折叠符号（本实验）—— 族伪影有两重，不是一重。

判据（先写下）：
 [0] 对照组 g = x+0.3x^3（同胚）、x~N(0,1)：任务 Gaussianizer 的 I(g;z) 应 == 容量，
     identity 的 I(g;z) 应 == 它的 E79 值（<容量，因为 g(x) 非高斯）。折叠不该带来额外增益。
 [1] 决定性组 g = x^2：报 I(g;任务Gaussianizer)、I(g;identity)、I(g;最优奇三次)、容量。
     判据：任务Gaussianizer == 容量 > identity 与最优奇三次 => tier (iv) 的一般陈述**否**，
           且增益幅度 = identity 丢掉的 fiber 熵 I(x;z|g)。
 [2] 报丢掉的 bit 与 1 bit（符号熵）的关系：小 V 时应接近 1 bit，大 V 时应趋 0。
"""
import sys
sys.stdout.reconfigure(encoding='utf-8')
import numpy as np
from scipy.stats import norm

ln2 = np.log(2.0)
NX, CHUNK = 6001, 200
print('== E81：任务版容量  sup_phi I(g(x); phi(x)+v) = 1/2 log2(1+P/V)，最优元 = 任务 Gaussianizer ==')


def hz_of(phi, V, lo, hi, nz):
    """输出密度熵 h(phi(x)+v)（nats）；x ~ N(0,1) 在 [lo,hi] 上离散。"""
    xs = np.linspace(lo, hi, NX)
    px = norm.pdf(xs) * (xs[1] - xs[0])
    ph = np.asarray(phi(xs), float)
    zs = np.linspace(min(ph.min(), -6.0) - 7 * np.sqrt(V),
                     max(ph.max(), 6.0) + 7 * np.sqrt(V), nz)
    pz = np.zeros(nz)
    s = np.sqrt(2 * np.pi * V)
    for i0 in range(0, NX, CHUNK):
        sl = slice(i0, min(i0 + CHUNK, NX))
        pz += px[sl] @ np.exp(-0.5 * (zs[None, :] - ph[sl, None]) ** 2 / V) / s
    dz = zs[1] - zs[0]
    return -np.trapezoid(pz * np.log(np.maximum(pz, 1e-300)), zs), float(pz.sum() * dz)


def hz_given_g(phi, V, nc=600, nz=3001):
    """g=x^2 的 fiber 熵：E_c h( 1/2 N(phi(c),V) + 1/2 N(phi(-c),V) )，c=|x| 服从折叠正态。"""
    # c 的密度 = 2*phi(c) 在 (0,8]；用 Gauss-Legendre
    u, w = np.polynomial.legendre.leggauss(nc)
    c = 4.0 * (u + 1.0)
    wc = 4.0 * w * 2 * norm.pdf(c)                     # 折叠正态在 (0,8) 上的权重
    mp = np.asarray(phi(c), float)
    mm = np.asarray(phi(-c), float)                    # 一般 != -mp：第一次跑硬编码了奇对称
    zlo = max(np.abs(mp).max(), np.abs(mm).max()) + 7 * np.sqrt(V)
    zg = np.linspace(-zlo, zlo, nz)
    a = np.sqrt(2 * np.pi * V)
    e = 0.0
    for i0 in range(0, nc, 50):
        sl = slice(i0, min(i0 + 50, nc))
        p = 0.5 * (np.exp(-0.5 * (zg[None, :] - mp[sl, None]) ** 2 / V)
                   + np.exp(-0.5 * (zg[None, :] - mm[sl, None]) ** 2 / V)) / a
        p = np.maximum(p, 1e-300)
        hc = -np.trapezoid(p * np.log(p), zg)
        e += float(np.sum(wc[sl] * hc))
    return e, float(np.sum(wc))


def power(phi, lo=-8.0, hi=8.0, n=200001):
    t = np.linspace(lo, hi, n)
    p = norm.pdf(t)
    return float(np.trapezoid(np.asarray(phi(t), float) ** 2 * p, t))


P = 1.0
print('\n[0] 对照：g = x + 0.3 x^3（同胚），x ~ N(0,1)。同胚时 F_{g(x)}(g(x)) = F_x(x)，')
print('    所以任务 Gaussianizer 对任何源都退化为"源的" Gaussianizer（对高斯源 = 线性）=> 折叠无额外增益')
eps = 0.3
gg = lambda x: x + eps * x ** 3
# 任务 Gaussianizer：phi*(x) = sqrt(P) * Phi^{-1}( F_{g(x)}(g(x)) )，用 g(x) 的分位数反演求其 CDF
qs = np.linspace(1e-7, 1 - 1e-7, 20001)
yq = gg(norm.ppf(qs))                                 # g 单调 => 这是 g(x) 的分位数函数
gz_raw = lambda x, yq=yq, qs=qs: norm.ppf(
    np.clip(np.interp(gg(x), yq, qs, left=qs[0], right=qs[-1]), 1e-11, 1 - 1e-11))
sgz = np.sqrt(P / power(gz_raw))
gz = lambda x, f=gz_raw, s=sgz: s * f(x)
gid = lambda x: x * np.sqrt(P / power(lambda t: t))
_tc = np.linspace(-3, 3, 200)
print('    自检：高斯源 + 同胚任务 => phi* 应是恒等（P=1）；max|phi*(x)-x| = %.2e，斜率 %.6f'
      % (float(np.abs(gz(_tc) - _tc).max()), float(np.polyfit(_tc, gz(_tc), 1)[0])))
for V in (0.5, 2.0):
    nz = 20001
    hz_id, _ = hz_of(gid, V, -8, 8, nz)
    hz_gz, _ = hz_of(gz, V, -8, 8, nz)
    hv = 0.5 * np.log(2 * np.pi * np.e * V)
    cap = (0.5 * np.log(2 * np.pi * np.e * (P + V)) - hv) / ln2
    I_id = (hz_id - hv) / ln2
    I_gz = (hz_gz - hv) / ln2                                # 同胚 => I(g;z)=I(x;z)
    print('    V=%.1f | 容量=%.4f | identity=%.4f | 任务Gaussianizer(Var=%.4f)=%.4f'
          '  | 增益 %+.4f bit（同胚对照：应 ~0）' % (V, cap, I_id, power(gz), I_gz, I_gz - I_id))

print('\n[1] 决定性：g = x^2（2-to-1）。I(g;z) = h(z) - E_c h(z|c)，c=|x|')
# 任务 Gaussianizer：phi*(x) = sqrt(P) Phi^{-1}(2Phi(|x|)-1)，偶函数，故意折叠
fold_raw = lambda x: norm.ppf(np.clip(2 * norm.cdf(np.abs(x)) - 1.0, 1e-11, 1 - 1e-11))
sfold = np.sqrt(P / power(fold_raw))
fold = lambda x, s=sfold: s * fold_raw(x)
gid2 = lambda x: x * np.sqrt(P / power(lambda t: t))
print('    phi* 是偶函数，Var(phi*)=%.6f（目标 %.4f）；它丢掉符号 => I(x;z|g)=0'
      % (power(fold), P))
# PIT 直检：phi*(x) 应精确 ~ N(0,P)。尖点(x->0 处 dphi/dx 发散)会让数值卷积偏 ~5e-3 bit，
# 所以把解析容量与经验分位数放在一起，作为"数值 vs 理论"的定标。
_s = np.linspace(-8, 8, 200001)
_u = np.sort(fold(_s))
_q = np.linspace(0.002, 0.998, 499)
_pit = float(np.abs(np.quantile(_u, _q) - np.sqrt(P) * norm.ppf(_q)).max())
print('    PIT 直检：phi* 经验分位 vs N(0,1) 分位，最大偏差 %.3e（>0 解释数值 h(z) 略超容量）' % _pit)
# 奇三次族（C 的族）：无法折叠 —— 逐 q 求 I(g;z)
def I_task(phi, V, nz=40001):
    hz, _ = hz_of(phi, V, -8, 8, nz)
    hg, wsum = hz_given_g(phi, V)
    return (hz - hg) / ln2, (hz - 0.5 * np.log(2 * np.pi * np.e * V)) / ln2, wsum


for V in (0.1, 0.5, 2.0):
    cap = 0.5 * np.log2(1 + P / V)
    It_id, Ix_id, _ = I_task(gid2, V)
    It_f, Ix_f, _ = I_task(fold, V)
    best = (-9e9, 0.0, -9e9)
    for q in np.linspace(-0.6, 0.6, 13):
        u = lambda t, q=q: t + q * t ** 3
        c = np.sqrt(P / power(u))
        it, ix, _ = I_task(lambda t, c=c, q=q: c * (t + q * t ** 3), V, nz=20001)
        if it > best[0]:
            best = (it, q, ix)
    print('    V=%.1f | 容量=%.4f | 任务Gaussianizer: I(g;z)=%.4f (I(x;z)=%.4f) '
          '| identity: I(g;z)=%.4f (I(x;z)=%.4f) | 最优奇三次 q*=%+.2f: I(g;z)=%.4f'
          % (V, cap, It_f, Ix_f, It_id, Ix_id, best[1], best[0]))
    print('        [1] 非线性(折叠)传感器增益 = I(g;z)_fold - I(g;z)_id = %+.4f bit'
          '        vs C 的 (iv)："没有非线性映射能提高率"' % (It_f - It_id,))
    print('        [2] identity 丢掉的 fiber 熵 I(x;z|g) = %+.4f bit（符号熵上限 1 bit；V 越小越接近）'
          % (Ix_id - It_id,))
print('\n判据回顾：任务 Gaussianizer 达到容量、奇三次族达不到 => (iv) 作为一般陈述为假；')
print('          但对同胚任务（他们全部数值实验所用的 g=x+eps x^3）折叠不带来增益 => 他们的数字不必改，')
print('          要改的是 (iv) 的**量词范围**（"我们测的族内" 而不是 "没有非线性映射"）。')
