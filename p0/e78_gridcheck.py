r"""E78 = 审 C 的 c56（两阶段受控观测舞台）的率估计器：分辨率能不能伪造"非线性增益"。

背景（不是我的舞台，是我按同一公式独立复算它的度量）：C 的 c56 用
    phi(x)=c(x+q x^3), z=phi(x)+v, v~N(0,V),  x~N(0,1)
    I = (h(z) - 1/2 log(2 pi e V))/log 2      [加性高斯噪声下 = I(x;z)，这条我认可]
其中 h(z) 用梯形法在 **z 网格 = linspace(phi.min()-6sqrt(V), phi.max()+6sqrt(V), 41)** 上算
（c56 第 12 行 zg0n=41；pushforward 第 15-19 行；第 27-28 行取熵）。
疑点：节点间距 dz = span/(n-1) 随 |q| 爆炸（c=1,q=0.2 时 phi(7)=75.6 => span~160 => dz~4），
而 p(z) 的变化尺度是 sqrt(V)~0.7。**被比的两臂（q=0 与 q!=0）落在不同分辨率区间**，
于是"同一 cost 下非线性类 J_info 更低"这条判据可以被离散化偏差单独满足。
c45/c53/c54 的率都是同一个式子（build() 第 17-18 行 ng=1201），所以这条不只打在 c56 上。

判据（先写下，跑完不许改）：
 [0] 收敛参考：固定 (c,q,V)，报 n_z = 41/201/1001/6001/20001 的 I 序列。
     若 n_z=41 与收敛值差 >= 0.1 bit（C 报告的增益量级），[1] 才有意义。
 [1] 偏差量级：扫 C 的搜索盒 c in [0.2,4], q in {0} 或 {±0.1,±0.25,0.45,0.6}, V in [0.05,0.5,2]，
     报 bias = I(41) - I(收敛) 的中位/最大/符号分布，线性臂与非线性臂分开。
     判据：若两臂的 bias 中位数之差与 C 报告的"非线性增益"同量级 => 该判据不可用作证据。
 [2] 同功率对照：把 c 归一到 Var(phi)=1（匹配输出功率，c=1/sqrt(1+6q+15q^2)），
     用收敛网格算 I(q) 曲线。这才是"非线性感知省不省率"的可比量。
     预期（我写的先在这）：匹配功率下 I 随 |q| 单调下降，且 q=0 最大 —— 若如此，
     c56 里非线性臂的"省率"必须来自 cost 侧（任务形状），而不是来自率本身。
 [3] x 网格截断：c56 的 x 网格固定 [-7,7]，而 x1 的条件中心 a x0 + b u0(z0) 里
     u0=a0+a1 z0+a3 z0^3 随 z0 支撑（受 |q| 驱动）爆炸；报几个代表参数下中心超出网格多远。
"""
import sys
sys.stdout.reconfigure(encoding='utf-8')
import numpy as np

print('== E78：审 c56 的率估计器（网格分辨率 / 同功率可比量 / 截断）==')
ln2 = np.log(2.0)
XS = None


def I_of(c, q, V, nz, nx=801, xspan=7.0, chunk=120):
    """独立实现 c56 的率：p(z)=sum_x p(x)dx N(z;phi(x),V)，I=(h(z)-0.5log(2pi e V))/log2。
    按 x 分块累加，避免 nx x nz 的稠密大矩阵。"""
    xs = np.linspace(-xspan, xspan, nx)
    dx = xs[1] - xs[0]
    px = np.exp(-0.5 * xs ** 2) / np.sqrt(2 * np.pi) * dx
    ph = c * (xs + q * xs ** 3)
    zs = np.linspace(ph.min() - 6 * np.sqrt(V), ph.max() + 6 * np.sqrt(V), nz)
    pz = np.zeros(nz)
    s = np.sqrt(2 * np.pi * V)
    for i0 in range(0, nx, chunk):
        sl = slice(i0, min(i0 + chunk, nx))
        # 只保留对 pz 有贡献的 x（|z-phi|<10 sigma 太宽会退化，直接全算但按块）
        pz += px[sl] @ (np.exp(-0.5 * (zs[None, :] - ph[sl, None]) ** 2 / V) / s)
    hz = -np.trapezoid(pz * np.log(np.maximum(pz, 1e-300)), zs)
    return float((hz - 0.5 * np.log(2 * np.pi * np.e * V)) / ln2), zs[1] - zs[0]


REF = 6001
print('\n[0] 收敛参考 (c=1, V=0.5)')
for q in (0.0, 0.2, 0.45):
    conv, _ = I_of(1.0, q, 0.5, 20001, nx=1201)
    row = []
    for nz in (41, 201, 1001, 6001):
        I, dz = I_of(1.0, q, 0.5, nz)
        row.append('%d:%+.4f(dz=%.2f)' % (nz, I - conv, dz))
    print('  q=%.2f 收敛值 I=%+.4f bit | 相对偏差 %s' % (q, conv, '  '.join(row)))

print('\n[1] C 搜索盒上的 bias = I(n_z=41) - I(n_z=%d)' % REF)
arm = {}
for label, qs in (('线性臂 q=0', (0.0,)),
                  ('非线性臂 |q|>0', (0.1, -0.1, 0.25, -0.25, 0.45, 0.6))):
    bs = []
    worst = (0.0, None)
    for c in (0.2, 0.6, 1.0, 2.0, 4.0):
        for q in qs:
            for V in (0.05, 0.5, 2.0):
                ic, dzc = I_of(c, q, V, 41)
                idn, _ = I_of(c, q, V, REF, nx=601)
                d = ic - idn
                bs.append(d)
                if abs(d) > abs(worst[0]):
                    worst = (d, (c, q, V, dzc))
    b = np.array(bs)
    arm[label] = b
    print('  %-16s n=%3d | bias 中位 %+.4f | 均值 %+.4f | |bias|max %.4f bit | 正 %d 负 %d'
          % (label, len(b), float(np.median(b)), float(b.mean()), float(np.max(np.abs(b))),
             int((b > 0).sum()), int((b < 0).sum())))
    d, (c, q, V, dz) = worst
    print('      最差: c=%.1f q=%.2f V=%.2f dz=%.2f -> bias %+.4f bit' % (c, q, V, dz, d))
diff = float(np.median(arm['非线性臂 |q|>0']) - np.median(arm['线性臂 q=0']))
print('  => 两臂 bias 中位之差 = %+.4f bit（C 报告的率增益量级 0.14~0.47 bit）' % diff)

print('\n[2] 同功率（Var phi=1）下的 I(|q|) —— 率侧的可比量')
for V in (0.1, 0.5, 2.0):
    row = []
    for q in (0.0, 0.1, 0.25, 0.45, 0.6):
        m2 = 1.0 + 6.0 * q + 15.0 * q * q
        if m2 <= 0:
            row.append('q=%.2f:无解' % q); continue
        I, _ = I_of(1.0 / np.sqrt(m2), q, V, REF, nx=601)
        row.append('q=%.2f:%.4f' % (q, I))
    print('  V=%.1f | ' % V + '  '.join(row))

print('\n[3] x 网格截断（c56: x in [-7,7]，N=301）')
for (c0, q0, V0, a1, a3) in ((1.0, 0.0, 0.5, 1.5, 0.0), (1.0, 0.2, 0.5, 1.5, 0.0),
                            (1.0, 0.2, 0.5, 1.5, 0.3), (2.0, 0.6, 0.5, 2.5, 0.6)):
    xs = np.linspace(-7, 7, 3001)
    ph = c0 * (xs + q0 * xs ** 3)
    zmax = ph.max() + 6 * np.sqrt(V0)
    u0 = abs(a1 * zmax + a3 * zmax ** 3)
    print('  c0=%.1f q0=%.2f a1=%.1f a3=%.2f | z0 支撑 ±%.1f -> |u0| 可达 %.1f '
          '-> x1 中心可达 ±%.1f，超出 ±7 %s' % (c0, q0, a1, a3, zmax, u0, 7.0 + u0,
                                                '**（质量被归一化抹掉）**'))
