r"""E58b = 修掉 E58 的两处口径问题后，把 C 的 §9.2 定理与 §9.5 的挑战一次做完。

E58 的 [2] 第二段"同 $D$ 对比"**作废**：我的二分区间取成 $[V_{\rm aff}\cdot10^{-3},\,3V_{\rm aff}]$，
而 $D_{\rm aff}\approx V\,\E[g'^2]=5.23V$，任务机制要 $V_{\rm task}\approx D_{\rm aff}$ 才能同 $D$，
上界 $3V_{\rm aff}$ 根本够不着 $\Rightarrow$ 二分撞到端点，打印的"差 0.0699"是端点值不是同 $D$ 差。
这里把区间改成以 $D$ 为目标、并在收敛后**核对** $D_{\rm task}$ 与 $D_{\rm aff}$ 的实测差。
第二处：小 $V$ 端 $R_{\rm task}-{\rm SLB}$ 变负（$R$ 不可能低于 SLB）$\Rightarrow$ 那是网格误差，
所以加一条 $n_x$ 加倍的收敛核对，并把定理检验限制在网格可分辨的 $V$ 上。
第三处：他的 $\varepsilon^4$ 系数用 sympy 精确对照（我手算 $2268$ vs 他印 $1964.25$），
并给出"该级数对任何 $\varepsilon>0$ 发散"的一行理由（矩 $(2k-1)!!$ 阶乘增长）。
"""
import sys
sys.stdout.reconfigure(encoding='utf-8')
import numpy as np
from math import log as ln, pi, e
LN2 = np.log(2.0)


def mech(eps, q, V, xmax=6.0, nx=0, hz_cap=0.002, kw_sig=7.0):
    """dens/dg 是**质量**（$\int=1$），kern 是**密度**样本并按 $\sum k\,h_z=1$ 修正截断，
    于是 $p=\sum_j \mathrm{dens}_j\,k(z_i-z_j)$ 直接是密度，$\sum p\,h_z=1$。"""
    if nx == 0:
        nx = int(max(20001, min(400001, np.ceil(2 * xmax / (np.sqrt(V) / 40.0)))))
        nx |= 1
    x = np.linspace(-xmax, xmax, nx)
    dx = x[1] - x[0]
    wx = np.exp(-0.5 * x ** 2) / np.sqrt(2 * pi) * dx          # 质量
    phi = x + q * x ** 3
    gx = x + eps * x ** 3
    sd = np.sqrt(V)
    hz = min(sd / 8.0, hz_cap)
    zlo, zhi = phi.min() - kw_sig * sd, phi.max() + kw_sig * sd
    ng = int(np.ceil((zhi - zlo) / hz)) + 1
    hz = (zhi - zlo) / (ng - 1)
    pos = (phi - zlo) / hz
    i0 = np.clip(np.floor(pos).astype(int), 0, ng - 2)
    fr = pos - i0
    dens = np.zeros(ng + 1); dg = np.zeros(ng + 1)
    np.add.at(dens, i0, wx * (1 - fr)); np.add.at(dens, i0 + 1, wx * fr)
    np.add.at(dg, i0, wx * gx * (1 - fr)); np.add.at(dg, i0 + 1, wx * gx * fr)
    dens, dg = dens[:ng], dg[:ng]
    kw = int(np.ceil(kw_sig * sd / hz))
    kk = np.arange(-kw, kw + 1)
    kern = np.exp(-0.5 * (kk * hz / sd) ** 2) / (sd * np.sqrt(2 * pi))
    kern = kern / (kern.sum() * hz)                            # Σk·hz = 1
    p = np.convolve(dens, kern, mode='same')                   # 密度
    m1 = np.convolve(dg, kern, mode='same') / np.maximum(p, 1e-300)
    mass = p.sum() * hz
    p = np.maximum(p / mass, 1e-300)
    Eg2 = float(np.sum(wx * gx ** 2))
    D = Eg2 - float(np.sum(p * m1 ** 2) * hz)
    hZ = -float(np.sum(p * np.log(p)) * hz)
    R = (hZ - 0.5 * np.log(2 * pi * e * V)) / LN2
    return D, R, dict(nx=nx, ng=ng, hz=hz, mass=mass)


def consts(eps, nx=400001, xmax=8.0):
    x = np.linspace(-xmax, xmax, nx); dx = x[1] - x[0]
    w = np.exp(-0.5 * x ** 2) / np.sqrt(2 * pi) * dx
    gp = 1 + 3 * eps * x ** 2
    EL = float(np.sum(w * np.log(gp))); AM = float(np.sum(w * gp ** 2))
    hX = 0.5 * np.log(2 * pi * e)
    return AM, EL, hX + EL, 0.5 * np.log2(AM) - EL / LN2      # AM, E ln g', h(g), ½log2 AM/GM


print(r'== E58b：定理本体（同 $D$ 口径修好）+ 级数系数 ==')
for eps in (0.3, 0.05, 0.01, 0.002):
    AM, EL, hG, gapT = consts(eps)
    print(r'\n[1] $\varepsilon=%.3g$：定理右边 $\tfrac12\log_2\mathrm{AM}/\mathrm{GM}=%.6f$ bit；'
          '$26\varepsilon^2=%.6f$（他的首项）；$E[g^{\prime2}]=%.6f$（解析 $1+6\varepsilon+27\varepsilon^2=%.6f$）'
          % (eps, gapT, 26 * eps ** 2, AM, 1 + 6 * eps + 27 * eps ** 2))
    Vs = np.geomspace(0.02, 2e-4, 7)
    print('    %-9s %-11s %-11s %-12s | %-11s %-11s %-12s' %
          ('V', 'D_aff', 'R_aff', 'R_aff−SLB', 'D_task', 'R_task', 'R_task−SLB'))
    for V in Vs:
        Da, Ra, d1 = mech(eps, 0.0, V)
        Dt, Rt, d2 = mech(eps, eps, V)
        SL = lambda D: (hG - 0.5 * np.log(2 * pi * e * D)) / LN2
        print('    %-9.4g %-11.6f %-11.6f %-12.6f | %-11.6f %-11.6f %-12.6f  (nx=%d/%d)' %
              (V, Da, Ra, Ra - SL(Da), Dt, Rt, Rt - SL(Dt), d1['nx'], d2['nx']))
    # 同 D 对比（以 D 为目标二分，并核对命中）
    print('    同 $D$ 的机制差（$R_{\rm aff}-R_{\rm task}$，$V_{\rm task}$ 由 $D$ 反解）：')
    for V in (0.02, 5e-3, 1e-3, 2e-4):
        Da, Ra, _ = mech(eps, 0.0, V)
        lo, hi = 1e-7, max(4 * Da, 2e-3)
        for _ in range(45):
            mid = np.sqrt(lo * hi)
            Dm, Rm, _ = mech(eps, eps, mid)
            if Dm > Da:
                hi = mid
            else:
                lo = mid
        Dm, Rm, _ = mech(eps, eps, mid)
        print('      $D=%.6f$：$R_{\rm aff}=%.5f$，$R_{\rm task}(同D)=%.5f$（$V_{\rm task}=%.5g$，'
              '实得 $D=%.6f$，失配 $%.1e$）差 $=%.5f$ vs 预言 $%.5f$'
              % (Da, Ra, Rm, mid, Dm, abs(Dm - Da) / Da, Ra - Rm, gapT))
    # 分辨率核对：nx 加倍
    V = 1e-3
    a = mech(eps, eps, V, nx=100001); b = mech(eps, eps, V, nx=400001)
    print('    分辨率核对 $V=10^{-3}$：$n_x$=1e5 → D=%.7f R=%.6f；$n_x$=4e5 → D=%.7f R=%.6f（差 %.1e/%.1e bit）'
          % (a[0], a[1], b[0], b[1], abs(a[0] - b[0]), abs(a[1] - b[1])))

print(r'\n[2] 他的 $\varepsilon$ 级数：精确系数对照 + 收敛半径')
import sympy as sp
x, ep = sp.symbols('x eps', positive=True)
gpr = 1 + 3 * ep * x ** 2
Eg2 = 1 + 6 * ep + 27 * ep ** 2                        # E[g'^2] 精确（矩已代）
# E[ln(1+3 eps x^2)] 的泰勒系数：sum_k (-1)^{k+1} 3^k E[x^{2k}]/k * eps^k，E[x^{2k}]=(2k-1)!!
import math
def EL_series(K):
    c = []
    for k in range(1, K + 1):
        mom = math.factorial(2 * k - 1) // (2 ** (k - 1) * math.factorial(k - 1))   # (2k-1)!!
        c.append((-1) ** (k + 1) * 3 ** k * mom / k)
    return c
def halflogAM_series(K):
    ser = sp.series(sp.log(Eg2) / 2, ep, 0, K + 1).removeO()
    return [ser.coeff(ep, k) for k in range(1, K + 1)]
K = 8
ELc = EL_series(K); AHc = halflogAM_series(K)
print('    gap·ln2 的系数 $a_k$（$=\tfrac12[\varepsilon^k]\ln E[g^{\prime2}]-[\varepsilon^k]E\ln g^\prime$）：')
for k in range(1, K + 1):
    a = AHc[k - 1] - ELc[k - 1]
    print('      $k=%d$: %14s   (除以 ln2: %12.5f)' % (k, sp.nsimplify(a), float(a) / math.log(2)))
print('    C 印的：$a_2/\ln2=18/\ln2=25.945$、$a_3/\ln2=-180/\ln2$、$a_4/\ln2=1964.25/\ln2$')
print('    我手算 $a_4=2268$。上表 $k=4$ 就是仲裁。')
print('    发散性：$|[\varepsilon^k]E\ln g^\prime|=3^k(2k-1)!!/k\to\infty$ 对每个 $\varepsilon>0$ ⇒ 零收敛半径。')
for e0 in (0.01, 0.002):
    terms = [abs(AHc[k - 1] - ELc[k - 1]) * e0 ** k / math.log(2) if k >= 1 else 0 for k in range(1, K + 1)]
    kn = int(np.argmin(terms)) + 1
    full = float(sum((AHc[k - 1] - ELc[k - 1]) * e0 ** k for k in range(1, K + 1)) / math.log(2))
    print('    $\varepsilon=%.3g$：逐项 $\to$ %s；最优截断在第 %d 项，截断和 %.6f' %
          (e0, np.array2string(np.array(terms), precision=3), kn, full))
