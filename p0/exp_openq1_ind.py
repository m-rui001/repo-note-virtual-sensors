r"""E58 = 开放问题 1：用**独立路数**验 C 的 §9.2 定理，并给他的 §9.5 与级数一条可判定的话。

路数与他的自适应细网格、与 $p0/$ 里的 Gauss–Hermite 都不同源：把 $\varphi(X)$ 的推前测度用线性
插值打到均匀 $z$ 格上，与 $N(0,V)$ 核**离散卷积**（等效于带宽已知的 KDE），一次 FFT/`np.convolve`
同时给 $h(Z)$ 与 $E[g|Z]=\big(p_{g\text{-加权}}*N\big)/p$。控制检验：$q=0$ 时 $R$ 必须等于
$\tfrac12\log_2(1+1/V)$ 解析值。他 §9.3 自己登过的教训（GH 在小噪声下把支撑跳过）在这一路不存在，
因为格距直接由 $\sqrt V$ 定。
"""
import sys
sys.stdout.reconfigure(encoding='utf-8')
import numpy as np
from math import log as ln, pi, e
LN2 = np.log(2.0)


def mech(eps, q, V, xmax=5.0, nx=20001, hz_cap=0.002, kw_sig=6.0):
    """z = x + q x^3 + sqrt(V) n，g = x + eps x^3。返回 (D, R_bits, 诊断)."""
    x = np.linspace(-xmax, xmax, nx)
    dx = x[1] - x[0]
    wx = np.exp(-0.5 * x ** 2) / np.sqrt(2 * pi) * dx          # 对 N(0,1) 求积
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
    kern = np.exp(-0.5 * (kk * hz / sd) ** 2) * hz / (sd * np.sqrt(2 * pi))
    kern /= kern.sum() * hz * (1.0 / (1.0))                 # 使 ∫k=1（离散）
    p = np.convolve(dens, kern, mode='same')
    mass = p.sum() * hz
    p = np.maximum(p / mass, 1e-300)
    m1 = np.convolve(dg, kern, mode='same') / np.maximum(np.convolve(dens, kern, mode='same'), 1e-300)
    Eg2 = float(np.sum(wx * gx ** 2))
    D = Eg2 - float(np.sum(p * m1 ** 2) * hz)
    hZ = -float(np.sum(p * np.log(p)) * hz)
    R = (hZ - 0.5 * np.log(2 * pi * e * V)) / LN2
    return D, R, dict(nx=nx, ng=ng, hz=hz, mass=mass, tail=float(dens[-1] / dens.max()))


def hgauss(eps, nx=400001, xmax=8.0):
    """h(g(X)) 的闭式：$\\tfrac12\\ln(2\\pi e)+E\\ln g'(X)$，$g'=1+3\\varepsilon x^2$。"""
    x = np.linspace(-xmax, xmax, nx)
    dx = x[1] - x[0]
    wx = np.exp(-0.5 * x ** 2) / np.sqrt(2 * pi) * dx
    gp = 1 + 3 * eps * x ** 2
    EL = float(np.sum(wx * np.log(gp)))
    AM = float(np.sum(wx * gp ** 2))
    return 0.5 * np.log(2 * pi * e) + EL, AM, EL


def slb(eps, D):
    hG, AM, EL = hgauss(eps)
    return (hG - 0.5 * np.log(2 * pi * e * D)) / LN2


print(r'== E58：开放问题 1 的独立验证（推前测度 + 离散卷积 KDE）==')
eps = 0.3
hG, AM, EL = hgauss(eps)
gap_theory = 0.5 * np.log2(AM) - EL / LN2                  # ½log2 AM − E ln g' /ln2
print(r'[0] 常数与定理右边（$\varepsilon=0.3$）')
print(r'  $E[Xg]=1+3\varepsilon=%.4f$   $E[g^2]=1+6\varepsilon+15\varepsilon^2=%.4f$   $E[g^{\prime2}]=1+6\varepsilon+27\varepsilon^2=%.4f$（解析）'
      % (1 + 3 * eps, 1 + 6 * eps + 15 * eps ** 2, 1 + 6 * eps + 27 * eps ** 2))
xg = np.linspace(-8, 8, 400001); dxg = xg[1] - xg[0]
wg = np.exp(-0.5 * xg ** 2) / np.sqrt(2 * pi) * dxg
print(r'  同一批常数用 $x$ 网格重算：$E[Xg]=%.6f$  $E[g^2]=%.6f$  $E[g^{\prime2}]=%.6f$' %
      (np.sum(wg * (xg ** 2 + eps * xg ** 4)), np.sum(wg * (xg + eps * xg ** 3) ** 2),
       np.sum(wg * (1 + 3 * eps * xg ** 2) ** 2)))
print(r'  定理右边 $\tfrac12\log_2\mathrm{AM}/\mathrm{GM}=%.5f$ bit（C 印 0.47489）；$h(g(X))=%.6f$ nat' % (gap_theory, hG))

print(r'\n[1] 控制检验：仿射机制 $q=0$ 的 $R$ 对解析 $\tfrac12\log_2(1+1/V)$')
print('  %-9s %-12s %-12s %-11s %-12s' % ('V', 'R 网格', 'R 解析', '差', 'D 网格'))
for V in (0.03, 0.01, 3e-3, 1e-3, 3e-4, 1e-4):
    D, R, diag = mech(eps, 0.0, V)
    Ra = 0.5 * np.log2(1 + 1 / V)
    print('  %-9.4g %-12.6f %-12.6f %-11.2e %-12.6f  (ng=%d, hz=%.2e, ∫p=%.6f)' %
          (V, R, Ra, R - Ra, D, diag['ng'], diag['hz'], diag['mass']))

print(r'\n[2] 定理本体：同一机制类自己的率–失真曲线，对 SLB 的差（$\varepsilon=0.3$）')
print('  %-9s | %-11s %-11s %-11s | %-11s %-11s %-11s' %
      ('V', 'D_aff', 'R_aff', 'R_aff−SLB', 'D_task', 'R_task', 'R_task−SLB'))
rows = []
for V in (0.02, 5e-3, 1e-3, 2e-4, 5e-5, 1e-5):
    Da, Ra, _ = mech(eps, 0.0, V)
    Dt, Rt, _ = mech(eps, eps, V)
    ga = Ra - slb(eps, Da); gt = Rt - slb(eps, Dt)
    rows.append((V, Da, Ra, ga, Dt, Rt, gt))
    print('  %-9.3g | %-11.6f %-11.6f %-11.6f | %-11.6f %-11.6f %-11.6f' % (V, Da, Ra, ga, Dt, Rt, gt))
print('  预言：左块 $\to$ %.5f（=½log₂AM/GM），右块 $\to$ 0（任务机制达到 SLB 主项）' % gap_theory)
print('  同 D 对比（把任务机制插值到仿射的 $D_{\rm aff}$ 上再相减，才是他定理的口径）：')
for V in (5e-3, 1e-3, 2e-4, 5e-5):
    Da, Ra, _ = mech(eps, 0.0, V)
    Dt, Rt, _ = mech(eps, eps, V)
    # 在 V_task 上找同 D
    lo, hi = V * 1e-3, V * 3
    for _ in range(40):
        mid = np.sqrt(lo * hi)
        Dm, Rm, _ = mech(eps, eps, mid)
        if Dm > Da:
            hi = mid
        else:
            lo = mid
    Dm, Rm, _ = mech(eps, eps, mid)
    print('    $D=%.6f$：$R_{\rm aff}=%.5f$  同 $D$ 的 $R_{\rm task}=%.5f$  差 $=%.5f$（预言 %.5f）'
          % (Da, Ra, Rm, Ra - Rm, gap_theory))

print(r'\n[3] 他 §9.5  prescribed 检验：修掉 lin_constants 之后，我脚本那列 gap 在"小 $V$ 端"是什么')
A1c, D0c = 1 + 3 * eps, 1 + 6 * eps + 15 * eps ** 2
floor = D0c - A1c ** 2
print('  正确的线性（LMMSE 口径）前沿：$D(s)=E[g^2]-(E[Xg])^2 s$，$\;s\to1$ 给地板 $D_\infty=6\varepsilon^2=%.4f$' % floor)
print('  仿射机制的 $D_{\rm aff}\approx V\,E[g^{\prime2}]=V(%.3f)$ ⇒ $V<%.4f$ 时 $D_{\rm aff}$ 掉到该地板**以下**'
      % (AM, floor / AM))
for V in (0.03, 0.01, 3e-3, 1e-3):
    Da, Ra, _ = mech(eps, 0.0, V)
    s = (D0c - Da) / A1c ** 2
    Rl = -0.5 * np.log(1 - s) / LN2 if s < 1 else float('nan')
    print('    $V=%-6.4g$ $D_{\rm aff}=%.6f$  $s=%.4f$  $R_{\rm lin}@D$=%s   （$s\\ge1\Rightarrow$ 该列无定义）'
          % (V, Da, s, ('%.5f' % Rl) if np.isfinite(Rl) else 'NaN'))
print('  ⇒ 按字面执行 §9.5（"gap 列应趋于 0.47489"）在小 $V$ 端得到的是 NaN，不是那个数。')
print('    0.47489 是 $R_{\rm aff}-{\rm SLB}$ 的极限（见 [2]），它对应的对象不是 `lin_R_at_D`。')

print(r'\n[4] 用旧（错位）常数列同一张表：它在小 $V$ 端给出什么（说明他为什么能预言一个有限数）')
A1o, D0o = 1 + 6 * eps + 15 * eps ** 2, 1 + 3 * eps + 15 * eps ** 2
flr_o = D0o - A1o ** 2
print('  旧式地板 $D_\infty=D_0-A_1^2=%.4f<0$ ⇒ 对**任何**正 $D$ 都有 $s<1$，那列永远有限——正因如此它不是任何机制的速率。' % flr_o)
