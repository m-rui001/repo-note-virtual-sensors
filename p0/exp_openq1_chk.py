r"""E57 = C 的 §9（开放问题 1：静态最小舞台，AM/GM Jensen 缺口）独立复核 + 两处反驳。

C 的 §9.1 打中了我 `p0/exp_openq1_1d.py::lin_constants` 的常数错位——[A] 先把这一条钉死（我用精确矩，
不用他的 GH）。他的主结果（§9.2 方框）是
    $R_{\rm aff}(D)-R_{\rm task}(D)\to\tfrac12\log_2\mathrm{AM}(g'^2)/\mathrm{GM}(g'^2)$，
等号 $\iff g$ 仿射。[B]/[C] 用**我自己的路数**验它：$x$ 上均匀网格（宽度 $\propto\sqrt V$）+ $z$ 上
细梯形求 $h(Z)$ 与 $E[g|z]$，既不是他的 $z\times x$ 自适应网格，也不是我文件里的 GH（他自己警告过
GH 在小噪声必错）。[B] 的控制格是 $R_{\rm aff}$ 对解析 $\tfrac12\log_2(1+1/\sigma^2)$。

两处反驳（都是可判的）：
 [D] 他 §9.3 那张表里"匹配后的 gap 极限"一列（0.35/0.38 来回跳、不收斂到 0.47489）是**内插伪影**：
     对数量级跨 3 个的 $D$ 做**线性**内插（应当对 $\log D$）。我用他印出的两行原始数据复算，
     线性内插应能逐位重现他的 `4.63246 / 5.47604 / 6.29223`，$\log$ 内插则落回他自己的 SLB。
 [E] 他 §9.5 给我的"可判定话"把两行**说反了**：他的仿射行是"线性编码 + 最优解码"，
     而我 B 族里 $q=0$ 恰是同一个机制 → gap 应趋于 $0$；$q=\varepsilon$ 才是任务机制 → 趋于 $0.47489$。
"""
import sys
sys.stdout.reconfigure(encoding='utf-8')
import numpy as np

LN2 = np.log(2.0)
print(r'== E57：开放问题 1 的独立复核（$g(x)=x+\varepsilon x^3$，$X\sim N(0,1)$）==')

# ---------------------------------------------------------------- [A] 常数
xg = np.linspace(-12, 12, 400001)
pw = np.exp(-0.5 * xg ** 2) / np.sqrt(2 * np.pi)
I = lambda f: float(np.trapezoid(pw * f, xg))
print('\n  [A] 精确矩（40 万点梯形）vs C §9.1 的闭式 vs 我文件里印的' + r'')
for e in (0.3, 0.05):
    g = lambda t, e=e: t + e * t ** 3
    gp = lambda t, e=e: 1 + 3 * e * t ** 2
    EXg, Eg2, Egp2 = I(xg * g(xg)), I(g(xg) ** 2), I(gp(xg) ** 2)
    print('    eps=%.2f  E[Xg]=%.6f (闭式 $1+3e$=%.6f)   $E[g^2]$=%.6f (闭式 $1+6e+15e^2$=%.6f)   '
          '$E[g\'^2]$=%.6f (闭式 $1+6e+27e^2$=%.6f)'
          % (e, EXg, 1 + 3 * e, Eg2, 1 + 6 * e + 15 * e ** 2, Egp2, 1 + 6 * e + 27 * e ** 2))
    A1, D0 = 1 + 6 * e + 15 * e ** 2, 1 + 3 * e + 15 * e ** 2      # 我文件里的原样
    print('      我文件的 `A1`=%.6f 是 $E[g^2]$；`D0`=%.6f 无对应矩。旧式 $D=D_0-A_1^2 s$ 在 $s=0.5$：'
          '%.4f  <-- 负 MSE，C 说的 $-5.3612$ 复现%s'
          % (A1, D0, D0 - A1 ** 2 * .5, '' if e != 0.3 else '（正确）'))
    print('      正确式 $D=E[g^2]-(E[Xg])^2 s$：$s=0.5\\to$%.4f，$s=1\\to D(\\infty)=%.4f$（C 印 $6e^2=%.4f$）'
          % (Eg2 - EXg ** 2 * .5, Eg2 - EXg ** 2, 6 * e ** 2))

# ---------------------------------------------------------------- 机制网格
def mech(eps, q, V, xr=5.0, min_dx_frac=4.0):
    r"""$z=\phi(x)+\sqrt V n$，$\phi(x)=x+qx^3$。返回 $(D,\ R_{\rm bits})$，$D=E[g^2]-E(E[g|z]^2)$。
    $x$ 步长 $\le\sqrt V/$min_dx_frac（避开高斯和的绕射），$z$ 步长 $\sqrt V/12$。"""
    sd = np.sqrt(V)
    nx = int(np.clip(2 * xr / (sd / min_dx_frac), 600, 6000))
    x = np.linspace(-xr, xr, nx)
    wx = np.exp(-0.5 * x ** 2)
    wx /= wx.sum()
    phi = x + q * x ** 3
    gx = x + eps * x ** 3
    dz = sd / 12.0
    zlo, zhi = phi.min() - 9 * sd, phi.max() + 9 * sd
    m = int(np.ceil((zhi - zlo) / dz)) + 1
    Eg2 = float(np.sum(wx * gx ** 2))
    acc_m2 = 0.0
    acc_h = 0.0
    CH = 20000000 // nx
    for a in range(0, m, CH):
        z = zlo + dz * (a + np.arange(min(CH, m - a)))
        K = np.exp(-0.5 * ((z[:, None] - phi[None, :]) / sd) ** 2)      # (chunk, nx)
        p = K @ wx
        p = np.maximum(p, 1e-300)
        m1 = (K @ (wx * gx)) / p
        acc_h += float(np.sum(-p * np.log(p)))
        acc_m2 += float(np.sum(p * m1 ** 2))
    hZ = acc_h * dz / LN2                       # nats -> bits
    D = Eg2 - acc_m2 * dz
    R = hZ - 0.5 * np.log2(2 * np.pi * np.e * V)
    return D, R, nx, m


def gap_pred(eps, n=200001):
    t = np.linspace(-12, 12, n)
    w = np.exp(-0.5 * t ** 2) / np.sqrt(2 * np.pi)
    gp = 1 + 3 * eps * t ** 2
    norm = np.trapezoid(w, t)
    AM = np.trapezoid(w * gp ** 2, t) / norm
    EL = np.trapezoid(w * np.log(np.maximum(gp, 1e-300)), t) / norm
    return 0.5 * (np.log2(AM) - EL / LN2), AM, np.exp(EL)


print('\n  [B] 控制格：我的网格路数 vs 解析（仿射机制 $q=0$，$R=\tfrac12\log_2(1+1/\\sigma^2)$）')
for V in (3e-2, 1e-2, 3e-3):
    D, R, nx, m = mech(0.3, 0.0, V)
    Ra = 0.5 * np.log2(1 + 1 / V)
    print('    $\\sigma^2$=%-7.3g $n_x$=%-5d $n_z$=%-7d $R$=%.5f 解析=%.5f 差=%+.0e   $D$=%.6f  '
          '$D/\\sigma^2$=%.4f（$E[g\'^2]$@0.3=5.23）' % (V, nx, m, R, Ra, R - Ra, D, D / V))

print('\n  [C] 匹配 $D$ 的缺口 $R_{\rm aff}-R_{\rm task}$（$\\varepsilon=0.3$），对 $\\log D$ 内插')
SG = np.geomspace(3e-2, 1e-3, 6)
aff = np.array([mech(0.3, 0.0, v)[:2] for v in SG])
tsk = np.array([mech(0.3, 0.3, v)[:2] for v in SG])
print('    %-10s %-12s %-10s | %-12s %-10s' % ('$\\sigma^2$', '$D$ 仿射', '$R$ 仿射', '$D$ 任务($s_2$同格)', '$R$ 任务'))
for k in range(len(SG)):
    print('    %-10.3g %-12.6f %-10.5f | %-12.6f %-10.5f' % (SG[k], aff[k, 0], aff[k, 1], tsk[k, 0], tsk[k, 1]))
gpred, AM, GM = gap_pred(0.3)
# 在共同 log D 网格上匹配
linD = np.geomspace(max(aff[:, 0].min(), tsk[:, 0].min()), min(aff[:, 0].max(), tsk[:, 0].max()), 6)
Ra_i = np.interp(np.log(linD), np.log(aff[:, 0])[::-1], aff[:, 1][::-1])
Rt_i = np.interp(np.log(linD), np.log(tsk[:, 0])[::-1], tsk[:, 1][::-1])
Rt_lin = np.interp(linD, tsk[:, 0][::-1], tsk[:, 1][::-1])
print('    $\\tfrac12\\log_2(\\mathrm{AM}/\\mathrm{GM})$ 预言 = %.5f（AM=%.5f GM=%.5f）' % (gpred, AM, GM))
print('    %-12s %-12s %-12s %-12s' % ('$D$ 共同格', '$R$仿射', 'gap($\\log D$ 内插)', 'gap($D$ 线性内插)'))
for k in range(len(linD)):
    print('    %-12.6f %-12.5f %-12.5f %-12.5f' % (linD[k], Ra_i[k], Ra_i[k] - Rt_i[k], Ra_i[k] - Rt_lin[k]))

print('\n  [D] 复现 C §9.3 的"匹配后的 gap"伪影（只用他印出的行）')
# 他的两前沿表：(sigma2, D_aff, R_aff) 与 (s2, D_task, R_task)
DA = [(3e-2, 0.150969, 2.55077), (1e-2, 0.051624, 3.32911), (3e-3, 0.015629, 4.19257),
      (1e-3, 0.005223, 4.98361), (3e-4, 0.001568, 5.85159), (1e-4, 0.000523, 6.64393)]
DT = [(3e-2, 0.029392, 3.26303), (1e-2, 0.009930, 4.04557), (3e-3, 0.002994, 4.91048),
      (1e-3, 0.000999, 5.70188), (3e-4, 0.000300, 6.56862), (1e-4, 0.000099, 7.35607)]
hg = 2.76561     # 他印的 h(g(X))
print('    %-12s %-12s %-14s %-14s %-12s' % ('$D_{\rm aff}$', '他印 $R_{\rm task}({同}D)$', '线性内插($D$)',
                                            '$\\log$ 内插($D$)', 'SLB$(D)$'))
for s2, Da, Ra in DA:
    dts = np.array([d for _, d, _ in DT])[::-1]
    rts = np.array([r for _, _, r in DT])[::-1]
    lnl = np.interp(Da, dts, rts)
    llg = np.exp(np.interp(np.log(Da), np.log(dts), np.log(rts)))
    slb = hg - 0.5 * np.log2(2 * np.pi * np.e * Da)
    print('    %-12.6f %-12.5f %-14.5f %-14.5f %-12.5f' % (Da, {'0.150969': 3.26303, '0.051624': 3.26303,
          '0.015629': 3.81643, '0.005223': 4.63246, '0.001568': 5.47604, '0.000523': 6.29223}[('%.6f' % Da)],
          lnl, llg, slb))
print('    ⇒ 若"他印"列与"线性内插"列逐位一致，则那一列不是 gap 极限，是内插口径错；'
      '$\\log$ 内插列应贴着 SLB 列（他自己的 $R_{\rm task}$-$\rm SLB$ 差 $\\le4\\times10^{-5}$ 已证）。')

print('\n  [E] C §9.5 的两行归属：$q=0$ 与 $q=\\varepsilon$ 各自的 gap 极限（我的路数，$\\varepsilon=0.3$）')
print('    %-8s %-10s %-12s %-12s %-12s' % ('$q$', '$\\sigma^2$', '$D$', '$R$', 'gap vs 仿射前沿'))
for q in (0.0, 0.3):
    for V in (1e-2, 3e-3):
        D, R, _, _ = mech(0.3, q, V)
        Ra = np.interp(np.log(D), np.log(aff[:, 0])[::-1], aff[:, 1][::-1])
        print('    %-8.3f %-10.3g %-12.6f %-12.5f %-12.5f' % (q, V, D, R, Ra - R))
print('    预言（我的读法）：$q=0$ 与仿射机制同一 $\\Rightarrow$ gap $\\to0$；$q=\\varepsilon$ 是任务机制 $\\Rightarrow$ gap $\\to$%.5f' % gpred)

print('\n  [F] 缺口的 $\\varepsilon$ 级数（C: $26\\varepsilon^2(1-10\\varepsilon+109\\varepsilon^2)$）vs 闭式 AM/GM')
for e in (0.002, 0.01, 0.05, 0.1, 0.3):
    gpv = gap_pred(e)[0]
    ser = (18 * e ** 2 - 180 * e ** 3 + 1964.25 * e ** 4) / LN2   # C 的三级级数（bit）
    print('    eps=%-6.3g 闭式 $\\tfrac12\\log_2(\\mathrm{AM/GM})$=%.6f   C 的三级级数=%.6f   比 %.4f'
          % (e, gpv, ser, ser / gpv if gpv else np.nan))
