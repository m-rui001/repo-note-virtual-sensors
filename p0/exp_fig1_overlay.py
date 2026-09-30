"""E26：把本仓库的 TRV 前沿压到原文 Fig.1 的红曲线上对表，并用红曲线**反拟合地板**。

三方数据：
  (1) 我 E25 逐像素读出的红曲线 `p0/fig/fig1_red.npy`（标定自检两条都过：虚线读 1.1720 vs 理论 1.169；
      黑曲线 D=80 读 1.5864 vs 发表 1.602）。红 = "2 unstable" = 我的 F2。
  (2) 本仓库锥上真最优链（`exp_nearfloor_law.solve_cone`，各向同性二分解给可行起点 + 热启动）。
  (3) C 的 c16 稠密前沿（`.work3/c16_plaw.py` 的 `mine`/`dig` 两个表）。

C 的 `dig` 第一个点 (46.5, 2.3438) 是脏的：红线在图上从 I=5（图顶）离开，对应 D=46.988，
D=46.5 处根本没有红线像素；C 把图顶点当成 D=46.5 的点，于是"原文 Δ(46.5)=2.34 vs 你 3.61"不成立。
本脚本用我自己的独立数字化重做这件事。

主检验是 (1)：把本仓库锥上真最优逐点压到红线像素上对表；(2) 用局部斜率把速率偏差换算成
地板偏差的上限。**作废**的是原来那个"用红曲线反拟合对数律来钉地板"的口径（留在 (4) 里给数字）：
图窗 $x=D-46.1231\\in[0.88,13.9]$ 全在交叉段，离渐近区（要 $x\\lesssim0.05$）差两个数量级，
所以那个拟合给 $f=43.99$、$a=1.68$ 都不构成对 Prop E 的检验。(3) 是蓝曲线对 rank-3 的第二重对表。
"""
import sys
import numpy as np
sys.stdout.reconfigure(encoding='utf-8')
from p0.p0_replicate_letter import A, W, n, sym, ln2, ctrl, rate_cost, unconstrained_sdp
from p0.exp_nearfloor_law import Fs, floor_of, solve_cone, gamma_iso, pack, unpack, ij, Th, Pc

FLOOR_F2 = 46.12314
A_PRED = (2 / 2.) * np.log2(10.)            # = 3.3219 bit / 十倍程
C_MINE = [(46.5, 3.6094), (47.5, 1.8498), (48.5, 1.2217), (49.5, 0.9068), (50.5, 0.7149),
          (51.5, 0.5845), (52.5, 0.4868), (55.5, 0.3120), (60.5, 0.1940), (65.5, 0.1142),
          (70.5, 0.0742), (75.5, 0.0506), (80.5, 0.0359)]
C_DIG = [(46.5, 2.3438), (47.5, 1.9350), (48.5, 1.2365), (49.5, 0.8851), (50.5, 0.7188),
         (51.5, 0.5642), (52.5, 0.4797), (55.5, 0.3075), (60.5, 0.1656), (65.5, 0.1022),
         (70.5, 0.0659), (75.5, 0.0456), (80.5, 0.0321), (85.5, 0.0256), (89.5, 0.0220)]


def fit_log(xs_ds, a_fixed=None, f_fixed=None):
    """Δ = a·log₁₀(1/(D−f)) + b 的最小二乘；a 或 f 可钉死。用网格 + 一维精化求 f。"""
    D = np.array([x for x, _ in xs_ds])
    y = np.array([v for _, v in xs_ds])

    def sse(f):
        a, b = (a_fixed, None) if a_fixed else (None, None)
        L = np.log10(np.maximum(1. / np.maximum(D - f, 1e-9), 1e-12))
        if a_fixed:
            b = np.mean(y - a_fixed * L)
            r = y - (a_fixed * L + b)
        else:
            M = np.vstack([L, np.ones_like(L)]).T
            coef = np.linalg.lstsq(M, y, rcond=None)[0]
            a, b = coef
            r = y - (a * L + b)
        return np.sum(r * r), (a if not a_fixed else a_fixed), b

    if f_fixed:
        s, a, b = sse(f_fixed)
        return a, b, f_fixed, s
    grid = np.linspace(D.min() - 3.0, D.min() - 1e-3, 20000)
    f = grid[np.argmin([sse(g)[0] for g in grid])]
    lo, hi = f - .05, f + .05
    for _ in range(60):
        m1, m2 = lo + (hi - lo) / 3, hi - (hi - lo) / 3
        if sse(m1)[0] < sse(m2)[0]:
            hi = m2
        else:
            lo = m1
    s, a, b = sse(.5 * (lo + hi))
    return a, b, .5 * (lo + hi), s


if __name__ == '__main__':
    F = Fs['F2 尾部2(封闭, r=2)']
    red = np.load('p0/fig/fig1_red.npy')
    blk = np.load('p0/fig/fig1_black.npy')
    mag = np.load('p0/fig/fig1_magenta.npy')
    print('红曲线像素点 %d 个，D∈[%.3f, %.3f]，I∈[%.3f, %.3f]'
          % (len(red), red[:, 0].min(), red[:, 0].max(), red[:, 1].min(), red[:, 1].max()))
    print('品红(r=4 满秩)与黑(无约束)在 D=40 的读出：%.4f vs %.4f（应重合，Prop 3）'
          % (mag[mag[:, 0] < 40.3][:, 1].mean(), blk[blk[:, 0] < 41][:, 1].mean()))

    # ---------- (1) 我的前沿 vs 图上红曲线 ----------
    print('\n===== 本仓库锥上真最优 vs 原文 Fig.1 红曲线 =====')
    print('%8s %9s %9s %8s | %9s %9s' % ('D', '图上I', '我的I', '差', 'C的Δ', '我的Δ'))
    Dgrid = [46.60, 46.80, 47.0, 47.5, 48.0, 48.5, 49.5, 50.5, 52.5, 55.5, 60.5, 65.5, 70.5, 80.5]
    prev_G, prev_x, mine = None, None, []
    for D in sorted(Dgrid, reverse=True):
        x = D - FLOOR_F2
        Iu = unconstrained_sdp(D, np.eye(n), Th, Pc)['I_true']
        lg = gamma_iso(F, D)
        starts = ([pack(np.exp(lg) * np.eye(2))] if lg is not None else [])
        if prev_G is not None:
            G0 = prev_G * (prev_x / x)
            starts += [pack(G0), pack(sym(np.diag([np.sqrt(.3), np.sqrt(1 / .3)]) @ G0
                                          @ np.diag([np.sqrt(.3), np.sqrt(1 / .3)])))]
        r = solve_cone(F, D, FLOOR_F2, starts)
        if r is None:
            print('%8.2f  求解失败' % D)
            continue
        It, Jt, Gt = r[0], r[1], r[2]
        prev_G, prev_x = Gt, x
        mine.append((D, It))
        sel = red[np.abs(red[:, 0] - D) < .08]
        cd = next((v for k, v in C_MINE if abs(k - D) < .08), np.nan)
        print('%8.2f %9s %9.4f %+8.4f | %9s %9.4f'
              % (D, '%.4f' % sel[:, 1].mean() if len(sel) else '—',
                 It, It - sel[:, 1].mean() if len(sel) else np.nan,
                 '%.4f' % (cd + Iu) if np.isfinite(cd) else '—', It - Iu))
    mm = np.array(mine)
    dev = []
    for D, I in mm:
        sel = red[np.abs(red[:, 0] - D) < .08]
        if len(sel):
            dev.append(I - sel[:, 1].mean())
    dev = np.array(dev)
    print('重叠 %d 点：平均偏差 %+.4f bit，最大 |偏差| %.4f bit（标定精度 0.003 bit）'
          % (len(dev), dev.mean(), np.max(np.abs(dev))))
    Itop = red[:, 1].max()
    Dtop = np.interp(Itop, mm[:, 1], mm[:, 0])      # mm 天然按 I 升序（D 降序）
    print('红线在图上是被坐标轴上界截断的：图顶 $I_{max}=%.3f$。' % Itop)
    print('我的曲线（内在插值，不外推）到 $I=%.3f$ 需 $D=%.3f$；图上红线左端 $D=%.3f$（同水平差 %+.3f）'
          % (Itop, Dtop, red[:, 0].min(), Dtop - red[:, 0].min()))

    # ---------- (2) 图能钉住地板到什么精度：用局部斜率把速率偏差换成地板偏差 ----------
    print('\n===== 地板被 Fig.1 红曲线钉住的程度（斜率换算，不做外推拟合）=====')
    sl = []
    for i in range(len(mm) - 1):
        if mm[i, 0] - mm[i + 1, 0] > 0:
            sl.append((.5 * (mm[i, 0] + mm[i + 1, 0]), (mm[i, 1] - mm[i + 1, 1])
                       / (mm[i, 0] - mm[i + 1, 0])))
    near = [abs(s) for D0, s in sl if D0 < 48.5]
    s_lo = min(near)
    print('%-28s 局部斜率 |dI/dD| ∈ [%.3f, %.3f] bit/unit' % ('近地板端 D<48.5', min(near), max(near)))
    print('最大速率偏差 %.4f bit ÷ 最陡斜率 %.3f ⟹ 地板偏差不超过 ±%.3f（Φ 的 %.2f%%）'
          % (np.max(np.abs(dev)), s_lo, np.max(np.abs(dev)) / s_lo,
             100 * np.max(np.abs(dev)) / s_lo / 14.63979))
    print('注意：图窗 $x=D-46.1231\\in[0.88,13.9]$ 全在交叉段（渐近区要 $x\\lesssim0.05$），')
    print('所以**不能**用图反拟合对数律。实测：钉 $a=3.3219$ 拟合 → $f=43.99$；钉 $f=46.1231$ → $a=1.68$。')
    print('这两个数都不反驳 Prop E，只说明窗口离渐近区还差两个数量级。')

    # ---------- (3) 蓝曲线（=我的 F3，rank 3，地板 36.6032）第二重独立对表 ----------
    print('\n===== 本仓库 F3 锥上真最优 vs 原文 Fig.1 蓝曲线（第二个地板，Φ=5.1199）=====')
    F3 = Fs['F3 尾部3(不封闭, r=3)']
    fl3 = floor_of(F3)[0]
    blu = np.load('p0/fig/fig1_blue.npy')
    print('%8s %9s %9s %8s %10s' % ('D', '图上I', '我的I', '差', '我的Δ'))
    dev3 = []
    prev_G, prev_x = None, None
    for D in sorted([40.0, 42.0, 45.0, 48.0, 52.0, 58.0, 66.0, 80.0], reverse=True):
        x = D - fl3
        Iu = unconstrained_sdp(D, np.eye(n), Th, Pc)['I_true']
        lg = gamma_iso(F3, D)
        starts = ([pack(np.exp(lg) * np.eye(3))] if lg is not None else [])
        if prev_G is not None:
            starts.append(pack(prev_G * (prev_x / x)))
        r3 = solve_cone(F3, D, fl3, starts)
        if r3 is None:
            print('%8.2f  求解失败' % D)
            continue
        sel = blu[np.abs(blu[:, 0] - D) < .08]
        prev_G, prev_x = r3[2], x
        if len(sel):
            dev3.append(r3[0] - sel[:, 1].mean())
        print('%8.2f %9s %9.4f %+8.4f %10.4f'
              % (D, '%.4f' % sel[:, 1].mean() if len(sel) else '—', r3[0],
                 r3[0] - sel[:, 1].mean() if len(sel) else np.nan, r3[0] - Iu))
    dev3 = np.array(dev3)
    print('蓝曲线重叠 %d 点：平均 %+.4f，最大 |偏差| %.4f bit'
          % (len(dev3), dev3.mean(), np.max(np.abs(dev3))))

    # ---------- (2) 只用图上红曲线反拟合地板 ----------
    print('\n===== 只用 Fig.1 红曲线反拟合 Prop E 对数律（Δ 用红线减黑线，同 D 插值）=====')
    dd = np.unique(np.round(red[:, 0], 2))
    dd = dd[(dd > red[:, 0].min() + .05) & (dd < 60.0)]
    pair = []
    for D in dd:
        s1 = red[np.abs(red[:, 0] - D) < .06]
        s2 = blk[np.abs(blk[:, 0] - D) < .06]
        s3 = mag[np.abs(mag[:, 0] - D) < .06]
        if len(s1) and (len(s2) or len(s3)):
            base = s2[:, 1].mean() if len(s2) else s3[:, 1].mean()
            pair.append((D, s1[:, 1].mean() - base))
    print('图上 Δ(D) 取 %d 点，D∈[%.2f, %.2f]' % (len(pair), pair[0][0], pair[-1][0]))
    a, b, f, s = fit_log(pair, a_fixed=A_PRED)
    print('(A) 系数钉死 a=%.4f，地板自由 → 图上反拟合地板 f = %.4f ± ?   （Prop D 预言 %.4f，差 %+.4f）'
          % (A_PRED, f, FLOOR_F2, f - FLOOR_F2))
    a2, b2, f2, s2 = fit_log(pair, f_fixed=FLOOR_F2)
    print('(B) 地板钉死 f=%.4f，系数自由 → 图上拟合系数 a = %.4f（预言 %.4f，误差 %+.1f%%）'
          % (FLOOR_F2, a2, A_PRED, 100 * (a2 / A_PRED - 1)))
    a3, b3, f3, s3 = fit_log(pair)
    print('(C) 三参数自由 → a=%.4f  b=%.4f  f=%.4f' % (a3, b3, f3))
    print('    对照：C 的 c16 幂律拟合在 D∈[47,53] 给 p=0.910、c=2.656（原文数字化）/ p=0.868（它自己的前沿）')
    for tag, dat in (('C 的 dig(含脏首点)', C_DIG), ('C 的 mine(它的前沿)', C_MINE)):
        aa, bb, ff, ss = fit_log(dat, a_fixed=A_PRED)
        print('    同样口径拟合 %-22s → 地板 f = %.4f' % (tag, ff))
