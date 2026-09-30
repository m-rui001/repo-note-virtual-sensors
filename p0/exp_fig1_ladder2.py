r"""E39b：图 1 的红曲线在拐点左侧**被图例框吃掉**——把这件事量化，再验阶梯的族切换。

E39 的两个发现需要这一步：
 (a) 黑/品红像素与我的解析 $I_{unc}$ 在 $D\in[41,89]$ 上对到 $\pm0.01\sim0.02$ bit，所以图的轴标定
     不是误差源；红曲线与我的可行点的偏差是**求解器**的。
 (b) 拐点预言在 $D=47.2075$，只有 7.9 列像素，而最左 4 列的红线顶端钉在 $I\simeq4.926$
     （图例框下沿）——E39 用逐列最长段中点做回归，把这 4 列当成了曲线上的点，于是"拐点前斜率 0.43"
     是遮挡伪影，不是数据。这里改成：只取未被图例框切断的列，并给出遮挡边界。

然后做真正有内容的一件事：**族切换的地点在图窗内**。E39 的竞赛表显示 $x\gtrsim1.7$ 时最优族从
$\{0,1\}$ 换成单通道 $\{0\}$，而 $\mathrm{gap}_{\{1\}}=74.28$，也就是说原文 F2 那条曲线画的整个预算区间
里，弱通道（$C$ 的小本征方向）根本不该传输。这可以在数值最优的 $\Gamma^\ast$ 谱上直接看。
"""
import sys
import numpy as np
sys.stdout.reconfigure(encoding='utf-8')
from PIL import Image
import p0.p0_replicate_letter as base
import p0.exp_plants_oos as E35
import p0.exp_ladder as X37
import p0.exp_C_exact as X36
from p0.exp_fig1_digitize import frame, legend_box, to_data, IMG
from p0.exp_nearfloor_law import Fs, solve_cone, gamma_iso, pack, unpack
from p0.p0_replicate_letter import n, unconstrained_sdp

ln2 = np.log(2.0)
F2 = Fs['F2 尾部2(封闭, r=2)']
Pl = E35.Pl(base.A, base.B, base.W, np.eye(4))
Th, Pc = Pl.Th, Pl.Pc


def unc(D):
    return unconstrained_sdp(D, np.eye(n), Th, Pc)['I_true']


def red_columns():
    """逐列返回红线的连续段：(D, 段中点 I, 段长, 顶端 I, 是否贴图例框)。"""
    a = np.asarray(Image.open(IMG).convert('RGB'))
    box = frame(a)
    leg = legend_box(a, box)
    r0, r1, c0, c1 = box
    R, G, B = a[..., 0].astype(int), a[..., 1].astype(int), a[..., 2].astype(int)
    m = (R > 150) & (G < 110) & (B < 110)
    mg = 3
    m[:r0 + mg, :] = m[r1 - mg:, :] = False
    m[:, :c0 + mg] = m[:, c1 - mg:] = False
    lr0, lr1, lc0, lc1 = leg
    me0, me1 = lr0 - 4, lr1 + 5          # 被掩掉的行区间 [me0, me1)
    m[me0:me1, lc0 - 4:lc1 + 5] = False
    rows = []
    for c in range(c0 + mg, c1 - mg):
        rr = np.where(m[r0 + mg:r1 - mg, c])[0]
        if not len(rr):
            continue
        inleg = (lc0 - 4) <= c <= (lc1 + 5)
        brk = np.where(np.diff(rr) > 3)[0]
        for sg in np.split(rr, brk + 1):
            top, bot = sg.min() + r0 + mg, sg.max() + r0 + mg
            Itop, Ibot = to_data(0, top, box)[1], to_data(0, bot, box)[1]
            cut = (top == r0 + mg) or (bot >= r1 - mg - 1) or \
                  (inleg and (top == me1 or bot == me0 - 1))
            rows.append((to_data(c, 0, box)[0], .5 * (Itop + Ibot), len(sg), Itop, Ibot, cut))
    return np.array([r[0] for r in rows]), np.array([r[1] for r in rows]), \
        np.array([r[2] for r in rows]), np.array([r[3] for r in rows]), \
        np.array([r[5] for r in rows]), leg, box


def slope_fit(D, I, mask, FLOOR, tag):
    s = mask
    if s.sum() < 4:
        print('    %-30s 可用列 %d，太少' % (tag, s.sum())); return
    uu = np.log2(1.0 / (D[s] - FLOOR))
    dd = I[s]
    M = np.vstack([uu, np.ones_like(uu)]).T
    coef, *_ = np.linalg.lstsq(M, dd, rcond=None)
    resid = dd - M @ coef
    sd = np.sqrt(np.sum(resid ** 2) / (s.sum() - 2))
    cov = sd ** 2 * np.linalg.inv(M.T @ M)
    print('    %-30s 列 %3d  斜率 %+.4f $\pm$ %.4f   截距 %+.3f  残差 $\sigma$=%.4f bit  $u$ 跨度 %.3f'
          % (tag, s.sum(), coef[0], np.sqrt(cov[0, 0]), coef[1], sd, uu.max() - uu.min()))
    return coef, np.sqrt(cov[0, 0])


if __name__ == '__main__':
    fx = X37.floor_exact(Pl, F2, None)
    phi0 = fx['闭式']
    C = X37.C_of(Pl, F2)
    fams, e, U, rho = X37.fams_of(Pl, F2, C, phi0)
    FLOOR = phi0 + float(np.trace(Pl.W @ Pl.Pc))
    g1 = min(f['gap'] for f in fams if f['a'] == 1)
    ffull = [f for f in fams if f['a'] == rho][0]
    f0 = [f for f in fams if f['name'] == '0'][0]
    print(r'  地板 $D_{\min}=%.6f$  $\mathrm{rank}C=%d$  $\lambda(C)=%s$  '
          r'$\mathrm{gap}_{\{0\}}=%.6f$  $\mathrm{gap}_{\{1\}}=%.4f$'
          % (FLOOR, rho, np.array2string(e[:rho], precision=4), g1,
             [f['gap'] for f in fams if f['name'] == '1'][0]))

    # ---------- [4] 遮挡与拐点前斜率 ----------
    D, I, ns, Itop, cut, leg, box = red_columns()
    print('\n===== [4] 红曲线逐列读（图例框 %s，图框 %s） =====' % (leg, box))
    print('  列数 %d，其中被图例框/图框切断 %d 列' % (len(D), int(cut.sum())))
    o = np.argsort(D)
    D, I, ns, Itop, cut = D[o], I[o], ns[o], Itop[o], cut[o]
    if cut.sum():
        print('  被切断的列 $D$ 范围：%.4f..%.4f（%d 列）' % (D[cut].min(), D[cut].max(), cut.sum()))
    print('  红线最左可见 $D=%.4f$，顶端 $I=%.4f$；前 12 列：'
          % (D.min(), Itop.max()) + '  '.join('%.3f/%.3f/n%d/%s'
                                              % (D[i], I[i], ns[i], '切' if cut[i] else 'ok')
                                              for i in range(12)))
    u_of = lambda d: np.log2(1.0 / (d - FLOOR))
    print('  预言拐点 $D=%.4f$（$x=%.4f$），它左边共 %d 列，其中可用 %d 列。'
          % (FLOOR + g1, g1, (D < FLOOR + g1).sum(), ((D < FLOOR + g1) & ~cut).sum()))
    print('  斜率 $\mathrm d\Delta/\mathrm d\log_2(1/x)$（未切断列；拐点前段 $\{0,1\}$ 族的预言 $\approx0.85$）：')
    ugrid = np.concatenate([np.arange(46.90, 48.50, .05), np.arange(48.6, 60.1, .2)])
    uvals = np.array([unc(q) for q in ugrid])
    uf = lambda dq: float(np.interp(dq, ugrid, uvals))
    DEL = I - np.interp(D, ugrid, uvals)             # 图上的 $\Delta(D)$（用我的解析 $I_{unc}$ 还原）
    print('  我的 $I_{unc}$ 网格 %d 点，最近地板处 $I_{unc}(47.0)=%.4f$，$\mathrm d I_{unc}/\mathrm dD=%.4f$ bit/单位'
          % (len(ugrid), uf(47.0), (uf(47.05) - uf(46.95)) / .1))
    capd = D.min() + 3 * np.median(np.diff(D))
    slope_fit(D, DEL, (D >= D.min()) & (D < capd), FLOOR, '线端帽 3 列（不参与回归）')
    slope_fit(D, DEL, (D >= capd) & (D < FLOOR + g1), FLOOR, '拐点前 $x<\mathrm{gap}$（去帽）')
    slope_fit(D, DEL, (D >= FLOOR + g1) & (D < FLOOR + 1.7), FLOOR, '拐点后 $x\in[1.08,1.7]$')
    slope_fit(D, DEL, (D >= FLOOR + 1.7) & (D < FLOOR + 3.0), FLOOR, '$x\in[1.7,3]$（$\{0\}$ 族）')
    slope_fit(D, DEL, (D >= FLOOR + 3.0) & (D < 60.0), FLOOR, '$x\in[3,14]$')

    # ---------- [5] 族切换地点 + 数值最优的谱 ----------
    print('\n===== [5] 族切换与 $\Gamma^\ast$ 的秩：原文 F2 曲线画的区间里弱通道值多少？ =====')
    from p0.exp_fig1_ladder import envelope          # 复用一维二分
    xs = np.concatenate([np.arange(.82, 1.7, .06), np.arange(1.7, 6.01, .25)])

    def ray(fm, x):
        xA = x - fm['gap']
        if xA <= 30 * fm['res']:
            return None
        lo, hi = -15.0, 3.0
        for _ in range(46):
            mid = .5 * (lo + hi)
            try:
                xF = X37.pv(Pl, fm['F'], 10 ** mid * np.diag(1.0 / fm['w']), phi0)[0]
            except Exception:
                return None
            if xF - fm['gap'] < xA:
                lo = mid
            else:
                hi = mid
        return X37.pv(Pl, fm['F'], 10 ** (.5 * (lo + hi)) * np.diag(1.0 / fm['w']), phi0)[1]
    last = None
    cross = []
    for x in xs:
        a_, b_ = ray(ffull, x), ray(f0, x)
        if a_ is None or b_ is None:
            continue
        who = '0+1' if a_ <= b_ else '0'
        if last is not None and who != last:
            cross.append(x)
        print('  %8.3f %10.4f %10.4f %10.4f %10s %+10.4f' % (x, x + FLOOR, a_, b_, who, abs(a_ - b_)))
        last = who
    print('  切换点（$E$ 口径）$x^*\approx%s$ → $D^*\approx%s$' %
          (', '.join('%.2f' % q for q in cross), ', '.join('%.2f' % (q + FLOOR) for q in cross)))
    print('\n  数值锥最优 $\Gamma^\ast$ 在 $C$ 本征基下的谱（比值大 = 弱通道实质关闭）：')
    print('  %8s %8s %10s %10s %12s %12s %12s %10s' %
          ('D', 'x', 'I_cone', 'E(min)', '$\lambda(V)$ 两方向', '比值', '锥最优-红', '红像素'))
    prev_G, prev_x = None, None
    Ci = U[:, :rho] @ np.diag(1.0 / e[:rho]) @ U[:, :rho].T
    for Dq in [47.0, 47.2075, 47.5, 47.8, 48.1231, 48.5, 49.5, 51.0, 55.0, 60.0]:
        x = Dq - FLOOR
        lg = gamma_iso(F2, Dq)
        starts = ([pack(np.exp(lg) * np.eye(2))] if lg is not None else [])
        if prev_G is not None:
            starts += [pack(prev_G * (prev_x / x))]
        r_ = solve_cone(F2, Dq, FLOOR, starts)
        if r_ is None:
            print('  %8.4f  失败' % Dq); continue
        It, Jt, Gt = r_[0], r_[1], r_[2]
        prev_G, prev_x = Gt, x
        lam = np.linalg.eigvalsh(np.linalg.pinv(Gt))[::-1]      # V 的本征值
        # 投到 C 的本征基，看两个方向各自的噪声方差
        Vc = U[:, :rho].T @ np.linalg.pinv(Gt) @ U[:, :rho]
        vc = np.linalg.eigvalsh((Vc + Vc.T) / 2.)[::-1]
        env = envelope(phi0, fams, x)
        sel = D[np.abs(D - Dq) < .03]
        selp = I[np.abs(D - Dq) < .03]
        cutp = cut[np.abs(D - Dq) < .03]
        rd = ('%.4f' % selp[~cutp].mean()) if np.any(~cutp) else '遮挡'
        print('  %8.4f %8.4f %10.4f %10.4f  [%9.3g,%9.3g] r=%7.1f %10s %10s'
              % (Dq, x, It, env[0] if env else float('nan'), vc[0], vc[-1],
                 vc[0] / max(vc[-1], 1e-30),
                 '%+.4f' % (It - float(rd)) if rd != '遮挡' else '—', rd))
