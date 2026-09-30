r"""E39：把阶梯的第一格**钉到原文图 1 的横轴上**，逐像素检验。

锚点植物就是原文的 $n=4$ 算例（E35 里 `PT=Pl(base.A, base.B, base.W, I)`），所以 E37 的
$\mathrm{gap}_A$ 已经是图 1 横轴（cost bound $D$）的单位：
$$\text{预言拐点 } D=\Phi_0^{F}+\operatorname{tr}(W P_c)+\mathrm{gap}_{\{0\}}.$$
拐点以下只有全支撑族可行，斜率应当是 $\operatorname{rank}C/2=1$ bit/十倍程（$\log_2$ 口径）。

三件事，按检验强度排：
 (1) **参数自由的公开曲线对表**：我的精确 $I_{unc}(D)$（Tanaka 式(18) max-det SDP）对图上黑曲线
     （1393 个像素）与品红曲线（满秩 $F$，Prop 3 说它应当与黑曲线重合）。这两个是逐点等式，
     不涉及任何渐近近似，等于把图 1 的轴标定与作者的无约束解一起复核。
 (2) **预言的梯级斜率**：$\Delta(D)=I_{TRV}(D)-I_{unc}(D)$ 的局部斜率（对 $\log_2(1/x)$）在拐点
     两侧各自应当是多少；用族射线上界 $E(D)=\min_A I_A(x)$ 给预测，用红像素给测量，
     误差棒来自像素离散。
 (3) **射线在拐点处是否真的紧**：同一 $D$ 上跑 `solve_cone`（数值最优），看 $E(D)$ 与它的差；
     若差 $\ll$ 像素噪声，(2) 的预测就是真曲线的预测。
"""
import sys
import numpy as np
sys.stdout.reconfigure(encoding='utf-8')
import p0.p0_replicate_letter as base
import p0.exp_plants_oos as E35
import p0.exp_ladder as X37
import p0.exp_C_exact as X36
from p0.exp_nearfloor_law import Fs, solve_cone, gamma_iso, pack, unpack, ij
from p0.p0_replicate_letter import n, sym, unconstrained_sdp

ln2 = np.log(2.0)
F2 = Fs['F2 尾部2(封闭, r=2)']
Pl = E35.Pl(base.A, base.B, base.W, np.eye(4))


def prep():
    fx = X37.floor_exact(Pl, F2, None)
    phi0 = fx['闭式']
    C = X37.C_of(Pl, F2)
    fams, e, U, rho = X37.fams_of(Pl, F2, C, phi0)
    unc0 = float(np.trace(Pl.W @ Pl.Pc))
    print(r'  $\Phi_0$ 闭式 $=%.12f$（残差 %.1e，%d 步）  $\operatorname{tr}(WP_c)=%.6f$  '
          r'地板 $D_{\min}=%.6f$' % (phi0, fx['res'], fx['iters'], unc0, phi0 + unc0))
    print(r'  $\operatorname{rank}C=%d$  $\lambda(C)=%s$' % (rho, np.array2string(e[:rho], precision=5)))
    for f in fams:
        print('    族 %-6s $|A|=%d$  $\mathrm{gap}=%+.6f$  $\Phi_0(A)=%.8f$  $w=%s$'
              % (f['name'], f['a'], f['gap'], f['phi0'], np.array2string(f['w'], precision=4)))
    return phi0, C, fams, rho, e, unc0


def ray_I(Pl, f, x, phi0):
    """族射线上一维二分：求 $t$ 使 $x(t)=\mathrm{gap}+x_A$，返回该 feasible 点的率 $I$。"""
    xA = x - f['gap']
    if xA <= 30 * f['res']:
        return None, None
    lo, hi = -15.0, 3.0
    for _ in range(46):
        mid = .5 * (lo + hi)
        try:
            xF = X37.pv(Pl, f['F'], 10 ** mid * np.diag(1.0 / f['w']), phi0)[0]
        except Exception:
            return None, None
        if xF - f['gap'] < xA:
            lo = mid
        else:
            hi = mid
    t = 10 ** (.5 * (lo + hi))
    xg, I, _P, _Pt, rr = X37.pv(Pl, f['F'], t * np.diag(1.0 / f['w']), phi0)
    return I, abs(xg - x)


def envelope(phi0, fams, x):
    cand = []
    for f in fams:
        I, err = ray_I(Pl, f, x, phi0)
        if I is not None and err < .02 * x + 1e-13:
            cand.append((I, f['name'], f['a']))
    if not cand:
        return None
    I, nm, a = min(cand)
    return I, nm, a, sorted((round(c[0], 6), c[1]) for c in cand)


if __name__ == '__main__':
    phi0, C, fams, rho, elam, unc0 = prep()
    FLOOR = phi0 + unc0
    red = np.load('p0/fig/fig1_red.npy')
    blk_all = np.load('p0/fig/fig1_black.npy')
    mag = np.load('p0/fig/fig1_magenta.npy')
    # 黑数组里混着那条水平虚线（$\sum\log|\lambda_u|$ 的标定线）：449 个像素全在 $I=1.1703$。
    dash = np.abs(blk_all[:, 1] - 1.1703) < .004
    blk = blk_all[~dash]
    print(r'  虚线（标定自检）：%d 个像素，读出 $I=%.4f$，理论 $\sum\log|\lambda_u|=%.4f$，差 %+.4f bit'
          % (dash.sum(), blk_all[dash, 1].mean(), base.plant_facts()[2],
             blk_all[dash, 1].mean() - base.plant_facts()[2]))

    _ucache = {}

    def unc(D):
        key = round(D, 9)
        if key not in _ucache:
            _ucache[key] = unconstrained_sdp(D, np.eye(n), Pl.Th, Pl.Pc)['I_true']
        return _ucache[key]

    def unc_seg(lo, hi):
        """段内稠密网格 + 线性插值（近地板段网格 0.02，长段自动放粗；插值误差 $\ll$ 像素噪声）。"""
        st = max(.02, (hi - lo) / 40.)
        g = np.arange(lo, hi + 1e-9, st)
        return lambda D: float(np.interp(D, g, [unc(q) for q in g]))

    # ---------- (1) 无约束率对表：我的 SDP vs 图上黑/品红像素 ----------
    print('\n===== [1] $I_{unc}(D)$ 的闭式对表（无渐近近似参与） =====')
    Dgrid = np.concatenate([np.arange(46.95, 48.01, .05), np.arange(48.1, 50.01, .2),
                            np.arange(51., 60.1, 1.), np.arange(62., 90.1, 4.)])
    rows = []
    for D in Dgrid:
        Iu = unc(D)
        sb = blk[np.abs(blk[:, 0] - D) < .06]
        sm = mag[np.abs(mag[:, 0] - D) < .06]
        rows.append((D, Iu, len(sb), sb[:, 1].mean() if len(sb) else np.nan,
                     np.nanstd(sb) if len(sb) > 1 else np.nan,
                     len(sm), sm[:, 1].mean() if len(sm) else np.nan))
    for r_ in rows[:8] + rows[16:20] + rows[::9]:
        print('  $D=%6.2f$  我的 $I_{unc}=%.4f$  黑像素 %3d 均值 %7s  散 $\\pm$%5s  品红 %3d %7s'
              % (r_[0], r_[1], r_[2], '%.4f' % r_[3] if np.isfinite(r_[3]) else '—',
                 '%.4f' % r_[4] if np.isfinite(r_[4]) else '—', r_[5],
                 '%.4f' % r_[6] if np.isfinite(r_[6]) else '—'))
    d_b = np.array([r_[3] - r_[1] for r_ in rows if r_[2] >= 2])
    print('  黑曲线 $-$ 我的 SDP：%d 格，均值 %+.4f，标准差 %.4f，最大 $|$%.4f$（$\le$ 半个像素线宽即算通过）'
          % (len(d_b), d_b.mean(), d_b.std(), np.abs(d_b).max()))
    d_m = np.array([r_[6] - r_[1] for r_ in rows if r_[5] >= 2])
    if len(d_m):
        print('  品红（满秩 $F$，Prop 3 应当 $=$ 黑）$-$ 我的 $I_{unc}$：%d 格，均值 %+.4f，最大 $|$%.4f$'
              % (len(d_m), d_m.mean(), np.abs(d_m).max()))

    # ---------- (3) 射线在拐点处是否紧 ----------
    print('\n===== [3] 族射线上界 $E(D)=\min_A I_A$ 与数值锥最优对表（检验预测的松紧） =====')
    unc_f = unc
    prev_G, prev_x = None, None
    print('  %8s %9s %10s %10s %10s %10s %s' %
          ('D', 'x', 'E(D)界', 'solve_cone', '界-数值', '最优族', '$|A|$'))
    cone_ok = []
    for D in [47.0, 47.1, 47.2, 47.3, 47.5, 48.0, 48.5, 49.5, 50.5]:
        x = D - FLOOR
        env = envelope(phi0, fams, x)
        lg = gamma_iso(F2, D)
        starts = ([pack(np.exp(lg) * np.eye(2))] if lg is not None else [])
        if prev_G is not None:
            starts += [pack(prev_G * (prev_x / x))]
        r_ = solve_cone(F2, D, FLOOR, starts)
        if env is None or r_ is None:
            print('  %8.2f %9.4f   ---' % (D, x)); continue
        It, Jt, Gt = r_[0], r_[1], r_[2]
        prev_G, prev_x = Gt, x
        cone_ok.append((D, env[0], It))
        print('  %8.2f %9.4f %10.4f %10.4f %+10.4f %10s %s'
              % (D, x, env[0], It, env[0] - It, env[1], env[2]))

    # ---------- (2) 局部斜率：预言 vs 像素 ----------
    print('\n===== [2] $\Delta(D)=I_{TRV}-I_{unc}$ 对 $\log_2(1/x)$ 的局部斜率 =====')
    print('  预言：$x<\\mathrm{gap}_{\{0\}}=%.4f$ 只有全支撑族可行 $\Rightarrow$ 斜率 $=%.1f$；'
          '拐点以上 $\{0\}$ 族（斜率 $0.5$）入场，包络应当下折。'
          % (min(f['gap'] for f in fams if f['a'] == 1) if any(f['a'] == 1 for f in fams) else float('nan'),
             rho / 2.))
    print('  拐点在图上的横轴位置 $D=%.4f$；红线左端 $D=%.4f$，即拐点之前只有 %.4f 个横轴单位 $\approx$ %d 列像素。'
          % (FLOOR + min(f['gap'] for f in fams if f['a'] == 1), red[:, 0].min(),
             FLOOR + min(f['gap'] for f in fams if f['a'] == 1) - red[:, 0].min(),
             (red[:, 0] < FLOOR + min(f['gap'] for f in fams if f['a'] == 1)).sum()))
    # 预言斜率：直接在 $u=\log_2(1/x)$ 上数值微分 $E(D)-I_{unc}(D)$
    pred = []
    for xq in (0.88, 0.95, 1.02, 1.0844, 1.15, 1.3, 1.6, 2.0, 2.6, 3.5, 5., 8., 13.):
        h = .012 * xq
        e_lo, e_hi = envelope(phi0, fams, xq - h), envelope(phi0, fams, xq + h)
        if e_lo is None or e_hi is None:
            continue
        d_lo, d_hi = e_lo[0] - unc_f(xq - h + FLOOR), e_hi[0] - unc_f(xq + h + FLOOR)
        pred.append((xq, (d_lo - d_hi) / np.log2((xq + h) / (xq - h)), e_hi[1], e_hi[2]))
    print('  预言斜率 $\mathrm d\Delta/\mathrm du$（$u=\log_2(1/x)$，$x=D-D_{\min}$）：')
    for xq, s_, nm, a_ in pred:
        print('    $x=%6.3f$（$D=%7.4f$）  %+.4f   最优族 %s $|A|=%d$' % (xq, xq + FLOOR, s_, nm, a_))
    # 测量：红像素 - 我的 Iunc，在拐点两侧回归
    u_of = lambda D: np.log2(1.0 / max(D - FLOOR, 1e-9))
    for tag, lo_, hi_ in (('拐点前 $x\in[0.86,1.084]$', red[:, 0].min(), FLOOR + 1.0844),
                          ('拐点后 $x\in[1.084,1.6]$', FLOOR + 1.0844, FLOOR + 1.6),
                          ('$x\in[1.6,3]$', FLOOR + 1.6, FLOOR + 3.0),
                          ('$x\in[3,13.7]$（图窗其余）', FLOOR + 3.0, 59.9)):
        uf = unc_seg(lo_, hi_)
        s = red[(red[:, 0] >= lo_) & (red[:, 0] < hi_)]
        if len(s) < 4:
            print('    %-28s 像素 %d 个，太少' % (tag, len(s))); continue
        uu = np.array([u_of(d) for d in s[:, 0]])
        dd = np.array([s[i, 1] - uf(s[i, 0]) for i in range(len(s))])
        M = np.vstack([uu, np.ones_like(uu)]).T
        coef, res, *_ = np.linalg.lstsq(M, dd, rcond=None)
        resid = dd - M @ coef
        sd = np.sqrt(np.sum(resid ** 2) / max(len(dd) - 2, 1))
        cov = sd ** 2 * np.linalg.inv(M.T @ M)
        print('    %-28s 像素 %4d  斜率 %+.4f $\pm$ %.4f bit/dec  截距 %+.3f  残差 $\sigma$=%.4f bit'
              % (tag, len(s), coef[0], np.sqrt(cov[0, 0]), coef[1], sd))
