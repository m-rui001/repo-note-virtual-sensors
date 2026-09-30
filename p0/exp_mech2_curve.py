r"""E44：机制二的见证植物上，"免费删"到底免费到哪一步——地板、近地板斜率、还是整条曲线。

E43 只证明了删掉 $\Theta$ 看不见的通道**不抬地板**（$\mathrm{gap}=0$，$\Phi_0(F)=\Phi_0(F_{-d})$）。
这句话的强度还不清楚，因为三种可能都跟它相容：

  (a) 两条 $E(D)$ 曲线在渐近段重合，只在大的 $x$ 分开；
  (b) 曲线处处重合（那"免费"就是 $\Delta_{\rm TRV}$ 层面的恒等，比门票强得多）；
  (c) 曲线处处分开，近地板斜率也不同——那门票为 0 只是个孤立的代数巧合。

判据现成：E23 的对数律给 $\Delta(D)=\frac{r_{\rm act}}2\log_2\frac1x+O(1)$，E24/B18 把系数从
$\operatorname{rank}F$ 改成了 $\operatorname{rank}C$。机制植物的 $\operatorname{rank}\Theta=1$，
第二条通道只往 $\ker\Theta$ 里买精度，所以预测是 $r_{\rm act}(F)=r_{\rm act}(F_{-d})=1$：
**同地板、同斜率、近段重合**。这条预测可证伪——只要 $\operatorname{rank}C(F)$ 真等于 $\operatorname{rank}C(F_{-d})$
而两曲线在 $x\to0$ 仍分开，我自己的斜率理论就漏了东西。

第二个看点是耦合：$F=\{e_1^\top,e_3^\top\}$ 允许把信道噪声做成相关的（$V$ 非对角），
于是第二条通道能变相给 $x_1$ 买精度——这是 $F_{-d}$ 一族根本做不到的动作。若它真的换来率，
分离点会出现在中等 $x$ 而不是 $x\to0$。所以本轮把两族曲线在 $x$ 从 $10^{-6}$ 到 $10^{0}$ 全程对拉，
而不是只看端点。

口径全部闭式（奇异 DARE + 抛光）或逐点带证书；下界腿用 $Pl.sdp(D)$（Tanaka 式(18) 的 max-det 松弛，
即作者黑线 $I_{unc}$ 的解析版，E42 已把它对到像素 $\pm0.02$）。
"""
import sys
import numpy as np
sys.stdout.reconfigure(encoding='utf-8')
from scipy.optimize import minimize
from scipy.linalg import solve_discrete_are
import p0.exp_plants_oos as E35
import p0.exp_ladder as X37
import p0.exp_C_exact as X36
from p0.exp_nearfloor_law import pack, unpack

sym = X36.sym
ln2 = np.log(2.0)
RT = 1e-9


def plant(c13=0.0, au=.6, b3=0.0, wu=1.0):
    A = np.array([[1.25, .30, c13], [0, .85, 0], [0, 0, au]])
    B = np.array([[1.], [0.], [b3]])
    return E35.Pl(A, B, np.diag([1.0, .7, wu]), np.eye(3), np.array([[1.0]]))


F  = np.array([[1., 0, 0], [0, 0, 1]])      # 行 1 = 冗余通道 $e_3^\top$（$\Theta$ 看不见）
Fm = np.delete(F, 1, axis=0)
Pl = plant()
UNC = Pl.UNC
wT = np.linalg.eigvalsh(sym(Pl.Th))
print('== E44 ==  机制二植物上的全程 $E(D)$ 对拉  $n=%d$  $\\lambda(\\Theta)$=%s  $\\operatorname{rank}\\Theta$=%d'
      % (Pl.nn, ', '.join('%.3g' % v for v in wT), int((wT > 1e-9 * max(wT.max(), 1e-30)).sum())))

fam = {}
for tag, Fx in [('F 完整($r$=2, 含冗余通道)', F), ('F 删冗余($r$=1)', Fm)]:
    ok, lam = Pl.detectable(Fx)
    fl = X37.floor_exact(Pl, Fx, None)
    phi0 = fl['闭式']
    C = X37.C_of(Pl, Fx)
    rho, elam, U = X36.rank_C(C)
    fam[tag] = dict(F=Fx, r=Fx.shape[0], phi0=phi0, floor=phi0 + UNC, C=C, rho=rho, elam=elam)
    print('\n  %-24s 可检测=%s  $\\Phi_0$=%.9f  地板 $D_{\\min}$=%.9f  DARE 残差 %.1e'
          % (tag, '是' if ok else '否', phi0, phi0 + UNC, fl['res']))
    print('     $\\lambda(C)$=%s  $\\operatorname{rank}C$=%d（预测：两族都应是 1，即 $C$ 的第二个本征方向是 $e_3$）'
          % (', '.join('%.6g' % v for v in elam), rho))

gf = fam['F 完整($r$=2, 含冗余通道)']
gd = fam['F 删冗余($r$=1)']
print('\n  门票复核：$\\Phi_0(F)-\\Phi_0(F_{-d})=%.3e$（E43 的 $\\mathrm{gap}=0$ 在这里应精确复现）'
      % (gf['phi0'] - gd['phi0']))
print('  两族地板之差 $=%.3e$ = 门票；非 0 说明我的抛光不够，本轮全部读数作废。')


def xI(Fx, V, phi0):
    return Pl.xI(Fx, sym(V), phi0)


def solve_iso(Fx, phi0, xt, lo=-16., hi=14.):
    r"""各向同性分支 $V=vI$：$x(v)$ 单调减，二分出 $v$。返回 $(I,x,v)$ 或 None。"""
    r = Fx.shape[0]
    f = lambda lg: xI(Fx, np.exp(lg) * np.eye(r), phi0)[0] - xt
    a, b = f(lo), f(hi)
    if a * b > 0:
        return None
    for _ in range(200):
        m = .5 * (lo + hi)
        if f(m) * a <= 0:
            hi = m
        else:
            lo = m
    lg = .5 * (lo + hi)
    xv, I, _ = xI(Fx, np.exp(lg) * np.eye(r), phi0)
    return I, xv, np.exp(lg)


def solve_cone_local(Fx, phi0, xt, starts, budget=3000):
    r"""$V=L L^\top$（Cholesky，对角取 log），罚函数 Nelder-Mead 再 SLSQP 硬约束。返回 $(I,x,V)$。"""
    r = Fx.shape[0]

    def IJ(p):
        try:
            x, I, _ = xI(Fx, unpack(p), phi0)
        except Exception:
            return 1e6, 1e6
        if not (np.isfinite(I) and np.isfinite(x)):
            return 1e6, 1e6
        return I, x

    best = None
    for p0 in starts:
        pen = lambda p: (lambda I, x: I + 2000. * max(0., x - xt) ** 2 + 200. * max(0., x - xt))(*IJ(p))
        r1 = minimize(pen, p0, method='Nelder-Mead',
                      options=dict(maxiter=budget, xatol=1e-11, fatol=1e-13))
        r2 = minimize(lambda p: IJ(p)[0], r1.x, method='SLSQP',
                      constraints=[dict(type='ineq', fun=lambda p: xt - IJ(p)[1])],
                      options=dict(maxiter=400, ftol=1e-13))
        for cand in (r1.x, r2.x):
            I, x = IJ(cand)
            if x <= xt * (1 + 1e-9) and x > -1e-9 and (best is None or I < best[0]):
                best = (I, x, unpack(cand))
    return best


print('\n[2] 两族曲线全程对拉（$x=D-D_{\\min}$，用**同一块地板**，因为门票是 0）')
print('  %10s %12s %12s %12s %12s %12s %10s %10s' %
      ('$x$', '$E_{\\rm full}$', '$E_{\\rm del}$', '$E_{\\rm full}-E_{\\rm del}$',
       r'$E_{\rm iso,full}$', r'$E_{\rm iso,del}$', r'$\rho(V^*)$', r'$\lambda(V^*)$'))
xs = np.logspace(-6, 0., 13)
cur = {'full': [], 'del': []}
prev = {}
for x in xs:
    row = []
    iso = {}
    for key, fx, phi0 in [('full', gf['F'], gf['phi0']), ('del', gd['F'], gd['phi0'])]:
        r = fx.shape[0]
        got_iso = solve_iso(fx, phi0, x)
        iso[key] = got_iso[0] if got_iso else float('nan')
        if got_iso is None:
            row.append((float('nan'), float('nan'), np.eye(r)))
            continue
        st = [pack(np.eye(r) / max(got_iso[2], 1e-30))]
        if key in prev:
            st.append(pack(prev[key]))
        if r == 1:
            got = (got_iso[0], got_iso[1], np.array([[got_iso[2]]]))   # r=1：各向同性就是全体
        else:
            got = solve_cone_local(fx, phi0, x, st)
            if got is None:
                got = (got_iso[0], got_iso[1], np.eye(r) / max(got_iso[2], 1e-30))
        prev[key] = got[2]
        row.append(got)
    Ifull, xfull, Vfull = row[0]
    Idel, xdel, Vdel = row[1]
    rho = Vfull[0, 1] / np.sqrt(Vfull[0, 0] * Vfull[1, 1]) if Vfull.shape == (2, 2) else float('nan')
    lv = ', '.join('%.3g' % v for v in np.linalg.eigvalsh(Vfull)) if Vfull.shape == (2, 2) else '—'
    cur['full'].append((x, Ifull, xfull))
    cur['del'].append((x, Idel, xdel))
    print('  %10.3e %12.6f %12.6f %12.6f %12.6f %12.6f %10.3f %10s'
          % (x, Ifull, Idel, Ifull - Idel, iso['full'], iso['del'], rho, lv))

print('\n[3] 近地板斜率：两族各拟合 $E=a\\log_2(1/x)+b$，比较 $a$ 与 $\\operatorname{rank}C/2$')
for key, fx, rho in [('full', gf['F'], gf['rho']), ('del', gd['F'], gd['rho'])]:
    pts = [p for p in cur[key] if p[0] < 0.05 and np.isfinite(p[1])]
    lx = np.log2(1. / np.array([p[0] for p in pts]))
    Iv = np.array([p[1] for p in pts])
    a, b = np.polyfit(lx, Iv, 1)
    print(r'  %-5s 窗口 $x\in[%.1e,%.2f]$ 共 %d 点：斜率 $a$=%.4f   $\operatorname{rank}C/2$=%.4f   截距 $b$=%.4f   拟合残差最大 %.3e'
          % (key, min(p[0] for p in pts), max(p[0] for p in pts), len(pts), a, rho / 2., b,
             np.max(np.abs(Iv - (a * lx + b)))))
print(r'  两族斜率一致 $\Rightarrow$ 第二条通道在渐近段完全 inactive，"免费"至少到 (a) 级；')
print(r'  斜率不一致 $\Rightarrow$ 我的 $\operatorname{rank}C$ 理论漏掉了 $\ker\Theta$ 通道之间的相关噪声耦合。')

print('\n[4] 分离点定位：$E_{\\rm full}-E_{\\rm del}$ 随 $x$ 的符号')
d = [(p[0], p[1] - q[1]) for p, q in zip(cur['full'], cur['del'])]
neg = [q for q in d if q[1] < -1e-9]
pos = [q for q in d if q[1] > 1e-9]
print('  负（冗余通道确实省率）：%s' % (', '.join('$x$=%.2e(%.4f)' % (a, b) for a, b in neg) if neg else '无'))
print('  正（冗余通道反而更贵，即优化器没跑赢）：%s' % (', '.join('$x$=%.2e(%.4f)' % (a, b) for a, b in pos) if pos else '无'))
print('  $|E_{\\rm full}-E_{\\rm del}|$ 最大 %.5f，最小 %.2e' % (max(abs(b) for _, b in d), min(abs(b) for _, b in d)))

print('\n[5] 下界腿与两族的表示惩罚 $\\Delta(D)=E_{\\rm TRV}(D)-I_{unc}(D)$')
print('  %10s %12s %10s %10s %10s %10s' % ('$D$', '$I_{\\rm sdp}$', '$\\Delta_{\\rm full}$', '$\\Delta_{\\rm del}$', '证书残差', '状态'))
for (x, Ifull, xf), (x2, Idel, xd) in zip(cur['full'], cur['del']):
    D = x + gf['floor']
    Isdp = Pl.sdp(D)
    if Isdp is None:
        print('  %10.4f  SDP 无解' % D)
        continue
    print('  %10.4f %12.6f %10.6f %10.6f %10.1e  %s'
          % (D, Isdp, Ifull - Isdp, Idel - Isdp, max(xf - x, xd - x), 'ok'))
print('  $\\Delta$ 用的是同一块地板，所以两族 $\\Delta$ 之差就是 [2] 的列 4；')
print('  近地板处两族 $\\Delta$ 应同步发散（同系数），中等 $x$ 才可能分开。')

print('\n[6] 两族重合到底是定理还是我的优化器偷懒')
print(r'  $E_{\rm full}\le E_{\rm del}$ 恒成立（把第二条通道的噪声推到无穷就是 $F_{-d}$ 一族），')
print(r'  所以 [2] 的列 4 只可能是 $0$ 或负。要判"真重合"，看放大后的锥里还有没有一阶下降方向。')
rng = np.random.default_rng(20260929)
print(r'  基线一阶条件：$\min I$ s.t. $x\le xt$ 的 KKT 是 $\partial_j I+\mu\,\partial_j x=0$，'
      r'$\mu=-\partial_1 I/\partial_1 x$。')
print(r'  三个参数 $\theta=(\log V_{11},\log V_{22},\rho)$ 在放大族里都是内部自由的，所以残差必须全为 0；')
print(r'  哪个残差为负，就是存在一阶下降方向（我的解不最优）。')
print('  %9s %11s %11s %10s | %11s %11s %11s %11s || %11s %11s %s' %
      ('$x$', '$I_{\\rm del}$', '$I_{\\rm full}$', '差', r'$\mu$',
       r'$\partial_1 I$', r'$\partial_2 I$', r'$\partial_\rho I$',
       r'$\partial_2 x$', r'$\partial_\rho x$', '一阶残差 → 判定'))
for xq in [1e-4, 1e-2, 1e-1]:
    iso1 = solve_iso(gd['F'], gd['phi0'], xq)
    v1 = iso1[2]
    full = solve_cone_local(gf['F'], gf['phi0'], xq,
                            [pack(np.diag([v1, 1e8])), pack(np.eye(2) / max(v1, 1e-30))])
    V2 = full[2]
    a1, a2 = V2[0, 0], V2[1, 1]
    rho0 = V2[0, 1] / np.sqrt(a1 * a2)

    def IOfTh(th):
        p1, p2 = np.exp(th[0]), np.exp(th[1])
        rr = np.tanh(th[2])
        Vq = np.array([[p1, rr * np.sqrt(p1 * p2)], [rr * np.sqrt(p1 * p2), p2]])
        xv, Iq, _ = xI(gf['F'], Vq, gf['phi0'])
        return Iq, xv

    th0 = np.array([np.log(a1), np.log(a2), np.arctanh(rho0)])
    g = np.zeros((2, 3))
    h = 1e-5
    for j in range(3):
        tp, tm = th0.copy(), th0.copy()
        tp[j] += h
        tm[j] -= h
        fp, fm = IOfTh(tp), IOfTh(tm)
        g[0, j] = (fp[0] - fm[0]) / (2 * h)
        g[1, j] = (fp[1] - fm[1]) / (2 * h)
    mu = -g[0, 0] / g[1, 0]
    res = [g[0, j] + mu * g[1, j] for j in range(3)]
    verdict = '驻点' if max(abs(v) for v in res[1:]) < 1e-3 * abs(g[0, 0]) else '非驻点'
    print(r'  %9.2e %11.5f %11.5f %10.2e | %11.3e %11.3e %11.3e %11.3e || %11.3e %11.3e  残差 %s → %s'
          % (xq, iso1[0], full[0], full[0] - iso1[0], mu,
             g[0, 0], g[0, 1], g[0, 2], g[1, 1], g[1, 2],
             ', '.join('%.2e' % v for v in res), verdict))
    # 暴力对照：相关噪声全扫，看有没有任何可行点把率打到 $E_{\rm del}$ 以下
    best = None
    for _ in range(4000):
        a1 = 10 ** rng.uniform(np.log10(v1) - 2, np.log10(v1) + 4)
        a2 = 10 ** rng.uniform(-4, 14)
        rr = rng.uniform(-.999, .999)
        Vq = np.array([[a1, rr * np.sqrt(a1 * a2)], [rr * np.sqrt(a1 * a2), a2]])
        if np.linalg.eigvalsh(Vq)[0] <= 1e-12 * a2:
            continue
        xv, Iq, _ = xI(gf['F'], Vq, gf['phi0'])
        if xv <= xq and (best is None or Iq < best[0]):
            best = (Iq, xv, a1, a2, rr)
    if best is None:
        print('     随机 4000 点没有一个可行（$x\\le xt$）——扫描范围要重设，本轮不据此下结论')
    else:
        print('     随机 4000 点里最好 $I$=%.6f（$x$=%.3e，$V=%s$，$\\rho$=%.3f）  vs $E_{\\rm del}$=%.6f  差 $=%.2e$'
              % (best[0], best[1], '%.2e/%.2e' % (best[2], best[3]), best[4], iso1[0],
                 best[0] - iso1[0]))
print(r'  读法：$I_{\rm full}-I_{\rm del}\ge0$ 恒成立（放大族只会更便宜），所以关键在残差符号。')
print(r'  残差 $\approx0$ $\Rightarrow$ 求解器找到的是放大锥里的驻点，"整条曲线免费"是实打实的实验事实（情形 (b)）；')
print(r'  残差 $<0$ $\Rightarrow$ 有方向能一阶降率，那是我的优化器偷懒，[2] 的重合读数要打折。')
print(r'  两种情形都不改门票结论（$\mathrm{gap}=0$ 是闭式，跟优化器无关），改的是"免费"的深度。')
print(r'  真要写成命题，得给"$A$ 块对角 $+\Theta$ 支撑在任务块上 $\Rightarrow$ 第二条通道在任何 $x$ 都不活跃"')
print(r'  的单调性论证——本轮只到数值驻点这一级。')

