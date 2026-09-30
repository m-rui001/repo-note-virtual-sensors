r"""E40：把"这段窗口在不在垂直区"变成一个可代入的判据（#19 欠的门票），并给出台阶合并的精确条件。

设定全是闭式，不需要优化器、不需要 SDP、不需要渐近展开：
  $\Phi_0(F)=\operatorname{tr}(\Theta P_0(F))$，$P_0(F)=\tilde P_0-\tilde P_0F^\top(F\tilde P_0F^\top)^{-1}F\tilde P_0$，
  $\tilde P_0=\mathrm{ARE}(A^\top,F^\top,W,0)$ 抛光到残差 $\sim10^{-16}$。
删掉一个传输方向 $d$（即只传 $(I-dd^\top)$ 那部分）得 $F_{-d}=(I-dd^\top)F$，门槛
$\mathrm{gap}_d=\Phi_0(F_{-d})-\Phi_0(F)=\operatorname{tr}(\Theta\,\Delta P_d)$，$\Delta P_d=P_0(F_{-d})-P_0(F)\succeq0$。
门票：全支撑族还能待在自己的渐近段里，前提是 $x<\min_d\mathrm{gap}_d$。

交付三件事：
 (1) **精确 iff**：$\mathrm{gap}_d=0\iff\operatorname{ran}(\Delta P_d)\subseteq\ker\Theta$
     （$\Delta P_d\succeq0$，所以 $\operatorname{tr}(\Theta\Delta P_d)=0$ 与 $\Theta^{1/2}\Delta P_d\Theta^{1/2}=0$ 等价）。
     这把 #19 的实测"把 $W$ 转到被砍通道的不可观测子空间上台阶就合并"换成一个可检验的代数条件，
     并且明确它跟"噪声被归零"是两回事。
 (2) **两界分解**：$\lambda_{\min}(\Theta)\operatorname{tr}(\Delta P_d)\le\mathrm{gap}_d\le\lambda_{\max}(\Theta)\operatorname{tr}(\Delta P_d)$。
     任务侧只剩一个落在 $[\lambda_{\min},\lambda_{\max}]$ 里的因子，几何部分 $\operatorname{tr}(\Delta P_d)$
     与 $\Theta$ 无关（两个零噪声后验的迹之差）。于是"$\mathrm{gap}$ 能不能用 $\lambda_i(C)$ 当代理"
     被 $\kappa(\Theta)$ 与 $\mathrm{gap}_d/\lambda_d(C)$ 同时限制住，后者本轮直接报数。
 (3) **便宜对照**：不解 $C$、不取本征基，删 $F$ 的行得 $\mathrm{gap}^{\rm row}$。
     若与精确门槛同量级，外人一行就能筛；差得远就必须老实解 $C$。

$\theta$ 网格从 #19 的 7 个点加密到 25 个：固定 $W$ 的谱 $\Sigma=\operatorname{diag}(1,10^{-3},10^{-6},10^{-9})$，
只在 $(e_0,e_1)$ 平面转，$W(\theta)=R\Sigma R^\top$。
"""
import sys
import numpy as np
sys.stdout.reconfigure(encoding='utf-8')
from scipy.linalg import solve_discrete_are
import p0.p0_replicate_letter as base
import p0.exp_plants_oos as E35
import p0.exp_ladder as X37
import p0.exp_C_exact as X36
from p0.exp_nearfloor_law import Fs

sym = X36.sym
RT = 1e-9


def floor_pair(Pl, F):
    r"""$(\Phi_0,P_0,\text{残差},\text{迭代数})$，$P_0$ 为 $V\to0$ 的后验极限。只用闭式。"""
    r = F.shape[0]
    Pt0 = sym(solve_discrete_are(Pl.A.T, F.T, Pl.W, np.zeros((r, r))))
    Pt0, Pk, res, k = X37.polish_Pt0(Pl, F, Pt0)
    return float(np.trace(Pl.Th @ Pk)), Pk, res, k


def gaps_of(Pl, F, drop='eig', C=None):
    """所有"删一个方向"的门槛。drop='eig' 在 $C$ 的本征基里删，'row' 直接删 $F$ 的行。"""
    r = F.shape[0]
    phi0F, P0F, resF, it = floor_pair(Pl, F)
    if drop == 'eig':
        rho, elam, U = X36.rank_C(C if C is not None else X37.C_of(Pl, F))
        dirs = [U[:, i] for i in range(rho)]
        tag = ['$-u_%d$（$\\lambda_C$=%.4g）' % (i, elam[i]) for i in range(rho)]
        lam = elam[:rho]
    else:
        dirs = None
        tag = ['删行 %d' % i for i in range(r)]
        lam = None
    out = []
    for i in range(len(tag)):
        if dirs is None:
            Feff = np.delete(F, i, axis=0)
        else:
            # 取 $d$ 正交补的标准正交基（投影 $I-dd^\top$ 本身奇异，不能直接当传感器）
            wv, Vv = np.linalg.eigh(np.eye(r) - np.outer(dirs[i], dirs[i]))
            Feff = Vv[:, wv > .5].T @ F
        det_ok, det_lam = Pl.detectable(Feff)
        if not det_ok:
            # 奇异 DARE 无稳定解 $\Rightarrow\Phi_0(F_{-d})=+\infty$：这个族在任何预算下都不可行
            out.append(dict(tag=tag[i], gap=float('inf'), tr=float('inf'), ev=np.array([np.nan]),
                            neg=0., lam=(lam[i] if lam is not None else None),
                            rank_dP=None, rank_T=None, dimT=None, resA=None, undet=det_lam))
            continue
        phiA, P0A, resA, _ = floor_pair(Pl, Feff)
        dP = sym(P0A - P0F)
        ev = np.linalg.eigvalsh(dP)
        emax = max(ev[-1], 1e-30)
        wT, QT = np.linalg.eigh(sym(Pl.Th))
        sup = wT > RT * max(wT[-1], 1e-30)
        S = sym(QT[:, sup].T @ dP @ QT[:, sup])
        wS = np.linalg.eigvalsh(S)
        out.append(dict(tag=tag[i], gap=phiA - phi0F, tr=float(np.trace(dP)), ev=ev,
                        neg=float(min(ev[0], 0.)), lam=(lam[i] if lam is not None else None),
                        rank_dP=int(np.sum(ev > RT * emax)),
                        rank_T=int(np.sum(wS > RT * max(wS[-1], 1e-30))),
                        dimT=int(sup.sum()), resA=resA))
    return phi0F, out, resF


def report(nm, Pl, F, with_C=True):
    wT = np.linalg.eigvalsh(sym(Pl.Th))
    print('\n===== %s  n=%d r=%d' % (nm, Pl.nn, F.shape[0]))
    print(r'  $\lambda(\Theta)=[%.4g,%.4g]$  $\kappa(\Theta)=%.3g$  $\operatorname{rank}\Theta$=%d  $\|W\|=%.3g$'
          % (wT[0], wT[-1], wT[-1] / max(wT[0], 1e-30),
             int(np.sum(wT > 1e-9 * wT[-1])), np.linalg.norm(Pl.W, 2)))
    modes = [('eig', '精确（删 $C$ 的本征方向）'), ('row', '便宜（删 $F$ 的行）')] if with_C \
        else [('row', '便宜（删 $F$ 的行）')]
    first = None
    for mode, lab in modes:
        phi0F, g, resF = gaps_of(Pl, F, drop=mode)
        if mode == 'eig':
            first = [q['gap'] for q in g]
        print(r'  -- %s   全支撑地板 $\Phi_0=%.8f$（残差 %.1e）' % (lab, phi0F, resF))
        print('     %-30s %11s %11s %10s %8s %8s %s' %
              ('删除', 'gap', 'tr(dP)', 'gap/tr', 'dP 秩', 'Θ 上秩', '两界 $[\\lambda_{\\min},\\lambda_{\\max}]\\operatorname{tr}$'))
        for q in g:
            if not np.isfinite(q['gap']):
                print('     %-30s %11s %11s %10s %8s %8s %s'
                      % (q['tag'], '+$\\infty$', '—', '—', '—', '—',
                         '不可检测（PBH 在 $|\\lambda|=%.4g$ 失败）$\\Rightarrow$该族在任何预算下都不可行'
                         % abs(q['undet'])))
                continue
            rat = q['gap'] / q['tr'] if q['tr'] > 1e-30 else float('nan')
            print('     %-30s %+11.6f %11.6f %10.4f %8d %8s %s%s'
                  % (q['tag'], q['gap'], q['tr'], rat, q['rank_dP'],
                     '%d/%d' % (q['rank_T'], q['dimT']),
                     '[%.4f,%.4f]' % (wT[0] * q['tr'], wT[-1] * q['tr']),
                     '' if q['neg'] > -1e-12 else '  ! dP 有负本征值 %.1e' % q['neg'])
                  + ('' if q['lam'] is None else '   gap/λ_C=%.3f' % (q['gap'] / q['lam']))
                  + ('   台阶合并（$\\operatorname{ran}\\Delta P\\subseteq\\ker\\Theta$）'
                     if q['rank_T'] == 0 else ''))
        fin = [q['gap'] for q in g if np.isfinite(q['gap'])]
        xs = [v for v in fin if v > 1e-9]
        print('     本口径的门票 $x_{\\min}=%s$   分类：%d 个 $+\\infty$（不可检测）/ %d 个有限正 / %d 个合并（$0$）'
              % ('%.6f' % min(xs) if xs else ('不存在：所有台阶合并' if fin else '不存在：所有删除都不可检测'),
                 sum(1 for q in g if not np.isfinite(q['gap'])), len(xs), len(fin) - len(xs)))
    return phi0F, g


if __name__ == '__main__':
    for nm, pl, F in E35.plants():
        if nm.find('锚点') >= 0:
            report(nm.replace('$', ''), pl, F)
    # 单输入 PL4：$\operatorname{rank}\Theta=1$，$\ker\Theta$ 真大，正是 (1) 该亮起来的地方
    for nm, pl, F in E35.plants():
        if nm.find('PL4 单输入 $r=3$') >= 0:
            print('\n' + r'  ↓ 下面这个算例 $\operatorname{rank}\Theta=1$，(1) 的 iff 应当把"任务看不见"和"后验不动"分开')
            report(nm.replace('$', ''), pl, F, with_C=False)

    # ---------- $\theta$ 网格 7 → 25 ----------
    print('\n\n' + r"""===== $\theta$ 扫描（锚点植物，$W(\theta)=R\Sigma R^\top$，$\Sigma=\mathrm{diag}(1,10^{-3},10^{-6},10^{-9})$，只转 $(e_0,e_1)$） =====""")
    F2 = Fs['F2 尾部2(封闭, r=2)']
    Pl0 = E35.Pl(base.A, base.B, base.W, np.eye(4))
    Sig = np.diag([1., 1e-3, 1e-6, 1e-9])
    print('  %6s %10s %11s %11s %11s %11s %9s %9s %s' %
          ('deg', r'$\Phi_0$', 'gap(行0)', 'gap(行1)', 'tr dP0', 'tr dP1',
           'gap0/tr0', 'gap1/tr1',
           r'$\lambda(C)$  $\kappa(C)$  $\kappa(\Theta)$   门票$^{\rm eig}$  gap/λ_C'))
    rows = []
    for deg in np.linspace(0, 90, 25):
        th = np.deg2rad(deg)
        R = np.eye(4)
        c, s = np.cos(th), np.sin(th)
        R[0, 0] = R[1, 1] = c
        R[0, 1], R[1, 0] = -s, s
        Pl = E35.Pl(Pl0.A, Pl0.B, sym(R @ Sig @ R.T), np.eye(4))
        C = X37.C_of(Pl, F2)
        rho, elam, _U = X36.rank_C(C)
        phi0F, ge, _ = gaps_of(Pl, F2, drop='eig', C=C)
        _p, g, _r = gaps_of(Pl, F2, drop='row')
        wT = np.linalg.eigvalsh(sym(Pl.Th))
        rat = [q['gap'] / q['tr'] if np.isfinite(q['tr']) and q['tr'] > 1e-30 else float('nan')
               for q in g]
        tke = min([q['gap'] for q in ge if np.isfinite(q['gap'])], default=float('inf'))
        gl = [q['gap'] / q['lam'] for q in ge if np.isfinite(q['gap']) and q['lam']]
        print('  %6.1f %12.6g %+11.5f %+11.5f %11.5f %11.5f %9.4f %9.4f  [%.4g,%.4g] %.3g %.3g  '
              r'门票$^{\rm eig}$=%.4f  gap/λ_C=%s'
              % (deg, phi0F, g[0]['gap'], g[1]['gap'], g[0]['tr'], g[1]['tr'], rat[0], rat[1],
                 elam[0], elam[1], elam[0] / max(elam[1], 1e-30),
                 wT[-1] / max(wT[0], 1e-30), tke,
                 ' '.join('%.3f' % v for v in gl)))
        rows.append((deg, [q['gap'] for q in g], [q['tr'] for q in g], elam[:rho], wT,
                     [q['rank_T'] for q in g], [q['gap'] for q in ge], tke, Pl.UNC, phi0F))
    g0 = np.array([q[1][0] for q in rows])
    g1 = np.array([q[1][1] for q in rows])
    tk = np.minimum(g0, g1)
    tke = np.array([q[7] for q in rows])
    tkrow = np.array([min(v for v in q[1] if np.isfinite(v)) for q in rows])
    print('  门票（精确口径，$C$ 本征基）跨 %.3e..%.3e（%.1f 倍）   门票（便宜口径，删行）跨 %.3e..%.3e'
          % (tke.min(), tke.max(), tke.max() / max(tke.min(), 1e-30), tkrow.min(), tkrow.max()))
    print('  两口径之比（便宜/精确）：%.2f..%.2f  ——  删行能不能当代理，看这一列'
          % ((tkrow / tke).min(), (tkrow / tke).max()))
    unc = np.array([q[8] for q in rows])
    ph = np.array([q[9] for q in rows])
    print(r'  同一个网格里三个量的分工：$\Phi_0$ 跨 %.6e..%.6e（相对 %.2f%%），'
          r'$\operatorname{tr}(WP_c)$ 跨 %.4f..%.4f（%.2f 倍，即 $D_{\min}$ 动 %.0f%%），'
          r'精确门票动 %.0f 倍'
          % (ph.min(), ph.max(), (ph.max() / ph.min() - 1) * 100,
             unc.min(), unc.max(), unc.max() / unc.min(), (unc.max() / unc.min() - 1) * 100,
             tke.max() / tke.min()))
    f0 = g0[np.isfinite(g0)]
    print('  gap(行0)：%d/25 个为 $+\\infty$（不可检测）' % int(np.sum(np.isinf(g0)))
          + ('   有限值跨 %.3e..%.3e' % (f0.min(), f0.max()) if len(f0) else '')
          + '   gap(行1) 跨 %.3e..%.3e（%.1f 倍）   门票（删行口径）跨 %.3e..%.3e（%.1f 倍）'
          % (g1.min(), g1.max(), g1.max() / max(g1.min(), 1e-30),
             tk.min(), tk.max(), tk.max() / max(tk.min(), 1e-30)))
    fin_tr = [q[2][k] for q in rows for k in (0, 1) if np.isfinite(q[2][k])]
    rt = [q[1][k] / q[2][k] for q in rows for k in (0, 1)
          if np.isfinite(q[1][k]) and np.isfinite(q[2][k]) and q[2][k] > 1e-12]
    print(r'  两界因子 $\mathrm{gap}/\operatorname{tr}(dP)$：跨 %.4f..%.4f（倍程 %.2f），'
          r'同期纯几何量 $\operatorname{tr}(\Delta P)$ 跨 %.3e..%.3e（倍程 %.0f）'
          % (min(rt), max(rt), max(rt) / min(rt), min(fin_tr), max(fin_tr),
             max(fin_tr) / min(fin_tr)))
    # 注意 $\Theta$ 是 $n\times n$（本例 $n=4$），$\lambda_{\max}$ 是最后一个而不是第二个
    okb = all(q[4][0] * q[2][k] - 1e-10 <= q[1][k] <= q[4][-1] * q[2][k] + 1e-10
              for q in rows for k in (0, 1)
              if np.isfinite(q[1][k]) and np.isfinite(q[2][k]))
    viol = [(q[0], k, q[1][k], q[2][k], q[4][0] * q[2][k], q[4][-1] * q[2][k])
            for q in rows for k in (0, 1)
            if np.isfinite(q[1][k]) and np.isfinite(q[2][k])
            and not (q[4][0] * q[2][k] - 1e-10 <= q[1][k] <= q[4][-1] * q[2][k] + 1e-10)]
    print('  两界 $\\lambda_{\\min}(\\Theta)\\operatorname{tr}(dP)\\le\\mathrm{gap}\\le\\lambda_{\\max}(\\Theta)\\operatorname{tr}(dP)$ '
          '逐点检验：%s（$\\Theta$ 随 $\\theta$ 一起变，所以用的是每个点自己的 $\\lambda(\\Theta)$）'
          % ('全部成立' if okb else '有违例'))
    for v in viol:
        print('     违例 $\\theta=%.1f$ 行%d  gap=%.6e  tr=%.6e  界=[%.6e,%.6e]' % v)
    merged = [q[0] for q in rows if any(np.isfinite(v) and v < 1e-9 for v in q[1])]
    print(r'  台阶合并（有限且 $<10^{-9}$）出现在：' + (', '.join('$\\theta=%.1f$' % v for v in merged) if merged else '这 25 个角度上一个都没有'))
    zero = [(q[0], q[5]) for q in rows if any(v == 0 for v in q[5] if v is not None)]
    print(r'  $\Theta^{1/2}\Delta P\Theta^{1/2}$ 出现零秩的点：%s'
          % (' '.join('%.1f°:%s' % (a, b) for a, b in zero) if zero else '无'))

    # ---------- (4) 方向扫描：给 (1) 的 iff 找见证 ----------
    print('\n\n===== (4) 方向扫描：$\\operatorname{rank}\\Theta=1$ 的植物上，是否存在删除方向使台阶恰好合并？ =====')
    for nm, Pl, F in E35.plants():
        if nm.find('PL4 单输入 $r=3$') < 0 and nm.find('PL3 $r=2$') < 0:
            continue
        r = F.shape[0]
        wT, QT = np.linalg.eigh(sym(Pl.Th))
        sup = wT > 1e-9 * wT[-1]
        print('  %s  r=%d  $\\Theta$ 的非零本征值 %d 个（$\\lambda_{\\max}=%.4g$），$\\dim\\ker\\Theta$=%d'
              % (nm.replace('$', ''), r, int(sup.sum()), wT[-1], int((~sup).sum())))
        _pF, P0F, _, _ = floor_pair(Pl, F)
        best, hits = None, []
        rng = np.random.default_rng(3)
        dirs = [q for q in rng.standard_normal((400, r))]
        for k in range(r):
            e = np.zeros(r); e[k] = 1.
            dirs.append(e)
        for d in dirs:
            d = d / np.linalg.norm(d)
            wv, Vv = np.linalg.eigh(np.eye(r) - np.outer(d, d))
            Feff = Vv[:, wv > .5].T @ F
            det_ok, _ = Pl.detectable(Feff)
            if not det_ok:
                continue
            phiA, P0A, _res, _ = floor_pair(Pl, Feff)
            dP = sym(P0A - P0F)
            S = sym(QT[:, sup].T @ dP @ QT[:, sup])
            wS = np.linalg.eigvalsh(S)
            gap = float(np.trace(Pl.Th @ dP))
            rk = int(np.sum(wS > 1e-9 * max(wS[-1], 1e-30))) if wS.size else 0
            if best is None or gap < best[0]:
                best = (gap, d, rk, float(wS[-1]) if wS.size else 0.)
            if rk == 0:
                hits.append((gap, d))
        print('    %d 个方向（400 随机 + %d 坐标）里最小的 $\\mathrm{gap}=%.3e$'
              r'（$\\Theta$-受限秩 %d，$\|\Theta^{1/2}\Delta P\Theta^{1/2}\|_{\max}=%.3e$）'
              % (len(dirs), r, best[0], best[2], best[3]))
        print('    合并（$\\Theta$-受限秩 0）的见证：%s'
              % ('、'.join('gap=%.3e, $d$=%s' % (h[0], np.array2string(h[1], precision=2))
                           for h in hits[:4]) if hits else '没有——这些方向上 $\\Delta P$ 从不整块落进 $\\ker\\Theta$'))

