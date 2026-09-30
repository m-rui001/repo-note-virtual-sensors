r"""E37 v2：垂直段不是一个斜率，是一副**阶梯**。$V\to0$ 的渐近展开要求**每个被使用的通道**自己的方差都 $\ll1$，
所以通道是一个一个进入对数发散的，斜率按 $1/2$ 递增，只有严格极限才到 $\operatorname{rank}C/2$。

E35 的两个读数都在这套机理下被解释掉：它的窗口是 $x\in[10^{-4},10^{-2}]$，而
$C$ 的条件数在出样算例里是 $10^{2}\!\sim\!10^{7}$，窗口里只有强通道进入了渐近段。
（v1 用 L-BFGS-B 找最优 $V$，数值梯度被 $x$ 的 $10^{-11}$  cancellation 噪声淹没，$\mu\gtrsim10^3$ 就停住；
这一版**不用优化器**：候选主动集是**降阶传感器**，一维二分即可精确求出。）

三件事：
 (1) 地板的闭式。$\Phi_0=\operatorname{tr}(\Theta P_0)$，$P_0$ 是 $V\to0$ 的极限后验
     （奇异 DARE $\tilde P=\mathrm{ARE}(A^\top,F^\top,W,0)$ 加 Woodbury）。三方对照：闭式 vs
     E35 的两点外推 vs $\Phi(t)-t\operatorname{tr}C$。若一致，$x$ 的绝对误差从 $10^{-11}$ 降到 $10^{-15}$，
     阶梯才能被读到 $x\sim10^{-12}$。
 (2) 主动集枚举。在 $C$ 的本征基下取 $U_A$，$F_{\rm eff}=U_A^\top F$（丢掉的通道 $v_i=\infty$ 即不传输），
     射线方向 $v_i=t/\lambda_i$。每个族的地板 $\Phi_0(A)\ge\Phi_0(F)$，差就是**丢通道的地板代价**
     $\mathrm{gap}_A$；族的自身坐标 $x_A=x-\mathrm{gap}_A$，渐近斜率应当是 $|A|/2$。
     同一预算 $x$ 下取 $\min_A$，得到真正的上界曲线——阶梯、平台、对数发散。
 (3) 水线填充判据。沿全支撑射线精确算 $\mathrm dI/\mathrm dv_i$（含 $\tilde P$ 对 $V$ 的依赖，
     用 E36 的 $Z\mapsto\mathrm dP$ 表示矩阵）和 $\mathrm dx/\mathrm dv_i=C_{ii}(V)$。
     内点最优的**充要**条件是所有比值 $R_i=(\mathrm dI/\mathrm dv_i)/(\mathrm dx/\mathrm dv_i)$ 相等。
     若它们在 $v_i\ll1$ 时才收敛、在 E35 窗口里分开——"$C^{-1}$ 方向渐近最优、有限 $x$ 不最优"就是被证明而非被假设。
"""
import sys
import numpy as np
sys.stdout.reconfigure(encoding='utf-8')
from scipy.linalg import solve_discrete_are
import p0.exp_plants_oos as E35
import p0.exp_C_exact as X36

ln2 = np.log(2.0)
sym = X36.sym
TS_C = (3e-3, 1e-3, 3e-4, 1e-4, 3e-5)


def pv(Pl, F, V, phi0, polish=600):
    """先验路线：$\tilde P=\mathrm{ARE}(A^\top,F^\top,W,V)$（scipy）+ 定点抛光，后验走 Woodbury。
    为什么不用信息型迭代：$\Pi=\tilde P^{-1}=F^\top V^{-1}F+O(1)$，其 $O(1)$ 部分只能分辨到
    $\epsilon\|F^\top V^{-1}F\|$，于是 $x=\operatorname{tr}(\Theta P)-\Phi_0$（真值是 $O(t)$）在
    $t\lesssim10^{-6}$ 就被彻底吃掉——锚点实测 $x/(\rho t)$ 在 $t=10^{-7}$ 就变成 $-4.4$。
    先验的两条 $O(1)$ 矩阵相减，绝对误差只有 $\epsilon\|\Theta\|\,\|\tilde P\|$。"""
    Pt = Pl.Pt(F, V)
    for k in range(polish):
        Pn = sym(Pt - Pt @ F.T @ np.linalg.solve(F @ Pt @ F.T + V, F @ Pt))
        new = sym(Pl.A @ Pn @ Pl.A.T + Pl.W)
        d = np.linalg.norm(new - Pt)
        Pt = new
        if d <= 1e-16 * max(np.linalg.norm(Pt), 1e-30):
            break
    Pn = sym(Pt - Pt @ F.T @ np.linalg.solve(F @ Pt @ F.T + V, F @ Pt))
    Sz = F @ Pt @ F.T + V
    _, l1 = np.linalg.slogdet(Sz)
    _, l2 = np.linalg.slogdet(V)
    res = np.linalg.norm(Pt - (Pl.A @ Pn @ Pl.A.T + Pl.W)) / np.linalg.norm(Pt)
    return (float(np.trace(Pl.Th @ Pn)) - phi0, .5 * (l1 - l2) / ln2, Pn, Pt, res)


def info(Pl, F, V, phi0, cap=6000):
    """信息型迭代（只用于交叉检验和 E36 的 $C$）：$\Pi=(A\Pi^{-1}A^\top+W)^{-1}+F^\top V^{-1}F$，
    上解 $\Pi_0=F^\top V^{-1}F+W^{-1}$ 单调下降收敛。"""
    A, W = Pl.A, Pl.W
    Vi = np.linalg.inv(V)
    RV = F.T @ Vi @ F
    Pi = RV + np.linalg.inv(W)
    for k in range(cap):
        Pn = np.linalg.inv(Pi)
        M = np.linalg.inv(A @ Pn @ A.T + W)
        new = sym(M + RV)
        d = np.linalg.norm(new - Pi)
        Pi = new
        if k > 60 and d <= 1e-16 * max(np.linalg.norm(M), 1.):
            break
    P = sym(np.linalg.inv(Pi))
    S = sym(A @ P @ A.T + W)
    rres = np.linalg.norm(Pi - np.linalg.inv(S) - RV) / np.linalg.norm(Pi)
    if rres > 1e-11:
        raise AssertionError('info-DARE 残差 %.2e (V=%.0e)' % (rres, V[0, 0]))
    Sz = F @ S @ F.T
    _, l1 = np.linalg.slogdet(Sz + V)
    _, l2 = np.linalg.slogdet(V)
    return float(np.trace(Pl.Th @ P)) - phi0, .5 * (l1 - l2) / ln2, P, S
def floor_exact(Pl, F, phiC):
    """$\\Phi_0=\\operatorname{tr}(\Theta P_0)$，$P_0$ 为 $V\\to0$ 极限后验；三方对照。"""
    r = F.shape[0]
    out = {}
    try:
        Pt0 = sym(solve_discrete_are(Pl.A.T, F.T, Pl.W, np.zeros((r, r))))
        raw = float(np.trace(Pl.Th @ (Pt0 - Pt0 @ F.T @ np.linalg.solve(F @ Pt0 @ F.T, F @ Pt0))))
        Pt0, Pk, res, k = polish_Pt0(Pl, F, Pt0)
        out['闭式'] = float(np.trace(Pl.Th @ Pk))
        out['raw'], out['res'], out['iters'] = raw, res, k
    except Exception as e:
        out['闭式'] = None
        out['闭式err'] = repr(e)
    return out


def polish_Pt0(Pl, F, Pt0, cap=4000):
    """奇异 DARE 解的定点抛光：scipy 在 $R=0$ 时给残差 $\sim10^{-9}$ 的解，对 $\Phi_0$ 是致命的
    ——锚点 $\Phi_0=14.64$，绝对误差 $10^{-8}$ 就把 $x$ 的可用下限钉死在 $10^{-7}$。"""
    for k in range(cap):
        Pk = sym(Pt0 - Pt0 @ F.T @ np.linalg.solve(F @ Pt0 @ F.T, F @ Pt0))
        new = sym(Pl.A @ Pk @ Pl.A.T + Pl.W)
        d = np.linalg.norm(new - Pt0)
        Pt0 = new
        if k > 5 and d <= 3e-17 * max(np.linalg.norm(Pt0), 1e-30):
            break
    Pk = sym(Pt0 - Pt0 @ F.T @ np.linalg.solve(F @ Pt0 @ F.T, F @ Pt0))
    res = np.linalg.norm(Pt0 - (Pl.A @ Pk @ Pl.A.T + Pl.W)) / np.linalg.norm(Pt0)
    return Pt0, Pk, res, k


def C_of(Pl, F):
    seq = X36.C_seq(Pl, F, TS_C)
    best = None
    for p in (1, 2):
        C0, a1, a2 = X36.extrap(seq, p)
        d = np.linalg.norm(a1 - a2) / max(np.linalg.norm(C0), 1e-30)
        if best is None or d < best[0]:
            best = (d, C0)
    return sym(best[1])


def fams_of(Pl, F, C, phi0F):
    """$C$ 的本征基 + 全部非空主动集。返回 (族列表, 本征值, $U$)。"""
    rho, e, U = X36.rank_C(C)
    out = []
    for bits in range(1, 1 << rho):
        A = [i for i in range(rho) if bits >> i & 1]
        if len(A) > 3:
            continue
        Feff = U[:, A].T @ F
        fx = floor_exact(Pl, Feff, None)
        phi0A = fx['闭式']
        if phi0A is None:                       # 降采样后 $(A,F_{\rm eff})$ 不可检测：地板不存在
            print('    族 %-8s 不可检测（奇异 DARE 无解），跳过：%s' % ('+'.join(map(str, A)), fx['闭式err']))
            continue
        w = e[A]                                # 射线方向 $v_i=t/\lambda_i$
        nrm = np.linalg.norm(Pl.Th) * np.linalg.norm(Pl.Pt(Feff, 1e-6 * np.eye(len(A))))
        out.append(dict(A=A, a=len(A), F=Feff, phi0=phi0A, gap=max(phi0A - phi0F, 0.), w=w,
                        res=6e-17 * nrm, name='+'.join(str(i) for i in A)))
    return out, e, U, rho


def curve(Pl, fm, ts, phi0):
    """族曲线：$x$(相对全地板) 与 $I$。"""
    rows = []
    for t in ts:
        V = t * np.diag(1.0 / fm['w'])
        try:
            xF, I, _P, _Pt, rr = pv(Pl, fm['F'], V, phi0)
        except Exception:
            continue
        rows.append((t, xF, I))
    return rows


def slope(ts_logx, Is):
    """$\\mathrm dI/\\mathrm d\\log_2(1/x)$：逐点中心差分（$x$ 对 $t$ 单调）。"""
    u = np.log2(1.0 / np.asarray(ts_logx))
    I = np.asarray(Is)
    return (I[2:] - I[:-2]) / (u[2:] - u[:-2]), u, I


def dP_op_p(Pl, F, V):
    r"""先验路线的 $Z\mapsto dP$，把 $t$ 的幂次全部吸收进 $P/t$。

    E36 的 $dP\_op$ 有两处量级陷阱：(i) $\Pi=P^{-1}$ 靠 4000 步定点迭代，
    $\|\Pi\|\sim t^{-1}$ 时残差判据 $\|\Pi-M-F^	op V^{-1}F\|/\|\Pi\|<10^{-10}$ 只约束到绝对误差 $\sim1$，
    而 $M$ 本身就是 $O(1)$——锚点（$A$ 有 $|\lambda|=1.71$ 的不稳定块，收敛因子 $
ho(B)^2	o1$）
    的 $M$ 因此是错的，$\operatorname{tr}(\Theta\,dP)$ 在 $t\le10^{-7}$ 直接给出 $-5	imes10^{2}$，
    而中心差分给正确的 $25.9$。(ii) $P\otimes P=O(t^2)$ 与 $J_G=O(t^{-2})$ 相乘，
    $10^{-20}	imes10^{20}$ 的绝对舍入误差被放大。

    修法：先验 $S=AP_pA^	op+W$ 由 pv 给出，$M=S^{-1}$ 不必迭代；写 $V=tW_0$，
    $F^	op V^{-1}=t^{-1}F^	op W_0^{-1}$、$P=t\,(P/t)$，于是
    $J=-(P/t\otimes P/t)\,(I-B\otimes B)^{-1}\,(F^	op W_0^{-1}\otimes F^	op W_0^{-1})$，全程 $O(1)$。
    """
    A, W, n = Pl.A, Pl.W, Pl.nn
    r = F.shape[0]
    t = float(np.trace(V)) / r
    W0 = V / t
    Fv0 = F.T @ np.linalg.inv(W0)
    _x, _I, Pp, Pt_, _res = pv(Pl, F, V, 0.0)
    M = np.linalg.inv(Pt_)
    Bm = M @ A @ Pp
    Pi = M + (Fv0 @ F) / t
    dev = float(np.linalg.norm(np.linalg.inv(Pi) - Pp) / np.linalg.norm(Pp))
    J = (np.kron(Pp / t, Pp / t) @ np.linalg.solve(np.eye(n * n) - np.kron(Bm, Bm),
                                                       np.kron(Fv0, Fv0)))
    return J, Pp, Bm, dev


def grads(Pl, F, V):
    """在对称基 $\{E_{ii},\,E_{ij}+E_{ji}\}$ 上给出 $(\mathrm dI,\mathrm dx)$，用于水线填充 KKT。"""
    rho = F.shape[0]
    J, P, _B, dev = dP_op_p(Pl, F, V)
    _x, _I, _Pn, S, _r = pv(Pl, F, V, 0.0)     # $S=	ilde P$（先验）；走先验路线，$info$ 在 $t\le10^{-7}$ 会崩
    K = np.linalg.inv(F @ S @ F.T + V)
    G = F.T @ K @ F
    Vi = np.linalg.inv(V)
    out = []
    for i in range(rho):
        for j in range(i, rho):
            Z = np.zeros((rho, rho))
            if i == j:
                Z[i, i] = 1.
            else:
                Z[i, j] = Z[j, i] = 1.
            dP = (J @ Z.reshape(-1)).reshape(P.shape)
            dS = Pl.A @ dP @ Pl.A.T
            dI = .5 / ln2 * (np.trace(G @ dS) + np.trace((K - Vi) @ Z))
            dx = float(np.trace(Pl.Th @ dP))
            out.append(((i, j), dI, dx))
    return out, dev


def dx_fd(Pl, F, V, phi0, h=2e-3):
    r"""dx/dV 的中心差分（对称基），用来顶替 Fréchet 算子。

    动机见 dP_op_p：dx = tr(Theta dP) 的解析式在坐标对齐的 plant 上到 t=1e-9 都对，
    但锚点和随机正交 F 在 t<=1e-7 给出 -1e2（真值 +26）——Stein 求逆
    (I-B kron B)^{-1} 把机器 epsilon 放大到与 dx 同量级。差分只有
    1e-12/(2h v_i) 的舍入地板，深 t 下反而更可信。两条 mu 并列打印，读者自己看裂口在哪。
    """
    rho = F.shape[0]
    out = []
    for i in range(rho):
        for j in range(i, rho):
            if i == j:
                d = np.zeros_like(V)
                d[i, i] = h * V[i, i]
                step = 2 * d[i, i]
            else:
                d = np.zeros_like(V)
                d[i, j] = d[j, i] = h * np.sqrt(V[i, i] * V[j, j])
                step = 2 * d[i, j]
            out.append(((i, j), (pv(Pl, F, V + d, phi0)[0] - pv(Pl, F, V - d, phi0)[0]) / step))
    return out


def report(name, Pl, F, targets=(2., 1., .5, 1e-1, 1e-2, 1e-3, 1e-5, 1e-7, 1e-9)):
    r = F.shape[0]
    print('\n' + '=' * 90 + '\n%s   n=%d r=%d  rank$\\Theta$=%d' %
          (name, Pl.nn, r, E35.rank_psd(Pl.Th)[0]))
    phi_ext = E35.floor_of(Pl, F)
    C = C_of(Pl, F)
    phi_sub = float(np.trace(Pl.Th @ Pl.P(F, 1e-7 * np.eye(r)))) - 1e-7 * np.trace(C)
    fx = floor_exact(Pl, F, C)
    print(r'  地板 $V\to0$：闭式(抛光) $\Phi_0=%s$  scipy 原始 $=%.12f$  两点外推(E35) $=%.12f$  '
          r'$\Phi(t)-t\operatorname{tr}C$ $=%.12f$'
          % (('%.12f' % fx['闭式']) if fx['闭式'] is not None else repr(fx.get('闭式err')),
             fx.get('raw', float('nan')), phi_ext[1], phi_sub))
    if fx['闭式'] is not None:
        print(r'  ↑ 抛光后奇异 DARE 残差 %.1e（%d 步）；闭式 vs 外推 %.2e，vs 减法 %.2e'
              % (fx['res'], fx['iters'], abs(fx['闭式'] - phi_ext[1]), abs(fx['闭式'] - phi_sub)))
    phi0 = fx['闭式'] if fx['闭式'] is not None else phi_ext[1]
    rho, elam, U = X36.rank_C(C)
    print(r'  $\operatorname{rank}C=%d$  $\lambda(C)=%s$  条件数 $%.2e$  旧猜测阈值 $\rho\lambda=%s$'
          % (rho, np.array2string(elam, precision=7), elam[0] / max(elam[rho - 1], 1e-30),
             np.array2string(rho * elam[:rho], precision=3)))
    chk = []
    Ci = U[:, :rho] @ np.diag(1.0 / elam[:rho]) @ U[:, :rho].T
    for t in (1e-3, 1e-5, 1e-7, 1e-9):
        try:
            xF = pv(Pl, F, t * Ci, phi0)[0]
            chk.append(r'$t=%.0e$:$x/(\rho t)=%.6f$' % (t, xF / (rho * t)))
        except Exception as e:
            chk.append('$t=%.0e$:$\times$%s' % (t, repr(e)[:20]))
    print(r'  地板自检（沿射线 $x$ 应当是 $\rho t$）：' + '  '.join(chk))

    fm, ee, uu, _ = fams_of(Pl, F, C, phi0)
    ffull = [f for f in fm if f['a'] == rho]
    # (3) 水线填充判据：$V=tC^{-1}$（在 $C$ 的本征基下对角）一阶条件
    if ffull:
        f0 = ffull[0]
        print(); print('  ' + r'[KKT] 射线 $v_i=t/\lambda_i$：内点最优的一阶条件是矩阵等式 $\nabla I=\mu\,\nabla x$。'
              r'$\mathrm dI$ 用闭式（对中心差分验过，5 位符合），$\mathrm dx$ 两条路：'
              r'Fréchet/Stein 解析式与 $\operatorname{tr}(\Theta P)$ 的差分。'
              r'残差 $(\nabla I-\mu\nabla x)/\nabla I$ 按差分 $\mathrm dx$ 算')
        print('    %8s %10s %13s %13s %13s %8s %-16s %-12s %s' % ('t', 'x', 'mu(解析dx)',
              'mu(差分dx)', '预言mu', '有效$|A|$', '对角残差', '非对角残差', '$v_i$'))
        for t in (1e-1, 1e-3, 1e-5, 1e-7, 1e-9):
            V = t * np.diag(1.0 / f0['w'])
            try:
                pvr = pv(Pl, f0['F'], V, phi0)
                xF, pvres = pvr[0], pvr[4]
                gl, dev = grads(Pl, f0['F'], V)
            except (AssertionError, np.linalg.LinAlgError):
                continue
            dd = [(dI, dx) for (ij, dI, dx) in gl if ij[0] == ij[1]]
            od = [(dI, dx) for (ij, dI, dx) in gl if ij[0] != ij[1]]
            mu = sum(dI for dI, _ in dd) / sum(dx for _, dx in dd)
            xmap = dict(dx_fd(Pl, f0['F'], V, phi0))
            ddf = [(dI, xmap[ij]) for ij, dI, _ in gl if ij[0] == ij[1]]
            od = [(dI, xmap[ij]) for ij, dI, _ in gl if ij[0] != ij[1]]
            mufd = sum(dI for dI, _ in ddf) / sum(e for _x, e in ddf)
            pred = -(f0['a'] / 2.) / (xF * ln2)
            rd = [(dI - mufd * e) / dI for dI, e in ddf]
            gn = np.sqrt(sum(dI ** 2 for dI, _ in ddf) + sum(dI ** 2 for dI, _ in od))
            ro = max([abs(dI - mufd * e) / gn for dI, e in od], default=0.)
            print('    %8.1e %10.2e %13.4e %13.4e %13.4e %8.3f %-16s %-12s %s'
                  % (t, xF, mu, mufd, pred, 2 * abs(mufd) * xF * ln2,
                     np.array2string(np.array(rd), precision=2, max_line_width=200),
                     '%.2e' % ro, ' '.join('%.1e' % q for q in np.diag(V))))
    # (2) 各族渐近斜率
    print(); print('  ' + r'[族斜率] 自身坐标 $x_A=x-\mathrm{gap}_A$，预言 $\mathrm dI/\mathrm d\log_2(1/x_A)=|A|/2$')
    ts = 10.0 ** np.linspace(-0.7, -13.0, 26)
    for f in fm:
        rows = curve(Pl, f, ts, phi0)
        xs = np.array([q[1] - f['gap'] for q in rows])
        Is = np.array([q[2] for q in rows])
        res = f['res']
        ok = xs > 30 * res
        if ok.sum() < 6:
            print('    %-8s 可分辨点不足（$x_A$ 全部 $\\le$ 地板噪声 %.1e）' % (f['name'], res))
            continue
        sl, u, Ia = slope(xs[ok], Is[ok])
        pick = []
        for want_u in (6.0, 14.0, 22.0):
            if len(sl) == 0:
                break
            k = int(min(np.argmin(np.abs(u[:-1] - want_u)), len(sl) - 1))
            pick.append('$x_A=%.0e$:%.4f' % (xs[ok][k], sl[k]))
        print('    %-8s $|A|=%d$  $\mathrm{gap}=%+.4e  预 %.1f   %s'
              % (f['name'], f['a'], f['gap'], f['a'] / 2., '   '.join(pick)))
    # (2b) 同预算竞赛：每族一维二分（$x_A$ 随 $t$ 单调增）
    print(r'\n  [同预算 $x$ 上的竞赛] 单元格 $=I_{\min}$ 在该族的射线上的值（$I_{\text{opt}}\le\min_A$）')
    print('    %-8s %10s   %s' % ('族', 'gap', ' '.join('$x=%7.0e$' % X for X in targets)))
    grid = {}
    for f in fm:
        cells = []
        for X in targets:
            xA = X - f['gap']
            if xA <= 30 * f['res']:
                cells.append(None); continue
            lo, hi = -15.0, 3.0
            bad = False
            for _ in range(48):
                mid = .5 * (lo + hi)
                try:
                    xF = pv(Pl, f['F'], 10 ** mid * np.diag(1.0 / f['w']), phi0)[0]
                except Exception:
                    bad = True; break
                if xF - f['gap'] < xA:                  # $x$ 随 $t$ 增：太小就把下界抬高
                    lo = mid
                else:
                    hi = mid
            if bad:
                cells.append(None); continue
            try:
                xg, I, _P, _Pt, rr = pv(Pl, f['F'], 10 ** (.5 * (lo + hi)) * np.diag(1.0 / f['w']), phi0)
                cells.append(I if abs(xg - X) < .05 * X + 1e-14 else None)
            except Exception:
                cells.append(None)
        grid[f['name']] = cells
        print('    %-8s %10.3e   %s' % (f['name'], f['gap'],
              ' '.join('    ---   ' if c is None else '$%9.4f$' % c for c in cells)))
    a_of = {f['name']: f['a'] for f in fm}
    bi = []
    for j, X in enumerate(targets):
        cand = [(c, k) for k, cs in grid.items() for c in [cs[j]] if c is not None]
        if cand:
            c, k = min(cand)
            bi.append((X, c, a_of[k], k))
            print('    最优 @ $x=%.0e$：族 %-8s $|A|=%d$  $I=%.4f$' % (X, k, a_of[k], c))
        else:
            bi.append((X, None, None, None))
    seg = []
    for j in range(len(bi) - 1):
        b, n = bi[j], bi[j + 1]
        if b[1] is not None and n[1] is not None:
            seg.append((b[0], n[0], (n[1] - b[1]) / np.log2(b[0] / n[0]), b[2], n[2]))
    print('  逐段斜率 $\mathrm dI/\mathrm d\log_2(1/x)$：' +
          '   '.join('$[%.0e,%.0e]$%.3f ($|A|$%s→%s)' % s for s in seg))
    return phi0, C, elam, rho


if __name__ == '__main__':
    from p0.exp_plants_oos import plants
    want = ('T-F2 尾部2（锚点）', 'PL1 $r=2$ 坐标', 'PL2 $r=3$ 坐标', 'PL3 $r=2$',
            'PL4 单输入 $r=2$', 'PL4 单输入 $r=3$')
    for nm, pl, F in plants():
        if nm not in want:
            continue
        try:
            report(nm.replace('$', ''), pl, F)
        except Exception as ex:
            import traceback
            print('\n===== %s 失败：%r' % (nm, ex))
            traceback.print_exc()
