"""E30：出样检验——把 E28/E29 的闭式渐近式用到从未调过的 $F$ 上。

$$\\Delta(x)=\\tfrac r2\\log_2\\tfrac1x+b_{pred},\\qquad
b_{pred}=\\tfrac r2\\log_2 r+\\tfrac12\\log_2\\bigl(\\det C\\cdot\\det(F\\tilde P_0F^\\top)\\bigr)-I_{unc}(D_{min}(\\mathcal K_F))$$

$F_2,F_3$ 是沿着 Schur 尾部挑的（我的整个推导都在那个基里做），所以它们过拟合的风险是真实的。
这里换三个完全不相干的 $F$：坐标选择 $\\{x_1,x_2\\}$、$\\{x_2,x_4\\}$，和一个随机正交行对。
每个都独立走一遍：地板（$\\Phi=\\operatorname{tr}(\\Theta P_0)$）$\\to C$（秩一 PSD 差商）$\\to b_{pred}\\to$
数值锥上最优（起点里放解析 $V^\\ast=(x/r)C^{-1}$）。要检验的是**残差随 $x\\to0$ 塌向 0**，
不是"$b$ 全局拟合得好"——$x$ 大时 $O(x)$ 项当然还在。

顺带把 $r=1$ 的可证伪预言算出来（系数 $0.5\\log_2 10=1.6609$ bit/十倍程），这是原文没有的格。
"""
import sys
import numpy as np
sys.stdout.reconfigure(encoding='utf-8')
from p0.p0_replicate_letter import A, W, n, sym, ln2, ctrl, unconstrained_sdp
from p0.exp_authoritative import noiseless
from p0.exp_nearfloor_law import Th, Pc, solve_cone, gamma_iso, pack, unpack
import p0.exp_cmatrix_intercept as e28

UNC = float(np.trace(W @ Pc))
vs = (1e-3, 1e-4, 1e-5, 1e-6)


def measure_C(F):
    r = F.shape[0]
    C = np.zeros((r, r))
    for i in range(r):
        e = np.zeros((r, 1)); e[i] = 1.
        C[i, i] = e28.directional(F, e @ e.T, vs)[0]
    for i in range(r):
        for j in range(i + 1, r):
            ei = np.zeros((r, 1)); ei[i] = 1.
            ej = np.zeros((r, 1)); ej[j] = 1.
            tp = e28.directional(F, (ei + ej) @ (ei + ej).T, vs)[0]
            tm = e28.directional(F, (ei - ej) @ (ei - ej).T, vs)[0]
            C[i, j] = C[j, i] = (tp - tm) / 4.
    return C


def case(name, F):
    r = F.shape[0]
    Pn, res, it = noiseless(F, A, W, iters=80000, tol=1e-15)
    e28.P0[0] = (Pn, sym(A @ Pn @ A.T + W))
    Pt0 = e28.P0[0][1]
    M0 = sym(F @ Pt0 @ F.T)
    evM0 = np.linalg.eigvalsh(M0)
    fl = UNC + float(np.trace(Th @ Pn))
    Iu_floor = unconstrained_sdp(fl, np.eye(n), Th, Pc)['I_true']
    C = measure_C(F)
    evC = np.linalg.eigvalsh(C)
    ok = bool(evM0[0] > 1e-10 and evC[0] > 1e-10)
    b_pred = (r / 2) * np.log2(r) + .5 * np.log2(np.linalg.det(C) * np.linalg.det(M0)) - Iu_floor
    a_pred = (r / 2.) * np.log2(10.)
    print('\n===== %s：rank=%d  可检测=%s  $\\Phi$=%.4f  地板=%.4f  $C\\succ0$=%s ====='
          % (name, r, ok, fl - UNC, fl, evC[0] > 1e-10))
    print('  $C$ 本征值 %s   $\\det C$=%.4f   $\\det(F\\tilde P_0F^\\top)$=%.4f'
          % (np.array2string(evC, precision=4), np.linalg.det(C), np.linalg.det(M0)))
    print('  预言 $\\Gamma^\\ast$ 本征值比 %.2f   $a$=%.4f  $b_{pred}$=%.4f'
          % (evC[-1] / evC[0], a_pred, b_pred))
    if not ok:
        print('  不可检测/奇异 $C$，渐近式不适用，跳过')
        return
    print('%9s %12s %12s %12s %12s %12s' % ('x', 'Δ数值最优', 'Δ渐近式', '残差', 'Δ解析V*', 'Γ*比'))
    Ci = np.linalg.inv(C)
    for x in (0.5, 0.2, 0.1, 0.05, 0.02, 0.01):
        V = (x / r) * Ci
        # 解析点的真实代价（用来把名义 x 换成实际 x）
        P, Pt = e28.filter_V(F, V + 1e-15 * np.eye(r))
        xa = float(np.trace(Th @ (P - Pn)))
        I_ana = 0.5 * np.log(np.linalg.det(F @ Pt @ F.T + V) / np.linalg.det(V)) / ln2
        starts = [pack(sym(V))]
        lg = gamma_iso(F, fl + x)
        if lg is not None:
            starts.append(pack(np.exp(lg) * np.eye(r)))
        rr = solve_cone(F, fl + x, fl, starts)
        Iu = unconstrained_sdp(fl + x, np.eye(n), Th, Pc)['I_true']
        if rr is None:
            print('%9.4f %12s %12s %12s %12.4f' % (x, '失败', a_pred * np.log10(1 / x) + b_pred,
                                                   '—', I_ana - Iu))
            continue
        Dasym = a_pred * np.log10(1 / x) + b_pred
        print('%9.4f %12.4f %12.4f %+12.4f %12.4f %12.2f'
              % (x, rr[0] - Iu, Dasym, (rr[0] - Iu) - Dasym, I_ana - Iu,
                 np.max(np.linalg.eigvalsh(np.linalg.inv(rr[2]))) /
                 np.min(np.linalg.eigvalsh(np.linalg.inv(rr[2])))))


if __name__ == '__main__':
    def coord(i, j):
        M = np.zeros((2, n)); M[0, i] = 1.; M[1, j] = 1.
        return M
    rng = np.random.default_rng(11)
    Q = np.linalg.qr(rng.standard_normal((n, 2)))[0].T
    case('F6 坐标 $\\{x_1,x_2\\}$', coord(0, 1))
    case('F7 坐标 $\\{x_2,x_4\\}$', coord(1, 3))
    case('F8 随机正交行对', Q)
    print('\n===== $r=1$ 的可证伪预言 =====')
    e1 = np.zeros((1, n)); e1[0, 0] = 1.
    case('F9 单通道 $y=x_1$', e1)
