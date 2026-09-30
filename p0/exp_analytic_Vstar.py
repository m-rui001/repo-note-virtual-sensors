"""E29：用 E28 测出的 $C$ 直接构造 Prop E 的解析最优 $V^\\ast=(x/r)C^{-1}$，
在**同一套常数**（$\\det C$、$\\det(F\\tilde P_0F^\\top)$、$I_{unc}(地板)$）下检验整条渐近式
    $\\Delta(x)=\\tfrac r2\\log_2(1/x)+b_{pred}+o(1)$
的斜率与截距，并拿解析点去修复 E23 在 $x<0.03$ 崩掉的数值最优链。

三件事：
 (1) 一维标定 $\\alpha$ 使真实代价恰好落在 $x$：$V=\\alpha C^{-1}$，$J(V)=D_{min}(\\mathcal K_F)+x$；
     用恒等式 $I_{dir}=\\tfrac12\\log_2\\det(F\\tilde PF^\\top+V)-\\tfrac12\\log_2\\det V$ 精确算率（避开大 $\\Gamma$ 求逆）。
 (2) 残差表 $\\Delta_{exact}(x)-[\\tfrac r2\\log_2(1/x)+b_{pred}]$：应当随 $x\\to0$ 单调塌到 0，
     这才说明 E23 全局拟合的 $b$ 吸收的是 $O(x)$ 项，而不是我的截距错了。
 (3) 把解析点与 E23 的 $\\Delta_{opt}$ 逐点比大小：解析点更低 ⟹ E23 深处的链不是最优（它退化到各向同性），
     顺便给 $x\\in\\{0.02,0.01,0.005,0.003,0.002,0.001\\}$ 重测真最优。
"""
import sys
import numpy as np
sys.stdout.reconfigure(encoding='utf-8')
from p0.p0_replicate_letter import A, W, n, sym, ln2, ctrl, unconstrained_sdp
from p0.exp_authoritative import noiseless
from p0.exp_nearfloor_law import Fs, floor_of, Th, Pc, solve_cone, gamma_iso, pack, unpack

UNC = float(np.trace(W @ Pc))
P0 = [None]
E23_F2 = {0.03: 7.0850, 0.05: 6.3521, 0.1: 5.3623, 0.2: 4.3827, 0.377: 3.5038,
          0.8: 2.5015, 1.4: 1.8074, 2.0: 1.4012,
          0.02: 8.3036, 0.01: 9.4923, 0.005: 10.5677, 0.003: 11.3637, 0.002: 11.9701}


def filter_V(F, V, iters=30000, tol=1e-16):
    P = P0[0][0].copy()
    for k in range(iters):
        Pt = sym(A @ P @ A.T + W)
        Pn = sym(Pt - Pt @ F.T @ np.linalg.solve(F @ Pt @ F.T + V, F @ Pt))
        if k > 300 and np.max(np.abs(Pn - P)) < tol:
            P = Pn
            break
        P = Pn
    return P, sym(A @ P @ A.T + W)


def I_of(F, V, Pt=None):
    if Pt is None:
        _, Pt = filter_V(F, V)
    return 0.5 * np.log(np.linalg.det(F @ Pt @ F.T + V) / np.linalg.det(V)) / ln2


def calibrate(F, C, x_target, floor):
    """求 $\\alpha$ 使 $V=\\alpha C^{-1}$ 的真实代价正好超地板 $x_target$。"""
    Ci = np.linalg.inv(C)

    def g(a):
        P, Pt = filter_V(F, a * Ci)
        return float(np.trace(Th @ (P - P0[0][0]))) - x_target, Pt

    lo, hi = 1e-14, 1e3
    for _ in range(55):
        mid = np.sqrt(lo * hi)
        d, _ = g(mid)
        if d > 0:
            hi = mid
        else:
            lo = mid
    a = np.sqrt(lo * hi)
    d, Pt = g(a)
    V = sym(a * Ci)
    return V, a, I_of(F, V, Pt), x_target + d


if __name__ == '__main__':
    for tag, F in Fs.items():
        r = F.shape[0]
        fl, phi, _ = floor_of(F)
        Pn, res, it = noiseless(F, A, W, iters=80000, tol=1e-15)
        P0[0] = (Pn, sym(A @ Pn @ A.T + W))
        M0 = sym(F @ P0[0][1] @ F.T)
        C = np.load('p0/fig/Cmat_%s.npy' % tag.split()[0])
        Iu_floor = unconstrained_sdp(fl, np.eye(n), Th, Pc)['I_true']
        a_pred = (r / 2.) * np.log2(10.)
        b_pred = (r / 2) * np.log2(r) + .5 * np.log2(np.linalg.det(C) * np.linalg.det(M0)) - Iu_floor
        print('\n===== %s   地板=%.4f   斜率系数 $a$=%.4f   截距 $b_{pred}$=%.4f ====='
              % (tag, fl, a_pred, b_pred))
        print('%9s %11s %11s %11s %11s %11s'
              % ('x', 'Δ解析(V*)', 'Δ渐近式', '残差', 'Δ_E23', 'E29−E23'))
        for x in (2.0, 1.0, 0.5, 0.2, 0.1, 0.05, 0.03, 0.02, 0.01, 0.005, 0.002, 0.001):
            V, al, Iex, xreal = calibrate(F, C, x, fl)
            Dasym = a_pred * np.log10(1. / xreal) + b_pred
            Iu = unconstrained_sdp(fl + xreal, np.eye(n), Th, Pc)['I_true']
            de23 = E23_F2.get(round(x, 4), np.nan) if r == 2 else np.nan
            print('%9.4f %11.4f %11.4f %+11.4f %11s %11s'
                  % (xreal, Iex - Iu, Dasym, (Iex - Iu) - Dasym,
                     '%.4f' % de23 if np.isfinite(de23) else '—',
                     '%+.4f' % (Iex - Iu - de23) if np.isfinite(de23) else '—'))
            if r == 2 and x <= 0.02:
                G = np.linalg.inv(V)
                rr = solve_cone(F, fl + x, fl, [pack(G)])
                if rr:
                    print('           以解析点为唯一起点的 $\\\\min_{\\mathcal K_F}I$：%.4f（残差 %+.4f）'
                          % (rr[0] - Iu, (rr[0] - Iu) - Dasym))
        print('  解析 $\\Gamma=V^{*-1}$ 的本征值比 vs 预言 $c_{max}/c_{min}$=%.2f'
              % (np.linalg.eigvalsh(C)[-1] / np.linalg.eigvalsh(C)[0]))
