"""E28：把 Prop E 的**截距**也算出来跟 E23 的拟合对表，并测出活动判据需要的 $C$ 矩阵。

推导（每步都是可数值验证的代数，不含未证断言）：
  (a) 行列式引理给恒等式 $\\det P=\\det\\tilde P\\,\\det V/\\det(F\\tilde PF^\\top+V)$，于是
      $I_{dir}=\\tfrac12\\log_2\\det(F\\tilde PF^\\top+V)-\\tfrac12\\log_2\\det V$。$V\\to0$ 时第一项
      $\\to\\tfrac12\\log_2\\det(F\\tilde P_0F^\\top)$（$\\tilde P_0\\succ0$ 由 $(A,F)$ 可检测性保证），只有 $-\\tfrac12\\log_2\\det V$ 发散。
  (b) $x:=J-J_{min}(\\mathcal K_F)=\\operatorname{tr}(CV)+O(|V|^2)$，$C_{ij}:=\\operatorname{tr}(\\Theta\\,\\partial P/\\partial V_{ij})|_{V=0}$。
      标度指数 $=1$（不是 $1/2$）由 E27 实测钉死。
  (c) $\\min-\\tfrac12\\log_2\\det V$ s.t. $\\operatorname{tr}(CV)\\le x$ 的驻点 $-V^{-1}+\\lambda C=0$ 强制
      $V=\\tfrac{x}{r}C^{-1}$（最优 $V$ 与 $C$ 同特征基），$\\lambda=r/x$，$\\det V=(x/r)^r/\\det C$。于是
$$\\Delta(x)=\\tfrac r2\\log_2\\tfrac1x+b_{pred}+o(1),\\qquad
  b_{pred}=\\tfrac r2\\log_2 r+\\tfrac12\\log_2\\bigl(\\det C\\cdot\\det(F\\tilde P_0F^\\top)\\bigr)-I_{unc}(D_{min}(\\mathcal K_F))$$
      且最优 $\\Gamma=V^{-1}$ 的本征值比 $=C$ 的本征值比。$C\\succ0$ 就是"$r$ 个通道全活动"的判据；
      若有 $k$ 个零本征值，$b$ 与斜率系数都要按 $k$ 改写（$C$ 只在活动子空间上取行列式）。

三件事：秩一 PSD 方向差商测 $C$；$\\operatorname{tr}C$ 对 E27 的各向同性极限 $\\lim\\gamma x$ 做**加法性**交叉检验；
$b_{pred}$、$C$ 本征值比对表 E23 的 $b=3.1988/5.4208$ 与实测本征值比 $\\approx18$。
"""
import sys
import numpy as np
sys.stdout.reconfigure(encoding='utf-8')
from p0.p0_replicate_letter import A, W, n, sym, ln2, ctrl, unconstrained_sdp
from p0.exp_authoritative import noiseless
from p0.exp_nearfloor_law import Fs, floor_of, Th, Pc

UNC = float(np.trace(W @ Pc))
P0 = [None]         # 当前 F 的 (后验 P0, 预报 Pt0)；numpy 数组不能当键


def filter_V(F, V, iters=40000, tol=1e-15):
    """稳态滤波，噪声协方差 $V$（可奇异地小）。从 $V=0$ 的解热启动。"""
    P = P0[0][0].copy()
    for k in range(iters):
        Pt = sym(A @ P @ A.T + W)
        Pn = sym(Pt - Pt @ F.T @ np.linalg.solve(F @ Pt @ F.T + V, F @ Pt))
        if k > 300 and np.max(np.abs(Pn - P)) < tol:
            P = Pn
            break
        P = Pn
    return P, sym(A @ P @ A.T + W)


def x_of(F, V):
    """$x=\\operatorname{tr}(\\Theta P)-\\operatorname{tr}(\\Theta P_0)$，即超出地板的那部分代价。"""
    P, _ = filter_V(F, V)
    return float(np.trace(Th @ (P - P0[0][0])))


def directional(F, D, vs=(1e-3, 1e-4, 1e-5, 1e-6)):
    """沿 PSD 方向 $D$ 测 $\\operatorname{tr}(CD)$：$x(vD)=v\\operatorname{tr}(CD)+O(v^2)$，
       对 $(v,v^2)$ 做最小二乘，把 $v^2$ 项剥掉外推到 $v=0$。"""
    xs = np.array([x_of(F, v * D + 1e-13 * np.eye(D.shape[0])) for v in vs])
    M = np.vstack([np.array(vs), np.array(vs) ** 2]).T
    a, b = np.linalg.lstsq(M, xs, rcond=None)[0]
    r = xs - (a * np.array(vs) + b * np.array(vs) ** 2)
    return a, np.max(np.abs(xs - (a * np.array(vs) + b * np.array(vs) ** 2)) / xs), a * 1e-13


if __name__ == '__main__':
    for tag, F in Fs.items():
        r = F.shape[0]
        fl, phi, _ = floor_of(F)
        Pn, res, it = noiseless(F, A, W, iters=80000, tol=1e-15)
        P0[0] = (Pn, sym(A @ Pn @ A.T + W))
        Pt0, M0 = P0[0][1], sym(F @ P0[0][1] @ F.T)
        evM0 = np.linalg.eigvalsh(M0)
        print('\n===== %s  rank=%d  地板=%.4f =====' % (tag, r, fl))
        print('  $F\\tilde P_0F^\\top$ 本征值 %s（全正 ⟺ 可检测，$\det$=%.4f）'
              % (np.array2string(evM0, precision=4), np.linalg.det(M0)))

        C = np.zeros((r, r))
        rep = []
        for i in range(r):
            e = np.zeros((r, 1)); e[i] = 1.
            D = e @ e.T
            a, rel, _ = directional(F, D)
            C[i, i] = a
            rep.append(('C%d%d' % (i + 1, i + 1), a, rel))
        for i in range(r):
            for j in range(i + 1, r):
                ei = np.zeros((r, 1)); ei[i] = 1.
                ej = np.zeros((r, 1)); ej[j] = 1.
                tp, rp, _ = directional(F, (ei + ej) @ (ei + ej).T)
                tm, rm, _ = directional(F, (ei - ej) @ (ei - ej).T)
                C[i, j] = C[j, i] = (tp - tm) / 4.
                rep.append(('C%d%d' % (i + 1, j + 1), (tp - tm) / 4., max(rp, rm)))
        print('  差商测到的 $C$ 元素（相对残差 = $v^2$ 项在最小 $v$ 处的占比）：')
        for k, v, rr in rep:
            print('    %-5s = %10.5f   残差 %.2e' % (k, v, rr))
        ev = np.linalg.eigvalsh(C)
        trC = float(np.trace(C))
        print('  $C$ 本征值 %s   $\\operatorname{tr}C$=%.4f   $\\det C$=%.5f'
              % (np.array2string(ev, precision=4), trC, np.linalg.det(C)))

        print('  加法性交叉检验：各向同性支 $\\lim\\gamma x$ 应 $=\\operatorname{tr}C$')
        for lg in (11., 13.):
            g = np.exp(lg)
            xl = x_of(F, (1. / g) * np.eye(r))
            print('    $\\gamma$=%.2e：$\\gamma x$=%.4f  vs $\\operatorname{tr}C$=%.4f（差 %+.2f%%）'
                  % (g, g * xl, trC, 100 * (g * xl / trC - 1)))

        Iu = unconstrained_sdp(fl, np.eye(n), Th, Pc)['I_true']
        b_pred = (r / 2) * np.log2(r) + .5 * np.log2(np.linalg.det(C) * np.linalg.det(M0)) - Iu
        b_meas = 3.1988 if r == 2 else 5.4208
        print('  截距：$I_{unc}(地板)$=%.4f，$b_{pred}$=%.4f  vs E23 拟合 $b$=%.4f → 差 %+.4f bit'
              % (Iu, b_pred, b_meas, b_pred - b_meas))
        print('  最优 $\\Gamma$ 本征值比：预言 $c_{max}/c_{min}$=%.2f（E23 收敛段实测 17.7–18.6）'
              % (ev[-1] / ev[0] if ev[0] > 1e-12 else np.inf))
        print('  活动判据 $C\\succ0$：%s（最小本征值 %.5f）' % (bool(ev[0] > 1e-9), ev[0]))
        np.save('p0/fig/Cmat_%s.npy' % tag.split()[0], C)
