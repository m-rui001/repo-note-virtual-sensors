r"""E35：出样**换植物**——三项渐近式此前只在 Tanaka 那一个 4 阶算例上验过；顺带正面攻 $C\succ0$。

两件事：
 (a) 精度升级。$C$ 用 Richardson 双点（$v=10^{-7},10^{-6}$）测，地板用 $\Phi(t)=\operatorname{tr}(\Theta P(tI))$
     的常数$+$线性外推，都不再依赖"迭代到定点"。上一轮卡在 $\kappa$ 上的 $\delta=8\times10^{-6}$
     其实就是 $C$ 的相对误差；这次应当掉到 $10^{-9}$。
 (b) 结构断言。单调性白给 $C\succeq0$（噪声变大代价变大），但 $C\succ0$ 是猜想，而且它**等价于**
     垂直斜率是 $r/2$。若 $\Theta$ 有核（单输入时必然，$\operatorname{rank}\Theta\le m$），某些测量方向
     再吵也不影响代价 $\Rightarrow C$ 奇异 $\Rightarrow$ 最优把那些通道的 $V$ 推到 $\infty$（不传输），
     斜率塌成 $\operatorname{rank}(C)/2$，截距里的 $\det C$ 换成非零本征值之积。
     脚本检出降秩时自动收缩到 $C$ 的支集 $F_{\rm eff}=U_1^{\top}F$ 再跑一遍——被检验的是这个加强版。

四个新植物：PL1 复不稳定对（$r=1$ 即可检测，新格）、PL2 $n=5$ 三个不稳定模态、
PL3 临界 Jordan 块（$\sum\log|\lambda_u|=0$，检验垂直式是否依赖"严格不稳定"，并附一个不可检测反例）、
PL4 单输入（$\operatorname{rank}\Theta=1$，主动造降秩）。
"""
import sys
import numpy as np
sys.stdout.reconfigure(encoding='utf-8')
from scipy.linalg import solve_discrete_are
import cvxpy as cp

ln2 = np.log(2.0)
RT = 1e-7


def sym(M):
    return .5 * (M + M.T)


def rank_psd(M):
    e = np.linalg.eigvalsh(sym(M))
    return int(np.sum(e > RT * max(e.max(), 1e-12))), e


class Pl:
    def __init__(self, A, B, W, Q, R=None):
        self.A, self.B, self.W = A, B, sym(W)
        self.Q = sym(Q)
        self.nn = A.shape[0]
        self.R = np.eye(B.shape[1]) if R is None else sym(R)
        self.Pc = solve_discrete_are(A, B, self.Q, self.R)
        K = np.linalg.solve(self.R + B.T @ self.Pc @ B, B.T @ self.Pc @ A)
        self.Th = sym(K.T @ (self.R + B.T @ self.Pc @ B) @ K)
        self.UNC = float(np.trace(self.W @ self.Pc))

    def detectable(self, F):
        for lam in np.linalg.eigvals(self.A):
            if abs(lam) < 1 - 1e-12:
                continue
            M = np.vstack([self.A - lam * np.eye(self.nn), F])
            if np.linalg.matrix_rank(M, tol=1e-8) < self.nn:
                return False, lam
        return True, None

    def Pt(self, F, V):
        return sym(solve_discrete_are(self.A.T, F.T, self.W, sym(V)))

    def P(self, F, V, Pt=None):
        if Pt is None:
            Pt = self.Pt(F, V)
        return sym(Pt - Pt @ F.T @ np.linalg.solve(F @ Pt @ F.T + V, F @ Pt))

    def xI(self, F, V, Phi0):
        Pt = self.Pt(F, V)
        x = float(np.trace(self.Th @ self.P(F, V, Pt))) - Phi0
        _, l1 = np.linalg.slogdet(F @ Pt @ F.T + V)
        _, l2 = np.linalg.slogdet(V)
        return x, .5 * (l1 - l2) / ln2, Pt

    def sdp(self, D):
        """Tanaka 式(18) max-det SDP；返回 recovered $P$ 上的真实率 $\tfrac12\log\det(\tilde P/P)$。"""
        N = self.nn
        P = cp.Variable((N, N), symmetric=True)
        Pi = cp.Variable((N, N), symmetric=True)
        APT = self.A @ P @ self.A.T + self.W
        cons = [Pi >> 0, P >> 0, cp.trace(self.Th @ P) + self.UNC <= D, P << APT,
                cp.bmat([[P - Pi, P @ self.A.T], [self.A @ P, APT]]) >> 0]
        prob = cp.Problem(cp.Minimize(-.5 * cp.log_det(Pi) / ln2
                                      + .5 * np.log(np.linalg.det(self.W)) / ln2), cons)
        prob.solve(solver=cp.CLARABEL, gp=False)
        if P.value is None:
            return None
        Pv = sym(P.value)
        Pt = self.A @ Pv @ self.A.T + self.W
        return float(.5 * np.log(np.linalg.det(Pt) / np.linalg.det(Pv)) / ln2)


def floor_of(Pl, F):
    I = np.eye(F.shape[0])
    f = lambda t: float(np.trace(Pl.Th @ Pl.P(F, t * I)))
    a1, a2 = f(1e-6), f(1e-7)
    phi0 = a2 + (a2 - a1) * 1e-7 / (1e-6 - 1e-7)
    return phi0 + Pl.UNC, phi0


def measure_C(Pl, F, Phi0, t1=1e-7, t2=1e-6):
    r = F.shape[0]

    def lin(D):
        a1 = Pl.xI(F, t1 * D, Phi0)[0] / t1
        a2 = Pl.xI(F, t2 * D, Phi0)[0] / t2
        return (t2 * a1 - t1 * a2) / (t2 - t1)
    C = np.zeros((r, r))
    for i in range(r):
        e = np.zeros((r, 1)); e[i] = 1.
        C[i, i] = lin(e @ e.T)
    for i in range(r):
        for j in range(i + 1, r):
            ei = np.zeros((r, 1)); ei[i] = 1.
            ej = np.zeros((r, 1)); ej[j] = 1.
            C[i, j] = C[j, i] = (lin((ei + ej) @ (ei + ej).T) - lin((ei - ej) @ (ei - ej).T)) / 4.
    return sym(C)


def three_term(Pl, F, Phi0, fl, C):
    r = F.shape[0]
    Ci = np.linalg.inv(C)
    M0 = sym(F @ Pl.Pt(F, 1e-13 * np.eye(r)) @ F.T)
    trm = float(np.trace(np.linalg.solve(M0, Ci)))
    TS = np.logspace(-6, -3.4, 9)
    G = np.array([Pl.xI(F, t * Ci, Phi0)[0] - r * t for t in TS])
    dlt, b = np.linalg.lstsq(np.vstack([TS, TS ** 2]).T, G, rcond=None)[0]
    trD = [float(np.trace(np.linalg.solve(M0, sym(F @ Pl.Pt(F, t * Ci) @ F.T - M0) / t)))
           for t in (1e-4, 1e-5, 1e-6)][-1]
    X, IV = [], []
    for t in np.logspace(-4.3, -1.6, 9):
        x, I, _ = Pl.xI(F, t * Ci, Phi0)
        X.append(x); IV.append(I)
    X, IV = np.array(X), np.array(IV)
    cf = np.polyfit(X, IV + (r / 2) * np.log2(X), 2)
    Mfree = np.vstack([np.log2(1 / X), np.ones_like(X), X]).T
    cf3 = np.linalg.lstsq(Mfree, IV, rcond=None)[0]
    k_pred = (b + trm + trD) / (2 * r * ln2)
    print('  $\\delta$=%.1e  $b=q(B)$=%.5f  $\operatorname{tr}(M_0^{-1}C^{-1})$=%.5f  '
          '$\operatorname{tr}(M_0^{-1}D_{\mathbb M})$=%.5f' % (dlt, b, trm, trD))
    print('  斜率：自由拟合 %.6f  vs  $r/2$=%.6f   差 %+.1e' % (cf3[0], r / 2., cf3[0] - r / 2.))
    Iu = Pl.sdp(fl)
    if Iu is not None:
        bp = (r / 2) * np.log2(r) + .5 * np.log2(np.linalg.det(C) * np.linalg.det(M0)) - Iu
        print('  截距：$B_{测}$=%.5f   $\tfrac r2\log_2 r+\tfrac12\log_2(\det C\det M_0)$=%.5f   '
              '$I_{unc}($地板$)$=%.5f   $B_{测}-b_{pred}-I_{unc}$=%+.5f'
              % (cf3[1], bp + Iu, Iu, cf3[1] - bp - Iu))
    h = 5e-4
    su, sd = Pl.sdp(fl + h), Pl.sdp(fl - h)
    sl = (su - sd) / (2 * h) if su is not None and sd is not None else float('nan')
    print('  曲率：$\kappa\'_{测}$=%.5f  $\kappa\'_{预}$=%.5f  比值=%.4f   $|I\'_{unc}|$=%.5f  '
          '$\kappa_{总预}$=%.5f  $\zeta$=%.2e'
          % (cf[-2], k_pred, k_pred / cf[-2] if cf[-2] else float('nan'), abs(sl),
             k_pred + abs(sl), cf[0]))


def case(name, Pl, F, note='', depth=0):
    r = F.shape[0]
    det_ok, lam = Pl.detectable(F)
    ev = np.linalg.eigvals(Pl.A)
    su = float(np.sum(np.log(np.abs(ev[np.abs(ev) >= 1]))) / ln2)
    fl, phi0 = floor_of(Pl, F)
    C = measure_C(Pl, F, phi0)
    rr, evC = rank_psd(C)
    print('\n===== %s  n=%d r=%d  地板=%.6f  $\Phi$=%.4f  Nair–Evans=%.4f  '
          '$\operatorname{rank}\Theta$=%d  %s' % (name, Pl.nn, r, fl, phi0, su,
                                                  rank_psd(Pl.Th)[0], note))
    print('  可检测=%s  $|\lambda(A)|$=%s  $C$ 本征值=%s'
          % (det_ok, np.array2string(np.sort(np.abs(ev))[::-1], precision=4),
             np.array2string(evC, precision=7)))
    if not det_ok:
        print('  ↑ PBH 在 $\lambda=%s$ 失败：按原文 Thm.2，$\mathcal K_F$ 上的下确界不该在这里取到' % lam)
    if rr == 0:
        print('  $C\equiv0$：代价对该传感器的全部噪声不敏感 → 垂直渐近线不存在')
        return
    if rr < r:
        vec, U = np.linalg.eigh(C)
        Feff = U[:, -rr:].T @ F
        print('  **$C$ 降秩** $\operatorname{rank}C=%d<r=%d$：预言斜率 $(\operatorname{rank}C)/2=%.1f$，'
              '而非 $r/2=%.1f$；收缩到支集重跑' % (rr, r, rr / 2., r / 2.))
        if depth < 1:
            case(name + '｜支集', Pl, Feff, '（收缩后）', depth + 1)
        return
    three_term(Pl, F, phi0, fl, C)


def plants(seed=7):
    rng = np.random.default_rng(seed)
    rho, th = 1.18, 0.55
    A1 = np.array([[rho * np.cos(th), -rho * np.sin(th), 0.],
                   [rho * np.sin(th), rho * np.cos(th), 0.],
                   [0., 0., .42]])
    B1 = rng.standard_normal((3, 2)) * .8
    W1 = sym(rng.standard_normal((3, 3))) + 1.6 * np.eye(3)
    P1 = Pl(A1, B1, W1, np.eye(3))
    A2 = np.zeros((5, 5))
    A2[:2, :2] = 1.12 * np.array([[np.cos(1.1), -np.sin(1.1)], [np.sin(1.1), np.cos(1.1)]])
    A2[2, 2] = 1.31; A2[3, 3] = .61; A2[4, 4] = -.33
    A2 = A2 + np.triu(rng.standard_normal((5, 5)) * .18, 1)
    B2 = rng.standard_normal((5, 2)) * .8
    W2 = sym(rng.standard_normal((5, 5))) + 2.2 * np.eye(5)
    P2 = Pl(A2, B2, W2, np.eye(5))
    A3 = np.array([[1., 1., 0.], [0., 1., 0.], [0., 0., .5]])
    B3 = np.array([[0.], [1.], [.4]])
    W3 = sym(rng.standard_normal((3, 3))) + 1.3 * np.eye(3)
    P3 = Pl(A3, B3, W3, np.eye(3), np.eye(1))
    A4 = np.array([[1.25, .3, .1], [-.2, .7, .4], [.1, .2, -.55]])
    B4 = np.array([[1.], [.15], [1.2]])
    W4 = sym(rng.standard_normal((3, 3))) + 1.5 * np.eye(3)
    P4 = Pl(A4, B4, W4, np.eye(3), np.eye(1))
    import p0.p0_replicate_letter as base
    from p0.exp_nearfloor_law import Fs
    PT = Pl(base.A, base.B, base.W, np.eye(4))
    I = np.eye
    cases = [
        ('T-F2 尾部2（锚点）', PT, Fs['F2 尾部2(封闭, r=2)']),
        ('PL1 $r=1$ 坐标（复对只需一通道）', P1, np.array([[1., 0., 0.]])),
        ('PL1 $r=2$ 坐标', P1, np.array([[1., 0., 0.], [0., 1., 0.]])),
        ('PL1 $r=2$ 随机正交', P1, np.linalg.qr(rng.standard_normal((3, 2)))[0].T),
        ('PL2 $r=3$ 坐标', P2, I(5)[:3]),
        ('PL2 $r=3$ 随机正交', P2, np.linalg.qr(rng.standard_normal((5, 3)))[0].T),
        ('PL3 $r=1$ 观测 $x_1$（临界 Jordan）', P3, np.array([[1., 0., 0.]])),
        ('PL3 $r=2$', P3, np.array([[1., 0., 0.], [0., 1., 0.]])),
        ('PL3 $r=1$ 观测 $x_2$（应不可检测）', P3, np.array([[0., 1., 0.]])),
        ('PL4 单输入 $r=2$', P4, I(3)[:2]),
        ('PL4 单输入 $r=3$', P4, I(3))]
    return cases


if __name__ == '__main__':
    import p0.p0_replicate_letter as base
    from p0.exp_nearfloor_law import Fs
    print('锚点自检：Tanaka 植物 + Schur 尾部 2 —— 地板应回到 46.123137，$b$ 回到 $-0.840$，'
          '$\delta$ 应从 $8\times10^{-6}$ 掉到 $10^{-9}$')
    for nm, pl, F in plants():
        try:
            case(nm, pl, F)
        except Exception as e:
            print('\n===== %s  失败：%r' % (nm, e))
