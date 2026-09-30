r"""E36：把 $C$ 从"差分测出来的"变成"解出来的"——$C$ 的秩决定垂直斜率，必须先把它算准。

E35 暴露的问题：$C$ 用两点 Richardson（$t=10^{-7},10^{-6}$）测，噪声底 $\sim10^{-5}$（绝对），
于是 PL4 $r=2$ 的假本征值 $2.35\times10^{-5}$ 混过了 $10^{-7}$ 的相对阈值，整条案子的
斜率、截距、$\kappa$ 全错（比值 23）。靠调阈值糊过去是可耻的，正确做法是**别用差分**。

$C$ 有闭式。走信息型 Riccati：$\Pi=\tilde P^{-1}$ 满足
$$\Pi=F^\top V^{-1}F+\bigl(A\Pi^{-1}A^\top+W\bigr)^{-1},\qquad
  P=(\Pi+F^\top V^{-1}F)^{-1}.$$
对 $V$ 求 Frechet 导数：$G(Z)=-F^\top V^{-1}ZV^{-1}F$，
$$d\Pi=\mathcal L[d\Pi]+G(Z),\qquad \mathcal L[U]=BUB^\top,\quad B=M A\Pi^{-1},\ M=(A\Pi^{-1}A^\top+W)^{-1},$$
一个 Stein 方程（$\mathcal I-\mathcal L$ 可逆 $\Leftrightarrow\rho(\mathcal L)<1$，即闭环稳定），
然后 $dP=-P\,(d\Pi+G)\,P$。整个映射 $Z\mapsto dP$ 是一个 $r^2\times n^2$ 矩阵，$C$ 就是它与
$\Theta$ 的配对。机器精度，没有差分噪声。

副产品是一个引理，不是数值：**$C\succeq0$**。$\mathcal L$ 是完全正映射且 $\rho(\mathcal L)<1$，
故 $(\mathcal I-\mathcal L)^{-1}=\sum_k\mathcal L^k$ 也保序；$Z\succeq0\Rightarrow G(Z)\preceq0
\Rightarrow d\Pi\preceq0$，于是 $\Pi+F^\top V^{-1}F$ 随 $V$ 单调减、$P$ 单调增：
$$\langle u,Cu\rangle=\operatorname{tr}\bigl(\Theta\,dP[uu^\top]\bigr)\ge0
  \qquad(\text{两个半定矩阵的配对}).$$
进而 $\ker C=\{u:\Theta^{1/2}dP[uu^\top]\Theta^{1/2}=0\}$——配对为零的**充要**条件，逐点可判定。
垂直斜率于是是 $\operatorname{rank}(C)/2$；"斜率只由 $\operatorname{rank}F$ 决定"的写法作废。
"""
import sys
import numpy as np
sys.stdout.reconfigure(encoding='utf-8')
import p0.exp_plants_oos as E35

ln2 = np.log(2.0)
sym = E35.sym
Pl = E35.Pl
floor_of = E35.floor_of
measure_C = E35.measure_C


def post_Pi(Pl, F, V, iters=4000):
    """信息型 Riccati 迭代到定点：$\\Pi\leftarrow(A\\Pi^{-1}A^\\top+W)^{-1}+F^\\top V^{-1}F$。
    不用 scipy 的对偶解，避免 $V\to0$ 时它的条件数问题；同时给两者之差作交叉检验。"""
    A, W = Pl.A, Pl.W
    Vi = np.linalg.inv(V)
    RV = F.T @ Vi @ F
    Pi = np.linalg.inv(Pl.Pt(F, V))         # 用 DARE 解起步
    for _ in range(iters):
        Pn = np.linalg.inv(Pi)
        Pi = sym(np.linalg.inv(A @ Pn @ A.T + W) + RV)
    P = sym(np.linalg.inv(Pi))
    sc = Pl.P(F, V)
    return Pi, P, float(np.linalg.norm(P - sc) / np.linalg.norm(sc))


def dP_op(Pl, F, V, diag=None):
    """$Z\\mapsto dP$ 的表示矩阵 $(n^2,r^2)$：$\\Pi=P^{-1}$ 满足 info-DARE
    $\\Pi=(A\\Pi^{-1}A^\\top+W)^{-1}+F^\\top V^{-1}F$，$d\\Pi-\\mathcal L[d\\Pi]=G(Z)$，$dP=-P\\,d\\Pi\\,P$。"""
    A, W, n = Pl.A, Pl.W, Pl.nn
    r = F.shape[0]
    Vi = np.linalg.inv(V)
    Pi, P, dev = post_Pi(Pl, F, V)
    M = np.linalg.inv(A @ P @ A.T + W)                    # $=S^{-1}$，先验的逆
    RV = F.T @ Vi @ F
    res = np.linalg.norm(Pi - M - RV) / np.linalg.norm(Pi)
    assert res < 1e-10, 'info-DARE 残差过大 %.2e' % res
    if diag is not None:
        Pt_ = Pl.Pt(F, V)
        res2 = np.linalg.norm(Pt_ - (A @ P @ A.T + W)) / np.linalg.norm(Pt_)
        diag.append((V[0, 0], res, res2, dev))
    Bm = M @ A @ P / 1.0                                  # $=M A \Pi^{-1}$
    Fv = F.T @ Vi                                         # n x r
    J_G = -np.kron(Fv, Fv)                                # vec(G) = J_G vec(Z)
    J_dPi = np.linalg.solve(np.eye(n * n) - np.kron(Bm, Bm), J_G)
    return -(np.kron(P, P) @ J_dPi), P, Bm


def pair(Pl, J, i, j):
    """$C_{ij}=\tfrac12\operatorname{tr}(\Theta\,dP[e_ie_j^\\top+e_je_i^\\top])$，$C_{ii}=\operatorname{tr}(\Theta\,dP[e_ie_i^\\top])$。
    极化必须用 $e_ie_j^\\top+e_je_i^\\top$，用 $e_ie_i^\\top+e_je_j^\\top$ 会把非对角元全测错（只留下对角）。"""
    n = Pl.nn
    r = int(round(np.sqrt(J.shape[1])))
    Z = np.zeros((r, r))
    if i == j:
        Z[i, i] = 1.
    else:
        Z[i, j] = Z[j, i] = 1.
    dPm = (J @ Z.reshape(-1)).reshape(n, n)
    v = float(np.trace(Pl.Th @ dPm).real)
    return v if i == j else v / 2.


def C_seq(Pl, F, ts, diag=None):
    out = []
    for t in ts:
        J, P, Bm = dP_op(Pl, F, t * np.eye(F.shape[0]), diag)
        r = F.shape[0]
        C = np.zeros((r, r))
        for i in range(r):
            for j in range(i, r):
                C[i, j] = C[j, i] = pair(Pl, J, i, j)
        out.append((t, sym(C), np.max(np.abs(np.linalg.eigvals(Bm)))))
    return out


def extrap(seq, p=1):
    """对 $t$ 做 $O(t^p)$ 外推到 0（末三点 Richardson）。"""
    (t1, C1, _), (t2, C2, _), (t3, C3, _) = seq[-3:]
    a1 = (t2 ** p * C1 - t1 ** p * C2) / (t2 ** p - t1 ** p)
    a2 = (t3 ** p * C2 - t2 ** p * C3) / (t3 ** p - t2 ** p)
    return .5 * (a1 + a2), a1, a2


def report_C(Pl, F, phi0, tag, ts=(3e-3, 1e-3, 3e-4, 1e-4, 3e-5)):
    seq = C_seq(Pl, F, ts)
    r = F.shape[0]
    V0 = ts[-1] * np.eye(r)
    _, P0, Bm = dP_op(Pl, F, V0)
    Pt0 = Pl.Pt(F, V0)
    Acl = Pl.A - Pl.A @ Pt0 @ F.T @ np.linalg.solve(F @ Pt0 @ F.T + V0, F)
    print('\n  [%s] $\\rho(B)=%.4f$ vs 滤波闭环 $\max|\lambda|=%.4f$（应当一致，Stein 可逆性的凭据）'
          % (tag, np.max(np.abs(np.linalg.eigvals(Bm))), np.max(np.abs(np.linalg.eigvals(Acl)))))
    for t, C, _ in seq:
        print('    $t=%.0e$  eig = %s' % (t, np.array2string(np.linalg.eigvalsh(C)[::-1],
                                                             precision=8)))
    Cfd = measure_C(Pl, F, phi0)
    print('    差分法的 $C$（E35 用的）= %s' % np.array2string(np.linalg.eigvalsh(Cfd)[::-1], precision=8))
    best = None
    for p in (1, 2):
        C0, a1, a2 = extrap(seq, p)
        d = np.linalg.norm(a1 - a2) / max(np.linalg.norm(C0), 1e-30)
        if best is None or d < best[0]:
            best = (d, p, C0)
    d, p, C0 = best
    print('    外推 $t\\to0$（$O(t^{%d})$，两分支差 %.1e 相对）eig = %s'
          % (p, d, np.array2string(np.linalg.eigvalsh(C0)[::-1], precision=10)))
    # 极化自检：随机的 $u$ 上比 $u^\top Cu$ 与 $\operatorname{tr}(\Theta dP[uu^\top])$
    rng = np.random.default_rng(11)
    worst = 0.
    for _ in range(4):
        u = rng.standard_normal(r)
        u /= np.linalg.norm(u)
        Jj, Pp, _b = dP_op(Pl, F, V0)
        dPm = (Jj @ np.outer(u, u).reshape(-1)).reshape(Pp.shape)
        worst = max(worst, abs(float(np.trace(Pl.Th @ dPm)) - u @ C0 @ u) / abs(u @ C0 @ u))
    print('    极化自检（$u^\top Cu$ vs 直接配对，4 个随机方向）最大相对偏差 %.2e' % worst)
    assert worst < 1e-2, '极化自检失败：$C$ 的非对角元装配有错'
    return C0, Cfd


def rank_C(C, rtol=1e-8):
    e, U = np.linalg.eigh(sym(C))
    e, U = e[::-1], U[:, ::-1]
    rho = int(np.sum(e > rtol * e[0]))
    return rho, e, U


def _sqrtm_psd(M):
    w, V = np.linalg.eigh(sym(M))
    return V @ np.diag(np.sqrt(np.clip(w, 0, None))) @ V.T


if __name__ == '__main__':
    from p0.exp_plants_oos import plants
    TS = (3e-3, 1e-3, 3e-4, 1e-4, 3e-5)
    rows = []
    for nm, pl, F in plants():
        r = F.shape[0]
        tag = nm.replace('$', '')
        print('\n===== %s   n=%d r=%d  rank$\\Theta$=%d' %
              (tag, pl.nn, r, E35.rank_psd(pl.Th)[0]))
        try:
            det_ok, lam = pl.detectable(F)
            if not det_ok:
                print('  不可检测（PBH 在 $\\lambda=%s$ 失败），跳过' % lam)
                continue
            phi0 = floor_of(pl, F)[1]
            C, Cfd = report_C(pl, F, phi0, tag, TS)
            rho, e, U = rank_C(C)
            print('  $\\operatorname{rank}C=%d$（相对阈 $10^{-8}$）  本征值 = %s'
                  % (rho, np.array2string(e, precision=9)))
            print('  半定性：$\\lambda_{min}=%.3e$，$\\operatorname{sym}$ 残差 %.2e，'
                  '差分/解析 $|C_{fd}-C|/|C|=%.2e$'
                  % (e[-1], np.linalg.norm(C - C.T) / np.linalg.norm(C),
                     np.linalg.norm(Cfd - C) / np.linalg.norm(C)))
            rows.append((tag, r, E35.rank_psd(pl.Th)[0], rho, e))
            if rho < r:
                print('  —— 核方向的机理：$\Theta^{1/2}dP[uu^\top]\Theta^{1/2}=0$ 还是 $dP=0$ ——')
                t0 = TS[-1]
                J, P, _ = dP_op(pl, F, t0 * np.eye(r))
                scale = np.linalg.norm(P) / t0          # $P=O(t)$，自然尺度
                for k in range(rho, r):
                    u = U[:, k]
                    dPm = sym((J @ np.outer(u, u).reshape(-1)).reshape(P.shape))
                    Th2 = _sqrtm_psd(pl.Th)
                    print(r'    $u_{%d}$: $|dP[uu^\top]|/(P/t)=%.3e$   '
                          r'$|\Theta^{1/2}dP\Theta^{1/2}|/|\Theta|=%.3e$   $u^\top Cu=%.2e$'
                          % (k, np.linalg.norm(dPm) / scale,
                             np.linalg.norm(Th2 @ dPm @ Th2) / np.linalg.norm(pl.Th),
                             u @ C @ u))
        except Exception as ex:
            print('  失败：%r' % ex)
    print('\n\n===== 汇总：$\\operatorname{rank}C$ 与 $r$、$\\operatorname{rank}\Theta$ 的关系')
    print('%-38s %3s %3s %3s   %s' % ('案', 'r', 'Θ秩', 'C秩', '本征值'))
    for tag, r, th_, rho, e in rows:
        print('%-38s %3d %3d %3d   %s' % (tag[:38], r, th_, rho,
                                          np.array2string(e[:r], precision=5)))
    r1 = sum(1 for _, r, t_, rho, _ in rows if rho > t_)
    print('  $\\operatorname{rank}C>\\operatorname{rank}\\Theta$ 的案子数：%d / %d  '
          '（若 $>0$，E35 里写的 $\\operatorname{rank}C\\le\\operatorname{rank}\\Theta$ 作废）'
          % (r1, len(rows)))

