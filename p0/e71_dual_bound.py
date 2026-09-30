r"""E71 = 前沿的**认证下界**：拉格朗日（信息正则）腿。修掉了 E70 的方向错误。

== 0. E70 错在哪（先写清楚，免得下一轮又顺着错的方向算）==
前沿的正确口径是**速率预算**：
    D_V(I) := inf{ D(S) : range S subsetof V,  I(S) <= I }.                (budget)
代价 D(S) 在 Löewin 序下单调不增（信息越多后验越小、定点越小），I(S) 单调不降，
所以最优一定落在边界 I(S)=I 上，曲线随 I 下降 —— 这才是图和 Table I 里的东西。

E70 的**代码**写的不是 I(S) >= I（这句是我上一轮的记忆，已被 e70b 否证）。
E70b 实测：E70 的切线约束 logdet M >= tan(G0) - 2 I0 ln2 蕴含 I <= I0，因为
    logdet G - logdet M <= logdet G - tan(G0) + 2 I0 ln2 <= 2 I0 ln2
（logdet 凹 ⇒ tan 是全局上界）。所以它是预算集的**子集 = CCP 内逼近**，
它的 min 是 D_V(I0) 的**上界**方向，给不出 Remark 2 要的下界；等号只在 G = G0 处成立，
故一个 I=I0 的物理设计对 E70 可行 **当且仅当** G0 恰取它自己的 G（e70b [T1] 两平面都验：
G0=I 时 slack=-42.5/-41.8 不可行，G0=G_phys 时 slack=-0.000000 可行）。
本文件 §2 的坍缩命题仍然成立，但它否定的是**另一条**腿（把集合放大成 0<=T<=I 盒子 + 下界口径
"至少 I0 比特"）——那条腿才坍缩到墙。两条腿各自死在不同方向上，别混用。

== 1. 拉格朗日腿（正确的方向，代价最小化 + 信息罚项）==
对任意 mu >= 0 和任意可行设计（I(S) <= I0）：
    D(S) = [D(S) + mu I(S)] - mu I(S)  >=  inf_designs[D + mu I] - mu I0.
所以
    L(mu) := inf_M [ D(M) + mu I(M) ] - mu I0     <=  D_V(I0)      （对每个 mu 都成立）
sup_mu L(mu) 就是认证下界；mu=0 给出墙。零对偶间隙 <=> 前沿在 I0 处有斜率 -mu 的支撑线
（即前沿关于 I 凸）。

注意符号：罚项是 **+ mu I**（信息要付钱），不是我上一轮推的 "- mu I"。
带 "- mu I" 时因为 D 有下界而 I 无上界，inf 恒为 -inf，对偶完全无用 —— 那条推论是
我自己算错的，作废。**这条腿在文献里已有名字**：Tanaka–Sandberg, arXiv:1503.01848
的 "information-regularized LQG"（二次代价 + 信息量的加权和，SDP 求解，线性传感假设下）；
Tanaka–Mohajerin Esfahani–Mitter, arXiv:1510.04214（最小 directed information 的 SDP）；
Tanaka et al., arXiv:1604.01227（prefix-free 码率的 SDP **下界**，给定 LQG 性能）。
它们都是"最小化速率"方向，所以凸性是站在它们那边的；我们是"最小化代价 + 速率预算"，
所以必须走罚项这条腿才不坍缩。

== 2. 坍缩命题（集合松弛腿的否定结果，可证）==
冻结先验 X，令 G := Z^t X Z，Y := (M^{-1}+G)^{-1}，T := G^{1/2} Y G^{1/2}，
Q := Z^t(X-K)Z = 后验在 V 上的限制。则
    D = j_c + tr(Theta X) - tr(Gamma T),  Gamma := G^{-1/2}(Z^t X Theta X Z)G^{-1/2},
    I = -1/2 log det(I - T),              0 <= T <= I,
而"至少 I0 比特"= det(I-T) <= 2^{-2 I0}。T=I 满足该约束（它对应 M=infinity，率无穷大），
又因为 Gamma >= 0，线性泛函 tr(Gamma T) 在整个盒子 0<=T<=I 上的最大值就在 T=I 取到。
=> **任何**介于可行集与盒子之间的松弛（凸包也算）的最大值都是 tr(Gamma) = 墙值，与 I0 无关。
下界只能靠"给信息定价"，不能靠"放大可行集"。

== 3. 实现：把 +mu I 里的先验项用可证常数下界替掉
X = A P A^t + W 且 P >= 0 => X >= W => G = Z^t X Z >= Z^t W Z > 0（W 是过程噪声协方差，正定）
=> logdet G >= logdet(Z^t W Z)。于是
    obj_mu(X,K) = j_c + tr(Theta(X-K)) + (mu/2ln2)[ logdet(Z^t W Z) - logdet(Z^t(X-K)Z) ]
在 (X,K) 上**凸**（仿射 + 仿射 - logdet(仿射)），约束
    X - A X A^t + A K A^t = W（Riccati 的线性化）,  X>=0, K>=0, X-K>=0
是凸的（仿射切片 ∩ 半正定锥）。丢掉 rank K <= r 与 range K ⊆ range(XZ) 得到的是**包含全部
物理设计**的松弛 => 它的最小值是 inf[D+mu I] 的合法下界 => L(mu) 合法。
代价：用常数 logdet(Z^tWZ) 替掉 logdet G 会漏掉 (mu/2ln2)(logdet G - logdet Z^tWZ) >= 0 的松弛量；
松多少在 [2] 里按解处实测报告。

== 判据（写在前）==
 [A] 口径退化演示：红平面 I=3 的"至少"口径最优应当 = 墙，而不是 Table I 的 50.6023。
     打印沿最优形状把尺度推到 1e13 的代价序列，应当单调降到墙。=> 证明 note.tex 第 102--103
     行印的 `I(S)\\ge I` 与全表数字不相容（数字是边界值，定义印反了）。
 [B] **合法性硬门**：L(mu) 绝不能超过任何已知可达值。已知：无约束 42.0427 (I=3)、
     红锥 50.6023 (I=3)、蓝 Powell 45.4537 (I=3)。任一行越界 => 实现有 bug，下界作废。
 [C] 验收线（板上的 §26.4）：蓝平面 I0=3 上 sup_mu L(mu) > 43.80 => Remark 2 从"相对上界的
     排除"升级成绝对不可行；<= 43.711 => 蓝线在该平面内可能可行，Remark 2 再降一档。
"""
import sys, time
sys.stdout.reconfigure(encoding='utf-8')
import numpy as np
import cvxpy as cp
from scipy.linalg import schur

np.set_printoptions(precision=4, suppress=True, linewidth=170)
_src = open('p0/exp_c_audit.py', encoding='utf-8').read().split("print(r'== E53")[0]
_ns = {'__name__': 'p'}
exec(compile(_src, 'p0/exp_c_audit.py[preamble]', 'exec'), _ns)
A, W, TH, JC, n, sym = _ns['A'], _ns['W'], _ns['TH'], _ns['JC'], _ns['n'], _ns['sym']
ln2 = np.log(2.0)
print('== E71：信息正则（拉格朗日）下界 ==')


def JI(S, tol=1e-12, maxit=20000):
    Pt = W.copy()
    for _ in range(maxit):
        Pm = sym(np.linalg.inv(np.linalg.inv(Pt) + S))
        Pn = sym(A @ Pm @ A.T + W)
        if not np.all(np.isfinite(Pn)) or np.max(np.abs(Pn)) > 1e14:
            return np.inf, np.inf
        if np.max(np.abs(Pn - Pt)) < tol * max(1.0, np.max(np.abs(Pn))):
            Pt = Pn; break
        Pt = Pn
    else:
        return np.inf, np.inf
    return (float(np.trace(TH @ Pm)) + JC,
            float(0.5 * (np.linalg.slogdet(Pt)[1] - np.linalg.slogdet(Pm)[1]) / ln2))


Ts, U = schur(A, output='real', sort=lambda a: abs(a) < 1.0)[:2]
PLANE = {'free_R4': U[:, :], 'blue_tail3': U[:, n - 3:], 'red_tail2': U[:, n - 2:]}
KNOWN = {'free_R4': 42.0427, 'blue_tail3': 45.4537, 'red_tail2': 50.6023}


def sdp_pen(Z, mu, solver='CLARABEL'):
    """inf over the relaxation of [D + mu*I]，先验项用 logdet(Z'WZ) 常数下界。返回 (lb, X, K)。"""
    r = Z.shape[1]
    X = cp.Variable((n, n), symmetric=True)
    K = cp.Variable((n, n), symmetric=True)
    Pm = X - K
    Q = Z.T @ Pm @ Z
    cons = [X - A @ X @ A.T + A @ K @ A.T == W, X >> 0, K >> 0, Pm >> 0, Q >> 0]
    c0 = np.linalg.slogdet(sym(Z.T @ W @ Z))[1]
    obj = JC + cp.trace(TH @ Pm) + (mu / (2.0 * ln2)) * (c0 - cp.log_det(Q))
    p = cp.Problem(cp.Minimize(obj), cons)
    p.solve(solver=solver, verbose=False)
    if p.status not in ('optimal', 'optimal_inaccurate'):
        return np.nan, None, None, p.status
    return float(p.value), np.array(X.value), np.array(K.value), p.status


def true_pen(Z, X, K):
    """在松弛解处算**真实**的 D 和 I（含 logdet G 项），用于报告松弛漏掉的量。"""
    Pm = sym(X - K)
    G = sym(Z.T @ X @ Z)
    Q = sym(Z.T @ Pm @ Z)
    D = JC + float(np.trace(TH @ Pm))
    I = 0.5 * (np.linalg.slogdet(G)[1] - np.linalg.slogdet(Q)[1]) / ln2
    res = np.max(np.abs(X - A @ X @ A.T + A @ K @ A.T - W))
    rankK = int(np.linalg.matrix_rank(sym(K), tol=1e-6))
    return D, I, float(res), rankK


def J_exact(S):
    """精确 Riccati 路线（sdare），不受定点迭代收敛/溢出影响；E68 的墙就用它算的。"""
    from scipy.linalg import solve_discrete_are as sdare
    ev, EV = np.linalg.eigh(sym(S))
    keep = ev > max(1e-11 * ev.max(), 1e-300)
    if not keep.any():
        return np.inf, np.inf
    F = (EV[:, keep] * np.sqrt(ev[keep])).T
    k = F.shape[0]
    try:
        Pt = sym(sdare(A.T, F.T, sym(W), np.eye(k)))
    except Exception:
        return np.inf, np.inf
    FPt = F @ Pt
    try:
        Rinv = np.linalg.solve(FPt @ F.T + np.eye(k), FPt)
    except Exception:
        return np.inf, np.inf
    P = sym(Pt - FPt.T @ Rinv)
    ld = np.linalg.slogdet(Pt)[1] - np.linalg.slogdet(P)[1]
    return JC + float(np.trace(TH @ P)), 0.5 * ld / ln2


I0 = 3.0
MUS = [0.0, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 9.0, 14.0, 22.0]

print('\n[A] 口径退化：note.tex 印的 I(S)>=I 会把最优推到墙上（沿红平面胜者形状放大尺度）')
Z2 = PLANE['red_tail2']
Rv = np.array([[np.cos(np.radians(128.0)), -np.sin(np.radians(128.0))],
               [np.sin(np.radians(128.0)), np.cos(np.radians(128.0))]])
S1 = sym(Z2 @ sym(Rv @ np.diag([1.0, 0.1]) @ Rv.T) @ Z2.T)
for v in (1e-1, 1e1, 1e3, 1e6, 1e10, 1e16):
    D, I = J_exact(v * S1)
    print('   尺度 v=%8.0e: D=%10.5f  I=%8.4f' % (v, D, I))
print('   E68 的红平面墙 = 46.1231，蓝平面墙 = 36.6032。')
print('   => 沿射线 I 单调增、D 单调降到墙。所以 inf{D: I>=3} = 墙 46.1231，'
      '不是 Table I 的 50.6023；表里的数是边界值 I(S)=3（= 预算口径 I(S)<=3 的最优）。')

print('\n[B]/[C] 拉格朗日下界 L(mu) = lb(mu) - mu*I0，I0=%.2f' % I0)
print('   硬门：L 不得超过已知可达值 %s' % KNOWN)
t0 = time.time()
res = {}
for name, Z in PLANE.items():
    rows = []
    for mu in MUS:
        lb, Xs, Ks, st = sdp_pen(Z, mu)
        if not np.isfinite(lb):
            rows.append((mu, np.nan, np.nan, np.nan, st))
            continue
        D, I, res_, rk = true_pen(Z, Xs, Ks)
        rows.append((mu, lb, lb - mu * I0, D, (I, res_, rk, st)))
    Ls = [(r[0], r[2]) for r in rows if np.isfinite(r[2])]
    mu_b, L_b = max(Ls, key=lambda a: a[1])
    res[name] = (mu_b, L_b, rows)
    bad = [L for _, L in Ls if L > KNOWN[name] + 1e-6]
    print('  %s (r=%d): sup L = %10.5f 在 mu*=%.2f | 已知可达 %.4f | 越界行数 %d %s'
          % (name, Z.shape[1], L_b, mu_b, KNOWN[name], len(bad),
             '**非法：实现有 bug**' if bad else '合法'))
    for r in rows:
        mu, lb, L, D, ex = r
        if not np.isfinite(lb):
            print('     mu=%6.2f  求解失败 %s' % (mu, ex))
            continue
        I, res_, rk, st = ex
        print('     mu=%6.2f  lb=%9.4f  L=%9.4f  | 解处真实 D=%9.4f I=%7.4f ARE残差=%.1e rankK=%d'
              % (mu, lb, L, D, I, res_, rk))

print('\n[结论判据]')
mu_b, L_b, _ = res['blue_tail3']
print('  蓝平面 I0=3：sup L = %.5f（mu*=%.2f）' % (L_b, mu_b))
print('   vs 43.711+0.09=43.80 验收线 -> %s' % ('Remark 2 升级为绝对不可行' if L_b > 43.80 else
      ('落在带内，符号仍不可读' if L_b > 43.711 else '蓝线在该平面内可能可行，Remark 2 再降一档')))
print('  与 Powell 上界 45.4537 的间隙 = %.4f（%.1f%% 的区间宽度）'
      % (45.4537 - L_b, 100.0 * (45.4537 - L_b) / 45.4537))
print('  耗时 %.0f s' % (time.time() - t0))
