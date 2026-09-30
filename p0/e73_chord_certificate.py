r"""E73 = **认证**下界：割线外逼近（chord / outer approximation）把 E71 的先验项补成合法凸下界。

问题回顾（E71/E72）：
  前沿 g_V(I0) = inf{ D(M) : range M subsetof V, I(M) <= I0 }（预算口径，见 note.tex 修好的定义）。
  对每个 mu>=0:  g_V(I0) >= p_V(mu) - mu*I0,  p_V(mu) := inf_M [D(M) + mu I(M)]。   (E71 §1)
  E71 的凸代理把先验项 (mu/2)logdet G(X) 用常数 logdet(Z'WZ) 替掉，损失正比于 mu：
  在 mu*~8（E72 定位的对偶最优）处丢掉了 20~30 个单位，所以 sup L 只有 33.0，够不到验收线 43.80。
  E72 说真对偶间隙几乎为零（自由平面 sup L_hat=42.021 vs 前沿 42.043），
  所以**瓶颈完全是那个常数替身**，而它是可以补的：

割线下界。X = A P A^t + W 且 P>=0 => X >= W => G = Z^t X Z 的每个特征值 >= l := lam_min(Z^t W Z) > 0。
再用"真最小点的目标值 <= Pbar"（Pbar 是 p(mu) 的任一上界，见下）挤出上界：
  D + mu I = j_c + tr(Theta P) + (mu/2ln2)[logdet G - logdet Q] >= j_c + (mu/2ln2) logdet G
  （因为 tr(Theta P)>=0，且 Q <= G 使括号 >=0）
  => 在最小点 logdet G <= b := (2 ln2/mu)(Pbar - j_c)。
若 Theta ≻ 0，还有 tr(Theta P) <= Pbar - j_c => lam_max(P) <= (Pbar-j_c)/lam_min(Theta)
  => lam_max(G) <= ||Z^t A||^2 (Pbar-j_c)/lam_min(Theta) + lam_max(Z^t W Z) =: u2。
取 u := min(exp(b - (r-1) log l), u2)。则最小点的每个特征值都在 [l,u] 内，而 log 在 [l,u] 上的
**割线**从下方托住它（凹函数弦在图下方）：
  log lam >= log l + s (lam - l),  s := (log u - log l)/(u - l)
  => logdet G >= r log l + s (tr G - r l)   —— 右端是 X 的**线性**函数。
把它代回目标，得到凸问题（线性 + 线性 - logdet(仿射) 在 PD 锥上 = 凸）：
  lb(mu) := min_{X,K} j_c + tr(Theta(X-K)) + (mu/2ln2)[ r log l + s(tr(Z'XZ) - r l) - logdet(Z'(X-K)Z) ]
        s.t.  X - A X A' + A K A' = W,  X>=0, K>=0, X-K>=0.
合法性：真最小点落在特征值 ∈[l,u] 的那块集合上，割线在那块集合上下界成立，而 SDP 在**更大**的集合
上取 inf，只会更小 => lb(mu) <= p(mu) => **L(mu) := lb(mu) - eps - mu I0 <= g_V(I0)**。
eps 是求解器容差余量（极小锥求解器给的是目标值的上界，不是精确下确界），这里显式扣掉并打印。
注意：**不能**往 SDP 里加 "Z'XZ <= u I" 这种约束——加约束会抬高 inf，破坏合法性。

Pbar 的合法取法（不需要任何搜索）：前沿设计本身。p(mu) <= g(I0) + mu*I0，因为我手上
"I0 处可达"的那个设计就是罚问题的一个可行点。于是
  自由平面：Pbar = 42.0427 + 3 mu；蓝平面：Pbar = 45.4537 + 3 mu。
（E72 的多重起点 p_hat 若干净，可取两者更小。）

判据（写在前）：
 [0] 硬门：L(mu) 绝不能超过已知可达值（自由 42.0427、蓝 45.4537、红 50.6023）。越界=推导或实现错。
 [1] 自由平面 mu=8：应当 lb 接近 p(8)=66.02，即 L 接近 42.0。若 L(8) >= 41.5，
     说明这条腿**能**认证，E71 的弱点确实只是常数替身。
 [2] 蓝平面验收线：sup_mu L(mu) > 43.80 => Remark 2 升级为绝对不可行（认证版）；
     43.711 < L <= 43.80 => 仍在读图带宽内，符号不可读；<= 43.711 => Remark 2 降档。
 [3] Pbar 敏感性：把 Pbar 乘 1.02 / 1.10 重算，看认证对"上界有多松"的敏感度。
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
I0 = 3.0
MARGIN = 1e-3          # 求解器容差余量，显式扣掉
print('== E73：割线外逼近的认证下界 ==')
print('   lam_min(Theta)=%.4e  lam_max(Theta)=%.4f  sigma(A)=%.4f'
      % (np.min(np.linalg.eigvalsh(TH)), np.max(np.linalg.eigvalsh(TH)),
         np.max(np.abs(np.linalg.eigvals(A)))))

Ts, U = schur(A, output='real', sort=lambda a: abs(a) < 1.0)[:2]
PLANE = {'free_R4': U[:, :], 'blue_tail3': U[:, n - 3:], 'red_tail2': U[:, n - 2:]}
# Pbar := 前沿可达值 + mu*I0（合法上界）
REACH = {'free_R4': 42.0427, 'blue_tail3': 45.4537, 'red_tail2': 50.6023}


def bounds(Z, mu, Pbar):
    r = Z.shape[1]
    Zw = sym(Z.T @ W @ Z)
    l = float(np.min(np.linalg.eigvalsh(Zw)))
    b = (2.0 * ln2 / mu) * (Pbar - JC)
    u1 = np.exp(b - (r - 1) * np.log(l))
    lt = float(np.min(np.linalg.eigvalsh(TH)))
    u2 = float(np.linalg.norm(Z.T @ A) ** 2 * (Pbar - JC) / lt
               + np.max(np.linalg.eigvalsh(Zw))) if lt > 1e-12 else np.inf
    u = min(u1, u2)
    if not (np.isfinite(u) and u > l * 1.000001):
        return None
    return r, l, u, u1, u2, (np.log(u) - np.log(l)) / (u - l)


def chord_lb(Z, mu, Pbar, solver='CLARABEL'):
    bb = bounds(Z, mu, Pbar)
    if bb is None:
        return np.nan, None, bb
    r, l, u, u1, u2, s = bb
    X = cp.Variable((n, n), symmetric=True)
    K = cp.Variable((n, n), symmetric=True)
    Pm = X - K
    Q = Z.T @ Pm @ Z
    cons = [X - A @ X @ A.T + A @ K @ A.T == W, X >> 0, K >> 0, Pm >> 0, Q >> 0]
    lin = r * np.log(l) + s * (cp.trace(Z.T @ X @ Z) - r * l)
    obj = JC + cp.trace(TH @ Pm) + (mu / (2.0 * ln2)) * (lin - cp.log_det(Q))
    p = cp.Problem(cp.Minimize(obj), cons)
    p.solve(solver=solver, verbose=False)
    if p.status not in ('optimal', 'optimal_inaccurate'):
        return np.nan, None, bb
    Xv = np.array(X.value)
    G = sym(Z.T @ Xv @ Z)
    return float(p.value), dict(lamG=np.linalg.eigvalsh(G), logdetG=float(np.linalg.slogdet(G)[1]),
                                P=np.array(Xv - K.value), are_res=float(np.max(np.abs(
                                    Xv - A @ Xv @ A.T + A @ np.array(K.value) @ A.T - W)))), bb


MUS = [2.0, 4.0, 6.0, 8.0, 11.0, 15.0]
t0 = time.time()
print('\n[0]/[1] L(mu) = lb(mu) - %.0e - mu*I0，Pbar = 可达值 + mu*I0' % MARGIN)
sup = {}
for name, Z in PLANE.items():
    print('  %s (r=%d, 可达上界 %.4f)' % (name, Z.shape[1], REACH[name]))
    best = (-np.inf, None)
    for mu in MUS:
        Pbar = REACH[name] + mu * I0
        lb, info, bb = chord_lb(Z, mu, Pbar)
        if not np.isfinite(lb):
            print('     mu=%5.1f  不可用（界太松 bb=%s）' % (mu, bb is None))
            continue
        L = lb - MARGIN - mu * I0
        if L > best[0]:
            best = (L, mu)
        viol = L > REACH[name] + 1e-6
        if bb is None:
            print('     mu=%5.1f  lb=%9.4f L=%9.4f  界构造失败' % (mu, lb, L))
            continue
        r, l, u, u1, u2, s = bb
        print('     mu=%5.1f  lb=%9.4f L=%9.4f %s | l=%.3e u=%.3e (u1=%.3e u2=%.3e) s=%.3e'
              ' | 解处 logdetG=%.3f ARE残差=%.1e'
              % (mu, lb, L, '**越界**' if viol else '', l, u, u1, u2, s,
                 info['logdetG'], info['are_res']))
    sup[name] = best
    print('     sup L = %.5f 在 mu=%s' % (best[0], best[1]))

print('\n[2] 验收线')
for name in ('free_R4', 'blue_tail3'):
    L, mu = sup[name]
    print('  %-11s sup L = %9.5f (mu=%s) | vs 43.711 %s | vs 43.80 %s'
          % (name, L, mu, '超过' if L > 43.711 else '未超过',
             '超过=>Remark 2 认证升级' if L > 43.80 else '未超过'))
print('\n[3] Pbar 敏感性（蓝平面 mu=8）')
Z = PLANE['blue_tail3']
for f in (1.0, 1.02, 1.10, 1.30):
    Pbar = f * (REACH['blue_tail3'] + 8.0 * I0)
    lb, info, bb = chord_lb(Z, 8.0, Pbar)
    print('   Pbar x%.2f = %8.4f -> lb=%9.4f  L=%9.4f' % (f, Pbar, lb, lb - MARGIN - 24.0))
print('  耗时 %.0f s' % (time.time() - t0))
