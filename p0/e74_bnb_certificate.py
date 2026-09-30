r"""E74 = **认证下界（分支定置版）**：对先验特征值做对数均匀分格，每格一条割线，全局界 = 各格 SDP 的极小。

E73 的教训（写下来）：一条覆盖 [l,u]=[0.30, 8e3] 的割线，斜率 s~1.3e-3，在下端误差 ~1.2 nat/特征值，
乘 (mu/2ln2)~5.8 就是 20+ 单位的松弛 —— sup L 只有 23.7，比 E71 的常数替身（33.0）还差。
根本原因：最小化一个**凹**项 (+mu/2 logdet G)，单条割线在整段区间上不可能又合法又紧。

修法：凸下界 + 逐格细化 = 标准分支定置。
 1) 全局 caps（可证，只切除不可能含最小点的区域，因此加约束合法）：
    X = A P A^t + W, P>=0 => X >= W => lam_i(G) >= l := lam_min(Z^t W Z)。
    罚目标在最小点 <= Pbar（Pbar = 任一已知可达设计的罚值，见下）且
        D + mu I >= j_c + (mu/2ln2) logdet G   （tr(Theta P)>=0，Q<=G）
    => logdet G <= b := (2 ln2/mu)(Pbar - j_c)；配合 lam_i >= l => lam_i <= u := exp(b-(r-1)log l)。
    => 最小点的每个特征值都在 [l,u]。
 2) 把 [l,u] 对数均分成 N 格 [a_k, b_k]。第 k 格 = {a_k I <= G <= b_k I}（LMI，凸）。
    凹函数 log 在 [a_k,b_k] 上的**弦**从下方托住它：
        log lam >= log a_k + s_k (lam - a_k),  s_k := (log b_k - log a_k)/(b_k - a_k)
        => logdet G >= r log a_k + s_k (tr G - r a_k)      （tr G = tr(Z^t X Z)，X 的线性函数）
 3) 每格解一个凸 SDP：
        min j_c + tr(Theta (X-K)) + (mu/2ln2)[ r log a_k + s_k(tr(Z^t X Z) - r a_k) - logdet(Z^t(X-K)Z) ]
        s.t. X - A X A^t + A K A^t = W, X>=0, K>=0, X-K>=0, a_k I <= Z^t X Z <= b_k I, Z^t(X-K)Z >= 0
    该格的 SDP 极小 <= 该格内真极小（弦是下界 + (X,K) 松弛丢掉了 rank<=r 与 range K subsetof range(XZ)，
    两者都只会让 inf 变小），所以
        p(mu) >= min_k SDP_k        （并集覆盖 [l,u]^r，最小点必在某格）
        L(mu) := min_k SDP_k - eps - mu I0  <=  g_V(I0).        【认证】
 4) 弦误差上界：格宽 w nat 时 log 与弦的最大偏差 <= w^2/8（每特征值），总松弛 <= r w^2/8 * mu/(2ln2)。
    mu=8, r=3, N=40 => w=0.25 => 松弛 ~0.13 单位。所以剩下的是 (X,K) 秩松弛的损失 —— 那正是这个证书要量的东西。

Pbar 的合法取法（零搜索）：前沿上已知的可达设计就是罚问题的可行点 => p(mu) <= 可达值 + mu*I0。
判据：
 [0] 硬门：L(mu) <= 已知可达值（自由 42.0427 / 蓝 45.4537 / 红 50.6023），越界=实现或推导错。
 [1] 自由平面（r=n，(X,K) 松弛可证**精确**）：L(mu) 应逼近 E72 的 p(8)-24 = 42.02。
     到了 => 整套证书机制正确，蓝平面差值就是秩松弛损失。
 [2] 蓝平面验收线：sup_mu L(mu) > 43.80 => Remark 2 升级为认证不可行。
 [3] N 敏感性：N=10/20/40/80，看界随分格的收敛（应当单调上升，因为每格弦更紧？注意 min over cells
     不是单调的，但整体上应收紧）。
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
EPS_SOLVER = 2e-3
print('== E74：分格割线的认证下界（分支定置）==')

Ts, U = schur(A, output='real', sort=lambda a: abs(a) < 1.0)[:2]
PLANE = {'free_R4': U[:, :], 'blue_tail3': U[:, n - 3:], 'red_tail2': U[:, n - 2:]}
REACH = {'free_R4': 42.0427, 'blue_tail3': 45.4537, 'red_tail2': 50.6023}
MUS = [2.0, 4.0, 6.0, 8.0, 11.0, 15.0]
NS = [10, 20, 40, 80]


def cell_sdp(Z, mu, a, b):
    """单格凸 SDP。返回 (值, 状态)。"""
    r = Z.shape[1]
    s = (np.log(b) - np.log(a)) / (b - a)
    X = cp.Variable((n, n), symmetric=True)
    K = cp.Variable((n, n), symmetric=True)
    Pm = X - K
    G = Z.T @ X @ Z
    Q = Z.T @ Pm @ Z
    cons = [X - A @ X @ A.T + A @ K @ A.T == W, X >> 0, K >> 0, Pm >> 0, Q >> 0,
            G >> a * np.eye(r), b * np.eye(r) - G >> 0]
    lin = r * np.log(a) + s * (cp.trace(G) - r * a)
    obj = JC + cp.trace(TH @ Pm) + (mu / (2.0 * ln2)) * (lin - cp.log_det(Q))
    p = cp.Problem(cp.Minimize(obj), cons)
    try:
        p.solve(solver='CLARABEL', verbose=False)
    except Exception as e:
        return np.inf, 'exc:%s' % type(e).__name__
    if p.status in ('infeasible', 'infeasible_inaccurate'):
        return np.inf, p.status
    if p.status not in ('optimal', 'optimal_inaccurate'):
        return np.nan, p.status
    return float(p.value), p.status


def cert_lb(Z, mu, Pbar, N):
    r = Z.shape[1]
    Zw = sym(Z.T @ W @ Z)
    l = float(np.min(np.linalg.eigvalsh(Zw)))
    b_bound = (2.0 * ln2 / mu) * (Pbar - JC)
    u = np.exp(b_bound - (r - 1) * np.log(l))
    if not (np.isfinite(u) and u > l * 1.00001):
        return np.nan, None
    edges = np.exp(np.linspace(np.log(l), np.log(u), N + 1))
    best, argk, nbad = (np.inf, None, 0), None, 0
    for k in range(N):
        v, st = cell_sdp(Z, mu, edges[k], edges[k + 1])
        if not np.isfinite(v):
            nbad += 1
            continue
        if v < best[0]:
            best, argk = (v, k), (edges[k], edges[k + 1], st)
    return best[0], dict(u=u, l=l, N=N, k=argk, nbad=nbad,
                         w=(np.log(u) - np.log(l)) / N,
                         chord_slack=r * ((np.log(u) - np.log(l)) / N) ** 2 / 8.0 * mu / (2 * ln2))


t0 = time.time()
out = {}
for name, Z in PLANE.items():
    r = Z.shape[1]
    print('\n%s (r=%d, 可达上界 %.4f)' % (name, r, REACH[name]))
    for N in NS:
        rows, best = [], (-np.inf, None)
        for mu in MUS:
            Pbar = REACH[name] + mu * I0
            lb, info = cert_lb(Z, mu, Pbar, N)
            if not np.isfinite(lb):
                rows.append((mu, np.nan, np.nan, info))
                continue
            L = lb - EPS_SOLVER - mu * I0
            rows.append((mu, lb, L, info))
            if L > best[0]:
                best = (L, mu)
        ok = all((not np.isfinite(L)) or L <= REACH[name] + 1e-6 for _, _, L, _ in rows)
        out[(name, N)] = best
        print('  N=%3d: sup L = %9.5f 在 mu=%-5s | 硬门 %s | 失败格 %d/%d | 弦松弛上界 ~%.3f'
              % (N, best[0], best[1], '合法' if ok else '**越界**',
                 sum(1 for _, _, _, i in rows if i and i.get('nbad', 0) > 0) if rows else 0, len(MUS),
                 np.nanmax([i['chord_slack'] for _, _, _, i in rows if i]) if any(
                     i for _, _, _, i in rows if i) else np.nan))
        for mu, lb, L, info in rows:
            if not np.isfinite(lb):
                print('     mu=%5.1f  不可用' % mu)
                continue
            print('     mu=%5.1f  lb=%9.4f  L=%9.4f  u=%.2e 最差格=[%.2f,%.2f] %s 废格%d'
                  % (mu, lb, L, info['u'], info['k'][0], info['k'][1], info['k'][2], info['nbad']))

print('\n[判据]')
for name in ('free_R4', 'blue_tail3'):
    for N in (20, 40, 80):
        L, mu = out[(name, N)]
        print('  %-11s N=%3d: sup L = %9.5f (mu=%s) | 验收线 43.80: %s'
              % (name, N, L, mu, '超过 => Remark 2 认证升级' if L > 43.80 else '未超过'))
print('  耗时 %.0f s' % (time.time() - t0))
