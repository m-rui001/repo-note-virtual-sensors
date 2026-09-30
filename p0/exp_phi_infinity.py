"""E17：Φ 什么时候是 +∞？——对 C 的"若 V^⊥ 含不稳定模态（即 ρ(ZᵀAZ)≥1）则 Φ=∞"这条事实做判决。

设置：传感预算固定为 rank n−1，不可测方向 = ker F = span(w)，w 为单位向量。
两条候选判据：
  (i) C 的压缩判据：ρ(ZᵀAZ) = |wᵀAw|（Z=w）。它把 A 压到不可测方向上。
  (ii) 经典能观性判据：V* = ker F 内最大的 A-不变子空间（Kalman 分解里的不可观子空间）；
      可检测 ⟺ ρ(A|_{V*}) < 1。
差别只在 ker F **不** A-不变时出现：此时 wᵀAw 是"把 A 砍掉一半"的伪造动力学，
而真实误差递推里，w 经 A 生成的分量落进**可测**方向、立刻被下一次完美量测抹掉。
于是 ρ(wᵀAw)≥1 完全可以配上一个有限的 Φ。

每个 w 同时报：|wᵀAw|、dim V*、V* 上的谱、定点迭代结果、Φ、以及 rank P∞ 与 max|F P∞|。
"""
import numpy as np
from scipy.linalg import schur
from p0.p0_replicate_letter import A, B, W, n, sym
from p0.exp_authoritative import noiseless, ctrl_at

np.set_printoptions(precision=4, suppress=True, linewidth=200)
T, U, sdim = schur(A, output='real', sort=lambda a: np.abs(a) < 1.0)
nu = n - sdim
Pc, K, Th = ctrl_at(A, B, np.eye(n), np.eye(B.shape[1]))
unc = float(np.trace(W @ Pc))
u = U  # 简写


def unobs_ker(z):
    """ker = span(z)（dim 1）时最大的 A-不变子空间：z 是特征向量则 dim 1，否则 0"""
    Az = A @ z
    c = float(z @ Az)
    r = np.linalg.norm(Az - c * z)
    return (1, c) if r < 1e-9 * max(1.0, np.linalg.norm(Az)) else (0, np.nan)


def run(tag, z):
    z = z / np.linalg.norm(z)
    Zc = np.eye(n) - np.outer(z, z)                       # 投到 span(z) 的正交补（rank 3）
    _, sv, Vt = np.linalg.svd(sym(Zc))
    Wb = Vt[:n - 1].T                                     # 前 3 个右奇异向量 = 补空间正交基（sv[2]≈1, sv[3]=0）
    Fm = Wb.T                                             # 行张成 = span(z)^⊥，ker Fm = span(z)
    comp = float(z.T @ A @ z)
    dv, ev = unobs_ker(z)
    try:
        P, res, it = noiseless(Fm, A, W, iters=40000, tol=1e-14)
    except np.linalg.LinAlgError:
        P, res, it = np.full((n, n), np.nan), np.nan, -1
    fin = np.isfinite(P).all() and np.max(np.abs(P)) < 1e10
    phi = float(np.trace(Th @ P)) if fin else np.inf
    rk = np.linalg.matrix_rank(P, tol=1e-8) if fin else -1
    print('  %-30s |wᵀAw|=%.4f  C 判据→%s   dimV*=%d%s   Φ=%s  rankP∞=%d  res=%.1e@%d' %
          (tag, abs(comp), '∞' if abs(comp) >= 1 else '有限', dv,
           '' if dv == 0 else ' (λ=%.3f%s)' % (ev, '，不稳定' if abs(ev) >= 1 else '，稳定'),
           ('+∞(发散)' if not fin else '%+9.5f (%+.1f%%)' % (phi, 100 * phi / unc)),
           rk, res, it))
    return phi if fin else np.inf


if __name__ == '__main__':
    print('D_min^unc = %.4f，|λ| 排序后 = %s' % (unc, np.round(np.abs(np.linalg.eigvals(T)).tolist(), 4)))
    print('\n---- 不可测方向 ker F = span(w)，传感 rank n-1 ----')
    print('  [特征向量的情形：两条判据应当一致]')
    run('w = u2（实不稳定 λ=1.7124）', u[:, 2])
    run('w = u3（实不稳定 λ=1.3127）', u[:, 3])
    run('w = u0（复稳定对内一员）', u[:, 0])
    print('  [关键情形：w 不是特征向量 ⟹ ker F 不 A-不变 ⟹ V*={0} ⟹ 可检测]')
    run('w = (u2+u3)/√2（两个不稳定混合）', u[:, 2] + u[:, 3])
    run('w = (u2+u0)/√2（不稳定混稳定）', u[:, 2] + u[:, 0])
    run('w = (u2−u3)/√2', u[:, 2] - u[:, 3])
    rs = np.random.default_rng(0)
    for i in range(3):
        run('w = 随机单位向量 #%d' % (i + 1), rs.standard_normal(n))
    print('\n---- 与 C 的 Φ_min(r) 对照：上面任何有限值都是 rank n-1 架构的一个上界 ----')
    print('  C 报 Φ_min(3)=0.1587 (+0.5%)。若上面出现更小的值，说明 C 的优化器还漏了方向。')
