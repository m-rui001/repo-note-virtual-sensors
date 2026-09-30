"""E17b：对 C 的散度判据做双向反证。

C 的规则（.work3/c13_joint.py 的 docstring）：
    "若 V^⊥ 含不稳定模态（等价 ρ(Z^⊤AZ) ≥ 1）则 Φ = ∞"
这里 V^⊥ = ker F，Z 是它的一组正交基。这条规则把 **ker F 本身**当成不变子空间来看，
但正确的对象是 **ker F 里最大的 A-不变子空间 V***（即经典不可观测子空间）：
    Φ < ∞  ⟺  (A, F) 可检测  ⟺  ρ(A|_{V*}) < 1
两者只在 ker F 恰好 A-不变（= 动态封闭）时重合。

本脚本用真实特征向量（不是 Schur 列；排序 Schur 的单列并非特征向量，
A u_j = Σ_{i≤j} u_i T[i,j]）构造三种 ker F：
  (1) span(不稳定特征向量)          —— 两条规则都该说 ∞（一致性对照）
  (2) span(不稳定特征向量, 随机方向) —— V* 含不稳定模态⇒Φ=∞，但压缩谱半径可能 <1⇒C 说有限
  (3) span(两个稳定 Schur 列)       —— 封闭、不可测部分全稳定⇒Φ 有限
"""
import sys
import numpy as np
sys.stdout.reconfigure(encoding='utf-8')
from p0.p0_replicate_letter import A, B, W, n, sym, ctrl, rate_cost
from p0.exp_authoritative import noiseless, ctrl_at
from scipy.linalg import schur

np.set_printoptions(precision=4, suppress=True)
Pc, K, Th = ctrl(np.eye(n))
unc = float(np.trace(W @ Pc))


def ker_to_F(Z):
    """给 ker F 的一组基 Z（n×m），返回行满秩的 (n-m)×n 量测算子 F，使 ker F = span(Z)。"""
    Z = np.linalg.qr(Z)[0] if Z.shape[1] else np.zeros((n, 0))
    Pp = sym(np.eye(n) - (Z @ Z.T if Z.shape[1] else 0.0))
    _, sv, Vt = np.linalg.svd(Pp)
    r = int(np.sum(sv > 1e-9 * sv[0]))
    return Vt[:r], r


def vstar(F):
    """ker F 里最大的 A-不变子空间（可观测性矩阵的零空间）+ 其限制谱半径。"""
    O = np.vstack([F @ np.linalg.matrix_power(A, i) for i in range(n)])
    _, sv, Vt = np.linalg.svd(O)
    d = n - int(np.sum(sv > 1e-9 * sv[0]))
    if d == 0:
        return 0, np.nan, np.nan
    Z = Vt[n - d:].T                                  # n×d，张成 V*
    resid = np.max(np.abs(A @ Z - Z @ (Z.T @ A @ Z)))  # A-不变性残差
    return d, float(np.max(np.abs(np.linalg.eigvals(Z.T @ A @ Z)))), resid


def compression(Z):
    return float(np.max(np.abs(np.linalg.eigvals(Z.T @ A @ Z))))


def report(tag, Z):
    F, r = ker_to_F(Z)
    d, rhoV, resid = vstar(F)
    Zk = np.linalg.svd(F)[2][r:].T                     # ker F 的正交基（C 的压缩判据用）
    try:
        P, res, it = noiseless(F, A, W, iters=60000, tol=1e-14)
        phi = float(np.trace(Th @ P))
        finite = np.all(np.isfinite(P)) and np.max(np.abs(P)) < 1e10
    except Exception:
        phi, finite, P = np.inf, False, None
    c_rule = '∞' if compression(Zk) >= 1.0 else '有限'
    mine = '∞' if (d and rhoV >= 1.0) else '有限'
    print('%-34s rank F=%d  dimV*=%d  ρ(A|V*)=%s  不变残差=%s'
          % (tag, r, d,
             ('%.4f' % rhoV) if d > 0 else '—',
             ('%.1e' % resid) if d > 0 else '—'))
    print('   C 的压缩判据 ρ(ZᵀAZ)=%.4f → 预测 %-4s ；我的可检测性判据 → 预测 %-4s ；实测 Φ = %s (%+.4f, %+.1f%%)'
          % (compression(Zk), c_rule, mine,
             ('%.5f' % phi) if finite else '发散',
             phi - unc if finite else np.inf, 100 * (phi - unc) / unc if finite else np.inf))
    return finite, phi


if __name__ == '__main__':
    print('D_min^unc = %.4f' % unc)
    w, V = np.linalg.eig(A)
    w = np.real_if_close(w, tol=1000)
    print('真实特征值 = %s' % np.array_str(np.real(w), precision=4))
    iu = int(np.argmax(np.abs(np.real(w)) - 1.0)) if np.any(np.abs(w) > 1) else 0
    iu = int(np.argmax(np.abs(w) - 1.0))
    vu = np.real(V[:, iu]); vu /= np.linalg.norm(vu)
    print('不稳定特征向量 λ=%.4f  |λ|>1: %s' % (np.real(w[iu]), abs(w[iu]) > 1))
    print('   校验 A v = λ v 残差 = %.2e' % np.max(np.abs(A @ vu - np.real(w[iu]) * vu)))

    T, U, sdim = schur(A, output='real', sort=lambda a: np.abs(a) < 1.0)
    iT = list(range(n - sdim, n)) if False else list(range(n - (n - sdim), n))
    # 排序后不稳定块在最后 n-sdim 列
    nu = n - sdim
    iUnstable = list(range(n - nu, n))
    iStable = list(range(sdim))

    print('\n--- (1) ker F = span(真实不稳定特征向量)，rank n-1 ---')
    report('ker F = span(v_u)', vu[:, None])

    print('\n--- (2) ker F = span(不稳定特征向量, 随机方向)，rank n-2 ---')
    rng = np.random.default_rng(3)
    for i in range(3):
        rr = rng.standard_normal(n)
        report('ker F = span(v_u, rnd%d)' % i, np.column_stack([vu, rr]))

    print('\n--- (2b) ker F = span(两个不稳定 Schur 列)（= 论文的尾部不稳定块）---')
    report('ker F = span(u_iT[0], u_iT[1])', U[:, iUnstable])

    print('\n--- (3) ker F = span(两个稳定 Schur 列)（封闭、全稳定不可测）---')
    report('ker F = span(u_stable, 2)', U[:, iStable[:2]] if sdim >= 2 else U[:, iStable])

    print('\n--- (4) ker F = span(不稳定特征向量 + 一个稳定 Schur 列) ---')
    report('ker F = span(v_u, u_s0)', np.column_stack([vu, U[:, iStable[0]]]))

    print('\n--- (5) 纯压缩对照：ker F 不含任何不变方向 ---')
    rng = np.random.default_rng(11)
    for i in range(2):
        zz = rng.standard_normal((n, 1))
        report('ker F = span(rnd%d)' % i, zz)
