"""E18：把 E17b 的两条结论钉死。

(A) 迭代数敏感性：case (2b)（ker F = span(两个不稳定 Schur 列)，C 判 ∞，我测得 Φ=+226.01）
    到底是**有限**还是**发散极慢**？对 iters 做 3 个量级，看 Φ 是否饱和。
    同时看 DARE 残差 ||A P Aᵀ + W − P||/||P||：真不动点 → 残差 → 0；慢发散 → 残差不降且 P 单调涨。
(B) C 的规则到底是"错"还是"保守"。对随机子空间族验证
        ρ(ZᵀAZ) ≥ ρ(A|_{V*})
    （Z = ker F 的正交基，V* = ker F 内最大 A-不变子空间）。
    若成立，则 C 的"ρ(ZᵀAZ)≥1 ⟹ Φ=∞"永远不会漏判真正的发散，只是会**过度宣告**；
    过度宣告的时机恰好是 ker F 不 A-不变（= 动态封闭被破坏），也就是我上一轮已证的两条判据重合之处。
"""
import sys
import numpy as np
sys.stdout.reconfigure(encoding='utf-8')
from p0.p0_replicate_letter import A, B, W, n, sym, ctrl
from p0.exp_authoritative import noiseless
from scipy.linalg import schur

Pc, K, Th = ctrl(np.eye(n))
unc = float(np.trace(W @ Pc))


def ker_to_F(Z):
    Z = np.linalg.qr(Z)[0] if Z.shape[1] else np.zeros((n, 0))
    Pp = sym(np.eye(n) - (Z @ Z.T if Z.shape[1] else 0.0))
    _, sv, Vt = np.linalg.svd(Pp)
    return Vt[:int(np.sum(sv > 1e-9 * sv[0]))], Z


def vstar(F):
    O = np.vstack([F @ np.linalg.matrix_power(A, i) for i in range(n)])
    _, sv, Vt = np.linalg.svd(O)
    d = n - int(np.sum(sv > 1e-9 * sv[0]))
    if d == 0:
        return 0, np.nan
    Zk = Vt[n - d:].T
    return d, float(np.max(np.abs(np.linalg.eigvals(Zk.T @ A @ Zk))))


if __name__ == '__main__':
    print('D_min^unc = %.4f\n' % unc)
    T, U, sdim = schur(A, output='real', sort=lambda a: np.abs(a) < 1.0)
    nu = n - sdim
    iU = list(range(n - nu, n))

    print('==== (A) case (2b)：ker F = span(不稳定 Schur 两列)，rank F=2 ====')
    Zb = np.ascontiguousarray(U[:, iU])
    Fb, _ = ker_to_F(Zb)
    print('%10s %14s %14s %22s' % ('iters', 'Φ=tr(ΘP)', '||P||_max', 'DARE 残差/||P||'))
    for it in (5000, 20000, 60000, 180000):
        try:
            P, _, _ = noiseless(Fb, A, W, iters=it, tol=1e-15)
            ok = np.all(np.isfinite(P))
            r = np.max(np.abs(A @ P @ A.T + W - P)) / max(np.max(np.abs(P)), 1e-300) if ok else np.nan
            print('%10d %14s %14s %22s' % (it, '%.5f' % np.trace(Th @ P) if ok else '发散',
                                           '%.3e' % np.max(np.abs(P)) if ok else '—',
                                           '%.3e' % r if np.isfinite(r) else '—'))
        except Exception as e:
            print('%10d  异常: %s' % (it, type(e).__name__))

    print('\n  对照：同样本族但 iters 递增时 Θ 权重是否饱和（若 P 沿 V* 无界则 Φ 必然继续涨）')
    print('  dim V* = %d, ρ(A|V*) = %s' % (vstar(Fb)[0],
          ('%.4f' % vstar(Fb)[1]) if vstar(Fb)[0] else '—'))

    print('\n==== (B) ρ(ZᵀAZ) ≥ ρ(A|V*) ？200 个随机子空间 ====')
    rng = np.random.default_rng(7)
    viol = 0
    n_inf_pred = n_over = 0
    for m in (1, 2, 3):
        for _ in range(200 // 3 + 1):
            Z = np.linalg.qr(rng.standard_normal((n, m)))[0]
            F, _ = ker_to_F(Z)
            d, rv = vstar(F)
            rc = float(np.max(np.abs(np.linalg.eigvals(sym(Z.T @ A @ Z)))))
            if d and rc < rv - 1e-9:
                viol += 1
            predC = rc >= 1.0
            truth = bool(d and rv >= 1.0)
            if predC:
                n_inf_pred += 1
            if predC and not truth:
                n_over += 1
    print('  违反 ρ(ZᵀAZ)≥ρ(A|V*) 的样本数 = %d / 200' % viol)
    print('  C 宣告 ∞ 的样本 = %d，其中 V* 确实含不稳定模态的 = %d → 过度宣告 %d'
          % (n_inf_pred, n_inf_pred - n_over, n_over))

    print('\n  再看 ker F **恰好 A-不变**的子空间族（动态封闭成立）：两判据应重合')
    cnt = 0
    for _ in range(4000):
        Z = np.linalg.qr(rng.standard_normal((n, 2)))[0]
        resid = np.max(np.abs(A @ Z - Z @ (Z.T @ A @ Z)))
        if resid < 1e-8:
            cnt += 1
            F, _ = ker_to_F(Z)
            d, rv = vstar(F)
            rc = float(np.max(np.abs(np.linalg.eigvals(sym(Z.T @ A @ Z)))))
            if d:
                assert abs(rc - rv) < 1e-6
    print('  找到 %d 个 2 维不变 ker F，全部满足 ρ(ZᵀAZ) = ρ(A|V*)（判据重合）' % cnt)

    print('\n  Θ 的零空间（Q=I 下 rank Θ = %d）：发散判据的精细版要看 Θ 是否看不见 V* 的不稳定模态'
          % np.linalg.matrix_rank(sym(Th)))
    wth, Vth = np.linalg.eigh(sym(Th))
    print('  Θ 特征值 = %s' % np.array_str(wth, precision=4))
    for name, vec in [('v_u(λ=-1.7124 特征向量)', None)]:
        w, V = np.linalg.eig(A)
        iu = int(np.argmax(np.abs(w) - 1.0))
        vu = np.real(V[:, iu]); vu /= np.linalg.norm(vu)
        print('  |Θ v_u| = %.4e   （>0 → Θ 看得见该不稳定模态 → Φ=∞ 成立）' % np.linalg.norm(Th @ vu))
