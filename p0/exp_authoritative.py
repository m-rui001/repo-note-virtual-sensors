"""E12（权威版）：只用 scipy.linalg.schur(sort=...) 的排序 Schur 形重建一切，
   彻底不碰 exp_trv_penalty.schur_task（它做朴素下标置换，行空间是否等于不稳定不变子空间未经证明）。
   同时把 E8/E9/E10 里互相矛盾的三件事一次判清：
     (a) (A,F) 的可观测性秩：若 ker F 是 A-不变的（前导 Schur 子空间），秩必 = rank F，V*≠0。
         ⟹ E8 报的 "dim V*=0 可观测=True" 是 tol 用错（相对 s[0] 太松）还是真的不变性破了？
     (b) 块结构预言：Q=FᵀF ⟹ Pc 的前导行列全零 ⟹ K 在余部上的增益=0 ⟹ Θ·Z=0（Z 张 ker F）
         ⟹ Φ(F;FᵀF) 必须精确为 0。E10 给出 +19.20 只能是 F 造错。
     (c) Φ(F;I) 的正确值，以及 D_min(K_F) 作为 Δ(D) 曲线左端点的可判决预测。
"""
import numpy as np
from scipy.linalg import schur, solve_discrete_are
from p0.p0_replicate_letter import A, B, W, n, sym, solve_

np.set_printoptions(precision=4, suppress=True, linewidth=200)
T, U, sdim = schur(A, output='real', sort=lambda a: np.abs(a) < 1.0)
nu = n - sdim                    # sdim = 前导(|λ|<1，即稳定)的个数；不稳定在尾部
iR, iT = list(range(n - nu)), list(range(n - nu, n))
F = U[:, iT].T                        # 行 = 尾部 nu 个 Schur 坐标（原文 §V 的 F=[0 I_ν]Ūᵀ）
print('|λ| =', np.round(np.abs(np.linalg.eigvals(A)).tolist(), 4))
print('sdim(不稳定块数) =', nu, ' rank F =', np.linalg.matrix_rank(F))
print('max|T[iT,iR]|（A_TR，封闭性）= %.3e' % np.max(np.abs(T[np.ix_(iT, iR)])))
print('max|T[iR,iT]|（A_RT）        = %.3e' % np.max(np.abs(T[np.ix_(iR, iT)])))

# ---- (a) 可观测性：三个容差下的秩
O = np.vstack([F @ np.linalg.matrix_power(A, k) for k in range(n)])
sv = np.linalg.svd(O, compute_uv=False)
print('\n可观测性矩阵奇异值 =', np.round(sv, 6), ' (s3/s0=%.2e)' % (sv[-1] / sv[0]))
for tol in (1e-8, 1e-10, 1e-12):
    print('  tol=%.0e → 观测秩=%d, dim V*=%d' % (tol, int((sv > tol * sv[0]).sum()),
                                                  n - int((sv > tol * sv[0]).sum())))
Zk = np.linalg.svd(O)[2][int((sv > 1e-10 * sv[0]).sum()):].T
print('  ker F 的一组正交基维数 =', np.linalg.matrix_rank(F.T, tol=1e-10) and
      (n - np.linalg.matrix_rank(F)), '；A-不变性残差 max|Zkᵀ A Z_perp|=%.3e' %
      np.max(np.abs(Zk.T @ A @ U[:, iT])))


def ctrl_at(Am, Bm, Qm, Rm):
    Pc = solve_discrete_are(Am, Bm, sym(Qm), sym(Rm))
    K = solve_(Rm + Bm.T @ Pc @ Bm, Bm.T @ Pc @ Am)
    return Pc, K, sym(K.T @ (Rm + Bm.T @ Pc @ Bm) @ K)


def noiseless(Fm, Am, Wm, iters=300000, tol=1e-15):
    r = Fm.shape[0]
    Pt = sym(Wm.copy()); res = np.inf
    for k in range(iters):
        M = Fm @ Pt @ Fm.T
        if np.min(np.linalg.eigvalsh(M)) < 1e-300:
            M = M + 1e-300 * np.eye(r)
        P = sym(Pt - Pt @ Fm.T @ np.linalg.solve(M, Fm @ Pt))
        Pn = sym(Am @ P @ Am.T + Wm)
        res = np.max(np.abs(Pn - Pt)); Pt = Pn
        if res < tol:
            break
    M = Fm @ Pt @ Fm.T
    return sym(Pt - Pt @ Fm.T @ np.linalg.solve(M, Fm @ Pt)), res, k + 1


def full(tag, Fm, Qm, Am=A, Wm=W, Bm=B):
    m = Bm.shape[1]
    try:
        Pc, K, Th = ctrl_at(Am, Bm, Qm, np.eye(m))
    except Exception as e:
        print('  %-30s LQR 无镇定解：%s' % (tag, str(e)[:60]))
        return None, None, None
    unc = float(np.trace(Wm @ Pc))
    P, res, it = noiseless(Fm, Am, Wm)
    if not np.isfinite(P).all():
        print('  %-30s 无噪声定点发散（Φ=+∞：买断信息也到不了任何代价）' % tag)
        return None, None, np.inf
    phi = float(np.trace(Th @ P))
    Zf = np.linalg.svd(Fm)[2][np.linalg.matrix_rank(Fm):].T     # ker F 的正交基
    tz = 0.0 if Zf.size == 0 else float(np.max(np.abs(Th @ Zf)))
    print('  %-30s Φ=%+10.5f (%+7.1f%%)  D_min^unc=%8.4f  D_min(K_F)=%9.4f' %
          (tag, phi, 100 * phi / unc, unc, unc + phi))
    print('      rankP∞=%d  max|F P∞|=%.2e  max|Θ Z_ker|=%.3e  res=%.1e@%d  '
          'λ(P∞)=%s' % (np.linalg.matrix_rank(P, tol=1e-9), np.max(np.abs(Fm @ P)),
                        tz, res, it,
                        np.round(np.sort(np.linalg.eigvalsh(P))[::-1], 4)))
    return phi, unc, unc + phi


if __name__ == '__main__':
    Qa = sym(F.T @ F)
    print('\n---- (b)(c) 两个代价下的地板')
    print(' F = 尾部 Schur 坐标（不稳定不变子空间），n=4, ν=2')
    full('Q = I（全状态代价）', F, np.eye(n))
    full('Q = FᵀF（任务对齐）', F, Qa)
    print('\n 对照：把不可测的稳定坐标也测（F 满秩）')
    full('Q = I, F=I', np.eye(n), np.eye(n))

    print('\n---- (a) 结论性：秩 1 的尾部坐标（只测 1 个不稳定 Schur 坐标，破 2×2 块）')
    F1 = U[:, [n - 1]].T
    print('  有限精度阶梯（Γ=γ·1）：J−D_min^unc 随 γ 是否收敛？')
    Pc, K, Th = ctrl_at(A, B, np.eye(n), np.eye(B.shape[1]))
    unc = float(np.trace(W @ Pc))
    for g in (1e2, 1e4, 1e6, 1e8, 1e10, 1e12):
        try:
            Pt = sym(solve_discrete_are(A.T, F1.T, W, sym(np.array([[1.0 / g]]))))
            Pp = sym(np.linalg.inv(np.linalg.inv(Pt) + F1.T @ np.array([[g]]) @ F1))
            print('    γ=%9.0e  tr(ΘP)=%+12.4f' % (g, float(np.trace(Th @ Pp))))
        except Exception as e:
            print('    γ=%9.0e  滤波 DARE 无有限解（%s）⟹ Φ=+∞' % (g, str(e)[:40]))
            break

    print('\n---- 三坐标（尾部 3 个：破坏前导复块的不变性 ⟹ A_TR≠0 ⟹ 不在 Thm.3 范围）')
    F3 = U[:, [1, 2, 3]].T
    full('F=尾3列, Q=I', F3, np.eye(n))
    full('F=尾3列, Q=FᵀF', F3, sym(F3.T @ F3))

    print('\n---- 命题检验：Φ=0 ⟺ Θ 在 ker F 上恒零（F=尾部2列不变子空间，Q 逐步加入余部权重）')
    Zc = U[:, iR]
    for eps in (0.0, 1e-4, 1e-2, 1e-1, 1.0):
        Qe = sym(Qa + eps * (Zc @ Zc.T))
        Pc, K, Th = ctrl_at(A, B, Qe, np.eye(B.shape[1]))
        P, res, it = noiseless(F, A, W)
        phi = float(np.trace(Th @ P))
        print('  ε=%-6g  max|Θ Z_ker|=%9.3e   Φ=%+10.5f  比值 Φ/max|ΘZ|=%8.2f'
              % (eps, np.max(np.abs(Th @ Zc)), phi, phi / max(np.max(np.abs(Th @ Zc)), 1e-300)))

    print('\n---- 与 E9/E10 的污染源对账：schur_task(2) 的行空间 vs 真正的不变子空间')
    from p0.exp_trv_penalty import schur_task
    Fb, _ = schur_task(2)
    Vt = U[:, iT]                                  # 真不变子空间的正交基
    Vb = np.linalg.qr(Fb.T)[0]                     # 朴素置换那版
    def inv_res(V):
        return float(np.max(np.abs((np.eye(n) - V @ V.T) @ A @ V)))
    sv = np.linalg.svd(Vt.T @ Vb, compute_uv=False)
    print('  不变性残差 ‖(I−VVᵀ)AV‖：真 Schur=%.3e   schur_task(2)=%.3e' % (inv_res(Vt), inv_res(Vb)))
    print('  两子空间主角度 cosσ =', np.round(sv, 6), ' → θ(度) =',
          np.round(np.degrees(np.arccos(np.clip(sv, -1, 1))), 2))
    print('  （这条同时给 C 的提案 1 一个现成的用法：主角度不是装饰，它是"子空间是否不变"的度量）')
    full('F=schur_task(2) 朴素置换, Q=I', Fb, np.eye(n))
    full('F=schur_task(2) 朴素置换, Q=FᵀF', Fb, sym(Fb.T @ Fb))
