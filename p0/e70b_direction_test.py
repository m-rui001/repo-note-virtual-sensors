r"""E70b = 判定 E70 的切线约束到底是**松弛（超集）**还是**内逼近（子集）**。

判据（先写下，再跑）：
  E70 的约束是  logdet M >= tan(G0) - 2 I0 ln2，M = Zᵗ(X-K)Z 是 V 上的后验，tan 是 logdet G 的切线。
  因为 logdet 凹，tan(G0) >= logdet G 全局成立，所以
      logdet G - logdet M <= logdet G - tan(G0) + 2 I0 ln2 <= 2 I0 ln2，
  即约束**蕴含** I <= I0：它是速率预算集的一个**子集**（CCP 内逼近），不是超集。
  ⇒ 它的 min 是 D_V(I0) 的**上界**方向，永远不可能是 Remark 2 要的认证下界。
  并且等号只在 G = G0 处成立 ⇒ 一个满足 I = I0 的物理设计可行 **当且仅当** 它的 G 恰为 G0。

  [T1] 取一个 I=3 的物理设计（标量 S=s·ZZᵗ，二分 s）。它在预算集的边界上。
       若它对 E70 的约束不可行（G0=I、3I），而只对 G0=G_phys 恰好取等，
       则上面两行的证明在数字上成立：E70 不是松弛。
  [T2] 用 G0=G_phys 解 E70：内逼近在切线点精确 ⇒ 期望 L = D_phys（而不是 42.0427 以下）。
       若 L 落在 [D_min, D_phys] 且从不跌破 42.0427，方向判定的实验部分完成。
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


def fixedpoint(S, tol=1e-13, cap=1e12):
    Pt = W.copy()
    conv = False
    for _ in range(40000):
        Pm = sym(np.linalg.inv(np.linalg.inv(Pt) + S))
        Pn = sym(A @ Pm @ A.T + W)
        if np.max(np.abs(Pn)) > cap:
            return None, None, np.inf, np.inf, False
        if np.max(np.abs(Pn - Pt)) < tol * max(1.0, np.max(np.abs(Pn))):
            conv = True
            break
        Pt = Pn
    X, P = Pt, Pm
    D = float(np.trace(TH @ P)) + JC
    I = 0.5 * (np.linalg.slogdet(X)[1] - np.linalg.slogdet(P)[1]) / ln2
    return X, X - P, D, I, conv


def scale_for_rate(Z, target, lo=1e-6, hi=1e4, it=120):
    """在平面 V=span(Z) 上取标量设计 S = s ZZᵗ，二分 s 使 I(S)=target。
    I(S) 在 s 上单调增；未收敛（信息不足以镇定）当作"信息太少"，抬 lo。"""
    P = Z @ Z.T
    for _ in range(it):
        mid = np.sqrt(lo * hi)
        _, _, _, I, conv = fixedpoint(mid * P)
        if (not conv) or I < target:
            lo = mid
        else:
            hi = mid
        if abs(np.log(hi / lo)) < 1e-11:
            break
    s = np.sqrt(lo * hi)
    X, K, D, I, conv = fixedpoint(s * P)
    return X, K, D, I, s, conv


def e70_set_value(Z, X, K, I0, G0):
    """在物理点上逐项评估 E70 的约束，并返回其切线残差（>0 才可行）。"""
    r = Z.shape[1]
    G = sym(Z.T @ X @ Z)
    M = sym(G - Z.T @ K @ Z)
    res_lin = np.max(np.abs(X - A @ X @ A.T + A @ K @ A.T - W))
    eigX = np.linalg.eigvalsh(sym(X))
    eigK = np.linalg.eigvalsh(sym(K))
    eigXM = np.linalg.eigvalsh(sym(X - K))
    eigM = np.linalg.eigvalsh(M)
    G0 = sym(G0)
    tan = np.linalg.slogdet(G0)[1] - r + np.trace(np.linalg.inv(G0) @ G)
    slack = np.linalg.slogdet(M)[1] - (tan - 2.0 * I0 * ln2)
    I_real = 0.5 * (np.linalg.slogdet(G)[1] - np.linalg.slogdet(M)[1]) / ln2
    eigK_full = np.linalg.eigvalsh(sym(K))
    rankK = int((eigK_full > 1e-7 * max(1.0, eigK_full.max())).sum())
    Q = np.eye(n) - Z @ Z.T
    Pm = sym(X - K)
    try:
        Sig = sym(np.linalg.inv(Pm) - np.linalg.inv(X))
    except np.linalg.LinAlgError:
        Sig = np.full((n, n), np.nan)
    sig_off = float(np.max(np.abs(Q @ Sig @ Q)))
    sig_cross = float(np.max(np.abs(Q @ Sig @ Z)))
    sig_min = float(np.linalg.eigvalsh(Sig).min())
    return dict(r_lin=float(res_lin), min_eig=dict(X=float(eigX.min()), K=float(eigK.min()),
                XmK=float(eigXM.min()), M=float(eigM.min())),
                lamG=[float(v) for v in np.linalg.eigvalsh(G)],
                I=float(I_real), tan=float(tan), logdetM=float(np.linalg.slogdet(M)[1]),
                slack_e70=float(slack), rankK=rankK,
                sig=(sig_off, sig_cross, sig_min),
                D=float(np.trace(TH @ (X - K)) + JC))


print('== E70b：方向判定（松弛 vs 内逼近）==')
Ts, U = schur(A, output='real', sort=lambda a: abs(a) < 1.0)[:2]
PLANES = {'free_R4': U[:, :], 'blue_tail3': U[:, n - 3:]}
REACH = {'free_R4': 42.0427, 'blue_tail3': 45.4537}
t0 = time.time()

for name, Z in PLANES.items():
    r = Z.shape[1]
    print('\n--- %s (r=%d) ---' % (name, r))
    X, K, D, I, s, conv = scale_for_rate(Z, 3.0)
    print('  [T1] 标量设计 s=%.6g：D=%.5f  I=%.6f  收敛=%s（I<=3 ⇒ 它是预算集里的合法物理设计）'
          % (s, D, I, conv))
    for gname, G0 in [('I', np.eye(r)), ('3I', 3.0 * np.eye(r)),
                      ('G_phys', sym(Z.T @ X @ Z))]:
        d = e70_set_value(Z, X, K, 3.0, G0)
        verdict = 'E70可行' if d['slack_e70'] >= -1e-9 else 'E70不可行'
        print('   G0=%-8s: slack(约束左-右)=%+12.6f  %s | logdetM=%.5f tan=%.5f | 率=%.5f D=%.5f'
              % (gname, d['slack_e70'], verdict, d['logdetM'], d['tan'], d['I'], d['D']))
        print('        逐项：ARE残差=%.2e  最小特征值 X=%.4f K=%.4f X-K=%.4f M=%.4f  lam(G)=%s'
              % (d['r_lin'], d['min_eig']['X'], d['min_eig']['K'], d['min_eig']['XmK'],
                 d['min_eig']['M'], np.round(d['lamG'], 4)))
    print('   预言：slack(G0=G_phys) = 2(I0-I)ln2 = %.6f；若 I=3 则应 → 0'
          % (2.0 * (3.0 - I) * ln2))
    d0 = e70_set_value(Z, X, K, 3.0, np.eye(r))
    print('   物理性核对 Σ:=P⁻¹-X⁻¹（物理 ⟺ Σ⪰0 且 V⊥块=交叉块=0）：'
          'V⊥块=%.2e 交叉块=%.2e 最小本征=%.4e rankK=%d'
          % (d0['sig'][0], d0['sig'][1], d0['sig'][2], d0['rankK']))

print('\n[T2] 切线点扫掠：内逼近的 min 只能 >= 已知可达值；若出现更低值则方向判定被推翻')


def sdp_e70(Z, I0, G0):
    r = Z.shape[1]
    X = cp.Variable((n, n), symmetric=True)
    K = cp.Variable((n, n), symmetric=True)
    IZ = Z.T @ X @ Z
    M = IZ - Z.T @ K @ Z
    G0 = sym(G0)
    tan = np.linalg.slogdet(G0)[1] - r + cp.trace(np.linalg.inv(G0) @ IZ)
    cons = [X - A @ X @ A.T + A @ K @ A.T == W, X >> 0, K >> 0, X - K >> 0, M >> 0,
            cp.log_det(M) >= tan - 2.0 * I0 * ln2]
    p = cp.Problem(cp.Minimize(JC + cp.trace(TH @ (X - K))), cons)
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        p.solve(solver='CLARABEL', verbose=False)
    return p, (None if X.value is None else (np.array(X.value), np.array(K.value)))


for name, Z in PLANES.items():
    X, K, D, I, s, conv = scale_for_rate(Z, 3.0)
    Gp = sym(Z.T @ X @ Z)
    vals = []
    for gname, G0 in [('I', np.eye(Gp.shape[0])), ('G_phys', Gp), ('G_phys/3', Gp / 3.0),
                      ('3*G_phys', 3.0 * Gp), ('10*G_phys', 10.0 * Gp),
                      ('G_phys/10', Gp / 10.0), ('G_phys/100', Gp / 100.0)]:
        p, sol = sdp_e70(Z, 3.0, G0)
        ok = p.status in ('optimal', 'optimal_inaccurate')
        if ok:
            vals.append(p.value)
        extra = ''
        if sol is not None:
            d = e70_set_value(Z, sol[0], sol[1], 3.0, G0)
            extra = (' | 解处 率=%.5f D=%.5f ARE残差=%.1e | Σ: V⊥块=%.3f 交叉=%.3f 最小本征=%.3f rankK=%d'
                     % (d['I'], d['D'], d['r_lin'], d['sig'][0], d['sig'][1], d['sig'][2], d['rankK']))
        print('  %-11s G0=%-11s: status=%-20s L=%s%s' % (
            name, gname, p.status,
            ('%10.5f' % p.value) if ok else '   --', extra))
    if vals:
        print('  %-11s min_gamma L = %.5f  vs 已知可达 %.5f  -> %s'
              % (name, min(vals), REACH[name],
                 '从未跌破（内逼近，方向判定成立）' if min(vals) >= REACH[name] - 1e-3
                 else '跌破可达值：它连上界都不是，只是无意义数'))
print('\n已知可达值（对照，E70 若为松弛则 min L 必须 <= 它们）：free_R4 42.0427, blue_tail3 45.4537；'
      '\nRemark 2 需要的认证下界方向是 L > 43.80（蓝平面）。耗时 %.1f s' % (time.time() - t0))
