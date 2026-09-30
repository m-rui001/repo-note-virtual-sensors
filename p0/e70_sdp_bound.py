r"""E70 = 锥前沿的**认证下界**（SDP + log-det），目标：把 Remark 2 的"相对上界排除"升级成绝对不可行。

== 推导（写下来，便于 C 复核）==
定点写成"先验 X = A P Aᵗ + W，后验 P"。设计变量是 V 上的信息 M ⪰ 0；令 Z 为 V 的正交基，
Y := (M⁻¹ + ZᵗXZ)⁻¹（V 上的后验协方差），G := ZᵗXZ。则
    P = X − (XZ) Y (ZᵗX),      X = A X Aᵗ + W − (AXZ) Y (ZᵗX A),
    I = −½ log det(I_r − G^{1/2} Y G^{1/2}) = ½[log det G − log det(G − ZᵗKZ)],
其中最后一个等号用了 ZᵗKZ = G Y G，K := (XZ)Y(ZᵗX)。

**线性化**：把 K 当独立变量，Riccati 变成**线性等式** X − AXAᵗ + AKAᵗ = W，代价是
j_c + tr(Θ(X−K))，率是 ½[logdet G − logdet(G − ZᵗKZ)]。

**松弛**：物理设计集是 {K = (XZ)Y(ZᵗX), Y ⪰ 0}；丢掉"K 必须这样分解"（即丢掉
rank K ≤ r 与 range K ⊆ range(XZ)）得到**包含全部物理设计**的凸集 ⇒ 最小值是真最小值的下界。
率约束 logdet G − logdet(G−ZᵗKZ) ≥ 2 I₀ ln2 是"凹 − 凸 ≥ 常数"，非凸；用
**凹函数 logdet G 的切线**（在任意 G₀ ≻ 0 处，切线是全局上界）替换：
    logdet(G−ZᵗKZ) − [logdet G₀ + tr(G₀⁻¹(G−G₀))] ≤ −2 I₀ ln2
这是凸约束，且因为切线 ≥ logdet G，原可行点仍可行 ⇒ **仍是松弛**。
每个 G₀ 给出一个有效下界 L(G₀)；对 G₀ 取最大（切线割平面式迭代）。

== 自检（判据先写下）==
 [0] V = R⁴（r=n）时松弛**是精确的**：K = XYX ⇔ Y = X⁻¹KX⁻¹ ⪰ 0 自动成立，
     所以 SDP 值应当等于已知的无约束锥最优 42.0427 (I=3)。不等就说明实现有 bug。
 [1] V = Schur-tail-3 平面、I₀=3：物理设计（Powell 的 45.4537）必须**可行**
     （否则松弛丢了约束），且下界 L 越接近 43.80 越好。
 [2] 验收线：L > 43.80 ⇒ Remark 2 升级为绝对不可行；L ≤ 43.711 ⇒ 蓝线在该平面内
     可能可行，Remark 2 必须再降一档措辞。
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
EPSR = 1e-11


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
PLANE = {'free_R4': U[:, :], 'schur_tail_3': U[:, n - 3:], 'schur_tail_2': U[:, n - 2:]}


def sdp_lower(Z, I0, G0=None, solver='CLARABEL'):
    """返回 (下界 L, X*, K*)。G0 = 切线展开点（r x r 正定）；None 时用 I_r。"""
    r = Z.shape[1]
    X = cp.Variable((n, n), symmetric=True)
    K = cp.Variable((n, n), symmetric=True)
    IZ = Z.T @ X @ Z
    IKZ = Z.T @ K @ Z
    M = IZ - IKZ
    cons = [X - A @ X @ A.T + A @ K @ A.T == W,
            X >> 0, K >> 0, X - K >> 0, M >> 0]
    if G0 is None:
        G0 = np.eye(r)
    G0 = sym(G0)
    G0i = np.linalg.inv(G0)
    # logdet G <= logdet G0 + tr(G0^-1 (G - G0))  (tangent = global upper bound of a concave f)
    # so the true feasible set implies  logdet(G - Z'KZ) >= tangent - 2 I0 ln2,
    # i.e.  logdet(M) <= tangent - 2 I0 ln2  is NOT the relaxation; we need the other side.
    tan = np.linalg.slogdet(G0)[1] - r + cp.trace(G0i @ IZ)
    cons += [cp.log_det(M) >= tan - 2.0 * I0 * ln2]
    obj = cp.Minimize(JC + cp.trace(TH @ (X - K)))
    p = cp.Problem(obj, cons)
    p.solve(solver=solver, verbose=False)
    if p.status not in ('optimal', 'optimal_inaccurate'):
        return np.nan, None, None
    return float(p.value), np.array(X.value), np.array(K.value)


def phys_XK(S):
    """把一个物理设计 S 变成 (X, K)，用于 [1] 的可行性自检。"""
    J, I = JI(S)
    if not np.isfinite(J):
        return None
    Pt = W.copy()
    for _ in range(20000):
        Pm = sym(np.linalg.inv(np.linalg.inv(Pt) + S))
        Pn = sym(A @ Pm @ A.T + W)
        if np.max(np.abs(Pn - Pt)) < 1e-12 * max(1.0, np.max(np.abs(Pn))):
            break
        Pt = Pn
    X = Pt
    P = Pm
    return X, X - P, J, I


def check_phys(Z, X, K, I0):
    G = sym(Z.T @ X @ Z)
    M = G - sym(Z.T @ K @ Z)
    r = np.linalg.matrix_rank(sym(K), tol=1e-7)
    I = 0.5 * (np.linalg.slogdet(G)[1] - np.linalg.slogdet(M)[1]) / ln2
    res = X - A @ X @ A.T + A @ K @ A.T - W
    return dict(rate=float(I), are_res=float(np.max(np.abs(res))),
                lam_G=[float(v) for v in np.linalg.eigvalsh(G)],
                lam_M=[float(v) for v in np.linalg.eigvalsh(M)],
                rankK=int(r), D=float(JC + np.trace(TH @ (X - K))))


print('== E70：锥前沿的 SDP 认证下界 ==')
print('\n[0] 自检：V = R^4 时松弛应当精确，对照无约束锥最优 42.0427 (I=3)')
t0 = time.time()
for I0 in (3.0,):
    for it, G0 in enumerate([None, sym(np.eye(4) * 3.0), None]):
        L, Xs, Ks = sdp_lower(PLANE['free_R4'], I0, G0)
        print('   I0=%.1f 切线点=%s: L=%10.5f  (%s, %.1f s)'
              % (I0, 'I' if G0 is None else 'G0@prev', L,
                 'optimal' if np.isfinite(L) else 'fail', time.time() - t0))
        if np.isfinite(L) and Xs is not None:
            c = check_phys(PLANE['free_R4'], Xs, Ks, I0)
            print('      SDP 解处：率=%.4f (>= %.1f?) ARE残差=%.2e rankK=%d D=%.5f'
                  % (c['rate'], I0, c['are_res'], c['rankK'], c['D']))
