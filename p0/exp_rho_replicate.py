"""E16：给第 7 条 §2 的两条加固。
   (a) ρ 的跨求解器稳定性：同一无约束 SDP 用 CLARABEL / ECOS / SCS 三个后端各解一次，
       若 ρ 在三个后端下一致非零，它就不是 CLARABEL 那条 "Solution may be inaccurate" 警告的伪影。
   (b) Δ̂ 的可重复性：紧约束求解换 4 个随机种子各跑 10 起点，报 Δ̂ 的散布。
   注意逻辑方向：minimization 里我返回的 I_TRV 是**上界**（可行点），所以 Δ̂>0 本身不证 Δ>0；
   真正证 Δ>0 需要 ρ>0 加上无约束最优唯一性。这个脚本只做"排除我的实现出错"这一层。
"""
import numpy as np
import cvxpy as cp
from scipy.optimize import minimize
from p0.p0_replicate_letter import A, B, W, n, sym, ln2, rate_cost, init_gamma
from p0.exp_authoritative import ctrl_at
from p0.exp_thm3_qhypothesis import F2, Z2, align_family
from p0.exp_epsilon_tight import tight_face, unpack3, ij

Pc0, K0, Th0 = ctrl_at(A, B, np.eye(n), np.eye(B.shape[1]))


def sdp_S(D, Q, Th, Pc, solver):
    """Tanaka 式(18) 的单后端求解，返回 S*=P^-1-Pt^-1 与求解器状态"""
    P = cp.Variable((n, n), symmetric=True)
    Pi = cp.Variable((n, n), symmetric=True)
    const = 0.5 * np.log(np.linalg.det(W)) / ln2
    cons = [Pi >> 0, P >> 0,
            cp.trace(Th @ P) + float(np.trace(W @ Pc)) <= D,
            P << A @ P @ A.T + W,
            cp.bmat([[P - Pi, P @ A.T], [A @ P, A @ P @ A.T + W]]) >> 0]
    prob = cp.Problem(cp.Minimize(-0.5 * cp.log_det(Pi) / ln2 + const), cons)
    kw = {}
    if solver == cp.SCS:
        kw = dict(eps=1e-11, max_iters=40000)
    elif solver == cp.ECOS:
        kw = dict(abstol=1e-11, reltol=1e-11, feastol=1e-11, maxit=400)
    prob.solve(solver=solver, gp=False, **kw)
    if P.value is None:
        return None
    Pv, Piv = sym(P.value), sym(Pi.value)
    Pt = A @ Pv @ A.T + W
    S = np.linalg.inv(Pv) - np.linalg.inv(Pt)
    I = 0.5 * np.log(np.linalg.det(Pt) / np.linalg.det(Pv)) / ln2
    r = float(np.linalg.norm(Z2.T @ S @ Z2) / np.linalg.norm(S))
    return dict(I=I, rho=r, status=prob.status,
                lead=float(np.linalg.norm(Z2.T @ S @ Z2, 'fro')))


if __name__ == '__main__':
    D0 = 80.0
    print('==== (a) ρ 跨求解器 ====')
    for eps in (0.0, 1e-4, 0.01, 0.1, 1.0):
        Q = align_family(F2, Z2, eps)
        Pc, K, Th = ctrl_at(A, B, Q, np.eye(B.shape[1]))
        outs = []
        for sv in (cp.CLARABEL, cp.ECOS, cp.SCS):
            try:
                o = sdp_S(D0, Q, Th, Pc, sv)
            except Exception as e:
                o = None
                print('   %s 抛异常 %s' % (str(sv)[:8], str(e)[:50]))
            outs.append((str(sv).split('.')[-1], o))
        s = '  ε=%-7g' % eps
        for name, o in outs:
            s += '  %s: ρ=%s I=%s [%s]' % (
                name[:7], '—' if o is None else '%.3e' % o['rho'],
                '—' if o is None else '%.6f' % o['I'],
                '—' if o is None else o['status'][:6])
        print(s)

    print('\n==== (b) Δ̂ 的可重复性（4 组随机起点 × 硬约束 J≤D(1+1e-9)）====')
    for eps in (0.0, 0.01, 0.1):
        Q = align_family(F2, Z2, eps)
        Pc, K, Th = ctrl_at(A, B, Q, np.eye(B.shape[1]))
        u = sdp_S(D0, Q, Th, Pc, cp.CLARABEL)
        vals = []
        for seed in range(4):
            rng = np.random.default_rng(seed)
            best = None
            cands = [np.eye(2), np.eye(2) * 20.0, np.diag([5.0, 0.5]), np.diag([0.5, 5.0])]
            cands += [sym(v @ v.T) + 0.05 * np.eye(2) for v in rng.standard_normal((6, 2, 2))]
            for G in cands:
                try:
                    L = np.linalg.cholesky(sym(G) + 1e-10 * np.eye(2))
                    p = np.array([L[0, 0], L[1, 0], L[1, 1]])
                except np.linalg.LinAlgError:
                    continue
                res = minimize(lambda q: ij(F2, q, Th, Pc)[0], p, method='SLSQP',
                               constraints=[dict(type='ineq', fun=lambda q: D0 - ij(F2, q, Th, Pc)[1])],
                               options=dict(maxiter=400, ftol=1e-13))
                I, J = ij(F2, res.x, Th, Pc)
                if J <= D0 * (1 + 1e-9) and (best is None or I < best):
                    best = I
            vals.append(best)
        vv = [v for v in vals if v is not None]
        print('  ε=%-7g  I_unc=%.6f   各起点最好 I_TRV: %s   Δ̂ 散布=%.2e bits' %
              (eps, u['I'], ' '.join('%.6f' % v for v in vv),
               (max(vv) - min(vv)) if vv else float('nan')))
