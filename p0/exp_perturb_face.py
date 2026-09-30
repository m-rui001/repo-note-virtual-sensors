"""E21：把"无约束最优唯一"这个前提换成可检验的量。

逻辑缺口（#8 §1、#9 §1）：$\hat\Delta>0$ 只说明我的求解器没找到无损解。真正证 $\Delta>0$ 需要
    (i) $\rho(D;\varepsilon)>0$，且 (ii) 无约束最优解集 $\mathrm{Opt}(D)$ 里**每一个**元素都满足 $\rho>0$。
原文 Prop.2 的 "$\Delta=0\iff\exists S^\ast\in\mathrm{Opt}(D)\cap\mathcal K_F$" 是存在量词，所以 (ii) 不能省。

不证唯一性，改测**扰动稳定性**：往目标加一个随机小项 $\eta\,\mathrm{tr}(\Sigma P)$（$\Sigma\succ0$，$\eta=10^{-3}/10^{-2}$），
这等价于在最优面上选点——若 $\mathrm{Opt}(D)$ 整个面上有 $\rho$ 任意小的点，则Generic 扰动会把解推向它，
$\rho$ 应当随 $\eta$ 掉到噪声底。若 $\rho$ 对 $\eta$ 稳健（保持在无扰值的同数量级），那么"$\exists S^\ast:\rho(S^\ast)=0$"
与数据不相容，$\Delta>0$ 就从"求解器没找到"升级成"整个最优面上都成立"。

变量与目标同 `unconstrained_sdp`（Tanaka 式18 的 $(P,\Pi)$ 提升形式，$-\tfrac12\log\det\Pi$ 严格凸）。
"""
import sys
import numpy as np
import cvxpy as cp
sys.stdout.reconfigure(encoding='utf-8')
from p0.p0_replicate_letter import A, W, n, sym, ln2, ctrl, const2
from scipy.linalg import schur

Pc, K, Th = ctrl(np.eye(n))
T, U, sdim = schur(A, output='real', sort=lambda a: abs(a) < 1.0)
nu = n - sdim
iT = list(range(n - nu, n))
F2 = U[:, iT].T
Z = np.ascontiguousarray(U[:, :sdim])           # ker F2（稳定不变块）
r2 = F2.shape[0]


def sdp_S(D, Qmat, Theta, Pcf, pert=None, eta=0.0):
    P = cp.Variable((n, n), symmetric=True)
    Pi = cp.Variable((n, n), symmetric=True)
    c = 0.5 * np.log(np.linalg.det(W)) / ln2
    obj = -0.5 * cp.log_det(Pi) / ln2 + c
    if pert is not None:
        obj = obj + eta * cp.trace(pert @ P) / ln2
    cons = [Pi >> 0, P >> 0, cp.trace(Theta @ P) + const2(Qmat, Pcf) <= D,
            P << A @ P @ A.T + W,
            cp.bmat([[P - Pi, P @ A.T], [A @ P, A @ P @ A.T + W]]) >> 0]
    prob = cp.Problem(cp.Minimize(obj), cons)
    prob.solve(solver=cp.CLARABEL, gp=False)
    if P.value is None:
        return None, None, prob.status
    Pv, Piv = sym(P.value), sym(Pi.value)
    Pt = sym(A @ Pv @ A.T + W)
    S = np.linalg.inv(Pv) - np.linalg.inv(Pt)
    I = 0.5 * np.log(np.linalg.det(Pt) / np.linalg.det(Pv)) / ln2
    return S, I, prob.status


def rho(S):
    return float(np.linalg.norm(Z.T @ S @ Z, 'fro') / max(np.linalg.norm(S, 'fro'), 1e-300))


if __name__ == '__main__':
    from p0.p0_replicate_letter import B
    from p0.exp_authoritative import ctrl_at
    rng = np.random.default_rng(5)
    D = 80.0
    print('F2 尾部2，D=%.0f；ρ = ‖ZᵀS*Z‖_F/‖S*‖_F，Z 张 ker F2' % D)
    for eps in (0.1, 1.0):
        Q = sym(F2.T @ F2 + eps * (Z @ Z.T))
        Pce, Ke, The = ctrl_at(A, B, Q, np.eye(n))
        S0, I0, st = sdp_S(D, Q, The, Pce, None, 0.0)
        r0 = rho(S0)
        print('\n  ε=%.2f  基线：%s  I=%.6f  ρ=%.4e' % (eps, st, I0, r0))
        print('%8s %10s %34s' % ('η', 'ρ/基线 范围', '六个随机 Σ 的 ρ'))
        for eta in (1e-3, 1e-2, 1e-1):
            vals = []
            for _ in range(6):
                L = rng.standard_normal((n, n))
                Sig = sym(L @ L.T) + 0.05 * np.eye(n)
                S, I, stt = sdp_S(D, Q, The, Pce, Sig, eta)
                vals.append(rho(S) if S is not None else np.nan)
            v = np.array(vals)
            print('%8.0e %10s %34s'
                  % (eta, '[%.2f, %.2f]×' % (np.nanmin(v) / r0, np.nanmax(v) / r0),
                     ' '.join('%.2e' % x for x in v)))

