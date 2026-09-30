#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
================================================================================
P0 判伪实验：同等信息预算下，"表示对齐的目标"是否改变决策代价？
================================================================================

问题
----
信息论表示学习（IB / 信息最大化 / 方差保留 / PCA）默认"保留方差大的方向"
（按 Sigma 的特征结构分配精度），因为那在 MSE / trace 意义下最优。
本实验检验：在**同等信息预算 C（bits/sample）**下，改按**决策敏感度矩阵 Pi**
的特征结构分配，无限时域 LQG 代价是否显著更低；差距是否随两者**失配角**增长。

理论骨架（本脚本逐条数值验证）
------------------------------
1) 确定性等价 + **每样本无记忆高斯测试信道**（"每样本信息预算"的标准理想化）
       x_t = x_hat_t + e_t,   x_hat_t ⟂ e_t,   e_t ~ N(0, Sigma_e) 白噪声
       u_t = -K x_hat_t
   闭环   x_{t+1} = T x_t + G e_t + w_t,     T = A - B K,  G = B K
   可得（本脚本用 20 万步蒙特卡洛验证）
       **J = tr(P W) + tr(Pi Sigma_e)**,   **Pi = K^T S K**,   S = R + B^T P B
   其中 P 为 LQR 的 cost-to-go，K 为 LQR 增益，W 为过程噪声协方差。
   Pi = K^T S K ⪰ 0 恒成立，且 rank(Pi) ≤ rank(B)：
   估计误差**只通过控制动作花钱**，某些方向的误差完全免费。
   （注：我最初猜的 Pi = P - K^T R K 是错的；蒙特卡洛只支持 Pi = K^T S K。）

2) 高斯率–失真（矩阵型/Loewner 型保真度）：
       可行集 {Sigma_e : 0 ⪯ Sigma_e ⪯ Sigma}
       率    R(Sigma_e) = 1/2 log2 det(Sigma)/det(Sigma_e)
   由测试信道 x = x_hat + e 且 x_hat ⟂ e 达成（与 2510.13025、经典率失真一致）。

3) 无约束最优误差形状（AM-GM）：
       min tr(Pi Sigma_e)  s.t.  det Sigma_e = det(Sigma) 2^{-2C}
       = d (det Pi det Sigma)^{1/d} 2^{-2C/d},   取等 ⟺ **Sigma_e ∝ Pi^{-1}**
   → 应按 **pi_i * sigma_i^2（敏感度 × 方差）** 分配，而不是按 sigma_i^2。
   **注意**：该形状常**不可达**（违反 Sigma_e ⪯ Sigma）。此时最优解落在边界上，
   表现为"整块放弃某些方向"。本脚本用参数化优化器求真实最优（见 eps_oracle）。

4) 基线（PCA / 能量 / 逆水填）：限制 Sigma_e 在 Sigma 特征基下对角。
   d=2 的闭式差距比（无截断时）：
       gap = sqrt(1 + sin^2(2 theta)/4 * (sqrt(rho) - 1/sqrt(rho))^2),  rho = pi_1/pi_2

判伪规则（预先写死）
-------------------
   (ii) 显著优于 (i) 且数值与公式同阶  -> 方向存活，进 P1
   (i) ≈ (ii)                          -> 方向被证伪，立即停手，写负结果

运行   cd p0 && python p0_experiment.py
输出   results/summary.txt, results/*.csv, results/*.png
================================================================================
"""
from __future__ import annotations

import csv
import math
import os

import numpy as np
from scipy.linalg import expm
from scipy.optimize import minimize

np.set_printoptions(precision=6, suppress=True, linewidth=170)
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "results")
os.makedirs(OUT, exist_ok=True)
LOG2 = math.log(2.0)

REPORT: list[str] = []


def say(*a):
    s = " ".join(str(x) for x in a)
    print(s)
    REPORT.append(s)


# ============================================================ 基础线性代数
def sym(M):
    return 0.5 * (M + M.T)


def logdet(M):
    s, ld = np.linalg.slogdet(sym(M))
    return ld if s > 0 else -np.inf


def rate_bits(Sigma, Sigma_e):
    """1/2 log2 det(Sigma)/det(Sigma_e)"""
    return 0.5 * (logdet(Sigma) - logdet(Sigma_e)) / LOG2


def min_eig(M):
    return float(np.linalg.eigvalsh(sym(M)).min())


def is_psd(M, tol=1e-9):
    return min_eig(M) >= -tol * max(1.0, float(np.abs(M).max()))


def psd_sqrt(M):
    w, V = np.linalg.eigh(sym(M))
    return V @ np.diag(np.sqrt(np.maximum(w, 0.0))) @ V.T


def loewner_le(A, B, tol=1e-7):
    return min_eig(sym(B - A)) >= -tol * max(1.0, float(np.abs(B).max()))


def dare(A, B, Q, R, iters=200000, tol=1e-15):
    P = Q.copy()
    for k in range(iters):
        S = R + B.T @ P @ B
        Pn = sym(Q + A.T @ P @ A - A.T @ P @ B @ np.linalg.solve(S, B.T @ P @ A))
        if np.max(np.abs(Pn - P)) < tol * max(1.0, np.max(np.abs(Pn))):
            return Pn, k
        P = Pn
    return P, iters


def lyap_d(Tmat, Qn):
    d = Tmat.shape[0]
    return np.linalg.solve(np.eye(d * d) - np.kron(Tmat, Tmat),
                           Qn.reshape(-1, order="F")).reshape((d, d), order="F")


def lqr(A, B, Q, R):
    P, _ = dare(A, B, Q, R)
    S = R + B.T @ P @ B
    K = np.linalg.solve(S, B.T @ P @ A)
    return K, P, S


def pi_matrix(K, S):
    """决策敏感度矩阵 Pi = K^T S K（保证 PSD，rank ≤ rank(B)）"""
    return sym(K.T @ S @ K)


def principal_angles_deg(U1, U2):
    s = np.linalg.svd(U1.T @ U2, compute_uv=False)
    return np.degrees(np.arccos(np.clip(s, -1.0, 1.0)))


def top_k_subspace(M, k):
    w, V = np.linalg.eigh(sym(M))
    return V[:, np.argsort(w)[::-1][:k]]


def dJ(Pi, Sigma_e):
    """ΔJ = tr(Pi Sigma_e)"""
    return float(np.trace(Pi @ Sigma_e))


# ============================================================ 误差协方差方案
def eps_analytic(Sigma, Pi, C):
    """AM-GM 无约束最优 Sigma_e = c Pi^{-1}；返回 (Sigma_e, 可行?, 说明)"""
    d = Sigma.shape[0]
    if not is_psd(Pi) or min_eig(Pi) < 1e-12 * max(1.0, np.abs(Pi).max()):
        return None, False, "Pi 奇异"
    c = math.exp((logdet(Pi) + logdet(Sigma)) / d) * 2.0 ** (-2.0 * C / d)
    Se = c * np.linalg.inv(Pi)
    return sym(Se), loewner_le(Se, Sigma), f"c={c:.6g}"


def amgm_bound(Sigma, Pi, C):
    """d (det Pi det Sigma)^{1/d} 2^{-2C/d}（可达当且仅当 c Pi^{-1} ⪯ Sigma）"""
    d = Sigma.shape[0]
    if logdet(Pi) == -np.inf:
        return 0.0
    return d * math.exp((logdet(Pi) + logdet(Sigma)) / d) * 2.0 ** (-2.0 * C / d)


def eps_oracle(Sigma, Pi, C, n_starts=6, seed=0, warm=None):
    """
    等率下的**真实最优**误差协方差：
        min tr(Pi Sigma_e)  s.t.  det Sigma_e = det(Sigma) 2^{-2C},  0 ⪯ Sigma_e ⪯ Sigma
    参数化保证两条约束**精确成立**：
        Sigma_e = Sigma^{1/2} expm(S) Sigma^{1/2},
        S = -U diag(s) U', U = expm(skew(v)), s = 2C ln2 * softmax(z)  (s_i>0, tr S = -2C ln2)
        => det Sigma_e = det Sigma * exp(tr S) = det Sigma 2^{-2C} 且 S ⪯ 0  => Sigma_e ⪯ Sigma
    """
    d = Sigma.shape[0]
    sq = psd_sqrt(Sigma)
    npair = d * (d - 1) // 2
    trg = 2.0 * C * LOG2

    def build(p):
        v, z = p[:npair], p[npair:]
        Sk = np.zeros((d, d))
        idx = 0
        for i in range(d):
            for j in range(i + 1, d):
                Sk[i, j], Sk[j, i] = v[idx], -v[idx]
                idx += 1
        U = expm(Sk)
        zz = np.concatenate([z, [0.0]])
        zz = zz - zz.max()
        w = np.exp(zz)
        s = trg * w / w.sum()
        S = -U @ np.diag(s) @ U.T
        return sym(sq @ expm(S) @ sq), S

    def obj(p):
        return float(np.trace(Pi @ build(p)[0]))

    rng = np.random.default_rng(seed)
    starts = [np.zeros(npair + d - 1)]
    if warm is not None:
        starts.insert(0, warm)
    for _ in range(n_starts):
        starts.append(rng.standard_normal(npair + d - 1))
    best, bx = np.inf, None
    for s0 in starts:
        r = minimize(obj, s0, method="Nelder-Mead",
                     options=dict(maxiter=6000, maxfev=12000, xatol=1e-11, fatol=1e-15))
        if r.fun < best:
            best, bx = r.fun, r.x
    for _ in range(2):
        r = minimize(obj, bx, method="Nelder-Mead",
                     options=dict(maxiter=20000, maxfev=40000, xatol=1e-13, fatol=1e-17))
        if r.fun < best:
            best, bx = r.fun, r.x
    Se, S = build(bx)
    return Se, best, bx


def eps_diag_in_sigma_basis(Sigma, Pi, C, mode="energy"):
    """
    Sigma_e 限制在 Sigma 特征基下对角：Sigma_e = V diag(D) V'。
    在 prod(D_i)=prod(sigma_i)2^{-2C}, 0<D_i<=sigma_i 下最小化 sum w_i D_i。
      mode='energy' -> w_i = 1          （能量/PCA/逆水填，MSE 最优）
      mode='pidiag' -> w_i = (V'PiV)_ii  （在 Sigma 基下仅重加权）
    """
    d = Sigma.shape[0]
    s, V = np.linalg.eigh(sym(Sigma))
    w = np.ones(d) if mode == "energy" else np.maximum(np.diag(V.T @ Pi @ V), 1e-300)
    active = list(range(d))
    D = s.copy()
    for _ in range(d + 1):
        if not active:
            break
        a = np.array(active)
        Kdet = float(np.prod(s[a])) * 2.0 ** (-2.0 * C)
        nu = (Kdet * float(np.prod(w[a]))) ** (1.0 / len(a))
        Da = nu / w[a]
        viol = Da > s[a]
        if not viol.any():
            D[a] = Da
            break
        bad = a[viol]
        D[bad] = s[bad]              # 该方向不传信息（误差=全方差，0 比特）
        active = [i for i in active if i not in set(bad.tolist())]
    return sym(V @ np.diag(D) @ V.T)


def eps_random(Sigma, Pi, C, rng, spread=0.35):
    """随机表示（下界参照）：Sigma_e = Sigma^{1/2} M Sigma^{1/2}, prod eig(M)=2^{-2C}, M ⪯ I"""
    d = Sigma.shape[0]
    Qm, _ = np.linalg.qr(rng.standard_normal((d, d)))
    base = -2.0 * C * LOG2 / d
    delta = rng.standard_normal(d) * spread * abs(base)
    delta -= delta.mean()
    m = np.exp(np.clip(base + delta, -700, 0.0))
    sq = psd_sqrt(Sigma)
    return sym(sq @ sym(Qm @ np.diag(m) @ Qm.T) @ sq)


# ============================================================ 校验 1：LQG 公式
def test_channel_sample(x, Sigma, Sigma_e, rng):
    """从 p(x_hat|x) 采样：E[x_hat|x]=(Sigma-Se)Sigma^{-1}x, cov=(Sigma-Se)Sigma^{-1}Se"""
    d = Sigma.shape[0]
    A0 = Sigma - Sigma_e
    L = psd_sqrt(sym(A0 @ np.linalg.solve(Sigma, Sigma_e)))
    return A0 @ np.linalg.solve(Sigma, x) + L @ rng.standard_normal(d)


def validate_formula(A, B, Q, R, W, Sigma_e, T_steps=150000, burn=3000, seed=0):
    K, P, S = lqr(A, B, Q, R)
    Pi = pi_matrix(K, S)
    T, G = A - B @ K, B @ K
    d = A.shape[0]
    Meff = sym(T @ Sigma_e @ G.T + G @ Sigma_e @ T.T + G @ Sigma_e @ G.T)
    Sig_x = lyap_d(T, Meff + W)
    J_lyap = float(np.trace((Q + K.T @ R @ K) @ Sig_x) - np.trace(K.T @ R @ K @ Sigma_e))
    J_formula = float(np.trace(P @ W)) + dJ(Pi, Sigma_e)
    rng = np.random.default_rng(seed)
    Sigma = lyap_d(T, W)
    x = psd_sqrt(Sigma) @ rng.standard_normal(d)
    acc = acc2 = 0.0
    n = 0
    for t in range(T_steps):
        xh = test_channel_sample(x, Sigma, Sigma_e, rng)
        e = x - xh
        u = -K @ xh
        c = float(x @ Q @ x + u @ R @ u)
        if t >= burn:
            acc += c
            acc2 += c * c
            n += 1
        x = T @ x + G @ e + rng.standard_normal(d) @ psd_sqrt(W).T
    J_mc = acc / n
    se = math.sqrt(max(acc2 / n - J_mc ** 2, 0.0) / n)
    return dict(J_mc=J_mc, se=se, J_formula=J_formula, J_lyap=J_lyap, n=n, Pi=Pi, P=P, K=K, S=S)


# ============================================================ 校验 2：d=2 闭式律
def gap_law_d2(rho, theta_deg):
    th = math.radians(theta_deg)
    return math.sqrt(1.0 + math.sin(2 * th) ** 2 / 4.0 * (math.sqrt(rho) - 1 / math.sqrt(rho)) ** 2)


def sweep_d2(C=5.0):
    rows = []
    Sigma = np.diag([4.0, 1.0])
    for rho in (4.0, 16.0, 64.0):
        for th in np.linspace(0, 90, 19):
            t = math.radians(th)
            Rm = np.array([[math.cos(t), -math.sin(t)], [math.sin(t), math.cos(t)]])
            Pi = Rm @ np.diag([rho, 1.0]) @ Rm.T
            Se_g, feas, _ = eps_analytic(Sigma, Pi, C)
            if not feas:
                continue
            Se_e = eps_diag_in_sigma_basis(Sigma, Pi, C, "energy")
            rows.append(dict(rho=rho, theta=th, gap_num=dJ(Pi, Se_e) / dJ(Pi, Se_g),
                             gap_closed=gap_law_d2(rho, th),
                             r_g=rate_bits(Sigma, Se_g), r_e=rate_bits(Sigma, Se_e)))
    return rows


# ============================================================ 物理系统
def plant():
    """4 维 4 输入（全驱动）线性系统；代价 Q 集中在方向 1,3。"""
    A = np.array([[0.95, 0.20, 0.00, 0.00],
                  [0.00, 0.90, 0.15, 0.00],
                  [0.00, 0.00, 0.85, 0.25],
                  [0.10, 0.00, 0.00, 0.92]])
    B = np.array([[1.0, 0.10, 0.00, 0.00],
                  [0.20, 1.0, 0.10, 0.00],
                  [0.00, 0.10, 1.0, 0.20],
                  [0.00, 0.00, 0.10, 1.0]])
    R = 0.5 * np.eye(4)
    Q = np.diag([10.0, 0.1, 10.0, 0.1])
    return A, B, Q, R


def rot4(phi_deg, i=0, j=1):
    Rm = np.eye(4)
    t = math.radians(phi_deg)
    Rm[i, i], Rm[i, j] = math.cos(t), -math.sin(t)
    Rm[j, i], Rm[j, j] = math.sin(t), math.cos(t)
    return Rm


def build_case(phi_deg, Q=None, B=None, Sigma0=None):
    """
    物理系统 + 指定"数据协方差"：
      T = A - B K（由 Q 定）固定；取目标 Sigma_target(phi)=R(phi) Sigma0 R(phi)'，
      W = PSD-clip(Sigma_target - T Sigma_target T')，再取 Sigma_actual = Lyap(T,W)。
    => 编码分布 Sigma_actual 与过程噪声 W 自洽，且其朝向随 phi 旋转。
    """
    A, B0, Q0, R = plant()
    B = B0 if B is None else B
    Q = Q0 if Q is None else Q
    K, P, S = lqr(A, B, Q, R)
    T, G = A - B @ K, B @ K
    Sigma0 = np.diag([4.0, 0.25, 4.0, 0.25]) if Sigma0 is None else Sigma0
    St = rot4(phi_deg) @ Sigma0 @ rot4(phi_deg).T
    W = sym(St - T @ St @ T.T)
    if not is_psd(W):
        W = W + (-min_eig(W) + 1e-8) * np.eye(A.shape[0])
    Sigma = lyap_d(T, W)
    Pi = pi_matrix(K, S)
    return dict(A=A, B=B, Q=Q, R=R, W=W, Sigma=Sigma, Pi=Pi, K=K, P=P, S=S, T=T, G=G)


def analyse(case, C, seed=0, warm=None, n_starts=6):
    Sigma, Pi = case["Sigma"], case["Pi"]
    d = Sigma.shape[0]
    Se_or, val_or, xvec = eps_oracle(Sigma, Pi, C, seed=seed, warm=warm, n_starts=n_starts)
    Se_an, feas, note = eps_analytic(Sigma, Pi, C)
    Se_e = eps_diag_in_sigma_basis(Sigma, Pi, C, "energy")
    Se_r = eps_diag_in_sigma_basis(Sigma, Pi, C, "pidiag")
    Se_rand = eps_random(Sigma, Pi, C, np.random.default_rng(seed + 991))
    k = max(1, d // 2)
    theta = principal_angles_deg(top_k_subspace(Sigma, k),
                                 top_k_subspace(Se_or, k)).max()
    eigPi = np.linalg.eigvalsh(sym(Pi))
    out = dict(
        C=C, theta_max=float(theta), pi_min=min_eig(Pi), pi_rank=int(np.linalg.matrix_rank(Pi, tol=1e-8)),
        pi_aniso=float(eigPi[-1] / max(eigPi[0], 1e-300)),
        feas=feas, dJ_oracle=dJ(Pi, Se_or), dJ_analytic=dJ(Pi, Se_an) if feas else np.nan,
        dJ_energy=dJ(Pi, Se_e), dJ_pidiag=dJ(Pi, Se_r), dJ_random=dJ(Pi, Se_rand),
        rate_oracle=rate_bits(Sigma, Se_or), rate_energy=rate_bits(Sigma, Se_e),
        rate_pidiag=rate_bits(Sigma, Se_r), rate_random=rate_bits(Sigma, Se_rand),
        amgm=amgm_bound(Sigma, Pi, C), J0=float(np.trace(case["P"] @ case["W"])),
    )
    for kk in ("energy", "pidiag", "random"):
        out["ratio_" + kk] = out["dJ_" + kk] / out["dJ_oracle"]
    out["feas_ok"] = (not out["feas"]) or loewner_le(Se_an, Sigma)
    out["oracle_feas_ok"] = loewner_le(Se_or, Sigma) and min_eig(Se_or) > -1e-10
    return out, xvec


# =================================================================================
def main():
    say("=" * 80)
    say("P0 判伪实验：同等信息预算下，表示对齐目标是否改变决策代价")
    say("=" * 80)
    A, B, Q, R = plant()
    K, P, S = lqr(A, B, Q, R)
    Pi = pi_matrix(K, S)

    # ------------------------------------------------------------ 校验 1
    say("\n[校验 1] 公式 J = tr(PW) + tr(Pi Sigma_e),  Pi = K^T S K   （蒙特卡洛 15 万步）")
    Sigma0 = lyap_d(A - B @ K, np.diag([0.02, 0.8, 0.02, 0.8]))
    say("   (1a) 与模拟闭环代价对比：")
    ok1 = True
    for frac, sd in ((0.03, 3), (0.15, 5), (0.35, 9)):
        r = validate_formula(A, B, Q, R, np.diag([0.02, 0.8, 0.02, 0.8]), frac * Sigma0, seed=sd)
        rel = (r["J_mc"] - r["J_formula"]) / r["J_formula"]
        z = rel * r["J_formula"] / r["se"]
        ok1 &= abs(z) < 4
        say(f"        Sigma_e={frac:.2f}*Sigma_x: MC={r['J_mc']:.5f}±{r['se']:.5f}  "
            f"公式={r['J_formula']:.5f}  Lyapunov={r['J_lyap']:.5f}  偏差={rel:+.2%} ({z:+.2f}σ)")
    say(f"        -> 公式与模拟在 {4}σ 内一致：{ok1}")
    say(f"   (1b) Pi 的结构：PSD={is_psd(Pi)}  min_eig={min_eig(Pi):.2e}  "
        f"rank(Pi)={np.linalg.matrix_rank(Pi, tol=1e-8)}  rank(B)={np.linalg.matrix_rank(B)}")
    say(f"        对照：tr(Pi)={np.trace(Pi):.4f} vs tr(P-K'RK)={np.trace(P - K.T @ R @ K):.4f}"
        f"（后者非 PSD 时无意义：min_eig={min_eig(P - K.T @ R @ K):.2e}）")
    say(f"    (1c) 闭环谱半径={float(np.abs(np.linalg.eigvals(A - B @ K)).max()):.4f} < 1，稳定")

    # ------------------------------------------------------------ 校验 2
    say("\n[校验 2] 等率最优解的三方交叉验证（解析 / 参数化优化器 / AM-GM 下界 / 约束）")
    ok2 = True
    for d, C in ((2, 6.0), (2, 1.5), (4, 3.0), (4, 6.0)):
        rng = np.random.default_rng(d * 10 + int(C))
        M = rng.standard_normal((d, d))
        Sigma = sym(M @ M.T) + np.eye(d)
        M2 = rng.standard_normal((d, d))
        Pi_ = sym(M2 @ M2.T) + 0.3 * np.eye(d)
        Se_or, val, _ = eps_oracle(Sigma, Pi_, C, seed=1, n_starts=8)
        bnd = amgm_bound(Sigma, Pi_, C)
        Se_an, feas, _ = eps_analytic(Sigma, Pi_, C)
        good = (rate_bits(Sigma, Se_or) - C) < 1e-8 and loewner_le(Se_or, Sigma) and val >= bnd - 1e-8
        if feas:
            good &= abs(val / dJ(Pi_, Se_an) - 1) < 1e-6
            tag = f"解析可行 -> 优化器与解析相对差={abs(val/dJ(Pi_,Se_an)-1):.1e}"
        else:
            tag = f"解析不可达 -> 优化器={val:.6f} > AM-GM 下界={bnd:.6f} (超出 {val/bnd-1:+.1%})"
        ok2 &= good
        say(f"    d={d}, C={C}: 率误差={abs(rate_bits(Sigma,Se_or)-C):.1e}, "
            f"Sigma_e⪯Sigma={loewner_le(Se_or,Sigma)}, {tag}")
    say(f"    -> 交叉验证通过：{ok2}")

    # ------------------------------------------------------------ 校验 3
    say("\n[校验 3] d=2 闭式差距律 gap=sqrt(1+sin^2(2θ)/4 (√ρ-1/√ρ)^2)   （C=5, 仅解析可达点）")
    rows_d2 = sweep_d2(C=5.0)
    errs = [abs(r["gap_num"] / r["gap_closed"] - 1) for r in rows_d2]
    say(f"    可用点 {len(rows_d2)} 个；数值/闭式最大相对偏差={max(errs):.2e}，均值={np.mean(errs):.2e}")
    say(f"    率一致性：max|r_g-C|={max(abs(r['r_g']-5.0) for r in rows_d2):.1e}, "
        f"max|r_e-C|={max(abs(r['r_e']-5.0) for r in rows_d2):.1e}")
    with open(os.path.join(OUT, "d2_gap_law.csv"), "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows_d2[0].keys()))
        w.writeheader()
        w.writerows(rows_d2)

    # ------------------------------------------------------------ 主实验
    say("\n[主实验] 真实 4 维控制系统（4 输入），数据协方差朝向扫描")
    say("  方案：(A) 等率最优（oracle）  (B) Sigma 基下重加权  (C) 能量/PCA 逆水填  (D) 随机表示")
    C_list = [2.0, 4.0, 6.0, 8.0]
    all_rows = []
    for C in C_list:
        case0 = build_case(0.0)
        _, warm = analyse(case0, C, seed=1, warm=None)
        rows = []
        for phi in np.linspace(0, 45, 16):
            case = build_case(float(phi))
            o, warm = analyse(case, C, seed=7, warm=warm)
            o["phi"] = float(phi)
            rows.append(o)
        all_rows += rows
        say(f"\n  --- C = {C:.1f} bits/sample ---")
        say("   phi  θ_max  ΔJ(oracle)  ΔJ(Σ重权)  ΔJ(PCA)   ΔJ(随机)  比PCA  比随机  比重权  解析可达")
        for o in rows:
            say(f"  {o['phi']:4.0f}  {o['theta_max']:5.1f}  {o['dJ_oracle']:.5f}   "
                f"{o['dJ_pidiag']:.5f}   {o['dJ_energy']:.5f}   {o['dJ_random']:.5f}  "
                f"{o['ratio_energy']:6.3f}  {o['ratio_random']:6.3f}  {o['ratio_pidiag']:6.3f}   {o['feas']}")
        med = float(np.median([o["ratio_energy"] for o in rows]))
        say(f"  中位：PCA 代价是 oracle 的 {med:.3f} 倍（{100*(med-1):+.1f}%）")
    with open(os.path.join(OUT, "sweep_phi.csv"), "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(all_rows[0].keys()))
        w.writeheader()
        w.writerows(all_rows)

    # 断言与统计
    say("\n  一致性/断言：")
    bad_rate = max(max(abs(o["rate_" + k] - o["C"]) for k in ("oracle", "energy", "pidiag", "random"))
                   for o in all_rows)
    say(f"    四方案的率与预算 C 的偏差最大值 = {bad_rate:.2e} bits  （等率公平性）")
    bad_feas = sum(1 for o in all_rows if not o["oracle_feas_ok"])
    say(f"    oracle 解满足 Sigma_e ⪯ Sigma 的违例数 = {bad_feas}")
    bad_bnd = max((o["amgm"] - o["dJ_oracle"]) / max(o["amgm"], 1e-12) for o in all_rows)
    say(f"    oracle ≥ AM-GM 下界的最大相对裕度 = {bad_bnd:.3f}（>0 表示下界不可达，符合预期）")
    order_ok = all(o["dJ_oracle"] <= o["dJ_pidiag"] + 1e-9 <= o["dJ_energy"] + 1e-7 for o in all_rows)
    say(f"    排序 ΔJ(oracle) ≤ ΔJ(Σ重权) ≤ ΔJ(PCA) 处处成立 = {order_ok}")

    # 低秩 B 情形：某些方向误差完全免费
    say("\n[附加实验] 低秩输入（rank(B)<d）：Pi=K'SK 低秩 -> 部分方向的误差完全免费")
    for n_in in (1, 2):
        A2, B2, Q2, R2 = plant()
        B2 = B2[:, :n_in]
        R2 = 0.5 * np.eye(n_in)
        K2, P2, S2 = lqr(A2, B2, Q2, R2)
        c2 = build_case(20.0, Q=Q2, B=B2, Sigma0=np.diag([4.0, 0.25, 4.0, 0.25]))
        C = 4.0
        o, _ = analyse(c2, C, seed=3)
        say(f"    输入数={n_in}: rank(Pi)={o['pi_rank']}/{len(c2['Sigma'])}  "
            f"θ_max={o['theta_max']:.1f}°  比PCA={o['ratio_energy']:.3f}  "
            f"比随机={o['ratio_random']:.3f}")

    # ------------------------------------------------------------ 图
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(1, 3, figsize=(16.5, 4.4))
        for C in C_list:
            rs = [o for o in all_rows if o["C"] == C]
            ax[0].plot([o["theta_max"] for o in rs], [o["ratio_energy"] for o in rs],
                       "o-", ms=3, label=f"PCA, C={C:.0f}b")
        ax[0].axhline(1.0, color="k", lw=1)
        ax[0].set_xlabel("principal angle theta_max (deg)")
        ax[0].set_ylabel("cost ratio  dJ(PCA)/dJ(oracle)")
        ax[0].set_title("P0-A: penalty of variance-aligned representation")
        ax[0].legend(fontsize=8)
        ax[0].grid(alpha=0.3)

        for C in C_list:
            rs = [o for o in all_rows if o["C"] == C]
            ax[1].plot([o["theta_max"] for o in rs], [o["ratio_random"] for o in rs],
                       "^--", ms=3, label=f"random, C={C:.0f}b")
            ax[1].plot([o["theta_max"] for o in rs], [o["ratio_pidiag"] for o in rs],
                       "s:", ms=3, label=f"reweighted, C={C:.0f}b")
        ax[1].axhline(1.0, color="k", lw=1)
        ax[1].set_xlabel("principal angle theta_max (deg)")
        ax[1].set_ylabel("cost ratio vs oracle")
        ax[1].set_title("P0-B: baselines and reweighting-only")
        ax[1].legend(fontsize=7)
        ax[1].grid(alpha=0.3)

        for rho in (4.0, 16.0, 64.0):
            rr = [r for r in rows_d2 if r["rho"] == rho]
            ax[2].plot([r["theta"] for r in rr], [r["gap_num"] for r in rr], "o", ms=3,
                       label=f"numeric rho={rho:g}")
            ax[2].plot([r["theta"] for r in rr], [r["gap_closed"] for r in rr], "-", lw=1)
        ax[2].axhline(1.0, color="k", lw=1)
        ax[2].set_xlabel("rotation theta (deg), d=2")
        ax[2].set_ylabel("cost ratio PCA/oracle")
        ax[2].set_title("P0-C: d=2 closed form (lines) vs numeric (dots)")
        ax[2].legend(fontsize=8)
        ax[2].grid(alpha=0.3)
        fig.tight_layout()
        fig.savefig(os.path.join(OUT, "p0_overview.png"), dpi=140)

        case = build_case(30.0)
        Se_or = eps_oracle(case["Sigma"], case["Pi"], 4.0, seed=1, n_starts=8)[0]
        Se_e = eps_diag_in_sigma_basis(case["Sigma"], case["Pi"], 4.0, "energy")
        fig2, ax2 = plt.subplots(1, 2, figsize=(11, 4.6))
        for M, nm, col in ((case["Sigma"], "Sigma (data)", "C0"),
                           (case["Pi"], "Pi (decision)", "C3")):
            w, V = np.linalg.eigh(sym(M))
            for i in range(len(w)):
                v = V[:, i] * math.sqrt(max(w[i], 0))
                ax2[0].arrow(0, 0, v[0], v[1], color=col, width=0.02,
                             head_width=0.12, alpha=0.8, length_includes_head=True)
        for M, nm, col in ((Se_or, "oracle error", "C2"), (Se_e, "PCA error", "C1")):
            w, V = np.linalg.eigh(sym(M))
            for i in range(len(w)):
                v = V[:, i] * math.sqrt(max(w[i], 0)) * 2
                ax2[1].arrow(0, 0, v[0], v[1], color=col, width=0.02,
                             head_width=0.12, alpha=0.8, length_includes_head=True)
        for a, ttl in ((ax2[0], "eigen-ellipses (first 2 coords): Sigma vs Pi"),
                       (ax2[1], "error ellipses: oracle vs PCA (x2 scale)")):
            a.set_aspect("equal")
            a.grid(alpha=0.3)
            a.set_title(ttl, fontsize=9)
        fig2.tight_layout()
        fig2.savefig(os.path.join(OUT, "p0_geometry.png"), dpi=140)
        say("\n  图：results/p0_overview.png, results/p0_geometry.png")
    except Exception as e:
        say(f"  [绘图跳过] {e}")

    # ------------------------------------------------------------ 判决
    ratios = np.array([o["ratio_energy"] for o in all_rows])
    thetas = np.array([o["theta_max"] for o in all_rows])
    med, worst = float(np.median(ratios)), float(np.max(ratios))
    corr = float(np.corrcoef(np.sin(np.radians(thetas)) ** 2, ratios)[0, 1])
    say("\n" + "=" * 80)
    say("判决")
    say("=" * 80)
    say(f"  能量/PCA 分配 vs 决策对齐分配（等率）：中位代价比 {med:.3f}（{100*(med-1):+.1f}%），"
        f"最坏 {worst:.3f}（{100*(worst-1):+.1f}%）")
    say(f"  corr(sin^2 θ_max, 代价比) = {corr:+.3f}      [θ_max 范围 {thetas.min():.1f}° ~ {thetas.max():.1f}°]")
    say(f"  仅重加权（不改基）方案的代价比中位 = {float(np.median([o['ratio_pidiag'] for o in all_rows])):.3f}")
    if med > 1.10 and corr > 0.4:
        v = ">>> 方向存活：对齐目标显著影响决策代价，且随失配角增长 -> 进入 P1"
    elif med > 1.03:
        v = ">>> 效应存在但中等（3%~10%）。需先把假设收窄（更强失配 / 更低预算 / 低秩输入）再判断"
    else:
        v = ">>> 方向被证伪：对齐无显著影响 -> 立即停手，写负结果"
    say("  " + v)
    say("\n  必须记住的建模前提（P1 要处理）：")
    say("   * 每样本**无记忆**高斯测试信道（无跨时间预测）。允许记忆信道时最优性需重做。")
    say("   * 编码分布取完美反馈稳态 Sigma_x^0（单点线性化，未做不动点闭合）。")
    say("   * 未包含 Kalman 滤波的结构（本实验的控制侧直接用 x_hat），故未使用控制对偶性。")
    say("   * oracle 需要求解带矩阵约束的小型优化问题，本身不是可部署算法；")
    say("     可部署版本是『按 Pi 加权分配』+『整块放弃最便宜方向』，属 P1 内容。")

    with open(os.path.join(OUT, "summary.txt"), "w", encoding="utf-8") as f:
        f.write("\n".join(REPORT))
    say("\n  报告：results/summary.txt")


if __name__ == "__main__":
    main()
