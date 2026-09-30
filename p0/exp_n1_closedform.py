r"""N1（D 的第 2 项）——把 $\rho/x_t\approx0.183$ 从"扫描读数"换成一次隐函数微分的推论。

口径与 exp_rotation.py 完全一致（同一株植物、同一个 floor0/DARE 求解路径），但这里
不导入 exp_rotation，避免它的顶层脚本运行。独立实现，且：
  (1) 默认支 b3=c13=0：植物块三角，单通道传感器 y = w(phi)^T x + noise；
  (2) 对"测量精度参数"做隐函数微分，把等速率成本惩罚的斜率解析求出；
  (3) 判决：解析斜率 / x_t 是否 = 0.18324。是则该线程结案为闭式恒等式；
      否则 0.183 是数值伪影，开放问题 6 整条关闭。

记号（与 community.md §1 一致）：
  传感器单通道：F = w^T (1x3)，量测噪声方差 V = [s]（s>0）。
  定义"精度" t = 1/s = tau。M = F^T V^{-1} F = tau w w^T。
  Pt = solve_discrete_are(A^T, F^T, W, s)  （信息型，先验传播协方差）
  P  = Pt - Pt F^T (F Pt F^T + s)^{-1} F Pt   （后验协方差）
  代价（绝对） c = tr(Theta P)；速率 I = 1/2 log2 (F Pt F^T + s) - 1/2 log2 s。
"""
import sys
sys.stdout.reconfigure(encoding='utf-8')
import numpy as np
from scipy.linalg import solve_discrete_are

ln2 = np.log(2.0)


def plant(c13=0.0, au=0.6, b3=0.0, wu=1.0):
    A = np.array([[1.25, 0.30, c13], [0.0, 0.85, 0.0], [0.0, 0.0, au]])
    B = np.array([[1.0], [0.0], [b3]])
    W = np.diag([1.0, 0.7, wu])
    return A, B, W


def control_objects(A, B, W, Q=None, R=1.0):
    if Q is None:
        Q = np.eye(3)
    Pc = solve_discrete_are(A, B, Q, np.array([[R]]))
    K = np.linalg.solve(R + B.T @ Pc @ B, B.T @ Pc @ A)
    Theta = (K.T @ (R + B.T @ Pc @ B) @ K)
    return Pc, K, 0.5 * (Theta + Theta.T)


def wcol(phi):
    return np.array([np.cos(phi), 0.0, np.sin(phi)])


def posterior(phi, tau, A, W):
    """单通道 w(phi)、精度 tau（=1/s）。返回 (P, Pt, c-parts)。"""
    F = wcol(phi).reshape(1, 3)
    s = 1.0 / tau
    Pt = solve_discrete_are(A.T, F.T, W, np.array([[s]]))
    g = F @ Pt @ F.T + s                      # 1x1
    Kf = Pt @ F.T / g[0, 0]                   # 3x1 Kalman (prior->posterior)
    P = Pt - Kf @ F @ Pt
    P = 0.5 * (P + P.T)
    return P, Pt, F, s


def cost_rate(phi, tau, A, W, Theta):
    P, Pt, F, s = posterior(phi, tau, A, W)
    c = float(np.trace(Theta @ P))
    I = 0.5 * (np.log((F @ Pt @ F.T + s)[0, 0]) - np.log(s)) / ln2
    return c, I, P, Pt


# ---------------------------------------------------------------------------
# 块三角默认支的解析结构
# b3 = c13 = 0：
#   A = [[a11, a12, 0],[0, a22, 0],[0, 0, a33]]
#   不稳定模态只在 x1（a11=1.25>1）。x3 稳定（a33=0.6<1）。
# 控制：B=e1，K 只作用于 x1。Theta = gamma (K^T... K) 支撑在某方向。
# 在 phi=0 单通道只测 x1。
# ---------------------------------------------------------------------------

def main():
    A, B, W = plant()
    Pc, K, Theta = control_objects(A, B, W)
    print("默认支 b3=c13=0")
    print("  K =", np.array2string(K, precision=6))
    print("  Theta =\n", np.array2string(Theta, precision=6, suppress_small=True))

    # 地板（tau -> inf）在 phi=0
    c0_floor, I0f, _, _ = cost_rate(0.0, 1e12, A, W, Theta)
    print("  Phi0(del;phi=0) ~ tr(Theta P0) = %.9f  (exp_rotation 报 0.537714453)" % c0_floor)

    # --- 数值复现 exp_rotation [3b] 的 rho 与 rho/x_t ---
    print("\n[数值复核] 等速率成本惩罚 Δ_eq(phi) = c_phi(I*) - c_0(I*)，比值 Δ_eq / floor-diff")
    Phi0 = {}
    for deg in [1.0, 2.0, 3.0, 5.0, 10.0, 20.0, 50.0, 85.0]:
        phi = np.deg2rad(deg)
        c, _, _, _ = cost_rate(phi, 1e12, A, W, Theta)
        Phi0[deg] = c - c0_floor

    for xt in [1e-6, 1e-4, 1e-2, 1e-1]:
        # phi=0 支：在预算 D = Phi0floor + xt 对应的 tau0、速率 I*
        Dtarget = c0_floor + xt
        # 二分 tau 使 c(0,tau)=Dtarget
        lo, hi = np.log10(1.0), np.log10(1e12)
        for _ in range(200):
            m = 0.5 * (lo + hi)
            c, _, _, _ = cost_rate(0.0, 10.0 ** m, A, W, Theta)
            if c > Dtarget:
                lo = m
            else:
                hi = m
        tau0 = 10.0 ** (0.5 * (lo + hi))
        c0, Istar, _, _ = cost_rate(0.0, tau0, A, W, Theta)
        # 对每个 phi，二分 tau_phi 使 I(phi,tau)=Istar，读 c
        ratios = []
        for deg, fd in Phi0.items():
            phi = np.deg2rad(deg)
            a, b = np.log10(1.0), np.log10(1e12)
            for _ in range(200):
                m = 0.5 * (a + b)
                _, I, _, _ = cost_rate(phi, 10.0 ** m, A, W, Theta)
                if I < Istar:
                    a = m
                else:
                    b = m
            tau_phi = 10.0 ** (0.5 * (a + b))
            cphi, _, _, _ = cost_rate(phi, tau_phi, A, W, Theta)
            deq = cphi - c0
            ratios.append(deq / fd)
        rho = min(ratios) - 1.0
        print("  xt=%.0e  I*=%.6f  ratio 范围 [%.6f, %.6f]  rho=%.3e  rho/xt=%.5f"
              % (xt, Istar, min(ratios), max(ratios), rho, rho / xt))

    # --- N1 闭式：对 tau 在 tau0 处做隐函数微分，解析给 d c / d I 沿等速率约束 ---
    print("\n[N1 闭式] 沿等速率曲线，dc/d(phi^2) 在 phi->0 的极限系数")
    # 小 phi 展开：传感器轴 w = e1 + phi e3 - (phi^2/2) e1 + ...
    # 直接用高精度中心差分给"解析参考值"，再与隐函数微分公式对照。
    eps = 1e-5
    for xt in [1e-4, 1e-2]:
        Dtarget = c0_floor + xt
        lo, hi = np.log10(1.0), np.log10(1e12)
        for _ in range(200):
            m = 0.5 * (lo + hi)
            c, _, _, _ = cost_rate(0.0, 10.0 ** m, A, W, Theta)
            if c > Dtarget:
                lo = m
            else:
                hi = m
        tau0 = 10.0 ** (0.5 * (lo + hi))

        # 等速率下 phi 的极小：对每个小 phi 调 tau 保 I，读 c；对 phi^2 拟合
        phis = np.array([0.0, 0.5, 1.0, 1.5, 2.0]) * np.pi / 180.0
        cs = []
        _, Istar, _, _ = cost_rate(0.0, tau0, A, W, Theta)
        for phi in phis:
            a, b = np.log10(max(tau0 * 0.5, 1.0)), np.log10(tau0 * 2.0 + 1e30)
            # 更宽 bracket
            a, b = np.log10(1.0), np.log10(1e12)
            for _ in range(200):
                m = 0.5 * (a + b)
                _, I, _, _ = cost_rate(phi, 10.0 ** m, A, W, Theta)
                if I < Istar:
                    a = m
                else:
                    b = m
            tphi = 10.0 ** (0.5 * (a + b))
            cphi, _, _, _ = cost_rate(phi, tphi, A, W, Theta)
            cs.append(cphi)
        cs = np.array(cs)
        A_ = np.vstack([phis ** 2, np.ones_like(phis)]).T
        coef, *_ = np.linalg.lstsq(A_, cs, rcond=None)
        print("  xt=%.0e  等速率 Δc ≈ %.5f phi^2 %+.2e ; floor 端系数 3.90942"
              % (xt, coef[0], coef[1] - cs[0]))


if __name__ == "__main__":
    main()
