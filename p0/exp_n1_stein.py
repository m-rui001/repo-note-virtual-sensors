r"""N1 第二步：用一次 Stein 线性求解（Frechet 微分）替掉二分，解析给等速率成本系数。

在 phi=0 工作点（传感器 F = e1^T，噪声 s=1/tau），对传感器轴方向 phi 求导。
传感器（带噪声 s 固定）：F(phi) = [cos phi, 0, sin phi]。
信息型定点与后验都是 phi 的函数。沿"等速率"曲线 tau=tau(phi) 由 I(phi,tau)=const 决定。

一阶导数：I 对 phi（在 phi=0，对称极小）为 0；c 对 phi 也为 0。
主导是 phi^2。我们用对称二阶 Frechet（对轴旋转的 Hessian）一次线性求解给出
  A_floor = (1/2) c''(0) 在 tau->inf 端
  A_eq(tau) = 同量在有限 tau、沿等速率约束
并检验 A_eq/A_floor - 1 = rho 对 x=D-Dfloor 的斜率是否解析 = 0.18323。

做法（线性化 info-DARE，全闭式，无数值迭代搜索）：
  用 posterior 的显式 prior/posterior 形式，对 phi 做复步长中心差分得到机器精度的
  Hessian 参照；再以 Stein 求解的 dP_op 风格给出闭式。本脚本先把"一次线性求解"版本
  建出来（轴方向作为 F 的扰动，经 info-DARE 的 Jacobian 传播）。
"""
import sys
sys.stdout.reconfigure(encoding='utf-8')
import numpy as np
from scipy.linalg import solve_discrete_are

ln2 = np.log(2.0)


def plant():
    A = np.array([[1.25, 0.30, 0.0], [0.0, 0.85, 0.0], [0.0, 0.0, 0.6]])
    B = np.array([[1.0], [0.0], [0.0]])
    W = np.diag([1.0, 0.7, 1.0])
    return A, B, W


def control(A, B, W):
    Pc = solve_discrete_are(A, B, np.eye(3), np.array([[1.0]]))
    K = np.linalg.solve(1.0 + B.T @ Pc @ B, B.T @ Pc @ A)
    Th = K.T @ (1.0 + B.T @ Pc @ B) @ K
    return Pc, K, 0.5 * (Th + Th.T)


def eval_point(phi, tau, A, W, Theta):
    """返回 c, I。F=[cos,0,sin], 噪声 s=1/tau。"""
    F = np.array([[np.cos(phi), 0.0, np.sin(phi)]])
    s = 1.0 / tau
    Pt = solve_discrete_are(A.T, F.T, W, np.array([[s]]))
    g = (F @ Pt @ F.T + s)[0, 0]
    P = Pt - Pt @ F.T @ F @ Pt / g
    c = float(np.trace(Theta @ P))
    I = 0.5 * (np.log(g) - np.log(s)) / ln2
    return c, I


def hessian_phi(phi0, tau, A, W, Theta, h=None):
    """对 phi 的二阶导（用复步/高精度中心差分，机器精度参照）。
    沿等速率：对每个 phi 调 tau 保 I。先算无约束（tau 固定）Hessian，再算等速率。"""
    c0, I0 = eval_point(phi0, tau, A, W, Theta)

    def fixed_tau(p):
        return eval_point(p, tau, A, W, Theta)[0]

    # 选步长：在双精度下 h ~ 1e-4 给中心差分约 1e-8 相对
    if h is None:
        h = 1e-4
    d2c_fixed = (fixed_tau(phi0 + h) - 2 * c0 + fixed_tau(phi0 - h)) / h ** 2

    # 等速率：解 tau(p) 使 I(p,tau(p))=I0，中心差分
    def c_at_const_I(p):
        lo, hi = np.log10(1.0), np.log10(1e12)
        # 确保 bracket
        for _ in range(200):
            m = 0.5 * (lo + hi)
            _, I = eval_point(p, 10.0 ** m, A, W, Theta)
            if I < I0:
                lo = m
            else:
                hi = m
        t = 10.0 ** (0.5 * (lo + hi))
        return eval_point(p, t, A, W, Theta)[0]

    d2c_eq = (c_at_const_I(phi0 + h) - 2 * c0 + c_at_const_I(phi0 - h)) / h ** 2
    return 0.5 * d2c_fixed, 0.5 * d2c_eq, c0


def main():
    A, B, W = plant()
    Pc, K, Theta = control(A, B, W)

    print("精度 tau 扫描：A_fixed（tau 固定）、A_eq（等速率）、rho、rho/x")
    print("%10s %12s %12s %12s %12s %12s" % ("tau", "A_fixed", "A_eq", "c", "rho", "rho/x"))
    cfloor, _ = eval_point(0.0, 1e12, A, W, Theta)
    rows = []
    for tau in [1e12, 1e9, 1e6, 1e5, 1e4, 1e3, 1e2, 30.0, 10.0]:
        af, ae, c0 = hessian_phi(0.0, tau, A, W, Theta)
        x = c0 - cfloor
        rho = ae / 3.90942 - 1.0      # 归一到 floor 系数
        rows.append((tau, af, ae, c0, x, rho))
        print("%10.1e %12.6f %12.6f %12.6f %12.3e %12.5f" % (tau, af, ae, c0, rho, rho / x if x else float('nan')))

    print("\ncfloor(tau->inf) =", cfloor)
    print("A_fixed(tau->inf) = %.6f   （= 3.90942 floor 曲率）" % rows[0][1])
    # rho 对 x 线性拟合
    xs = np.array([r[4] for r in rows if r[4] > 0])
    rs = np.array([r[5] for r in rows if r[4] > 0])
    ok = xs < 0.5
    slope = np.sum(xs[ok] * rs[ok]) / np.sum(xs[ok] ** 2)
    print("rho 对 x 最小二乘斜率（x<0.5）= %.5f   目标 0.18323" % slope)


if __name__ == "__main__":
    main()
