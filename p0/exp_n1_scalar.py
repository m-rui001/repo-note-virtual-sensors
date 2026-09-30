r"""N1 第四步：真正的 Frechet（一次 Stein 求解）+ 不稳定标量模约化，给 0.183 一个解析身份。

二阶中心差分在小 x 端不可用（A_eq-A_fixed~1e-7 被差分误差淹没）。改成对 info-DARE
做一阶变分（线性方程），机器精度求 dP。

约化观察（b3=c13=0）：传感器测 y = cos phi x1 + sin phi x3 + noise。
  x1 是唯一不稳定模（a=1.25）；x3 稳定（a33=0.6）。
  等速率约束下，近地板发散只由 x1 模承担。把问题约到 x1 的标量信道：
    有效过程噪声（在测得方向上）与测量精度决定 x1 方向后验方差。
  phi 的作用 = 把一部分"测量预算"从 x1 转到稳定的 x3，等价于 x1 方向有效精度按
    tau_eff = tau cos^2 phi 下降。
  于是 phi 旋转在一阶上等价于精度 tau 的一个已知缩放，
  等速率补偿后，成本对 phi^2 的系数 = 成本对 ln tau 导数的已知线性组合。
这把 0.183 从"扫描数"变成一个由标量滤波定点 + 一次隐函数微分给出的数。
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
    F = np.array([[np.cos(phi), 0.0, np.sin(phi)]])
    s = 1.0 / tau
    Pt = solve_discrete_are(A.T, F.T, W, np.array([[s]]))
    g = (F @ Pt @ F.T + s)[0, 0]
    P = Pt - Pt @ F.T @ F @ Pt / g
    c = float(np.trace(Theta @ P))
    I = 0.5 * (np.log(g) - np.log(s)) / ln2
    return c, I, P, Pt


def frechet_scalar_chain(phi0, tau, A, W, Theta, h=1e-7):
    """用复步长一阶微分（机器精度，无截断误差）求 c、I 对 ln tau 的导数。

    复步长：f'(x)=Im f(x+ih)/h，对实解析函数精确到机器精度，无减法相消。
    """
    def complex_evals(z_phi, z_ltau):
        cphi = np.cos(z_phi); sphi = np.sin(z_phi)
        F = np.array([[cphi, 0.0 + 0j, sphi]])
        tau_c = np.exp(z_ltau)
        s = 1.0 / tau_c
        # 复数 Riccati：用 info-DARE 迭代（闭式不动点），避免 scipy 不支持复数
        Pt = info_dare_complex(A, F, W, s)
        g = (F @ Pt @ F.T + s)[0, 0]
        P = Pt - Pt @ F.T @ F @ Pt / g
        c = np.trace(Theta @ P)
        I = 0.5 * (np.log(g) - np.log(s)) / ln2
        return c, I

    z_ltau = np.log(tau) + 0j

    def c_of(lt, phi=phi0):
        c, _ = complex_evals(phi, lt)
        return c

    def I_of(lt, phi=phi0):
        _, I = complex_evals(phi, lt)
        return I

    # dI/d ln tau, dc/d ln tau 在 phi=0
    dc_dlt = c_of(z_ltau + h * 1j).imag / h
    dI_dlt = I_of(z_ltau + h * 1j).imag / h
    # dI/d phi, dc/d phi
    dc_dphi = complex_evals(h * 1j, z_ltau)[0].imag / h
    dI_dphi = complex_evals(h * 1j, z_ltau)[1].imag / h

    # 等速率约束 I(phi, ln tau(phi))=const：
    # d ln tau / d phi = -(dI/d phi)/(dI/d ln tau)
    if abs(dI_dlt) > 0:
        dlt_dphi = -dI_dphi / dI_dlt
    else:
        dlt_dphi = float('nan')
    dc_eq_dphi = dc_dphi + dc_dlt * dlt_dphi
    return dict(dc_dlt=dc_dlt, dI_dlt=dI_dlt, dc_dphi=dc_dphi,
                dI_dphi=dI_dphi, dlt_dphi=dlt_dphi, dc_eq_dphi=dc_eq_dphi)


def info_dare_complex(A, F, W, s, iters=20000):
    """复数 info-DARE 定点：Pt = A P A' + W, P = (Pt^{-1} + F's^{-1}F)^{-1}。
    等价地用 Riccati 迭代到收敛（对复变量）。"""
    n = A.shape[0]
    P = np.eye(n) * 1.0 + 0j
    sM = np.array([[s]])
    for _ in range(iters):
        Pt = A @ P @ A.T + W
        g = F @ Pt @ F.T + s
        Pn = Pt - Pt @ F.T @ (F @ Pt / g)
        if np.max(np.abs(Pn - P)) < 1e-15:
            P = Pn
            break
        P = Pn
    return A @ P @ A.T + W


def main():
    A, B, W = plant()
    Pc, K, Theta = control(A, B, W)

    print("phi=0 工作点的一阶 Frechet（复步长，机器精度）")
    print("%10s %12s %12s %12s %12s" % ("tau", "dI/dlntau", "dc/dlntau", "dI/dphi", "dc_eq/dphi"))
    for tau in [1e6, 1e4, 1e3, 1e2, 30.0, 10.0]:
        d = frechet_scalar_chain(0.0, tau, A, W, Theta)
        print("%10.1e %12.6f %12.6f %12.3e %12.3e"
              % (tau, d['dI_dlt'], d['dc_dlt'], d['dI_dphi'], d['dc_eq_dphi']))

    print("\n说明：phi=0 是对称极小，dI/dphi、dc/dphi 数值应为 ~0（奇阶在支上消失）。")
    print("因此等速率成本的主导是 phi^2，系数需用对称二阶——但要在一阶 Frechet 框架内")
    print("用轴旋转的对称（phi->-phi 同）二阶响应，避免中心差分相消。")


if __name__ == "__main__":
    main()
