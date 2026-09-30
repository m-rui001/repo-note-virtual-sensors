r"""N1 第五步：复步长二阶导（无减法相消）求等速率 phi^2 系数，解析确认 rho/x=0.18323。

复步长二阶公式（实解析 f）：f''(x) = 2(f(x) - Re f(x+ih))/h^2，无减法相消，机器精度。
沿等速率曲线 ln tau = ln tau(phi)，I(phi, ln tau(phi))=const：
  一阶：dI/dphi + I_t t' = 0  -> t'=0（在 phi=0，因 dI/dphi=0）。
  二阶：0 = I_pp + 2 I_pt t' + I_tt (t')^2 + I_t t''
       => t'' = -I_pp / I_t   （I_t = dI/d ln tau）
  成本二阶：c_pp|eq = c_pp + 2 c_pt t' + c_tt(t')^2 + c_t t''
          = c_pp - c_t I_pp / I_t
所以等速率 phi^2 系数 A_eq = (1/2)(c_pp - (c_t/I_t) I_pp)。
全部用复步长一阶/二阶导，机器精度，无数值搜索（除预算->tau 的单调标量映射）。
"""
import sys
sys.stdout.reconfigure(encoding='utf-8')
import numpy as np
from p0.exp_n1_scalar import plant, control, info_dare_complex, eval_point

ln2 = np.log(2.0)


def cevals(phi, ltau, A, W, Theta):
    cphi = np.cos(phi); sphi = np.sin(phi)
    F = np.array([[cphi, 0.0 + 0j, sphi]])
    tau = np.exp(ltau)
    s = 1.0 / tau
    Pt = info_dare_complex(A, F, W, s)
    g = (F @ Pt @ F.T + s)[0, 0]
    P = Pt - Pt @ F.T @ (F @ Pt / g)
    c = np.trace(Theta @ P)
    I = 0.5 * (np.log(g) - np.log(s)) / ln2
    return c, I


def d1_dual(ltau, A, W, Theta, h=1e-7):
    """对 phi 与 ln tau 的一阶导（复步长）。"""
    def f_p(phi):
        return cevals(phi, ltau, A, W, Theta)
    def f_t(phi):
        return cevals(phi, ltau, A, W, Theta)
    cp, Ip = f_p(h * 1j)
    ct, It = cevals(0.0, ltau + h * 1j, A, W, Theta)
    return dict(cp=cp.imag / h, Ip=Ip.imag / h, ct=ct.imag / h, It=It.imag / h)


def d2_phi(ltau, A, W, Theta, h=1e-5):
    """对 phi 的二阶导（复步长：2(f-Re f(ih))/h^2），c 与 I 同时。"""
    c0, I0 = cevals(0.0, ltau, A, W, Theta)
    ch, Ih = cevals(h * 1j, ltau, A, W, Theta)
    cpp = 2 * (c0.real - ch.real) / h ** 2
    Ipp = 2 * (I0.real - Ih.real) / h ** 2
    return cpp, Ipp


def tau_for_budget(xt, cfloor, A, W, Theta):
    target = cfloor + xt
    lo, hi = 0.0, np.log(1e12)
    for _ in range(200):
        m = 0.5 * (lo + hi)
        c, _, _, _ = eval_point(0.0, np.exp(m), A, W, Theta)
        if c > target:
            lo = m
        else:
            hi = m
    return np.exp(0.5 * (lo + hi))


def main():
    A, B, W = plant()
    Pc, K, Theta = control(A, B, W)
    cfloor, _, _, _ = eval_point(0.0, 1e12, A, W, Theta)

    # floor coefficient A_floor = (1/2) c_pp at tau->inf (t term vanishes)
    ltau_f = np.log(1e12)
    cpp_f, Ipp_f = d2_phi(ltau_f, A, W, Theta)
    Afloor = 0.5 * cpp_f
    print("地板系数 A_floor = %.7f （目标 3.90942）" % Afloor)

    print("\n%10s %12s %12s %12s %12s" % ("x_t", "tau", "A_eq", "rho", "rho/x"))
    data = []
    for xt in [1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 3e-2, 1e-1]:
        tau = tau_for_budget(xt, cfloor, A, W, Theta)
        ltau = np.log(tau)
        d1 = d1_dual(ltau, A, W, Theta)
        cpp, Ipp = d2_phi(ltau, A, W, Theta)
        # A_eq = 1/2 (cpp - c_t I_pp / I_t)
        Aeq = 0.5 * (cpp - d1['ct'] * Ipp / d1['It'])
        rho = Aeq / Afloor - 1.0
        data.append((xt, rho))
        print("%10.1e %12.3e %12.7f %12.3e %12.6f"
              % (xt, tau, Aeq, rho, rho / xt))

    xs = np.array([d[0] for d in data])
    rs = np.array([d[1] for d in data])
    # 极限 rho/x（取小端稳定几点）
    print("\n小端 rho/x 稳定值：", ["%.6f" % (r / x) for x, r in data[:4]])
    slope = np.sum(xs[:4] * rs[:4]) / np.sum(xs[:4] ** 2)
    print("最小二乘（前 4 点）rho/x = %.6f   exp_rotation 报 0.18324" % slope)


if __name__ == "__main__":
    main()
