r"""N1 第三步：把 0.183 识别为极限常数，给它一个可复现的精确值与误差阶。

定义（无量纲）：
  A_eq(tau) = phi->0 时等速率成本 c(phi)|I=const 对 phi^2 的系数
  x(tau) = c(0,tau) - c(0,inf)
  R(x) = (A_eq(tau)/A_floor - 1) / x
我们已知 R -> 0.18323（x->0）。本脚本：
  (1) 加密 x（走 tau 的对数网格）外推 R(0)；
  (2) 量 R(x) 的展开阶（判断是 R=R0 + R1 x 还是别的）；
  (3) 报告 R0 的数值，与 exp_rotation 的 0.18324 对照。
"""
import sys
sys.stdout.reconfigure(encoding='utf-8')
import numpy as np
from p0.exp_n1_stein import plant, control, eval_point, hessian_phi

AFLOOR_COEF = 3.909423


def tau_for_budget(xt, cfloor, A, W, Theta):
    """二分 tau 使 c(0,tau)=cfloor+xt（单调标量；仅用于把预算 x 映射到工作点）。"""
    target = cfloor + xt
    lo, hi = np.log10(1.0), np.log10(1e12)
    for _ in range(200):
        m = 0.5 * (lo + hi)
        c, _ = eval_point(0.0, 10.0 ** m, A, W, Theta)
        if c > target:
            lo = m
        else:
            hi = m
    return 10.0 ** (0.5 * (lo + hi))


def R_of_x(xt, cfloor, A, W, Theta):
    tau = tau_for_budget(xt, cfloor, A, W, Theta)
    af, ae, c0 = hessian_phi(0.0, tau, A, W, Theta)
    R = (ae / AFLOOR_COEF - 1.0) / xt
    return R, af, ae, tau


def main():
    A, B, W = plant()
    Pc, K, Theta = control(A, B, W)
    cfloor, _ = eval_point(0.0, 1e12, A, W, Theta)

    print("%10s %12s %12s %12s" % ("x_t", "tau", "A_eq", "R(x)"))
    data = []
    for xt in [1e-6, 3e-6, 1e-5, 3e-5, 1e-4, 3e-4, 1e-3, 3e-3, 1e-2, 3e-2, 1e-1]:
        R, af, ae, tau = R_of_x(xt, cfloor, A, W, Theta)
        data.append((xt, R))
        print("%10.1e %12.3e %12.6f %12.7f" % (xt, tau, ae, R))

    # x->0 Richardson：x 在小端近似按 ~1/tau 走，取最左几点线性外推
    data.sort()
    xs = np.array([d[0] for d in data])
    Rs = np.array([d[1] for d in data])

    # 用最左 4 点，对 x 线性拟合
    k = 4
    p1 = np.polyfit(xs[:k], Rs[:k], 1)
    print("\n最左 %d 点线性外推 R(0) = %.7f  (斜率 %.4f)" % (k, p1[1], p1[0]))

    # 用 2、3、4 点看收敛
    for k in (2, 3, 4, 5):
        p = np.polyfit(xs[:k], Rs[:k], 1)
        print("  k=%d: R(0)=%.7f" % (k, p[1]))

    # 检验展开阶：相邻点 (R-R0)/x 是否稳定
    R0 = p1[1]
    print("\n残差 (R-R0)/x（若趋于常数则 R=R0+R1 x）：")
    for x, R in data[:7]:
        print("  x=%.2e  (R-R0)/x = %.4f" % (x, (R - R0) / x))


if __name__ == "__main__":
    main()
