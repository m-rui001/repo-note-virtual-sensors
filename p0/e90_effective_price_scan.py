r"""E90：§56 的新规则"首阈值丢弃稳定轴"到底是**计价律**还是**经验偏好**？

§56 的 6 点表里已经藏着矛盾（我自己先算出来再跑）：
  规则 (a)「丢稳定轴」命中 C1/C2/T1/T2/T3，T4 违反。
  若存在逐轴可分的"有效价格" $\pi_i=\lambda_i(\Theta)\phi(|\lambda_i(A)|)$，
  那么"稳定轴相对便宜程度" $r=\lambda_{\rm stab}/\lambda_{\rm unstab}$ 应当决定切换：
  $r$ 大（稳定轴贵）$\Rightarrow$ 丢稳定轴；$r$ 小 $\Rightarrow$ 丢不稳定的便宜轴。
  但实测 T1 的 $r=1.335/0.598=2.23$ **丢稳定轴**（服从 (a)），
  T4 的 $r=5.849/0.743=7.87$ **丢不稳定轴**（违反 (a)），
  T2 的 $r=5.849/0.598=9.78$ 又**丢稳定轴**。⇒ $r$ 随机的非单调，逐轴价格 $\phi$ 已经被 6 点否掉一半。

本脚本沿一维路径 $Q=(1,1,q_3)$（以及对照路径 $Q=(q_1,q_1,1)$）对每个 $\rho$ 扫 $q_3$，
记录每次首阈值丢的是稳定轴还是不稳定轴，看切换点是不是只由 $\rho$ 决定。

预先写下的判据（不许事后改）：
 (P) 每个 $\rho$ 的序列里切换**至多一次**，且切换时的 $r^*$（或 $\lambda_{\rm unstab}^*$ 等标量）
     只依赖 $\rho$、不依赖路径 $\Rightarrow$ 逐轴有效价格成立，$\phi$ 从 $r^*(\rho)$ 曲线读出；
 (K) 出现多次切换、或同一 $\rho$ 不同路径给不同切换点、或 $r^*$ 与 $\rho$ 无单调关系
     $\Rightarrow$ **逐轴可分性失败**。§56 的规则 (a) 立刻降级为"经验偏好"，
     我不再宣称任何"有效价格"，只保留"丢弃方向仍是 $\Theta$ 的本征向量"这条弱事实。
"""
import sys
import numpy as np

sys.stdout.reconfigure(encoding="utf-8")
import cvxpy as cp  # noqa: E402
from p0.e89_unstable_cheap_axis import Pl, solve, rk, sigma, TOL  # noqa: E402

LN2 = np.log(2.0)


def first_threshold(pl, UN, ST):
    """自适应括号 + 二分：谓词 rank==n 在 D<Dc 为真。返回 (Dc, 落零方向, rel, I)。"""
    n = pl.n
    prev = pl.jc + 1e-3
    r = solve(pl, prev)
    if r is None or rk(r) < n:
        return None
    hi = None
    for k in range(11):
        d = pl.jc + 1e-3 + 0.25 * (2.0 ** k)
        r = solve(pl, d)
        if r is None:
            continue
        if rk(r) < n:
            hi = d
            break
        prev = d
    if hi is None:
        return None
    lo = prev
    for _ in range(13):
        m = 0.5 * (lo + hi)
        r = solve(pl, m)
        if r is None:
            return None
        if rk(r) == n:
            lo = m
        else:
            hi = m
    dc = 0.5 * (lo + hi)
    r = solve(pl, dc - 5e-3, eps=1e-11, it=40000)
    if r is None:
        return None
    idx = n - 1
    v = r["V"][:, idx] / np.linalg.norm(r["V"][:, idx])
    return dc, v, r["ev"][idx] / max(r["ev"][0], 1e-300), r["I"]


RHOS = [1.02, 1.05, 1.15, 1.40, 2.00]
Q3S = [1.0, 1.6, 2.5, 4.0, 6.4, 10.0, 16.0, 25.0]

if __name__ == "__main__":
    print("== E90：沿 Q=(1,1,q3) 路径扫首阈值丢哪根轴（逐轴有效价格成立吗）==")
    for rho in RHOS:
        print("\n--- rho=%.2f  稳定化下界 %.4f bit ---" % (
            rho, 2 * np.log(rho) / LN2))
        print("   %6s %22s %8s %9s  %s" % ("q3", "lam(Theta) 升序", "r=S/U", "D_c", "丢弃轴 重叠^2 sigma_Theta"))
        for q3 in Q3S:
            pl = Pl(rho=rho, Q=[1., 1., q3])
            UN, ST = pl.subspaces()
            wT, UT = np.linalg.eigh(pl.Theta)
            pst = np.array([float(np.sum((ST.T @ UT[:, i]) ** 2)) for i in range(3)])
            pun = np.array([float(np.sum((UN.T @ UT[:, i]) ** 2)) for i in range(3)])
            stable = np.where(pst > 0.5)[0]
            unstab = np.where(pun > 0.5)[0]
            if len(stable) != 1:
                print("   %6.1f  稳定轴不唯一（p_stab=%s），跳过" % (
                    q3, np.array2string(pst, precision=2)))
                continue
            si = stable[0]
            ui = min(unstab, key=lambda i: wT[i]) if len(unstab) else None
            out = first_threshold(pl, UN, ST)
            if out is None:
                print("   %6.1f  没扫到降秩（$D_c$ 超出括号）" % q3)
                continue
            dc, v, rel, I = out
            ov = np.abs(UT.T @ v) ** 2
            j = int(ov.argmax())
            lab = "稳定" if j == si else ("不稳定" if j in tuple(unstab) else "混合")
            rr = wT[si] / wT[ui] if ui is not None else float("nan")
            print("   %6.1f %22s %8.2f %9.4f  %s轴 %d(%.4f) sig=%.2e rel=%.1e I=%.3f" % (
                q3, np.array2string(wT, precision=3), rr, dc, lab, j + 1,
                ov[j], sigma(v, pl.Theta), rel, I))
