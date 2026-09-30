r"""E87：E86b 的"落零方向≈Θ 本征向量、且按 λ(Θ) 升序被丢弃"是不是锚点偶然？

做法：换 $Q$（从而换 $P_c$、换 $\Theta=K^\top(R+B^\top P_cB)K$、换渐近线 $j_c$），
重新用二分定位三个秩跳变阈值，在每个阈值下方量落零方向对 $\Theta$ 本征轴的重叠。
锚点（$Q=I$）的已知结果用来做自检：
  $D_c$ = 32.3119 / 34.4139 / 56.6493，最大重叠^2 = 0.994 / 0.974 / 0.972，
  命中的 $\Theta$ 轴（按代价升序编号）= 1 / 2 / 3。
预言（若谱序结构是真的）：任意 $Q$ 下命中轴的编号随 $D_c$ 单调上升，且重叠接近 1。
"""
import sys
import numpy as np

sys.stdout.reconfigure(encoding="utf-8")
import cvxpy as cp  # noqa: E402
import p0.p0_replicate_letter as base  # noqa: E402

LN2 = np.log(2.0)
A, B, W, n = base.A, base.B, base.W, base.n
TOL = 1e-4          # 秩判定容差（方向本身对容差不敏感，见 §55 的方法学注）


def solve(D, Theta, jc, solver=cp.SCS, eps=1e-9, it=8000):
    P = cp.Variable((n, n), symmetric=True)
    Pi = cp.Variable((n, n), symmetric=True)
    const = 0.5 * np.log(np.linalg.det(W)) / LN2
    AP = A @ P
    Pt = A @ P @ A.T + W
    prob = cp.Problem(cp.Minimize(-0.5 * cp.log_det(Pi) / LN2 + const),
                      [Pi >> 0, P >> 0, cp.trace(Theta @ P) + jc <= D, P << Pt,
                       cp.bmat([[P - Pi, AP.T], [AP, Pt]]) >> 0])
    try:
        prob.solve(solver=solver, gp=False, eps=eps, max_iters=it)
    except Exception:
        return None
    if P.value is None:
        return None
    Pv = base.sym(P.value)
    Pt = A @ Pv @ A.T + W
    S = base.sym(np.linalg.inv(Pv) - np.linalg.inv(Pt))
    ev, V = np.linalg.eigh(S)
    return dict(ev=ev[::-1], V=V[:, ::-1],
                I=float(0.5 * np.log(np.linalg.det(Pt) / np.linalg.det(Pv)) / LN2))


def rk(r):
    return int(np.sum(r["ev"] / max(abs(r["ev"][0]), 1e-300) > TOL))


def scan(Q, Theta, jc, hi=None):
    """粗扫 + 二分，返回 [(Dc, 消失方向的 ev/V), ...] 按 Dc 升序。"""
    hi = hi or jc + 80.0
    Ds = np.arange(jc + 0.5, hi, 1.0)
    prev = None
    jumps = []
    for d in Ds:
        r = solve(d, Theta, jc)
        if r is None:
            continue
        k = rk(r)
        if prev is not None and k != prev[1]:
            jumps.append((prev[0], d, prev[1], k))
        prev = (d, k)
    out = []
    for lo, up, klo, khi in jumps:
        for _ in range(9):
            m = 0.5 * (lo + up)
            r = solve(m, Theta, jc, solver=cp.SCS)
            if r is None:
                break
            if rk(r) == khi:
                up = m
            else:
                lo = m
        else:
            dc = 0.5 * (lo + up)
            # 阈值下方 5e-3 处取落零方向（rank=klo 时最小非零本征值的方向）
            r = solve(dc - 5e-3, Theta, jc, solver=cp.SCS, eps=1e-11, it=40000)
            if r is not None and klo >= 1:
                idx = min(klo, n) - 1
                v = r["V"][:, idx] / np.linalg.norm(r["V"][:, idx])
                out.append((dc, klo, khi, idx, v, r["ev"][idx] / max(r["ev"][0], 1e-300)))
    return out


def sigma(v, M):
    """=0 iff v 是 M 的本征向量（与"最接近哪根轴"无关，是基无关的判据）。返回相对 ||M||_op 的值。"""
    mv = M @ v
    return float(np.sqrt(max(v @ (M @ mv) - (v @ mv) ** 2, 0.0))) / float(np.linalg.norm(M, 2))


CASES = [
    ("锚点 Q=I", np.eye(n)),
    ("Q=diag(4,1,1,1)", np.diag([4., 1., 1., 1.])),
    ("Q=diag(1,1,1,8)", np.diag([1., 1., 1., 8.])),
    ("Q=diag(1,3,0.2,5)", np.diag([1., 3., .2, 5.])),
]
rng = np.random.default_rng(11)
M = rng.standard_normal((n, n))
CASES.append(("Q 随机 SPD", base.sym(M @ M.T) + 0.5 * np.eye(n)))

if __name__ == "__main__":
    print("== E87：阈值落零方向 vs Θ 本征轴（换 Q 的泛化检验）==")
    for nm, Q in CASES:
        Pc, K, Theta = base.ctrl(Q)
        jc = float(np.trace(W @ Pc))
        wT, UT = np.linalg.eigh(Theta)          # 升序 = 控制代价从小到大
        res = scan(Q, Theta, jc)
        print("\n--- %s ---  j_c=%.4f  lam(Theta) 升序=%s" % (
            nm, jc, np.array2string(wT, precision=3)))
        if not res:
            print("    没扫到秩跳变（窗口或容差不对），跳过")
            continue
        for dc, klo, khi, idx, v, rel in res:
            ov = np.abs(UT.T @ v) ** 2
            j = int(ov.argmax())
            print("    D_c=%8.4f  rank %d->%d（D 减小时）  落零方向命中 Theta 第 %d 轴(升序) "
                  "重叠^2=%.4f  lam(Theta)=%.4f  rel=%.2e  全重叠=%s" % (
                      dc, khi, klo, j + 1, ov[j], wT[j], rel,
                      np.array2string(ov, precision=3)))
            print("           sigma/||.||_op :  Theta %.2e   Pc %.2e   W %.2e  "
                  "=> 分辨结论：%s" % (
                      sigma(v, Theta), sigma(v, Pc), sigma(v, W),
                      "Theta 明显最贴" if sigma(v, Theta) < 0.5 * min(sigma(v, Pc), sigma(v, W))
                      else "Theta 不独占"))
        order = [int((np.abs(UT.T @ v) ** 2).argmax()) for _, _, _, _, v, _ in res]
        print("    命中轴编号（按 D_c 升序）= %s   单调上升? %s" % (
            [o + 1 for o in order], all(o1 > o2 for o1, o2 in zip(order[1:], order[:-1]))))
