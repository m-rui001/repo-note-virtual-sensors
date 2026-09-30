"""E85：把 §2.4 的"秩跳变阈值"从容差读数改造成可复核的量。

E84 看到两件事：(i) 最小的那几个 λ_i(S*) 从不在有限 D 处落到 0，而是掉到求解器噪声底；
(ii) 但在噪声底之上，λ 随 D **近似线性**下降。若 λ(D) = c·(D_c − D)^α 且 α=1，
则 D_c 可用线性外推**无容差**地读出，§2.4 的三个数才有地位。

本脚本：
[A] 对三个窗口、两个求解器，拟合 log λ = log c + α·log(D_c − D)（网格搜索 D_c、OLS 出 α）；
[B] 同时拟合"线性版"（α 固定为 1）给 D_c 的置信邻域，并与 §2.4 的板面值对账；
[C] 把"地板读数"D_min(r)（rate→0）与"秩读数"（1e-6 容差）分开打印 ⇒ 直接量出 §2.4 表 r*(D) 的系统偏差。
"""
import sys
import numpy as np

sys.stdout.reconfigure(encoding="utf-8")
import cvxpy as cp  # noqa: E402
import p0.p0_replicate_letter as base  # noqa: E402

LN2 = np.log(2.0)
A, W, n = base.A, base.W, base.n
Pc, K, Theta = base.ctrl(np.eye(n))
JC = float(np.trace(W @ Pc))


def solve(D, solver=cp.CLARABEL, **kw):
    P = cp.Variable((n, n), symmetric=True)
    Pi = cp.Variable((n, n), symmetric=True)
    const = 0.5 * np.log(np.linalg.det(W)) / LN2
    obj = cp.Minimize(-0.5 * cp.log_det(Pi) / LN2 + const)
    AP = A @ P
    cons = [Pi >> 0, P >> 0, cp.trace(Theta @ P) + JC <= D, P << A @ P @ A.T + W,
            cp.bmat([[P - Pi, AP.T], [AP, A @ P @ A.T + W]]) >> 0]
    prob = cp.Problem(obj, cons)
    prob.solve(solver=solver, gp=False, **kw)
    if P.value is None:
        return None
    Pv = base.sym(P.value)
    Pt = A @ Pv @ A.T + W
    ev = np.linalg.eigvalsh(base.sym(np.linalg.inv(Pv) - np.linalg.inv(Pt)))[::-1]
    return dict(I=float(0.5 * np.log(np.linalg.det(Pt) / np.linalg.det(Pv)) / LN2),
                rel=ev / max(ev[0], 1e-300), st=prob.status)


def fit(win, idx, rows):
    """rows: list of (D, rel)。窗口内取噪声底以上的点，网格搜 D_c，OLS 出 α。"""
    pts = [(d, r[idx]) for d, r in rows if win[0] <= d <= win[1] and r[idx] > 3e-9]
    if len(pts) < 5:
        return None
    D = np.array([p[0] for p in pts])
    y = np.log(np.array([p[1] for p in pts]))
    best = None
    for Dc in np.arange(D.max() + 1e-4, D.max() + 4.0, 1e-4):
        x = np.log(Dc - D)
        b, a = np.polyfit(x, y, 1)          # y = a + b·x，b 即 α
        res = y - (a + b * x)
        sse = float(res @ res)
        if best is None or sse < best[0]:
            best = (sse, Dc, b, float(np.sqrt(sse / max(len(x) - 2, 1))))
    sse, Dc, alpha, rms = best
    x = np.log(Dc - D)
    r_lin = 1.0 - sse / float(((y - y.mean()) ** 2).sum())
    # 线性（α=1）版本：y = log c + log(Dc − D) ⇒ z = y − log(Dc−D) 应为常数
    lin = [np.exp(yk + np.log(Dc - dk)) for dk, yk in zip(D, y)]   # 反推 λ·(Dc−D)^{-1}... 用比值代替
    return dict(npts=len(pts), Dc=float(Dc), alpha=float(alpha), rms=rms,
                R2=float(r_lin), span=(float(D.min()), float(D.max())),
                lin_cv=float(np.std(lin) / np.mean(lin)), lin=lin, D=D, pts=pts)


BOARD = {3: 32.312, 2: 34.414, 1: 56.656}     # §2.4（来自 .work3/c9_thresholds.py，tau=1e-6）

if __name__ == "__main__":
    print("== E85：阈值的标度律拟合（α=1 ⇒ 有限 D_c 是真分岔；α→∞/无解 ⇒ 容差伪影）==")
    windows = ((31.7, 32.9, 3, "4->3"), (33.5, 35.2, 2, "3->2"), (54.5, 58.0, 1, "2->1"))
    for solver, tag, kw in ((cp.CLARABEL, "CLARABEL", {}), (cp.SCS, "SCS(eps=1e-10)", dict(eps=1e-10, max_iters=20000))):
        print("\n--- 求解器 %s ---" % tag)
        for lo, hi, idx, name in windows:
            Ds = np.round(np.arange(lo, hi + 1e-9, (hi - lo) / 34.0), 4)
            rows = []
            for d in Ds:
                try:
                    r = solve(d, solver, **kw)
                except Exception as e:  # noqa: BLE001
                    print("   !! %.4f %r" % (d, e))
                    continue
                if r is not None:
                    rows.append((d, r["rel"]))
            f = fit((lo, hi), idx, rows)
            if f is None:
                print("   %s：噪声底以上点数不足，无法拟合" % name)
                continue
            print("   %s  n=%d  窗口 D∈[%.3f,%.3f]  D_c(外推)=%.4f  α=%.4f  R2=%.6f  rms(log)=%.3f" % (
                name, f["npts"], f["span"][0], f["span"][1], f["Dc"], f["alpha"], f["R2"], f["rms"]))
            print("        线性残差 CV=%.3f（α 固定为 1 时的相对散布）  板面值 %.3f ⇒ 差 %+.3f" % (
                f["lin_cv"], BOARD[idx], f["Dc"] - BOARD[idx]))

    print("\n[C] 地板读数 vs 1e-6 容差读数（rate→0 ⇒ 秩的几何定义）")
    for idx, name, lo, hi in ((3, "4->3", 31.55, 31.75), (2, "3->2", 32.15, 32.30), (1, "2->1", 41.58, 41.78)):
        vals = []
        for d in np.round(np.arange(lo, hi, 0.02), 3):
            r = solve(float(d), cp.CLARABEL)
            if r:
                vals.append((d, r["I"], r["rel"][idx]))
        first = next((v for v in vals if v[1] < 1e-6), None)
        last = next((v for v in reversed(vals) if v[1] >= 1e-6), None)
        print("   %s：I_true<1e-6 的最小 D = %s；其 λ_rel 读数 = %s；上一个仍持秩的 D = %s" % (
            name,
            "%.3f" % first[0] if first else "窗内未达到",
            "%.2e" % first[2] if first else "--",
            "%.3f (I=%.2e, rel=%.2e)" % (last[0], last[1], last[2]) if last else "--"))
        print("        扫过的 rel 列：%s" % np.array2string(np.array([v[2] for v in vals]), precision=1))
