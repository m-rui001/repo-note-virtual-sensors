"""E84：§2.4 的"秩跳变阈值" 32.31/34.41/56.66 是真分岔，还是判定容差的伪影？

背景：板面 §2.4 引的三个 D_c 来自 C 的 `.work3/c9_thresholds.py`，其判据是
"S* 的第 idx 大特征值首次跌破 1e-6*max"（相对容差硬编码），且日志里 3->2 那一行
伴随 cvxpy 的 "Solution may be inaccurate" 警告。本脚本在**我自己车道**重做：
[A] 细扫 D，打印 S* 的全部相对谱 + I_true + 实际代价 J + 求解器状态；
[B] 同一批 D 换 SCS 复解，量小特征值的**跨求解器散布**；
[C] 从扫描曲线直接插值出 D_c(tau)，tau ∈ {1e-3,1e-4,1e-6,1e-8,1e-10} ⇒ 阈值的容差敏感性；
[D] 判据：λ 是否真在有限 D 处落到 0（相对值 <1e-9 并保持），还是只按幂律/指数渐近。
只读 C 的实现，不调用它的脚本。
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
AW_const = W


def solve(D, solver=cp.CLARABEL, **kw):
    P = cp.Variable((n, n), symmetric=True)
    Pi = cp.Variable((n, n), symmetric=True)
    const = 0.5 * np.log(np.linalg.det(W)) / LN2
    obj = cp.Minimize(-0.5 * cp.log_det(Pi) / LN2 + const)
    AP = A @ P
    cons = [Pi >> 0, P >> 0,
            cp.trace(Theta @ P) + JC <= D,
            P << A @ P @ A.T + W,
            cp.bmat([[P - Pi, AP.T], [AP, A @ P @ A.T + W]]) >> 0]
    prob = cp.Problem(obj, cons)
    prob.solve(solver=solver, gp=False, **kw)
    if P.value is None or Pi.value is None:
        return None
    Pv = base.sym(P.value)
    Piv = base.sym(Pi.value)
    Pt = A @ Pv @ A.T + W
    S = base.sym(np.linalg.inv(Pv) - np.linalg.inv(Pt))
    ev = np.linalg.eigvalsh(S)[::-1]
    rel = ev / max(ev[0], 1e-300)
    # 信息 DARE 是否取等：LMI 块的第二行给出的 Schur 补 vs P 的差
    dare_gap = float(np.max(np.linalg.eigvalsh(base.sym(A @ Piv @ A.T + W - Piv))))
    return dict(I=float(0.5 * np.log(np.linalg.det(Pt) / np.linalg.det(Pv)) / LN2),
                J=float(np.trace(Theta @ Pv) + JC), ev=ev, rel=rel,
                status=prob.status, dare_gap=dare_gap,
                minPi=float(np.min(np.linalg.eigvalsh(Piv))))


def scan(Ds, solver=cp.CLARABEL, tag=""):
    rows = []
    for D in Ds:
        try:
            r = solve(D, solver)
        except Exception as e:  # noqa: BLE001
            print("   !! D=%.3f 失败 %r" % (D, e))
            continue
        if r is None:
            continue
        r["D"] = D
        rows.append(r)
    print("  [%s] 求解 %d 点，状态分布 %s" % (tag, len(rows),
          dict((s, sum(1 for x in rows if x["status"] == s))
               for s in set(x["status"] for x in rows))))
    return rows


def crossing(rows, idx, tau):
    """从扫描曲线插值：rel[idx] 首次跌破 tau 的 D（rel 随 D 单调降）。"""
    prev = None
    for r in rows:
        cur = r["rel"][idx]
        if prev is not None and prev[1] >= tau:
            if cur < tau:
                if cur <= 0.0:
                    return r["D"]
                d0, v0 = prev
                d1, v1 = r["D"], cur
                if v0 <= 0:
                    return d1
                t = (np.log(tau) - np.log(v1)) / (np.log(v0) - np.log(v1))
                return d1 + t * (d0 - d1)
        if cur > 0:
            prev = (r["D"], cur)
    return float("nan")


if __name__ == "__main__":
    print("== E84：秩跳变阈值的容差/求解器敏感性 ==")
    print("  锚点：n=%d  eig(A)=%s  tr(WPc)=%.4f  lam(Theta)=%s" % (
        n, np.round(np.abs(np.linalg.eigvals(A)), 4), JC,
        np.round(np.linalg.eigvalsh(Theta), 4)))

    Ds = np.round(np.arange(31.8, 60.01, 0.2), 3)
    rows = scan(Ds, cp.CLARABEL, "CLARABEL")
    if not rows:
        raise SystemExit("无解")

    print("\n[A] 三个阈值邻域的相对谱 rel=lam_i/lam_max")
    hdr = "   %8s %10s %10s %10s %10s %9s %9s %8s" % (
        "D", "rel2", "rel3", "rel4", "I_true", "J", "dare_gap", "status")
    print(hdr)
    for lo, hi, idx, name in ((32.0, 32.9, 3, "4->3"), (33.9, 35.0, 2, "3->2"),
                              (55.4, 57.6, 1, "2->1")):
        print("   --- %s（C 板面值在窗内）---" % name)
        for r in rows:
            if lo <= r["D"] <= hi:
                print("   %8.2f %10.2e %10.2e %10.2e %10.4f %9.4f %8.1e %8s" % (
                    r["D"], r["rel"][idx], r["rel"][min(idx + 1, 3)],
                    r["rel"][3], r["I"], r["J"], r["dare_gap"],
                    "opt" if r["status"] == "optimal" else r["status"][:8]))

    print("\n[C] 阈值对判定容差的敏感性（同一批 CLARABEL 解，只换 tau）")
    print("   %6s %10s %10s %10s" % ("tau", "Dc(4->3)", "Dc(3->2)", "Dc(2->1)"))
    for tau in (1e-3, 1e-4, 1e-6, 1e-8, 1e-10):
        vals = [crossing(rows, 3, tau), crossing(rows, 2, tau), crossing(rows, 1, tau)]
        print("   %6.0e %10s %10s %10s" % (tau, *["%10.3f" % v if v == v else "%10s" % "--" for v in vals]))

    print("\n[D] 小特征值是否真落到 0：记录 rel 的最小读数与是否出现非正值")
    for idx, name in ((3, "第4大"), (2, "第3大"), (1, "第2大")):
        vals = np.array([r["rel"][idx] for r in rows])
        neg = int(np.sum(np.array([r["ev"][idx] for r in rows]) <= 0))
        print("   %s：min(rel)=%.3e @D=%.2f  非正读数 %d/%d  末端三点 rel=%s" % (
            name, vals.min(), rows[int(vals.argmin())]["D"], neg, len(rows),
            np.array2string(vals[-3:], precision=2)))

    print("\n[B] 跨求解器散布（SCS 复解同一批 D 的子集）")
    windows = ((32.0, 32.9, 3), (33.9, 35.0, 2), (55.4, 57.6, 1))
    sub = [(r["D"], idx) for r in rows for lo, hi, idx in windows if lo <= r["D"] <= hi][::2]
    print("   %8s %12s %12s %10s %10s" % ("D", "rel CLAR", "rel SCS", "I CLAR", "I SCS"))
    for D, idx in sub:
        rc = next((x for x in rows if x["D"] == D), None)
        try:
            rs = solve(D, cp.SCS, eps=1e-9, max_iters=8000)
        except Exception as e:  # noqa: BLE001
            print("   %8.2f  SCS 失败 %r" % (D, e))
            continue
        if rc is None or rs is None:
            continue
        print("   %8.2f %12.2e %12.2e %10.4f %10.4f" % (
            D, rc["rel"][idx], rs["rel"][idx], rc["I"], rs["I"]))
