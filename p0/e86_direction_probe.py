r"""E86：秩阈值 $D_c$ 有没有"谱式闭式"的可能？——查落零方向的几何。

§54 把 $D_c$ 做成容差无关的截距后，剩下的唯一问题是它能不能写成"$\lambda_i(\cdot)$ 追上某个水位"。
这种闭式**要求**消失方向是某个固定对称矩阵的本征向量（否则不存在"第 $i$ 个特征值"这种编号）。
本脚本直接把 $S^*(D)$ 最小本征方向的实/复投影量出来，与四个候选基对表：
  (a) $\Theta=K^\top(R+B^\top P_cB)K$ 的本征基（注水式闭式的天然候选）；
  (b) 控制 DARE 解 $P_c$ 的本征基；
  (c) $A$ 的不稳定不变子空间（实 Schur 不稳定块的实/虚部张成）；
  (d) 过程噪声 $W$ 的本征基。
判据：$\max_j|\langle v,u_j\rangle|^2$ 接近 1 且**随 $D$ 稳定** ⇒ 该方向是固定谱向量，谱闭式有可能；
若方向本身随 $D$ 旋转（相邻 $D$ 的 $|\langle v(D),v(D')\rangle|^2<1$ 明显），则"$\gamma_i$ 追上水位"这类陈述**没有指称对象**。
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


def solve(D, solver=cp.SCS, **kw):
    P = cp.Variable((n, n), symmetric=True)
    Pi = cp.Variable((n, n), symmetric=True)
    const = 0.5 * np.log(np.linalg.det(W)) / LN2
    obj = cp.Minimize(-0.5 * cp.log_det(Pi) / LN2 + const)
    AP = A @ P
    cons = [Pi >> 0, P >> 0, cp.trace(Theta @ P) + JC <= D, P << A @ P @ A.T + W,
            cp.bmat([[P - Pi, AP.T], [AP, A @ P @ A.T + W]]) >> 0]
    prob = cp.Problem(obj, cons)
    prob.solve(solver=solver, gp=False, eps=1e-11, max_iters=30000, **kw)
    if P.value is None:
        return None
    Pv = base.sym(P.value)
    Pt = A @ Pv @ A.T + W
    S = base.sym(np.linalg.inv(Pv) - np.linalg.inv(Pt))
    ev, V = np.linalg.eigh(S)
    return dict(ev=ev[::-1], V=V[:, ::-1],
                I=float(0.5 * np.log(np.linalg.det(Pt) / np.linalg.det(Pv)) / LN2))


def orth_basis(mats):
    """把一组（可能复的）向量拆成实部/虚部，SVD 取列空间正交基（秩由奇异值 >1e-9 判定）。"""
    cols = []
    for M in mats:
        for j in range(M.shape[1]):
            c = M[:, j]
            if np.abs(np.imag(c)).max() > 1e-9:
                cols += [np.real(c), np.imag(c)]
            else:
                cols.append(np.real(c))
    if not cols:
        return np.zeros((n, 0))
    Q = np.array(cols).T
    U, s, _ = np.linalg.svd(Q, full_matrices=False)
    return U[:, s > 1e-9 * max(s[0], 1e-30)]


def proj_sq(v, Q):
    if Q.shape[1] == 0:
        return float("nan")
    return float(np.sum((Q.T @ v) ** 2))


# 四个候选基
wT, UT = np.linalg.eigh(Theta)
wP, UP = np.linalg.eigh(Pc)
wW, UW = np.linalg.eigh(W)
evA = np.linalg.eigvals(A)
Uall = np.linalg.eig(A)[1]
Unst = Uall[:, np.abs(evA) > 1.0 + 1e-12]
Stb = Uall[:, np.abs(evA) <= 1.0 + 1e-12]
BASIS = [("Theta 全体", orth_basis([UT])), ("Pc 全体", orth_basis([UP])),
         ("W 全体", orth_basis([UW])), ("A 不稳定子空间", orth_basis([Unst])),
         ("A 稳定子空间", orth_basis([Stb]))]

WINDOWS = ((31.7, 32.3, 3, "4->3"), (33.5, 34.4, 2, "3->2"), (54.5, 56.6, 1, "2->1"))

if __name__ == "__main__":
    print("== E86：落零方向的几何（谱闭式有没有指称对象）==")
    print("  eig(A) 模 = %s   不稳定维数 = %d" % (
        np.array2string(np.sort(np.abs(evA))[::-1], precision=4), Unst.shape[1]))
    print("  lam(Theta) = %s   lam(Pc) = %s" % (
        np.array2string(wT, precision=4), np.array2string(wP, precision=4)))
    prev = {}
    for lo, hi, idx, name in WINDOWS:
        Ds = np.round(np.arange(lo, hi + 1e-9, (hi - lo) / 12.0), 4)
        print("\n--- 窗口 %s（消失的是第 %d 大的特征值）---" % (name, idx + 1))
        print("   %8s %10s %9s  %s" % ("D", "lam_rel", "lam_abs", "  ".join("%-14s" % b for b, _ in BASIS)))
        for d in Ds:
            r = solve(float(d))
            if r is None:
                continue
            rel = r["ev"] / max(abs(r["ev"][0]), 1e-300)
            if rel[idx] < 1e-4:
                continue
            v = r["V"][:, idx]
            v = v / np.linalg.norm(v)
            ps = [proj_sq(v, Q) for _, Q in BASIS]
            key = idx
            dot = "|v(D)·v(D_prev)|^2 = %.4f" % float(prev[key] @ v) ** 2 \
                if key in prev else "首点"
            prev[key] = v
            print("   %8.3f %10.3e %9.3e  %s   %s" % (
                d, rel[idx], r["ev"][idx],
                "  ".join("%-14.4f" % p for p in ps), dot))

    print("\n" + r'[A] 消失方向的分量（$|v_i|$，按坐标）与它最接近的 $\Theta$ 本征向量')
    for lo, hi, idx, name in WINDOWS:
        Ds = np.round(np.arange(lo, hi + 1e-9, (hi - lo) / 6.0), 4)
        for d in Ds[-1:]:
            r = solve(float(d))
            if r is None or (r["ev"] / max(abs(r["ev"][0]), 1e-300))[idx] < 1e-4:
                continue
            v = r["V"][:, idx] / np.linalg.norm(r["V"][:, idx])
            ov = np.abs(UT.T @ v) ** 2
            print("   %s D=%.3f  |v|=%s   对 Theta 本征向量的重叠^2=%s（最大 %.4f @lam=%.4f）" % (
                name, d, np.array2string(np.abs(v), precision=3),
                np.array2string(ov, precision=3), ov.max(), wT[int(ov.argmax())]))
