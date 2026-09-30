r"""E88：E87 的两个现象（首阈值命中 lam_min(Theta)；高阶阈值失序）有没有机制解释？

候选机制：一个状态方向被丢弃的"边际价值"由两件事决定——
  (a) 它在 LQR 代价里的权重 lam_i(Theta)；
  (b) 它是不是 A 的不稳定方向（不稳定方向无论如何要为稳定性付比特，不能被丢）。
若首阈值命中 lam_min(Theta) 的前提是"该轴恰好落在 A 的稳定子空间里"，那么
  对每个 Q 检查 Theta 各本征轴在 A 稳定/不稳定子空间上的投影，就能预测
  E87 里哪些 Q 会出现"高阶阈值失序 / 重叠掉到 0.6"。
判据（可事前写下）：记 $p_i=\lVert P_{stab}u_i\rVert^2$。
  若 lam_min(Theta) 的 $p$ 接近 1（便宜且稳定）=> 首阈值干净命中；
  若最便宜的两个轴一个稳定一个不稳定（p 分裂）=> 丢弃顺序不再按 lam 排序（E87 的 [1,3,2]、[1,3,3]）。
"""
import sys
import numpy as np

sys.stdout.reconfigure(encoding="utf-8")
import p0.p0_replicate_letter as base  # noqa: E402

A, W, n = base.A, base.W, base.n
evA = np.linalg.eigvals(A)
Uall = np.linalg.eig(A)[1]


def subspaces():
    """A 的实不变子空间基（复对取实部/虚部），按模是否 >1 分成不稳定/稳定。"""
    cols_un, cols_st = [], []
    used = np.zeros(n, bool)
    for i, lam in enumerate(evA):
        if used[i]:
            continue
        if abs(lam.imag) > 1e-10:
            j = int(np.argmin(np.abs(evA - np.conj(lam))))
            blk = [np.real(Uall[:, i]), np.imag(Uall[:, i])]
            used[i] = used[j] = True
        else:
            blk = [np.real(Uall[:, i])]
            used[i] = True
        (cols_un if abs(lam) > 1.0 else cols_st).extend(blk)

    def basis(cs):
        if not cs:
            return np.zeros((n, 0))
        Q = np.array(cs).T
        U, s, _ = np.linalg.svd(Q, full_matrices=False)
        return U[:, s > 1e-9 * max(s[0], 1e-30)]

    return basis(cols_un), basis(cols_st)


UN, ST = subspaces()


def proj(v, Q):
    return float(np.sum((Q.T @ v) ** 2)) if Q.shape[1] else 0.0


def oblique(v):
    """A 的稳定/不稳定不变子空间不正交（A 非对称），正交投影 p_stab+p_unstab 会超过 1，
    不能读成分配。改用斜投影分解 v = a + b (a∈稳定, b∈不稳定)，返回两者能量占比与交叉项。"""
    Vs = np.hstack([ST, UN]) if ST.shape[1] and UN.shape[1] else None
    if Vs is None:
        return (float("nan"),) * 4
    k = ST.shape[1]
    try:
        c = np.linalg.solve(Vs, np.real(v))
    except np.linalg.LinAlgError:
        return (float("nan"),) * 4
    a, b = Vs[:, :k] @ c[:k], Vs[:, k:] @ c[k:]
    na, nb = float(a @ a), float(b @ b)
    nv = float(v @ v)
    return na / nv, nb / nv, 2 * float(a @ b) / nv, float(np.linalg.cond(Vs))


CASES = [("锚点 Q=I", np.eye(n)),
         ("Q=diag(4,1,1,1)", np.diag([4., 1., 1., 1.])),
         ("Q=diag(1,1,1,8)", np.diag([1., 1., 1., 8.])),
         ("Q=diag(1,3,0.2,5)", np.diag([1., 3., .2, 5.]))]
rng = np.random.default_rng(11)
M = rng.standard_normal((n, n))
CASES.append(("Q 随机 SPD", base.sym(M @ M.T) + 0.5 * np.eye(n)))

E87 = {"锚点 Q=I": [0.9948, 0.9742, 0.9718],
       "Q=diag(4,1,1,1)": [0.9925, 0.9923, 0.9599],
       "Q=diag(1,1,1,8)": [0.9994, 0.5631, 0.6121],
       "Q=diag(1,3,0.2,5)": [0.9999, 0.6194],
       "Q 随机 SPD": [0.9763, 0.9717]}

if __name__ == "__main__":
    print("== E88：Theta 本征轴相对 A 的稳定/不稳定子空间的斜投影分解 ==")
    print("  dim 不稳定=%d  dim 稳定=%d  两子空间主角度 cos=%s" % (
        UN.shape[1], ST.shape[1],
        np.array2string(np.linalg.svd(ST.T @ UN, compute_uv=False), precision=3)))
    print("  （cos=奇异值；接近 1 说明两子空间几乎重合，分开'稳定/不稳定'在此植物里退化）")
    for nm, Q in CASES:
        Pc, K, Theta = base.ctrl(Q)
        wT, UT = np.linalg.eigh(Theta)
        print("\n--- %s （E87 命中重叠^2 = %s）---" % (
            nm, E87.get(nm, "未测")))
        for i in range(n):
            v = UT[:, i]
            fa, fb, fx, cd = oblique(v)
            print("    轴 %d  lam(Theta)=%8.4f   斜投影：稳定份额=%.4f 不稳定份额=%.4f "
                  "交叉=%.4f（和=1）  cond[S|U]=%.2f   正交 p_stab=%.4f p_unstab=%.4f" % (
                      i + 1, wT[i], fa, fb, fx, cd, proj(v, ST), proj(v, UN)))
        print("    机制检查：最便宜轴的不稳定份额=%.4f ；"
              "不稳定份额最大的轴编号=%s" % (
                  oblique(UT[:, 0])[1],
                  [i + 1 for i in range(n) if oblique(UT[:, i])[1] > 0.5]))
    print("\n== 自我更正 ==")
    print("  E88 第一版把正交投影 p_stab/p_unstab 当分配读，出现 p_stab+p_unstab=1.55 这种")
    print("  不可能的数（A 非对称时不变子空间不正交）。上表以'斜投影份额'为准，正交列只留作对照。")
