r"""E86b：把 E86 的"看着像不像本征向量"换成可判定的谱残差。

E86 的输出有两个问题，必须自己先挑出来：
  (i) 对 $\Theta$/$P_c$/$W$ 这种**完整正交基**做子空间投影恒等于 1，那一列是废的；
  (ii) 逐向量的 $\max_i|\langle v,u_i\rangle|^2$ 只告诉我们"最接近谁"，不告诉我们
       "是不是本征向量"。真正的判据是方差为零：
           $\sigma^2(v,M)=v^\top M^2v-(v^\top Mv)^2=0\iff v$ 是 $M$ 的本征向量。
本脚本沿三个窗口逼近阈值（$D_c$ 由 §54 的 E85 截距给定，不用容差），量：
  [A] 落零方向 $v(D)$ 对五个候选对称矩阵的 $\sigma(v,M)/\|M\|_{op}$；
  [B] 方向随 $D$ 的累积旋转 $|\langle v(D),v(D_{\rm far})\rangle|^2$；
  [C] 每个 $D$ 上 $S^*$ 值域对 $\Theta$ 本征轴的分配 $\sum_{\rm keep}(v_j^\top u_i)^2$
      —— 检验"丢的是不是 $\Theta$ 代价最小的方向"这条注水式预言。
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

# §54/E85 的容差无关截距（SCS 与 CLARABEL 一致到 1e-4）
DC = {"4->3": 32.3119, "3->2": 34.4139, "2->1": 56.6493}
IDX = {"4->3": 3, "3->2": 2, "2->1": 1}
FARD = {"4->3": 31.70, "3->2": 33.50, "2->1": 54.50}


def solve(D, solver=cp.SCS, eps=1e-11):
    P = cp.Variable((n, n), symmetric=True)
    Pi = cp.Variable((n, n), symmetric=True)
    const = 0.5 * np.log(np.linalg.det(W)) / LN2
    obj = cp.Minimize(-0.5 * cp.log_det(Pi) / LN2 + const)
    AP = A @ P
    Pt = A @ P @ A.T + W
    cons = [Pi >> 0, P >> 0, cp.trace(Theta @ P) + JC <= D, P << Pt,
            cp.bmat([[P - Pi, AP.T], [AP, Pt]]) >> 0]
    prob = cp.Problem(obj, cons)
    prob.solve(solver=solver, gp=False, eps=eps, max_iters=40000)
    if P.value is None:
        return None
    Pv = base.sym(P.value)
    Pt = A @ Pv @ A.T + W
    S = base.sym(np.linalg.inv(Pv) - np.linalg.inv(Pt))
    ev, V = np.linalg.eigh(S)
    return dict(ev=ev[::-1], V=V[:, ::-1], status=prob.status,
                I=float(0.5 * np.log(np.linalg.det(Pt) / np.linalg.det(Pv)) / LN2))


def sigma(v, M):
    """v 是单位实向量，M 实对称：sigma = sqrt(v'M^2 v - (v'Mv)^2)，=0 iff v 是本征向量。"""
    mv = M @ v
    return float(np.sqrt(max(v @ (M @ mv) - (v @ mv) ** 2, 0.0))), float(np.linalg.norm(M, 2))


wT, UT = np.linalg.eigh(Theta)          # 升序：代价从小到大
wP, UP = np.linalg.eigh(Pc)
wW, UW = np.linalg.eigh(W)
ATA = A.T @ A
symTP = base.sym(Theta @ Pc)
CAND = [("Theta", Theta), ("Pc", Pc), ("W", W), ("A^T A", ATA), ("sym(Theta Pc)", symTP)]

OFFS = [2e-2, 1e-2, 5e-3, 2e-3, 5e-4, 1e-4, 2e-5, 5e-6, 1e-6]

if __name__ == "__main__":
    print("== E86b：落零方向的谱残差与注水式分配 ==")
    print("  eig(A) 模 %s   lam(Theta) 升序 %s" % (
        np.array2string(np.sort(np.abs(np.linalg.eigvals(A)))[::-1], precision=4),
        np.array2string(wT, precision=4)))

    for name, dcv in DC.items():
        idx = IDX[name]
        print("\n--- 窗口 %s：$D_c=%.4f$，盯第 %d 大特征值 ---" % (name, dcv, idx + 1))
        print("   %9s %9s %8s  %s   %8s  %s" % (
            "D", "rel", "lam_abs", "  ".join("sig/%s" % c for c, _ in CAND),
            "|v.v_far|^2", "对 Theta 本征轴重叠^2（升序）"))
        v_far = None
        r0 = solve(FARD[name])
        if r0 is not None:
            v_far = r0["V"][:, idx] / np.linalg.norm(r0["V"][:, idx])
        for off in OFFS:
            d = dcv - off
            r = solve(d)
            if r is None:
                print("   %9.4f  不可行/无解（%s）" % (d, r))
                continue
            ev = r["ev"]
            rel = ev[idx] / max(abs(ev[0]), 1e-300)
            if rel < 1e-9:
                print("   %9.4f  rel=%.2e —— 掉进求解器噪声底，停" % (d, rel))
                break
            v = r["V"][:, idx] / np.linalg.norm(r["V"][:, idx])
            sg = "  ".join("%9.2e" % (sigma(v, M)[0] / sigma(v, M)[1]) for _, M in CAND)
            ov = np.abs(UT.T @ v) ** 2
            rot = "%.5f" % float((v_far @ v) ** 2) if v_far is not None else "  --  "
            print("   %9.4f %9.2e %8.2e  %s   %8s  %s" % (
                d, rel, ev[idx], sg, rot, np.array2string(ov, precision=3)))

        # [C] 见循环后的全局扫描

    print("\n--- [C] 保留子空间对 Theta 本征轴的分配（列按 Theta 代价升序）---")
    print("    注水式预言：$D$ 变小（信息变贵）时，保留的内容应越来越集中在代价最小的轴上。")
    for d in [31.0, 32.0, 32.3119, 33.0, 34.4139, 36.0, 40.0, 45.0, 56.6493, 60.0, 70.0, 80.0]:
        r = solve(d)
        if r is None:
            print("      D=%8.4f  无解" % d)
            continue
        ev = r["ev"]
        keep = ev / max(abs(ev[0]), 1e-300) > 1e-4
        V = r["V"]
        P_of = np.array([float(np.sum([(V[:, j] @ UT[i]) ** 2 for j in range(n) if keep[j]]))
                         for i in range(n)])
        print("      D=%8.4f rank=%d  I=%.4f  分配=%s" % (
            d, int(keep.sum()), r["I"], np.array2string(P_of, precision=3)))

    print("\n== 收尾：把 E86 的错也复述一遍，免得下轮再犯 ==")
    print("  E86 表里 Theta 全体/Pc 全体/W 全体三列恒为 1.0000（完整基的子空间投影=1，无信息量），")
    print("  能用的只有不稳定/稳定子空间两列与 [A] 的逐向量重叠。")
