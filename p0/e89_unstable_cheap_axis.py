r"""E89：§55(5) 预先注册的证伪实验——让 $\lambda_{\min}(\Theta)$ 落在**不稳定**方向上。

§55 的两条经验事实：首阈值（满秩 $\to$ 满秩$-1$）丢的方向总是 $\Theta$ 的 $\lambda_{\min}$ 轴（5/5），
而该轴在 E88 里 5/5 都"几乎全在 $A$ 的稳定子空间"（不稳定份额 $\le0.09$）。
⇒ "便宜"与"稳定"在这 5 个设计里**混在一起**，无法分辨谁在起作用。

E89 把它拆开：构造块对角植物，使稳定/不稳定子空间正交（不稳定 $=\operatorname{span}\{e_1,e_2\}$、
稳定 $=\operatorname{span}\{e_3\}$），再用 $Q$ 把 $\lambda_{\min}(\Theta)$ 轴推到不稳定平面里。

预先写下的判据（不许事后改）：
  记首阈值 $D_c^{*}$ 的落零方向 $v$，$s=\lVert P_{\rm unstab}u_{\min}\rVert^2$。
  (T) 若 $s>0.5$ 的组里仍有 $|\langle v,u_{\min}\rangle|^2>0.9$：
      ⇒ "便宜"是因，"稳定"不是；命题更强，且必须解释稳定性从哪来（猜想：LMI $P\prec APA^\top+W$
        与 $S\succeq0$ 自动禁止丢弃不稳定方向，于是那种设计的首阈值会**推迟**到更大的 $D$）。
  (F) 若 $s>0.5$ 的组里 $v$ 改贴到别的轴（尤其贴稳定轴）：
      ⇒ §55 的"首阈值 $=\lambda_{\min}$ 轴"降级为"便宜 $\wedge$ 稳定"的合取，立刻自我降级。
对照组：$s<0.2$（应与 §55 一致，命中 $\lambda_{\min}$ 轴），用来证明我的构造没把问题搞坏。
"""
import sys
import numpy as np

sys.stdout.reconfigure(encoding="utf-8")
import cvxpy as cp  # noqa: E402
from scipy.linalg import solve_discrete_are  # noqa: E402

LN2 = np.log(2.0)
TOL = 1e-4


def sym(M):
    return 0.5 * (M + M.T)


class Pl:
    """块对角植物：前两块是旋转放大（不稳定），第三块是收缩（稳定）。"""

    def __init__(self, rho=1.15, th=0.7, sig=0.55, b3=0.0, b13=0.0, w=(1., 1., 1.), Q=None):
        self.A = np.array([[rho * np.cos(th), -rho * np.sin(th), 0.],
                           [rho * np.sin(th), rho * np.cos(th), 0.],
                           [0., 0., sig]])
        # 必须 m=n=3：K 是 m x n，Θ=K'(R+B'PcB)K 只有在全输入时才满秩；
        # m<n 时 lam_min(Theta)=0（ker Theta 方向代价恒零），"便宜轴"退化成平凡情形。
        self.B = np.array([[1., 0., b13], [0., 1., 0.], [b3, 0., 1.]])
        self.W = sym(np.diag(w) + 0.15 * np.array([[0., .3, .2], [.3, 0., .1], [.2, .1, 0.]]))
        self.Q = np.diag(Q)
        self.n = 3
        self.m = 3
        self.Pc = solve_discrete_are(self.A, self.B, sym(self.Q), np.eye(3))
        R = np.eye(3)
        S = R + self.B.T @ self.Pc @ self.B
        self.K = np.linalg.solve(sym(S) + 1e-12 * np.eye(self.m), self.B.T @ self.Pc @ self.A)
        self.Theta = sym(self.K.T @ S @ self.K)
        self.jc = float(np.trace(self.W @ self.Pc))

    def subspaces(self):
        ev = np.linalg.eigvals(self.A)
        V = np.linalg.eig(self.A)[1]
        un, st = [], []
        used = np.zeros(3, bool)
        for i, lam in enumerate(ev):
            if used[i]:
                continue
            if abs(lam.imag) > 1e-10:
                j = int(np.argmin(np.abs(ev - np.conj(lam))))
                blk = [np.real(V[:, i]), np.imag(V[:, i])]
                used[i] = used[j] = True
            else:
                blk = [np.real(V[:, i])]
                used[i] = True
            (un if abs(lam) > 1.0 else st).extend(blk)

        def bs(cs):
            Qm = np.array(cs).T
            U, s, _ = np.linalg.svd(Qm, full_matrices=False)
            return U[:, s > 1e-9 * s[0]]

        return bs(un), bs(st)


def solve(pl, D, eps=1e-9, it=8000):
    n = pl.n
    A, W, Theta = pl.A, pl.W, pl.Theta
    P = cp.Variable((n, n), symmetric=True)
    Pi = cp.Variable((n, n), symmetric=True)
    const = 0.5 * np.log(np.linalg.det(W)) / LN2
    AP = A @ P
    Pt = A @ P @ A.T + W
    prob = cp.Problem(cp.Minimize(-0.5 * cp.log_det(Pi) / LN2 + const),
                      [Pi >> 0, P >> 0, cp.trace(Theta @ P) + pl.jc <= D, P << Pt,
                       cp.bmat([[P - Pi, AP.T], [AP, Pt]]) >> 0])
    try:
        prob.solve(solver=cp.SCS, gp=False, eps=eps, max_iters=it)
    except Exception:
        return None
    if P.value is None:
        return None
    Pv = sym(P.value)
    Pt = A @ Pv @ A.T + W
    S = sym(np.linalg.inv(Pv) - np.linalg.inv(Pt))
    ev, V = np.linalg.eigh(S)
    return dict(ev=ev[::-1], V=V[:, ::-1],
                I=float(0.5 * np.log(np.linalg.det(Pt) / np.linalg.det(Pv)) / LN2))


def rk(r):
    return int(np.sum(r["ev"] / max(abs(r["ev"][0]), 1e-300) > TOL))


def sigma(v, M):
    mv = M @ v
    return float(np.sqrt(max(v @ (M @ mv) - (v @ mv) ** 2, 0.0))) / float(np.linalg.norm(M, 2))


def first_threshold(pl, UN):
    """粗扫找满秩->满秩-1 的跳变（最小的那个 D_c），二分后取阈值下方 5e-3 的方向。"""
    lo0 = pl.jc + 0.3
    grid = np.arange(lo0, lo0 + 60.0, 1.0)
    prev = None
    for d in grid:
        r = solve(pl, d)
        if r is None:
            continue
        k = rk(r)
        if prev is not None and k != prev[1]:
            lo, up, klo, khi = prev[0], d, prev[1], k
            break
        prev = (d, k)
    else:
        return None
    for _ in range(9):
        m = 0.5 * (lo + up)
        r = solve(pl, m)
        if r is None:
            return None
        if rk(r) == khi:
            up = m
        else:
            lo = m
    dc = 0.5 * (lo + up)
    r = solve(pl, dc - 5e-3, eps=1e-11, it=40000)
    if r is None:
        return None
    idx = min(klo, pl.n) - 1
    v = r["V"][:, idx] / np.linalg.norm(r["V"][:, idx])
    wT, UT = np.linalg.eigh(pl.Theta)
    ov = np.abs(UT.T @ v) ** 2
    j = int(ov.argmax())
    return dict(dc=dc, klo=klo, khi=khi, v=v, ov=ov, j=j, wT=wT,
                sig=sigma(v, pl.Theta), sigP=sigma(v, pl.Pc), sigW=sigma(v, pl.W),
                un_share=float(np.sum((UN.T @ UT[:, 0]) ** 2)),
                rel=r["ev"][idx] / max(r["ev"][0], 1e-300), I=r["I"])


CAND = [
    # 事先筛好的设计：s = |lambda_min(Theta) 轴| 在不稳定平面里的份额；gap = (lam2-lam1)/lam1
    # gap>0.25 保证 lambda_min simple（否则"第 1 轴"不唯一，判据失去指称对象）
    ("C1 对照 s=0.00 gap=2.16 Q=(0.2,1,1)", dict(Q=[.2, 1., 1.])),
    ("C2 对照 s=0.00 gap=2.93 Q=(0.2,5,1)", dict(Q=[.2, 5., 1.])),
    ("T1 便宜轴全不稳定 s=1.00 gap=1.23 Q=(0.2,1,5)", dict(Q=[.2, 1., 5.])),
    ("T2 便宜轴全不稳定 s=1.00 gap=1.45 Q=(0.2,1,20)", dict(Q=[.2, 1., 20.])),
    ("T3 便宜轴全不稳定 s=1.00 gap=2.34 Q=(1,5,20)", dict(Q=[1., 5., 20.])),
    ("T4 便宜轴全不稳定 s=1.00 gap=6.87 Q=(0.2,5,20)", dict(Q=[.2, 5., 20.])),
]

if __name__ == "__main__":
    print("== E89：把'便宜'与'稳定'拆开——lambda_min(Theta) 落在不稳定方向时的首阈值 ==")
    for nm, kw in CAND:
        pl = Pl(**kw)
        UN, ST = pl.subspaces()
        wT, UT = np.linalg.eigh(pl.Theta)
        s = float(np.sum((UN.T @ UT[:, 0]) ** 2))
        tag = "TEST(s>0.5)" if s > 0.5 else ("CONTROL(s<0.2)" if s < 0.2 else "中间地带")
        print("\n--- %s ---  j_c=%.4f  lam(Theta)=%s  lambda_min 轴的不稳定份额 s=%.4f  [%s]" % (
            nm, pl.jc, np.array2string(wT, precision=3), s, tag))
        res = first_threshold(pl, UN)
        if res is None:
            print("    没找到满秩->降秩的跳变（窗口内秩没变过）")
            continue
        floor = float(np.sum(np.log(np.abs(np.linalg.eigvals(pl.A))[
            np.abs(np.linalg.eigvals(pl.A)) > 1.0])) / LN2)
        print(r"    首阈值 $D_c=%.4f$（rank %d->%d，$D$ 减小时）rel=%.2e   "
              r"该处 $I=%.4f$ bit，稳定化下界 $\sum\log_2|\lambda_u|=%.4f$ bit" % (
                  res["dc"], res["klo"], res["khi"], res["rel"], res["I"], floor))
        print("    落零方向：命中 Theta 第 %d 轴（升序）重叠^2=%.4f  全重叠=%s" % (
            res["j"] + 1, res["ov"][res["j"]], np.array2string(res["ov"], precision=3)))
        print("    sigma/||.||：Theta %.2e  Pc %.2e  W %.2e   I=%.4f" % (
            res["sig"], res["sigP"], res["sigW"], res["I"]))
        if s > 0.5:
            verdict = "(T) 便宜轴虽不稳定仍被丢弃" if res["j"] == 0 and res["ov"][0] > 0.9 \
                else "(F) 改贴别的轴 => 降级为'便宜 and 稳定'合取" if res["j"] != 0 \
                else "(F?) 命中 1 轴但重叠 %.3f <=0.9 => 也算降级" % res["ov"][0]
            print("    >>> 预注册判决 %s" % verdict)
