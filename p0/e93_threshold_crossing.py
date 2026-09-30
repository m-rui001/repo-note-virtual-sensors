r"""E93：E92 抓到 (C)（丢弃方向在 $q_1\approx0.30$ 连续旋转）之后的机制问题。

E92 的表（$\rho=1.15$，$Q=(q_1,5,20)$，稳定轴恒为第 2 轴、$\lambda_2=5.849$ 逐位不动）：

  q1=0.20 丢轴1(不稳定) ov=[.980,.002,.018]   Dc1=41.153  次跳变 44.1  (Δ=2.95)
  q1=0.25 丢轴1(不稳定) ov=[.968,.013,.020]   Dc1=42.876  次跳变 44.2  (Δ=1.33)
  q1=0.30 **混叠**     ov=[.343,.649,.008]   Dc1=43.934  次跳变 44.8  (Δ=0.87)  <- $D_c-0.30$ 处还是 [.776,.208]
  q1=0.40 丢轴2(稳定)  ov=[.022,.978,.001]   Dc1=44.287  次跳变 46.9  (Δ=2.61)
  q1=1.00 丢轴2(稳定)  ov=[.001,.999,.000]   Dc1=45.777  次跳变 68.1  (Δ=22.3)

⇒ 首阈值本身光滑单调（$41.15\to45.78$，切换处没有台阶，与 §54"值函数无奇性"一致），
但**首/次阈值的间距 $\Delta$ 在混叠点上收缩到最小**（$0.87$，两侧都更大）。
这给出一个机制假设：**"丢哪根轴"是候选阈值族 $\{D_c(v)\}$ 的 argmin 问题**，
规则 (c) 只是"稳定候选胜出"的那一侧；切换发生在两个候选阈值**交叉**处，
交叉点附近 SDP 近乎无差别，于是落零方向是两根轴的混合物（纪律 20 的条件数警告在此以
$S$ 的两根最小非零本征值并置的形式重现）。

判据（跑前写死，两支都可执行）：
  在 $q_1\in\{0.25,0.275,0.30,0.325,0.35,0.375,0.40\}$ 上二分出**前两个**跳变 $D_{c1}<D_{c2}$，
  记 $\Delta=D_{c2}-D_{c1}$，并在 $D_{c1}-5\mathrm e{-}3$ 精解取丢弃方向的完整重叠三元组
  （同时记录 $S$ 的最小/次小非零本征值之比 $\gamma$）。
  **(X) 交叉机制成立**：$\Delta$ 在某个 $q_1^*$ 处取严格最小（且 $\min\Delta<1.0$），
    且该点正是三元组最不平衡处（$\max_j\mathrm{ov}_j$ 最小）⇒ 我可以写
    "丢弃轴 $=\arg\min$ 候选阈值族；规则 (c) 是交叉的一侧"，并给出可检验的边界条件 $\Delta=0$。
  **(Y) 不成立**：$\Delta$ 的最小点与最不平衡点错开（$>0.05$ 的 $q_1$ 差），或 $\Delta$ 全程 $\ge2$
    ⇒ 混叠与次阈值无关，改判"$S$ 本征值并置导致的标签失效"，此线（哪个轴/为什么）终止，转去
    把 §54–§57 的成果写进 `note.tex` 并审计 C 的复现。
  硬性纪律：任一设计取不到两个跳变 ⇒ 打印"只取到一个"，该行不参与 $\Delta$ 排序。
"""
import sys
import numpy as np

sys.stdout.reconfigure(encoding="utf-8")
from p0.e89_unstable_cheap_axis import Pl, solve, rk, sigma, LN2  # noqa: E402

RHO = 1.15
Q2Q3 = (5.0, 20.0)
Q1S = [0.25, 0.275, 0.30, 0.325, 0.35, 0.375, 0.40]
FINE_EPS, FINE_IT = 1e-11, 40000
TOL_RANK = 1e-4


def jumps(pl, span=45.0, step=0.5):
    """粗扫列出所有 (下界, 上界, rank_hi, rank_lo) 跳变，按 D 升序。"""
    lo = pl.jc + 0.3
    prev, out = None, []
    d = lo
    while d <= lo + span:
        r = solve(pl, round(float(d), 4))
        if r is not None:
            k = rk(r)
            if prev is not None and k != prev[1]:
                out.append([prev[0], d, prev[1], k])
            prev = (d, k)
        d += step
    return out


def refine(pl, br):
    """把 [a,b] 内的跳变二分到中点。谓词：rank == br[2]（左侧秩）。"""
    a, b, klo = br[0], br[1], br[2]
    for _ in range(14):
        m = 0.5 * (a + b)
        r = solve(pl, m)
        if r is None:
            return None
        if rk(r) == klo:
            a = m
        else:
            b = m
    return 0.5 * (a + b)


def drop_vec(pl, UT, D):
    r = solve(pl, D, eps=FINE_EPS, it=FINE_IT)
    if r is None:
        return None
    idx = pl.n - 1
    v = r["V"][:, idx] / np.linalg.norm(r["V"][:, idx])
    ov = np.abs(UT.T @ v) ** 2
    ev = r["ev"]
    nz = ev[ev / max(abs(ev[0]), 1e-300) > TOL_RANK]
    gamma = float(nz[-1] / nz[0]) if len(nz) >= 2 else float("nan")
    return ov, sigma(v, pl.Theta), gamma, r["I"]


if __name__ == "__main__":
    print(r"== E93：候选阈值交叉机制，$\rho=%.2f$，$Q=(q_1,%.0f,%.0f)$ ==" % (RHO, *Q2Q3))
    print("   %6s %10s %10s %8s %8s %8s  %s" % ("q1", "Dc1", "Dc2", "delta", "gamma", "sig_Th", "丢弃方向重叠三元组"))
    rows = []
    for q1 in Q1S:
        pl = Pl(rho=RHO, Q=[q1, Q2Q3[0], Q2Q3[1]])
        wT, UT = np.linalg.eigh(pl.Theta)
        UN, ST = pl.subspaces()
        si = int(np.argmax([float(np.sum((ST.T @ UT[:, i]) ** 2)) for i in range(pl.n)]))
        js = jumps(pl)
        if len(js) < 2:
            print("   %6.3f  只取到 %d 个跳变（不参与 delta 排序）" % (q1, len(js)))
            continue
        d1, d2 = refine(pl, js[0]), refine(pl, js[1])
        if d1 is None or d2 is None:
            print("   %6.3f  二分失败" % q1)
            continue
        rr = drop_vec(pl, UT, d1 - 5e-3)
        if rr is None:
            print("   %6.3f  精解失败" % q1)
            continue
        ov, sg, gm, I = rr
        j = int(ov.argmax())
        print("   %6.3f %10.4f %10.4f %8.4f %8.1e %8.1e  %s  轴%d(%s) I=%.3f" % (
            q1, d1, d2, d2 - d1, gm, sg, np.array2string(ov, precision=3),
            j + 1, "稳定" if j == si else "不稳定", I))
        rows.append(dict(q1=q1, d1=d1, d2=d2, delta=d2 - d1, gamma=gm, mx=float(ov.max()),
                         sg=sg, lab="S" if j == si else "U"))

    if rows:
        dmin = min(rows, key=lambda t: t["delta"])
        omax = min(rows, key=lambda t: t["mx"])
        print("\n== 判决 ==")
        print("   delta 最小点 q1=%.3f (delta=%.4f)   max ov 最小点 q1=%.3f (max ov=%.4f)   错位=%.3f" % (
            dmin["q1"], dmin["delta"], omax["q1"], omax["mx"], abs(dmin["q1"] - omax["q1"])))
        shift = abs(dmin["q1"] - omax["q1"])
        if dmin["delta"] < 1.0 and shift <= 0.05:
            print(r"   (X) 交叉机制成立：Δ 严格最小且与最不平衡点重合（≤0.05）"
                  r" ⇒ 丢弃轴 = 候选阈值族的 argmin；规则 (c) 只是交叉的一侧。")
        else:
            print(r"   (Y) 交叉机制不成立（Δmin=%.3f / 错位=%.3f）⇒ 改判 S 本征值并置导致标签失效，"
                  r"此线终止。" % (dmin["delta"], shift))
        print("   序列（q1 升序）：" + "  ".join("%.3f:%s Δ%.2f ovmax%.3f" % (
            t["q1"], t["lab"], t["delta"], t["mx"]) for t in rows))
