r"""E92：把 T4->T3 那一步连续化——首阈值丢弃的轴到底是"切换"还是"旋转"？

§57 留下的事实：规则 (c)「首阈值丢稳定子空间里的轴」33/34，唯一反例 T4 $Q=(0.2,5,20)$
被 E91 细扫证实是**真**的（窗口内只有一个跳变，读数随 $D$ 单调缓变），而 T3 $Q=(1,5,20)$
与它**只差 $q_1$ 一个参数**、稳定轴价格逐位相同（$\lambda_2=5.849$），丢弃轴却相反。
⇒ 候选解释只剩两类，且它们对"沿路径怎么走"给出**不同**的可测预言：

  (B) 突变／相边界：存在 $q_1^*$，其上丢弃轴从"稳定"翻到"不稳定"，两边的方向各自钉在一根 $\Theta$ 轴上
      （$\max_j\mathrm{ov}_j\ge0.98$）。若边界由 $\lambda_1/\lambda_2$ 定价，则第二关：
      在另一条路径 $Q=(q_1,2,20)$ 上切换点应落在同一个 $\lambda_1/\lambda_2$（$\pm20\%$）。
  (C) 连续／旋转：某个 $q_1$ 上丢弃方向**同时**不贴两根轴（前两轴重叠都 $<0.95$）
      ⇒ "第几根轴"没有指称对象，本车道这条线（"丢哪个轴、为什么"）当场作废并转向。

判据在跑之前写死（不许事后改）：
  取 $\rho=1.15$，$q_1\in\{0.2,0.25,0.3,0.4,0.5,0.65,0.8,1.0\}$，$Q=(q_1,5,20)$。
  每个设计：粗扫 $[j_c+0.3,\,j_c+45]$ 步长 $0.5$ **数清秩跳个数**并确认二分取到的是第一个；
  二分首跳变（谓词 rank$=n$）；在 $D_c-5\mathrm e{-}3$ 与 $D_c-0.30$ 各做一次 $\epsilon=10^{-11}$ 精解，
  打印丢弃方向对**三根轴**的完整重叠三元组与 $\sigma_\Theta$。
  窗口内出现 $\ge2$ 个跳变 ⇒ 该行标"多跳变"；$\mathrm{gap}_{\rm drop}<0.05$ ⇒ 标"简并"（纪律 21/22）。
  若 (B)：报告 $q_1^*$ 与切换处的 $\lambda_1/\lambda_2$，并跑第二条路径复核；
  若 (C)：立刻自我降级并把 §57 读法 1 的"总是 $\Theta$ 本征向量"改成"近似"。
"""
import sys
import numpy as np

sys.stdout.reconfigure(encoding="utf-8")
from p0.e89_unstable_cheap_axis import Pl, solve, rk, sigma, LN2  # noqa: E402

RHO = 1.15
Q2Q3 = (5.0, 20.0)
Q1S = [0.2, 0.25, 0.3, 0.4, 0.5, 0.65, 0.8, 1.0]
FINE_EPS, FINE_IT = 1e-11, 40000


def screen(pl):
    """纯 numpy：Theta 谱、各轴稳定份额、被丢候选轴的简单性。"""
    wT, UT = np.linalg.eigh(pl.Theta)
    UN, ST = pl.subspaces()
    pst = np.array([float(np.sum((ST.T @ UT[:, i]) ** 2)) for i in range(pl.n)])
    si = int(np.argmax(pst))
    ui = int(np.argmin([wT[i] if pst[i] < 0.5 else np.inf for i in range(pl.n)]))
    gap = min(abs(wT[i] - wT[si]) / max(wT[si], 1e-12) for i in range(pl.n) if i != si)
    return wT, UT, si, ui, pst, gap


def first_threshold(pl):
    """粗扫数跳变 + 二分第一个跳变。返回 (Dc, jumps, n_jumps)。"""
    n, lo = pl.n, pl.jc + 0.3
    prev, jumps = None, []
    d = lo
    while d <= lo + 45.0:
        r = solve(pl, round(float(d), 4))
        if r is not None:
            k = rk(r)
            if prev is not None and k != prev[1]:
                jumps.append((prev[0], d, prev[1], k))
            prev = (d, k)
        d += 0.5
    if not jumps:
        return None
    a, b, klo, khi = jumps[0]
    for _ in range(13):
        m = 0.5 * (a + b)
        r = solve(pl, m)
        if r is None:
            return None
        if rk(r) == klo:
            a = m
        else:
            b = m
    return 0.5 * (a + b), jumps, len(jumps)


def readout(pl, UT, D):
    r = solve(pl, D, eps=FINE_EPS, it=FINE_IT)
    if r is None:
        return None
    idx = pl.n - 1
    v = r["V"][:, idx] / np.linalg.norm(r["V"][:, idx])
    ov = np.abs(UT.T @ v) ** 2
    return ov, sigma(v, pl.Theta), r["ev"][idx] / max(abs(r["ev"][0]), 1e-300), r["I"]


if __name__ == "__main__":
    print(r"== E92：$\rho=%.2f$，$Q=(q_1,%.0f,%.0f)$ 上首阈值丢弃轴的连续性/突变性 ==" % (RHO, *Q2Q3))
    print("   %5s %22s %10s %6s %8s %9s  %s" % ("q1", "lam(Theta) 升序", "各轴稳定份额",
                                                "lam1/lam2", "D_c", "跳变数", "丢弃方向（精解）"))
    rows = []
    for q1 in Q1S:
        pl = Pl(rho=RHO, Q=[q1, Q2Q3[0], Q2Q3[1]])
        wT, UT, si, ui, pst, gap = screen(pl)
        out = first_threshold(pl)
        if out is None:
            print("   %5.2f %22s  扫不到降秩" % (q1, np.array2string(wT, precision=3)))
            continue
        dc, jumps, nj = out
        parts = []
        for off in (5e-3, 0.30):
            rr = readout(pl, UT, dc - off)
            if rr is None:
                parts.append("精解失败")
                continue
            ov, sg, rel, I = rr
            j = int(ov.argmax())
            parts.append("D-%.3f: 轴%d(%s) ov=%s sig=%.1e rel=%.1e I=%.3f" % (
                off, j + 1, "稳定" if j == si else "不稳定",
                np.array2string(ov, precision=3), sg, rel, I))
            if off == 5e-3:
                rows.append((q1, float(wT[0] / wT[1]), si + 1, j + 1, float(ov[j]), nj, dc))
        deg = "  [简并]" if gap < 0.05 else ""
        multi = "  [多跳变: %s]" % ";".join("%.1f:%d->%d" % (b, kl, kh) for _, b, kl, kh in jumps[1:]) if nj > 1 else ""
        print("   %5.2f %22s %10s %8.3f %9.4f %9d  %s%s%s" % (
            q1, np.array2string(wT, precision=3), np.array2string(pst, precision=2),
            wT[0] / wT[1], dc, nj, " || ".join(parts), deg, multi))

    rs = [x for x in rows if x[5] == 1]
    print("\n== 判决（跑前写死）==")
    mx = [x[4] for x in rs]
    if any(x[4] < 0.95 for x in rs):
        print(r"   (C) 连续/旋转：有设计的 $\max_j\mathrm{ov}_j<0.95$（min=%.4f）" % min(mx))
        print(r"       => '第几根轴' 无指称对象，§57 读法 1 降为'近似'，本线终止并转向。")
    else:
        labs = ["S" if x[2] == x[3] else "U" for x in rs]
        sw = [k for k in range(1, len(labs)) if labs[k] != labs[k - 1]]
        print(r"   max ov 全场最小 = %.4f（判据 (B) 要求 $\ge0.98$）  标签序列（$q_1$ 升序）=%s  切换处=%s" % (
            min(mx), "".join(labs), [rs[k][0] for k in sw]))
        if len(sw) == 1:
            k = sw[0]
            print(r"   (B) 单相边界：$q_1^*\in(%.2f,%.2f)$，切换处 $\lambda_1/\lambda_2\in[%.4f,%.4f]$" % (
                rs[k - 1][0], rs[k][0], rs[k][1], rs[k - 1][1]))
            print(r"   第二关：在 $Q=(q_1,2,20)$ 上复核切换是否落在同一 $\lambda_1/\lambda_2$（$\pm20\%$）")
        elif not sw:
            print("   路径上没切换（全 %s）——说明这一段还在同一侧，需要扩窗口或换 q2" % labs[0])
        else:
            print("   切换次数 %d != 1：非单调，按 (C) 的精神处理（轴标签不可靠）" % len(sw))
