r"""E90b：E90 的两个自身缺陷修好之后再判"逐轴有效价格"。

E89/E90 自己写下的纪律被我自己违反了，先记账：
 ① E90 的路径 $Q=(1,1,q_3)$ 让 $\Theta$ 在不稳定平面里出**重根**（实测
    $[0.189,1.067,1.067]$、$[1.067,1.067,2.833]$……），违反纪律 21 的"目标特征值必须简单"。
    被丢的那根轴在重根情形下不唯一（$q_3\le4$ 时稳定轴恰好唯一、$q_3\ge10$ 时也唯一，
    所以标签还能读，但"最便宜的不稳定轴" $\lambda_U$ 在重根下是任取代表，$r$ 的值不可信）。
 ② 括号起点 $j_c+10^{-3}$ 在强不稳定植物上直接不可行（SCS 无解），
    导致 $\rho=2.0$ 整行、$\rho=1.4$ 的后半行全部"没扫到降秩"。这是实现失败，不是物理结论。

修法：路径取 $Q=(1,0.55,q_3)$（打破平面内各向异性），括号起点 $j_c+\max(0.05,0.02j_c)$，
并对每个设计报告 $\mathrm{gap}_{\rm drop}=$（被丢轴的本征值与最近邻之比），
$\mathrm{gap}<0.05$ 的行标"简并，忽略"。
"""
import sys
import numpy as np

sys.stdout.reconfigure(encoding="utf-8")
from p0.e89_unstable_cheap_axis import Pl, solve, rk, sigma, LN2  # noqa: E402

Q1Q2 = (1.0, 0.55)
RHOS = [1.02, 1.05, 1.15, 1.40, 2.00]
Q3S = [1.0, 2.0, 4.0, 8.0, 16.0, 32.0]


def first_threshold(pl):
    """自适应括号 + 二分（谓词 rank==n）。返回 (Dc, 落零方向, rel, I)。"""
    n = pl.n
    start = pl.jc + max(0.05, 0.02 * pl.jc)
    lo = start
    r = solve(pl, lo)
    if r is None:
        for k in range(4):
            lo = start * (2.0 ** k) + pl.jc
            r = solve(pl, lo)
            if r is not None:
                break
        if r is None:
            return None
    if rk(r) < n:
        return None
    hi = None
    for k in range(13):
        d = lo + 0.5 * (2.0 ** k)
        rr = solve(pl, d)
        if rr is None:
            continue
        if rk(rr) < n:
            hi = d
            break
        lo = d
    if hi is None:
        return None
    for _ in range(13):
        m = 0.5 * (lo + hi)
        rr = solve(pl, m)
        if rr is None:
            return None
        if rk(rr) == n:
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


if __name__ == "__main__":
    print("== E90b：非简并路径上的首阈值选择 ==")
    rows = []
    for rho in RHOS:
        print("\n--- rho=%.2f  下界 %.4f bit  Q=(1,0.55,q3) ---" % (rho, 2 * np.log(rho) / LN2))
        print("   %6s %24s %8s %10s  %s" % ("q3", "lam(Theta) 升序", "r=S/U", "D_c", "丢弃轴"))
        for q3 in Q3S:
            pl = Pl(rho=rho, Q=[Q1Q2[0], Q1Q2[1], q3])
            wT, UT = np.linalg.eigh(pl.Theta)
            UN, ST = pl.subspaces()
            pst = np.array([float(np.sum((ST.T @ UT[:, i]) ** 2)) for i in range(3)])
            si = int(np.argmax(pst))
            ui = int(np.argmin([wT[i] if pst[i] < 0.5 else np.inf for i in range(3)]))
            gap_drop = min(abs(wT[i] - wT[si]) / max(wT[si], 1e-12)
                           for i in range(3) if i != si)
            out = first_threshold(pl)
            if out is None:
                print("   %6.1f %24s  扫不到降秩（括号/可行性失败）" % (
                    q3, np.array2string(wT, precision=3)))
                continue
            dc, v, rel, I = out
            ov = np.abs(UT.T @ v) ** 2
            j = int(ov.argmax())
            lab = "稳定" if j == si else "不稳定"
            deg = "  [简并,忽略]" if gap_drop < 0.05 else ""
            print("   %6.1f %24s %8.2f %10.4f  %s轴 %d(重叠 %.4f, sig %.1e) rel=%.1e I=%.3f%s" % (
                q3, np.array2string(wT, precision=3), wT[si] / wT[ui], dc, lab, j + 1,
                ov[j], sigma(v, pl.Theta), rel, I, deg))
            if not deg:
                rows.append((rho, q3, wT[si] / wT[ui], lab, j + 1, si + 1, float(ov[j])))
        print("    （注：$r=\lambda_{\rm stab}/\lambda_{\rm unstab,\,min}$，只在这两个本征值都简单时才有意义）")

    print("\n== 判决（E90 里预先写死的 (P)/(K)）==")
    for rho in RHOS:
        sub = [x for x in rows if x[0] == rho]
        seq = ["S" if x[3] == "稳定" else "U" for x in sorted(sub, key=lambda t: t[2])]
        sw = sum(1 for a, b in zip(seq[1:], seq[:-1]) if a != b)
        print("   rho=%.2f  按 r 升序的丢弃类型 %s  切换次数=%d  切换处 r=%s" % (
            rho, "".join(seq), sw,
            [round(sub[k][2], 2) for k in range(1, len(sub))
             if sorted(sub, key=lambda t: t[2])[k][3] != sorted(sub, key=lambda t: t[2])[k - 1][3]]))
    allsw = sum(1 for rho in RHOS for _ in [0]) * 0
    rs = sorted(rows, key=lambda t: t[2])
    sw_all = sum(1 for a, b in zip(rs[1:], rs[:-1]) if a[3] != b[3])
    print("   全部 %d 个非简并设计按 $r$ 排序后的切换次数=%d（若逐轴价格成立应为 1）" % (len(rs), sw_all))
