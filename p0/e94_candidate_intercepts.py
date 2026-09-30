r"""E94：给 E93 的 (X) 配第二把尺——"候选阈值"用 §54 的无容差外推截距量，而不是只看跳变间距。

E93 落在 (X)：$\Delta=D_{c2}-D_{c1}$ 在身份切换处收缩到 $0.454$，与 $\max\mathrm{ov}$ 最小点错位 $0.025$。
但 (X) 的判据只是**相关性**，而且 $\Delta$ 量的是"已实现的两次降秩间距"，不是"两个反事实的第一丢弃间距"。
本脚本把两支候选方向的**外推截距**直接算出来：§54 已证 $\lambda_{\min}(S^*(D))$ 在阈值下方对 $D$ **线性**
（$\alpha=0.988\text{–}1.032$、$R^2\ge0.9999$、截距容差无关）。于是对 rank 3 分支的两根最小本征值
$\lambda_3(D)<\lambda_2(D)$ 各拟合 $\lambda=a(D-\hat D)$，得到
  $\hat D_3$ ＝ 第一支候选的阈值（应当 $\approx D_{c1}$，用来验方法自身），
  $\hat D_2$ ＝ **第二支候选在 rank-3 分支里的"反事实"阈值**，
  $\hat\Delta=\hat D_2-\hat D_3$。
若 $\hat\Delta$ 与实现的 $\Delta$ 一致，"丢弃轴 $=\arg\min_v D_c(v)$"就有两把独立的尺；
若不一致（例如 $\lambda_2$ 分支根本不外推到它的真实跳变），(X) 只剩跳变间距一把尺。

判据（跑前写死）：
  路径 A $Q=(q_1,5,20)$、路径 B $Q=(q_1,2,20)$，$q_1\in\{0.25,0.275,0.30,0.325,0.35,0.40\}$。
  **(T) 两把尺一致**：$|\hat\Delta-\Delta|\le0.5\,\Delta$ 且 $\hat\Delta$ 的最小点与 $\Delta$ 的最小点、
        与 $\max\mathrm{ov}$ 的最小点三者在同一路径上重合（网格内同点或相邻点），**两条路径都成立**
        ⇒ 可以写成相图：$\hat\Delta=0$ 的零集就是规则 (c) 的失效边界（本车道第一个可检验的边界条件）。
  **(V) 只有一把尺**：上式不成立，或路径 B 上没有翻转（$\hat\Delta$ 无最小点/标签全程同一侧）
        ⇒ 把 (X) 降为"间距型经验指标"，**终止**这条线，转 `note.tex` 定稿与对 C 的复现审计。
  硬性纪律：拟合窗口取 $D\in[D_{c1}-3.5,\,D_{c1}-1.0]$（离阈值 $\ge1.0$，避开纪律 20 的特征值噪声区）；
  $R^2<0.999$ 的行标"非线性"、不参与判决；采样点秩不足 3 ⇒ 弃点。
"""
import sys
import numpy as np

sys.stdout.reconfigure(encoding="utf-8")
from p0.e89_unstable_cheap_axis import Pl, solve, rk, sigma  # noqa: E402
from p0.e93_threshold_crossing import jumps, refine  # noqa: E402

RHO = 1.15
Q1S = [0.25, 0.275, 0.30, 0.325, 0.35, 0.40]
PATHS = [("A", 5.0), ("B", 2.0)]
DSHIFT = np.arange(1.0, 3.6, 0.5)


def branch(pl, UT, D):
    """返回 (lam2, lam3, v3, ov3, I)；rank 不足 3 时 None。"""
    r = solve(pl, D, eps=1e-10, it=20000)
    if r is None or rk(r) < pl.n:
        return None
    ev = r["ev"]
    v3 = r["V"][:, pl.n - 1] / np.linalg.norm(r["V"][:, pl.n - 1])
    return float(ev[1]), float(ev[2]), v3, np.abs(UT.T @ v3) ** 2, r["I"]


def intercept(xs, ys):
    """lambda = a*(D - Dhat) 的最小二乘：返回 Dhat, a, R^2。"""
    if len(xs) < 4:
        return None
    X = np.asarray(xs, float)
    Y = np.asarray(ys, float)
    a, c = np.polyfit(X, Y, 1)
    if abs(a) < 1e-14:
        return None
    yh = a * X + c
    ss = float(np.sum((Y - Y.mean()) ** 2))
    r2 = 1.0 - float(np.sum((Y - yh) ** 2)) / ss if ss > 0 else float("nan")
    return -c / a, a, r2


if __name__ == "__main__":
    print(r"== E94：候选阈值的外推截距 vs 实现的跳变间距 ==")
    allrows = {}
    for nm, q2 in PATHS:
        print(r"\n--- 路径 %s：$Q=(q_1,%.0f,20)$ ---" % (nm, q2))
        print("  %6s %9s %9s %8s %9s %9s %8s %7s  %s" % (
            "q1", "Dc1", "Dc2", "delta", "Dhat_3", "Dhat_2", "Dhat_delta", "R2_min", "ov@(Dc1-1.0) / 标签"))
        rows = []
        for q1 in Q1S:
            pl = Pl(rho=RHO, Q=[q1, q2, 20.0])
            wT, UT = np.linalg.eigh(pl.Theta)
            UN, ST = pl.subspaces()
            pst = [float(np.sum((ST.T @ UT[:, i]) ** 2)) for i in range(pl.n)]
            si = int(np.argmax(pst))
            js = jumps(pl, span=30.0)
            if len(js) < 2:
                print("  %6.3f  只取到 %d 个跳变（弃行）" % (q1, len(js)))
                continue
            d1, d2 = refine(pl, js[0]), refine(pl, js[1])
            if d1 is None or d2 is None:
                print("  %6.3f  二分失败" % q1)
                continue
            l2s, l3s, ds = [], [], []
            ov0 = sg0 = None
            lab = "?"
            for off in DSHIFT:
                b = branch(pl, UT, d1 - float(off))
                if b is None:
                    continue
                l2, l3, v3, ov, I = b
                l2s.append(l2)
                l3s.append(l3)
                ds.append(d1 - float(off))
                if abs(float(off) - 1.0) < 1e-9:
                    ov0, sg0 = ov, sigma(v3, pl.Theta)
                    lab = "S" if int(ov.argmax()) == si else "U"
            i3, i2 = intercept(ds, l3s), intercept(ds, l2s)
            if i3 is None or i2 is None:
                print("  %6.3f %9.4f %9.4f %8.4f   拟合点不足（%d 个）" % (q1, d1, d2, d2 - d1, len(ds)))
                continue
            r2m = min(i3[2], i2[2])
            flag = "" if r2m >= 0.999 else "  [非线性,弃行]"
            print("  %6.3f %9.4f %9.4f %8.4f %9.4f %9.4f %8.4f %7.4f  %s %s%s" % (
                q1, d1, d2, d2 - d1, i3[0], i2[0], i2[0] - i3[0], r2m,
                np.array2string(ov0, precision=3) if ov0 is not None else "NA", lab, flag))
            if not flag:
                rows.append(dict(q1=q1, delta=d2 - d1, dhat=i2[0] - i3[0], lab=lab,
                                 mx=float(ov0.max()) if ov0 is not None else 1.0))
        allrows[nm] = rows

    print("\n== 判决 ==")
    ok = True
    for nm, rows in allrows.items():
        if not rows:
            print("   路径 %s：没有可用行 => (V)" % nm)
            ok = False
            continue
        devs = [abs(t["dhat"] - t["delta"]) / max(t["delta"], 1e-9) for t in rows]
        md = min(rows, key=lambda t: t["delta"])
        mh = min(rows, key=lambda t: t["dhat"])
        mo = min(rows, key=lambda t: t["mx"])
        seq = "".join(t["lab"] for t in rows)
        flips = sum(1 for a, b in zip(seq[1:], seq[:-1]) if a != b)
        agree = (max(devs) <= 0.5) and (md["q1"] == mh["q1"] == mo["q1"] or
                                         abs(md["q1"] - mh["q1"]) <= 0.026 and abs(md["q1"] - mo["q1"]) <= 0.026)
        print("   路径 %s：标签序列 %s（翻转 %d 次）  max|Dhat-delta|/delta = %.3f  "
              "delta 最小 q1=%.3f / Dhat 最小 q1=%.3f / maxov 最小 q1=%.3f  => %s" % (
                  nm, seq, flips, max(devs), md["q1"], mh["q1"], mo["q1"], "一致" if agree else "不一致"))
        ok = ok and agree and flips >= 1
    print("\n   >>> %s" % ("(T) 两把尺一致且都在切换点取最小：可写相图 hat_delta=0 = 规则(c) 失效边界"
                           if ok else
                           "(V) 只有一把尺或路径 B 无翻转：(X) 降为间距型经验指标，此线终止，转 note.tex 定稿"))
