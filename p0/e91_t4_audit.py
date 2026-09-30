r"""E91：审计 E89 的唯一反例 T4——它是"首阈值丢了不稳定轴"，还是我的括号把两个跳变混成一个？

E90b 之后局势变了：28/28 个非简并设计的首阈值都丢**稳定**轴（$r$ 从 $0.02$ 到 $14.9$ 都不切换），
加上 E89 的 C1/C2/T1/T2/T3 与 §55 锚点的 5 个 $Q$，规则
  (c)「首阈值丢弃的方向总在 $A$ 的稳定子空间里」有 33/34 的支持，
唯一反例是 E89 的 T4（$Q=(0.2,5,20)$，丢的方向与全不稳定的第 1 轴重叠 $0.9805$，
而它当时 $\sigma_\Theta=1.2e{-}1$ 也是全场最大）。

T4 用的是 E89 的粗扫（步长 1.0）$+$ 谓词 `rk==khi` 的二分。若两个跳变落在同一步长里
（rank $3\to2\to1$ 之间不到 1.0），二分可能收敛到**第二**个跳变却按第一个解释，
于是"阈值下方 $5e{-}3$"取的 $\lambda_{\min}$ 方向就已经是前一个跳变丢掉的（重叠 $0.98$ 但不算数）。
本脚本在 $D\in[40.0,42.5]$ 以 $0.05$ 步长把秩与三个特征值全部打出来，直接看有几个跳变、各在哪。
"""
import sys
import numpy as np

sys.stdout.reconfigure(encoding="utf-8")
from p0.e89_unstable_cheap_axis import Pl, solve, rk, sigma  # noqa: E402

pl = Pl(Q=[.2, 5., 20.])
wT, UT = np.linalg.eigh(pl.Theta)
UN, ST = pl.subspaces()
pst = np.array([float(np.sum((ST.T @ UT[:, i]) ** 2)) for i in range(3)])
print("== E91：T4 细扫，$Q=(0.2,5,20)$，j_c=%.4f ==" % pl.jc)
print("   lam(Theta) 升序 = %s   各轴稳定份额 = %s" % (
    np.array2string(wT, precision=3), np.array2string(pst, precision=3)))
prev = None
for d in np.arange(40.0, 42.55, 0.05):
    r = solve(pl, round(float(d), 4), eps=1e-11, it=40000)
    if r is None:
        print("   D=%.2f  无解" % d)
        continue
    ev = r["ev"]
    rel = ev / max(abs(ev[0]), 1e-300)
    k = rk(r)
    v = r["V"][:, k - 1] / np.linalg.norm(r["V"][:, k - 1])
    ov = np.abs(UT.T @ v) ** 2
    lab = ["稳定" if i == int(np.argmax(pst)) else "不稳定" for i in [int(ov.argmax())]][0]
    mark = ""
    if prev is not None and k != prev:
        mark = " <== 秩跳变 %d->%d（D 减小方向）" % (prev, k)
    prev = k
    print("   D=%.2f rank=%d ev=%s  最小非零方向命中轴 %d (%s) 重叠=%.4f sig=%.1e%s" % (
        d, k, np.array2string(rel, precision=2), int(ov.argmax()) + 1, lab,
        ov.max(), sigma(v, pl.Theta), mark))
