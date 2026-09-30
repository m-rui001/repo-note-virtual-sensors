"""E24：Prop E 的系数公式在**第二个独立构型**上是否成立——顺便量原文 Fig.5 那条"垂直渐近线"。

文献核对（必须先做，因为它砍掉了我一部分新颖性）：1510.04214 §V [page 7] 原文——
    "The vertical asymptote D = Tr(WS) corresponds to the best achievable control performance
     when unrestricted amount of information about the state is available."
即**无约束**问题的 DI(D) 在 $D\downarrow D_{\min}^{unc}=\mathrm{tr}(WP_c)$ 处发散，作者已明说（垂直渐近线）。
所以"$D\to$ 地板时率发散"不是新观察。Prop E 若要站住，新东西只能在：
  (i) 把这条渐近线**搬位置**：任务受限下 $D_{\min}(\mathcal K_F)=D_{\min}^{unc}+\Phi(V;Q)$；
  (ii) 给出**系数** $\tfrac{\mathrm{rank}}{2}$（原文只画曲线，未给指数）；
  (iii) 由此得出 L-CSS 2026 Theorem 3 的失败模式是一段**不可行区间**而非"有损失"。

本脚本独立检验 (ii)：无约束构型 rank$=n=4$、$\Phi=0$，预言
    $\mathrm{DI}(D)=\tfrac{n}{2}\log_2\tfrac1x+O(1)$，$x=D-\mathrm{tr}(WP_c)$，斜率 $=6.6439$ bit/十倍程。
用 `unconstrained_sdp` 直接求最优（不是参数化近似），扫 $x$ 从 $10^{1}$ 到 $10^{-3}$。
另附与已发表锚点的同台对照：$D=33\to6.133$，$D=40\to3.266$，$D=80\to1.602$ bits/sample。
"""
import sys
import numpy as np
sys.stdout.reconfigure(encoding='utf-8')
from p0.p0_replicate_letter import A, W, n, sym, ln2, ctrl, rate_cost, unconstrained_sdp

Pc, K, Th = ctrl(np.eye(n))
FLOOR = float(np.trace(W @ Pc))
PUB = {33.0: 6.133, 40.0: 3.266, 80.0: 1.602}      # 1510.04214 §V 正文数值


def iso_branch(x_list):
    """F=I 的各向同性分支：给定 x 反解 γ，给出可行 (x, DI) 作为最优值的上界参照。"""
    Fi = np.eye(n)
    out = []
    for x in x_list:
        lo, hi = -3., 14.
        f = lambda lg: rate_cost(Fi, np.exp(lg) * np.eye(n), Th, Pc)[1] - (FLOOR + x)
        if f(lo) * f(hi) > 0:
            out.append((x, np.nan))
            continue
        for _ in range(120):
            mid = .5 * (lo + hi)
            if f(mid) * f(lo) <= 0:
                hi = mid
            else:
                lo = mid
        out.append((x, rate_cost(Fi, np.exp(.5 * (lo + hi)) * np.eye(n), Th, Pc)[0]))
    return out


if __name__ == '__main__':
    print('无约束地板 tr(WP_c)=%.5f   预言斜率 (n/2)log2(10)=%.4f bit/十倍程 (n=%d)'
          % (FLOOR, (n / 2.) * np.log2(10.), n))
    xs = [1e1, 5e0, 2e0, 1e0, 3e-1, 1e-1, 3e-2, 1e-2, 3e-3, 1e-3]
    rows, iso = [], dict(iso_branch(xs))
    print('\n%10s %10s %10s %10s' % ('x', 'DI_SDP', 'DI_iso', '局部指数'))
    prev = None
    for x in xs:
        D = FLOOR + x
        try:
            u = unconstrained_sdp(D, np.eye(n), Th, Pc)
            di = u['I_true']
        except Exception as e:
            print('%10.4g  求解异常 %s' % (x, e))
            continue
        if not np.isfinite(di):
            print('%10.4g  非有限，跳过' % x)
            continue
        loc = np.nan if prev is None else np.log(di / prev[1]) / np.log(x / prev[0])
        rows.append((x, di))
        print('%10.4g %10.4f %10s %10s'
              % (x, di, '%.4f' % iso[x] if np.isfinite(iso[x]) else '—',
                 '%+0.3f' % loc if np.isfinite(loc) else '—'))
        prev = (x, di)

    print('\n与 1510.04214 §V 已发表锚点同台对照（同一条无约束曲线）：')
    print('%8s %10s %10s %10s' % ('D', '原文 DI', '我的 DI', '差'))
    for D, ref in sorted(PUB.items()):
        try:
            mine = unconstrained_sdp(D, np.eye(n), Th, Pc)['I_true']
            print('%8.0f %10.3f %10.4f %+10.3f' % (D, ref, mine, mine - ref))
        except Exception as e:
            print('%8.0f  异常 %s' % (D, e))

    for lo_x in (1e-3, 1e-2, 3e-2):
        sel = [(x, d) for x, d in rows if lo_x <= x <= 1e0]
        if len(sel) >= 4:
            lx = np.log10(1. / np.array([a for a, _ in sel]))
            dd = np.array([b for _, b in sel])
            a, b = np.polyfit(lx, dd, 1)
            res = np.max(np.abs(dd - (a * lx + b)))
            print('\n窗口 x∈[%.0e,1]  %d 点：DI = %+.4f·log₁₀(1/x) %+.4f，最大残差 %.3f bit；'
                  '预言系数 %.4f，误差 %+.1f%%'
                  % (lo_x, len(sel), a, b, res, (n / 2.) * np.log2(10.),
                     100 * (a / ((n / 2.) * np.log2(10.)) - 1)))
            print('  同窗口幂律指数（log₂DI 对 log₁₀(1/x)）= %+.4f（→0 即对数律）'
                  % np.polyfit(lx, np.log2(dd), 1)[0])
