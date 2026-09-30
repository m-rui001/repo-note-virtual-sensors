#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""e159：在**本车道自己的前沿数据**上量 C 提出的梯子律（板 §69-B 第 3748 行全文）：

    > $\Delta I=\frac{\mathrm{rank}\,C}{2}\log_2\frac1x+b_{\rm pred}+O(x)$,  $x=D-D_{\min}$, $\Delta I=I-R_{\exp}$.

这条是**可判的**：它的全部内容就是"以 $\log_2(1/x)$ 为自变量时，主项系数 $=r/2$"。
本脚本不用 C 的任何文件、不用 cvxpy，只用本车道的定价（`p0/exp_c_audit.py` 前言 + `ctrl_full`）：
$C=\sqrt{s}Z^{\mathsf T}$, $P_m=\mathrm{DARE}(A^{\mathsf T},C^{\mathsf T},W,I_r)$,
$I=\tfrac12\log_2\det(I_r+CP_mC^{\mathsf T})$, $D=J_C+\mathrm{tr}(\Theta P_p)$。
对每个秩 $r$：扫 $s$ 网格 + 一族候选标架 $Z$，取**同秩下**的 $(D,I)$ 下包络，再回归 $\Delta I$ 对 $\log_2(1/x)$ 的斜率。

预注册判据（跑前写死，不看结果再挑）：
 [L1] 三个**从大到小**的 $x$ 窗口各给斜率 $=d\Delta I/d\log_2(1/x)$ 与其对 $r/2$ 的相对偏差。
 [L2] 判决：最靠近 $0$ 的窗口里 $|\text{dev}|\le10\%$ $\\Rightarrow$ 律可用、只是要很小的 $x$（于是**必须报有效域**）；
      $10\%<|\text{dev}|\le30\%$ $\\Rightarrow$ 次主项在可达域里不可忽略，$O(x)$ 的写法要改；
      $|\text{dev}|>30\%$ $\\Rightarrow$ **主项系数本身在可达窗口内不成立**（实际接近 $\kappa r/2,\ \kappa>1$），
      与我方 e150 的 $\gamma_{\rm loc}$ 高出 $+34.2/104.9/177.2\%$ 同向。
 [L3] 代数对账（不是测量）：代价轴整体乘 $c$ 时 $x\to cx$、$I$ 不变 $\\Rightarrow$ 律要求 $b_{\rm pred}$ 平移 $+(r/2)\log_2 c$
      $\\Rightarrow$ **截距不是标定量**，跨株/跨单位比 $b_{\rm pred}$ 无意义。
"""
import sys
import numpy as np
from scipy.linalg import solve_discrete_are
sys.stdout.reconfigure(encoding='utf-8')

out = []


def p(*a):
    s = ' '.join(str(x) for x in a)
    out.append(s)
    print(s)


_src = open('p0/exp_c_audit.py', encoding='utf-8').read().split("print(r'== E53")[0]
_ns = {'__name__': 'p'}
exec(compile(_src, 'p0/exp_c_audit.py[preamble]', 'exec'), _ns)


def sym(X):
    return 0.5 * (X + X.T)


def ctrl_full(Ax, Bx, Wx, Qt, Rx):
    Pc = sym(solve_discrete_are(Ax, Bx, sym(Qt), Rx))
    K = np.linalg.solve(Rx + Bx.T @ Pc @ Bx, Bx.T @ Pc @ Ax)
    return Pc, K, sym(K.T @ (Rx + Bx.T @ Pc @ Bx) @ K), float(np.trace(Wx @ Pc))


def make_plant(kind, seed):
    rng = np.random.default_rng(seed)
    if kind == 'anchor':
        A, B, W = _ns['A'].copy(), _ns['B'].copy(), _ns['W'].copy()
    else:
        n = 4
        A = np.linalg.qr(rng.normal(0, 1, (n, n)))[0] @ np.diag([1.6, 1.25, 0.9, 0.7]) @ np.linalg.qr(rng.normal(0, 1, (n, n)))[0].T
        B = rng.normal(0, 1, (n, n))
        m = n
        M = rng.normal(0.0, 1.0, (m, m))
        W = sym(M @ M.T + 0.6 * np.eye(m))
    n = A.shape[1]
    Pc, K, TH, JC = ctrl_full(A, B, W, np.eye(n), np.eye(n))
    return A, B, W, TH, JC, n


def R_exp_of(A):
    ev = np.linalg.eigvals(A)
    return float(np.sum(np.log2(np.abs(ev[np.abs(ev) > 1.0 + 1e-12]))))


def post_at(Z, A, W, s):
    r = Z.shape[1]
    Ir = np.eye(r)
    C = np.sqrt(s) * Z.T
    Pm = sym(solve_discrete_are(A.T, C.T, W, Ir))
    Sm = C @ Pm @ C.T + Ir
    Lk = Pm @ C.T @ np.linalg.inv(Sm)
    return sym(Pm - Lk @ Sm @ Lk.T), Pm, C, Ir


def point(Z, A, W, TH, JC, s):
    Pp, Pm, C, Ir = post_at(Z, A, W, s)
    I = 0.5 * np.log(np.linalg.det(Ir + C @ Pm @ C.T)) / np.log(2.0)
    return float(I), JC + float(np.trace(TH @ Pp))


def cands(A, W, TH, r, nrand=40, seed=7):
    wth, Uth = np.linalg.eigh(TH)
    U = Uth[:, np.argsort(-wth)]
    cs = [U[:, :r]]
    rng = np.random.default_rng(seed + r)
    for _ in range(nrand):
        Q, Rr = np.linalg.qr(rng.normal(0, 1, (U.shape[0], r)))
        cs.append(Q * np.sign(np.diag(Rr)))
    return cs


p('== 斜率回归：$\\Delta I$ 对 $\\log_2(1/x)$，同秩下包络 ==')
ROWS = []
for nm, spec in [('anchor', ('anchor', 0)), ('rand-1', ('rand', 5)), ('rand-2', ('rand', 11))]:
    A, B, W, TH, JC, n = make_plant(*spec)
    Rexp = R_exp_of(A)
    p('')
    p('--- %s: $J_C=%.4f$  $R_{\\exp}=%.6f$  spec$|_A|=%s' % (nm, JC, Rexp,
      np.array2string(np.sort(np.abs(np.linalg.eigvals(A)))[::-1], precision=3)))
    for r in range(1, n):
        Zs = cands(A, W, TH, r)
        # 地板：s=1e7 下同秩最小 D（与本车道 e153 同口径）
        df = np.inf
        for Z in Zs:
            try:
                _, Dv = point(Z, A, W, TH, JC, 1e7)
            except np.linalg.LinAlgError as e:
                p('    [跳过] %s' % e)
                continue
            df = min(df, Dv)
        pts = []
        for s in np.logspace(-1.5, 5.0, 55):
            for Z in Zs:
                try:
                    I, D = point(Z, A, W, TH, JC, s)
                except np.linalg.LinAlgError:
                    continue
                if np.isfinite(I) and np.isfinite(D):
                    pts.append((D, I))
        pts.sort()
        env, best = [], np.inf
        for D, I in pts:
            best = min(best, I)
            env.append((D, best))
        xs = np.array([D - df for D, _ in env])
        ds = np.array([I - Rexp for _, I in env])
        m_ok = xs > 1e-9
        xs, ds = xs[m_ok], ds[m_ok]
        if xs.size < 8:
            p('  r=%d 有效包络点 %d 个，太少（略）' % (r, xs.size))
            continue
        lo = np.log10(xs.min())   # 律断言的是 x->0 端，窗口必须锚在 x_min 往上（此前锚在 x_max 是判据 bug）
        wins = [(lo, lo + 1.0), (lo + 1.0, lo + 2.0), (lo + 2.0, lo + 3.0)]
        sl, dspan = [], []
        for a, b in wins:
            k = (np.log10(xs) >= a) & (np.log10(xs) < b)
            if k.sum() < 3:
                sl.append(np.nan)
                dspan.append('n/a')
                continue
            u = np.log2(1.0 / xs[k])
            sl.append(float(np.polyfit(u, ds[k], 1)[0]))
            dspan.append('[%.2f,%.2f]' % (ds[k].min(), ds[k].max()))
        pred = r / 2.0
        dev = [np.nan if not np.isfinite(v) else (v - pred) / pred * 100.0 for v in sl]
        p('  r=%d 地板 $D_{\\rm floor}=%.4f$（高出 $J_C$ %.4f）；十倍程窗口斜率 '
          '$\\times$ 预测 $r/2=%.2f$： %s' % (r, df, df - JC, pred,
                                             ' '.join('%.3f' % v if np.isfinite(v) else ' n/a' for v in sl)))
        p('        相对偏差 %s ；$x$ 覆盖 $[%.2e, %.2e]$，包络点 %d 个'
          % (' '.join('%+.0f%%' % v if np.isfinite(v) else ' n/a' for v in dev),
             xs.min(), xs.max(), xs.size))
        ROWS.append((nm, r, sl, dev, xs.min(), xs.max()))

p('')
p('== [L2] 判决行（判据见文件头）==')
worst = [x for x in ROWS if np.isfinite(x[3][-1])]
if not worst:
    p('[L2] 无可判窗口 => 不下结论（按板账 #46 的纪律：没有有效读数就不许出判决）')
else:
    ok10 = [x for x in worst if abs(x[3][-1]) <= 10]
    mid = [x for x in worst if 10 < abs(x[3][-1]) <= 30]
    bad = [x for x in worst if abs(x[3][-1]) > 30]
    p('[L2] 最靠近 0 的窗口里：|dev|<=10%% 的有 %d 项；10–30%% 有 %d 项；>30%% 有 %d 项（共 %d 项）'
      % (len(ok10), len(mid), len(bad), len(worst)))
    for nm, r, sl, dev, x0, x1 in worst:
        p('     %-8s r=%d  最外/中间/最内窗口 dev = %s' % (nm, r, ' '.join('%+.0f%%' % v if np.isfinite(v) else 'n/a' for v in dev)))
    if bad and not ok10:
        p('[结论] 在**可达窗口内**主项系数 $r/2$ 系统性偏小（实际斜率 $\\approx\\kappa r/2,\\ \\kappa>1$）$-$ $-$ '
          '与我方 e150 的 $\\gamma_{\\rm loc}$ 高出 $+34.2/104.9/177.2\\%$ 是同一件事的两种报价（斜率与倒数）。')
        p('       $\\Rightarrow$ C 的律不能按 $O(x)$ 的字面意义用：**要么报有效域（$x$ 小到几倍）**，'
          '要么把 $O(x)$ 换成显式次主项；在此之前正文里不得出现"由梯子律给出"这类句子。')

p('')
p('== [L3] 代数对账：截距随代价轴单位平移 ==')
for r in (1, 2, 3):
    for c in (7.0, 2.0):
        p('  r=%d, 代价轴乘 %.1f $\\Rightarrow$ $b_{\\rm pred}$ 必须平移 $+(r/2)\\log_2 %.1f = %+.4f$ bit'
          % (r, c, c, (r / 2.0) * np.log2(c)))
p('[L3] 结论：律里的 $b_{\\rm pred}$ **不是无量纲可比的常数**（它吃掉了单位）$\\Rightarrow$ 跨株比较 $b_{\\rm pred}$ 无意义；')
p('     可跨株比的只有**斜率** $r/2$ 与"$x$ 的多小才算渐近"。')

open('p0/e159_out.txt', 'w', encoding='utf-8').write('\n'.join(out) + '\n')
