#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""e159b：给 C 的梯子律**量出有效域**，并与本车道 e137b 的 $\gamma_{\rm loc}$ 表对账。

律（板 §69-B 第 3748 行）：$\Delta I=\frac{r}{2}\log_2\frac1x+b_{\rm pred}+O(x)$，$x=D-D_{\min}$，$\Delta I=I-R_{\exp}$。
e159 已证主项系数在最靠近 $0$ 的三个十倍程窗口 $|\text{dev}|\le1\%$（9/9 条梯级）。
本轮不再问"成不成立"，只问三件事：

 [L4] 有效域：从 $x_{\min}$ 往上一路扫十倍程窗口，每个窗口报斜率 $A$、对 $r/2$ 的偏差、点数；
      给出**最后一个仍满足 $|\text{dev}|\le10\%$ 的连续窗口上界** $x\le x_{\min}\cdot10^{j+1}$，
      并翻成 $\Delta I$ 的下界（小 $x$ = 大 $\Delta I$）。$O(x)$ 是"逐点有意义"还是只在域内有用，由这一行判。
 [L5] 口径对账：e137b 定义 $\gamma_{\rm loc}:=-\mathrm{d}\ln x/\mathrm{d}\Delta I$，
      故恒等式 $\gamma_{\rm loc}=\ln2/A$，且 $A=r/2\\Leftrightarrow\gamma_{\rm loc}=2\ln2/r$ —— C 的等价式成立。
      于是 e137b 的 $\gamma_{\rm loc}(1.00)=1.860/1.420/1.281$（每秩 6 个设计的**中位曲线**）
      与 e159 的**同秩下包络**必须在 $\Delta I=1.00$ 处可比：取 $\Delta I\in[0.85,1.15]$ 的包络点重拟斜率报 $\gamma$。
      两者相对差 $\le5\%$ ⇒ "包络 vs 中位曲线"不影响读数；否则须请 C 指明它的律是关于哪一条曲线的。
 [L6] 自我降级：若 $A$ 在整个 $[\,x_{\min},x_{\max}\,]$ 上都 $\approx r/2$，那"梯子律"就是高斯率失真函数
      在近地板端的一阶展开，**不含任何本站的动力学内容** $\\Rightarrow$ 不得作为"我们的发现"写进正文。
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
ln2 = np.log(2.0)


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


def point(Z, A, W, TH, JC, s):
    r = Z.shape[1]
    Ir = np.eye(r)
    C = np.sqrt(s) * Z.T
    Pm = sym(solve_discrete_are(A.T, C.T, W, Ir))
    Sm = C @ Pm @ C.T + Ir
    Lk = Pm @ C.T @ np.linalg.inv(Sm)
    Pp = sym(Pm - Lk @ Sm @ Lk.T)
    I = 0.5 * np.log(np.linalg.det(Ir + C @ Pm @ C.T)) / ln2
    return float(I), JC + float(np.trace(TH @ Pp))


def cands(TH, r, nrand=40, seed=7):
    wth, Uth = np.linalg.eigh(TH)
    U = Uth[:, np.argsort(-wth)]
    cs = [U[:, :r]]
    rng = np.random.default_rng(seed + r)
    for _ in range(nrand):
        Q, Rr = np.linalg.qr(rng.normal(0, 1, (U.shape[0], r)))
        cs.append(Q * np.sign(np.diag(Rr)))
    return cs


def envelope(A, W, TH, JC, r):
    Zs = cands(TH, r)
    df = np.inf
    for Z in Zs:
        try:
            _, Dv = point(Z, A, W, TH, JC, 1e7)
        except np.linalg.LinAlgError:
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
    ds = np.array([I - R_exp_of(A) for _, I in env])
    m = xs > 1e-9
    return xs[m], ds[m], df


def fit(u, v):
    if u.size < 3:
        return np.nan
    return float(np.polyfit(u, v, 1)[0])


ROWS = []
for nm, spec in [('anchor', ('anchor', 0)), ('rand-1', ('rand', 5)), ('rand-2', ('rand', 11))]:
    A, B, W, TH, JC, n = make_plant(*spec)
    p('')
    p('--- %s  $J_C=%.4f$  $R_{\\exp}=%.6f$' % (nm, JC, R_exp_of(A)))
    for r in range(1, n):
        xs, ds, df = envelope(A, W, TH, JC, r)
        if xs.size < 8:
            p('  r=%d 有效包络点太少（略）' % r)
            continue
        lo = np.log10(xs.min())
        hi = np.log10(xs.max())
        pred = r / 2.0
        p('  [r=%d] $x_{\\min}=%.3e$，$x$ 跨度 %.1f 个十倍程，包络点 %d' % (r, xs.min(), hi - lo, xs.size))
        devs, marks, dom = [], [], []
        for j in range(0, 13):
            a, b = lo + j, lo + j + 1
            if a >= hi:
                break
            k = (np.log10(xs) >= a) & (np.log10(xs) < b)
            if k.sum() < 3:
                marks.append('   .  '); devs.append(np.nan); dom.append(None)
                continue
            A_ = fit(np.log2(1.0 / xs[k]), ds[k])
            dv = (A_ - pred) / pred * 100.0
            devs.append(dv)
            marks.append('%+5.0f%%' % dv)
            dom.append((xs[k].max(), ds[k].min(), ds[k].max(), int(k.sum())))
            p('     窗口 $x\\in[%.1e,%.1e]$（$x_{\\min}\\cdot10^{%d\\sim%d}$） 点 %2d  $A=%.4f$  dev %+.1f%%  $\\gamma=\\ln2/A=%.3f$  $\\Delta I\\in[%.2f,%.2f]$'
              % (xs[k].min(), xs[k].max(), j, j + 1, k.sum(), A_, dv, ln2 / A_ if np.isfinite(A_) and A_ > 0 else np.nan,
                 ds[k].min(), ds[k].max()))
        j_ok = 0
        while j_ok < len(devs) and np.isfinite(devs[j_ok]) and abs(devs[j_ok]) <= 10:
            j_ok += 1
        if j_ok == 0:
            p('     [L4] 最内窗口就不满足 $\\le10\\%$ => 本轮该梯级不报有效域')
        else:
            xd, dinv, dmax, _ = dom[j_ok - 1]
            p('     [L4] $|\\text{dev}|\\le10%%$ 连续成立到第 %d 个十倍程 $\\Rightarrow$ 有效域 $x\\le %.2e$（$=x_{\\min}\\times10^{%.1f}$），对应 $\\Delta I\\ge%.2f$ bit'
              % (j_ok, xd, np.log10(xd / xs.min()), dinv))
        # [L5] 与 e137b 的 gamma_loc(1.00) 对账（只有 anchor 有那张表）
        kb = (ds >= 0.85) & (ds <= 1.15)
        g = np.nan
        if kb.sum() >= 3:
            A_ = fit(np.log2(1.0 / xs[kb]), ds[kb])
            g = ln2 / A_
        ref = {1: 1.860, 2: 1.420, 3: 1.281}.get(r) if nm == 'anchor' else None
        p('     [L5] $\\Delta I\\in[0.85,1.15]$ 包络点 %d 个 $\\Rightarrow$ $\\gamma_{\\rm loc}$(包络)=%s vs e137b $\gamma_{\\rm loc}$(中位曲线)=%s'
          % (int(kb.sum()), ('%.3f' % g) if np.isfinite(g) else 'n/a', ('%.3f' % ref) if ref else '-'))
        if ref and np.isfinite(g):
            p('          相对差 %+.1f%% $\\Rightarrow$ %s' % (100 * (g / ref - 1), '口径无关' if abs(g / ref - 1) <= 0.05 else '两条曲线给出的局域斜率不同，须指明律是关于哪一条的'))
        ROWS.append((nm, r, devs, j_ok))

p('')
p('== [L4/L5/L6] 汇总判决（判据见文件头）==')
p('  有效域（十倍程数，$|\\text{dev}|\\le10\\%$ 连续）：' + '  '.join(
    '%s r=%d:%d' % (nm, r, j) for nm, r, _, j in ROWS))
allc = [abs(d) for _, _, dv, _ in ROWS for d in dv if np.isfinite(d)]
p('  全部窗口 $|\\text{dev}|$ 的分位：中位 %.1f%%，90%% 分位 %.1f%%，最大 %.1f%%（有效窗口 %d 个）'
  % (np.median(allc), np.percentile(allc, 90), max(allc), len(allc)))
n_in = sum(1 for _, _, dv, _ in ROWS for d in dv if np.isfinite(d) and abs(d) <= 10)
p('  其中 $\\le10%%$ 的窗口 %d/%d $\\Rightarrow$ 律的**主项系数正确、$O(x)$ 的写法只在近地板端成立**'
  % (n_in, len(allc)))
p('[L6] 检查：若 anchor 的 $A$ 在全部秩上都 $\\approx r/2$（本表可见），则该律等于高斯 RDF 在近地板端的一阶展开，')
p('     与 e150 [W1] "常数 $2\\ln2/r$ 由反向水填算术给出" 同结论 $\\Rightarrow$ 正文不得把 $r/2$ 算作本站发现。')

open('p0/e159b_out.txt', 'w', encoding='utf-8').write('\n'.join(out) + '\n')
