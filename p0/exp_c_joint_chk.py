r"""E62 = 审 C 的 §11（联合设计：同一成本 J 下所需信息率 I 的对照表）。

他这张表要证明的是"$r$ 固定时该建哪几条有两个坐标：地板看 $\Phi_0$、率侧另算"，并且给原文 $F_2$
封了一个"率侧最优"的名号。我这边三个怀疑，逐个钉：

 (1) **评价口径**：他用 $P_m=(\tilde P^{-1}+\tau F^\top F)^{-1}$ 的信息矩阵形式迭代。
     我改成协方差/Kalman 增益形式（$V=I/\tau$，$K=\tilde PF^\top(F\tilde PF^\top+V)^{-1}$），
     全程不求 $\tilde P$ 的逆，逐位对照他的 $I,J$。
 (2) **设计空间维度**：$(I,J)$ 曲线只依赖 $G=F^\top F$（秩 $r$ 半正定，$rn-r(r-1)/2$ 维），
     不是只依赖 $\operatorname{span}(F)$。$\tau\to\infty$ 的地板确实只看支撑平面（尺度/形状盲），
     但**率**还吃行的相对增益与剪切。他把 7 维（$r{=}2$）问题画成 4 维平面问题，再拿 4 个候选比较。
     ⇒ 先验证地板的规范不变性，再验证改一行增益能在**同一下界**上白拿比特。
 (3) **"率侧最优"没有依据**：$F_2$ 从未在率侧被优化过。我直接在给定 $J$ 上对 $F$ 做无导数优化
     （240 个随机形状 + Powell），看 $J=80$ 那一列他印的 $1.986$ 是不是最优。

跑法：`python -X utf8 p0/exp_c_joint_chk.py > p0/e62_out.txt 2>&1`
"""
import sys
sys.stdout.reconfigure(encoding='utf-8')
import numpy as np
from scipy.linalg import schur
from scipy.optimize import minimize

np.set_printoptions(precision=5, suppress=True, linewidth=170)
_src = open('p0/exp_c_audit.py', encoding='utf-8').read().split("print(r'== E53")[0]
_ns = {'__name__': 'e53_preamble'}
exec(compile(_src, 'p0/exp_c_audit.py[preamble]', 'exec'), _ns)
A, W, TH, JC, n, sym = _ns['A'], _ns['W'], _ns['TH'], _ns['JC'], _ns['n'], _ns['sym']
phi_iter, phi_dare = _ns['phi_iter'], _ns['phi_dare']

print('== E62：C §11 联合前沿的独立复算 + 形状(增益分配)这一维 ==')
print('jc = %.4f   eig(A) = %s' % (JC, np.round(np.linalg.eigvals(A), 4)))

# ------------------------------------------------------------------ 候选 F（(r,n) 口径）
wT, VT = np.linalg.eigh(TH)
Vt = VT[:, ::-1].T
T2, U2 = schur(A, output='real', sort=lambda a: abs(a) < 1.0)[:2]
Sun = U2[:, 2:].T          # C 的"不稳 Schur 尾"：只有 2 行
best = {r: np.atleast_2d(np.load('.work3/c24_bestF_r%d.npy' % r)) for r in (1, 2, 3)}  # 文件已是 (r,n)

CANDS = {}
for r in (1, 2, 3):
    CANDS[(r, 'gdesc')] = best[r]
    CANDS[(r, 'Theta_top')] = Vt[:r]
    CANDS[(r, 'std_basis')] = np.eye(n)[:r]
CANDS[(2, 'Schur_tail(=原文F2)')] = Sun
CANDS[(3, 'Schur_tail[:3]')] = Sun[:3]     # 他 r=3 表里那一行就是这个切片
CANDS[(1, 'Schur_tail[0]')] = Sun[:1]

print('\n[A] 候选的行范数（尺度是率侧的自由度，行范数不齐=增益分配不同）')
for (r, nm), F in CANDS.items():
    print('  r=%d %-22s shape=%s 行范数=%s  ||F||_F=%.4f'
          % (r, nm, F.shape, np.round(np.linalg.norm(F, axis=1), 4), np.linalg.norm(F)))

# ------------------------------------------------------------------ 协方差形式的滤波定点
from scipy.linalg import solve_discrete_are


def filt_iter(F, tau, it=4000, tol=1e-12, cap=1e13):
    """纯迭代（不用 DARE 解器），只作为 filt 的备用/对照。"""
    F = np.atleast_2d(F); r = F.shape[0]
    V = np.eye(r) / tau
    Pt = W.copy(); P = Pt
    for k in range(it):
        FPt = F @ Pt
        S = F @ FPt.T + V
        K = np.linalg.solve(S.T, FPt).T
        P = sym(Pt - K @ FPt)
        Pn = sym(A @ P @ A.T + W)
        d = np.max(np.abs(Pn - Pt))
        Pt = Pn
        if d < tol * max(1.0, np.max(np.abs(Pt))):
            break
        if np.max(np.abs(Pt)) > cap:
            return None
    FPt = F @ Pt
    S = F @ FPt.T + V
    K = np.linalg.solve(S.T, FPt).T
    P = sym(Pt - K @ FPt)
    sign, ld = np.linalg.slogdet(np.eye(r) + tau * (F @ Pt @ F.T))
    if sign <= 0:
        return None
    return 0.5 * ld / np.log(2), JC + float(np.trace(TH @ P))


def filt(F, tau):
    """(I,J) @ 精度 tau（V=I/tau）。路线：滤波 Riccati = DARE(A^T, F^T, W, V) 的稳定解，
    与 C 的信息矩阵迭代、与我此前的定点/奇异 DARE 两条都不同；发散(不可检测)返回 None。"""
    F = np.atleast_2d(F); r = F.shape[0]
    V = np.eye(r) / tau
    try:
        Pt = sym(solve_discrete_are(A.T, F.T, sym(W), V))
    except Exception:
        return None
    if not np.all(np.isfinite(Pt)) or np.max(np.abs(Pt)) > 1e13:
        return None
    FPt = F @ Pt
    P = sym(Pt - FPt.T @ np.linalg.solve(FPt @ F.T + V, FPt))
    sign, ld = np.linalg.slogdet(np.eye(r) + tau * (F @ Pt @ F.T))
    if sign <= 0:
        return None
    return 0.5 * ld / np.log(2), JC + float(np.trace(TH @ P))


def rate_at(F, Jt, lo=-2.0, hi=13.0, steps=24):
    """给定总成本目标 Jt，解出所需 I（J 对 log tau 单调下降）。不可行返回 inf。"""
    a = filt(F, 10.0 ** hi)
    if a is None or a[1] > Jt:
        return np.inf
    b = filt(F, 10.0 ** lo)
    if b is not None and b[1] < Jt:
        return b[0]
    for _ in range(steps):
        m = 0.5 * (lo + hi)
        c = filt(F, 10.0 ** m)
        if c is None or c[1] > Jt:
            lo = m
        else:
            hi = m
    return filt(F, 10.0 ** hi)[0]

# ------------------------------------------------------------------ [B] 他的表逐项复算
print('\n[B] 地板：三条路线（我的定点 phi_iter / 奇异 DARE / 协方差滤波 tau=1e13）')
print('  %-26s %-9s %-9s %-9s %-9s' % ('候选', 'phi_iter', 'phi_dare', 'filt 1e13', '极差'))
floors = {}
for (r, nm), F in CANDS.items():
    p1 = phi_iter(np.atleast_2d(F))[0]
    p2 = phi_dare(np.atleast_2d(F))[0]
    a = filt(np.atleast_2d(F), 1e13)
    v3 = a[1] - JC if a else np.nan
    row = [p1, p2, v3]
    floors[(r, nm)] = p1
    print('  r=%d %-20s %-9.4f %-9.4f %-9.4f %-9.2e'
          % (r, nm, p1, p2, v3, max(row) - min(row)))

print('\n[B2] C 的 §11.2 表逐项复算：他印的值 / 我的协方差-DARE 路线 / 差')
HIS = {(2, 'gdesc'): {40: 3.920, 45: 3.368, 46.1231: 3.280, 50: 3.033, 60: 2.627, 80: 2.211},
       (2, 'Theta_top'): {40: 3.851, 45: 3.306, 46.1231: 3.219, 50: 2.975, 60: 2.577, 80: 2.169},
       (2, 'Schur_tail(=原文F2)'): {40: np.nan, 45: np.nan, 46.1231: np.nan, 50: 4.073, 60: 2.695, 80: 1.986},
       (2, 'std_basis'): {40: np.nan, 45: np.nan, 46.1231: np.nan, 50: np.nan, 60: 5.056, 80: 3.707},
       (3, 'gdesc'): {40: 4.650, 45: 3.979, 46.1231: 3.870, 50: 3.560, 60: 3.038, 80: 2.493},
       (3, 'Theta_top'): {40: 4.626, 45: 3.970, 46.1231: 3.863, 50: 3.560, 60: 3.049, 80: 2.515},
       (3, 'Schur_tail[:3]'): {40: np.nan, 45: np.nan, 46.1231: np.nan, 50: 4.073, 60: 2.695, 80: 1.986},
       (1, 'gdesc'): {40: np.nan, 45: 3.010, 46.1231: 2.817, 50: 2.417, 60: 1.971, 80: 1.644},
       (1, 'Theta_top'): {40: np.nan, 45: 3.664, 46.1231: 3.228, 50: 2.587, 60: 2.023, 80: 1.653},
       (1, 'Schur_tail[0]'): {40: np.nan, 45: np.nan, 46.1231: np.nan, 50: np.nan, 60: np.nan, 80: np.nan}}
HFL = {(2, 'gdesc'): 32.2185, (2, 'Theta_top'): 32.2751, (2, 'Schur_tail(=原文F2)'): 46.1235,
       (2, 'std_basis'): 50.1115, (3, 'gdesc'): 31.6420, (3, 'Theta_top'): 31.6736,
       (3, 'Schur_tail[:3]'): 46.1235, (1, 'gdesc'): 41.6527, (1, 'Theta_top'): 43.7057,
       (1, 'Schur_tail[0]'): 396.66}
TGT = [40.0, 45.0, 46.1231, 50.0, 60.0, 80.0]
worst = 0.0
for (r, nm) in HIS:
    F = np.atleast_2d(CANDS[(r, nm)])
    fl = JC + floors[(r, nm)]
    cells = []
    for t in TGT:
        v = rate_at(F, t); hv = HIS[(r, nm)][t]
        if np.isfinite(v) and np.isfinite(hv):
            worst = max(worst, abs(v - hv))
        cells.append('%6s/%-6s' % ('  --  ' if not np.isfinite(hv) else '%6.3f' % hv,
                                   'inf' if not np.isfinite(v) else '%6.3f' % v))
    print('  r=%d %-22s 地板 %8.4f/%-8.4f (他印 %.4f) |%s'
          % (r, nm, fl, fl, HFL[(r, nm)], ''.join('  ' + c for c in cells)))
print('  列口径：他值/我值，目标成本 J = %s' % TGT)
print('  最大符合偏差 = %.4f bit' % worst)

# ------------------------------------------------------------------ [C] 交叉点
print('\n[C] r=2：原文 F_2 vs Theta 前 2 的率曲线交叉位置（他报 J in (60,80)）')
F2, FT = np.atleast_2d(Sun), np.atleast_2d(Vt[:2])
for t in [46.2, 50, 55, 60, 65, 70, 75, 80, 90, 100, 120, 160]:
    a, b = rate_at(F2, t), rate_at(FT, t)
    print('  J=%-6.1f  F2 %7.3f   Theta前2 %7.3f   差 %+7.3f  %s'
          % (t, a, b, a - b, 'F2 省' if a < b else 'Theta 省'))
lo, hi = 60.0, 80.0
for _ in range(40):
    m = 0.5 * (lo + hi)
    if rate_at(F2, m) < rate_at(FT, m):
        hi = m
    else:
        lo = m
print('  二分交叉点 J* = %.4f（此处两者所需率相等，I*=%.4f bit）' % (m, rate_at(F2, m)))

# ------------------------------------------------------------------ [D] 给定 J 上真正优化 F
print('\n[D] 固定 J 上最小化 I：设计空间是 F^T F（r=2 时 7 维，模掉整体尺度还剩 6），不只是平面(4)')
print('    他拿 4 个候选就写"浅成本区取不稳 Schur 尾=率侧最优"——这里直接优化。')


def obj(x, Jt, r):
    F = x.reshape(r, n)
    if np.linalg.matrix_rank(F) < r:
        return 1e3
    v = rate_at(F, Jt)
    return 1e3 if not np.isfinite(v) else v


rng = np.random.default_rng(20260929)
import time
for r, Jt in [(2, 80.0), (2, 45.0), (2, 40.0), (1, 80.0)]:
    t0 = time.time()
    cands = [np.atleast_2d(CANDS[k]) for k in CANDS if k[0] == r]
    cands += [rng.standard_normal((r, n)) for _ in range(80)]
    vals = [(rate_at(F, Jt), i) for i, F in enumerate(cands)]
    vals = [(v, i) for v, i in vals if np.isfinite(v)]
    vals.sort()
    best_v, best_x = np.inf, None
    for v0, i in vals[:3]:
        res = minimize(obj, cands[i].ravel(), args=(Jt, r), method='Powell',
                       options=dict(xtol=1e-5, ftol=1e-8, maxfev=2500))
        if res.fun < best_v:
            best_v, best_x = float(res.fun), res.x.reshape(r, n)
    fl_opt = JC + phi_iter(np.atleast_2d(best_x))[0]
    print('  r=%d J=%-6.1f 随机/候选里最好 %7.3f -> Powell 后最优率 %7.3f  (他表 F2 %.3f / Theta %.3f / gdesc %.3f)'
          % (r, Jt, vals[0][0] if vals else np.nan, best_v,
             rate_at(F2, Jt) if r == 2 else np.nan,
             rate_at(FT, Jt) if r == 2 else np.nan,
             rate_at(np.atleast_2d(best[r]), Jt)))
    print('        优化解 地板 J=%.4f（率 %.3f 可行说明它不在地板上） 行范数=%s  ||F||_F=%.4f  [%.0fs]'
          % (fl_opt, rate_at(best_x, fl_opt + 1e-4), np.round(np.linalg.norm(best_x, axis=1), 4),
             np.linalg.norm(best_x), time.time() - t0))
    np.save('p0/e62_best_rate_F_r%d_J%g.npy' % (r, Jt), best_x)

# ------------------------------------------------------------------ [E] 只改增益能白拿多少
print('\n[E] 形状维的价值：固定 span(F) 不变，只调行增益/剪切，地板不动而率下降')
print('    额外维数 = r(r+1)/2 - 1（正定左乘模尺度）：r=1 为 0（无形状自由），r=2 为 2，r=3 为 5')
for (r, nm) in [(2, 'Schur_tail(=原文F2)'), (2, 'Theta_top'), (2, 'gdesc'), (3, 'gdesc'), (3, 'Theta_top')]:
    F = np.atleast_2d(CANDS[(r, nm)])
    for Jt in [40.0, 80.0]:
        base = rate_at(F, Jt)
        if not np.isfinite(base):
            print('  r=%d %-22s J=%-6.1f 不可行(%.4f 地板)' % (r, nm, Jt, JC + floors[(r, nm)]))
            continue
        # 只允许对角缩放（r 维，模整体尺度 = r-1 个真自由度）
        res = minimize(lambda u: rate_at(np.diag(np.exp(u)) @ F, Jt), np.zeros(r),
                       method='Nelder-Mead', options=dict(xatol=1e-6, fatol=1e-9, maxiter=1200))
        # 任意左乘 L（对称，模尺度：非对角 r(r-1)/2 个 + 前 r-1 个对角增益，最后一个对角固定为 1）
        n_off, nz = r * (r - 1) // 2, r * (r + 1) // 2 - 1

        def full(z):
            L = np.eye(r)
            k = 0
            for i in range(r):
                for j in range(i + 1, r):
                    L[i, j] = L[j, i] = z[k]; k += 1
            for i in range(r - 1):
                L[i, i] = np.exp(z[n_off + i])
            return rate_at(L @ F, Jt)
        res2 = minimize(full, np.zeros(nz), method='Nelder-Mead',
                        options=dict(xatol=1e-6, fatol=1e-9, maxiter=2500))
        print('  r=%d %-22s J=%-6.1f 基线 %6.3f -> 仅对角缩放 %6.3f (%+0.3f) -> 任意左乘 %6.3f (%+0.3f) bit'
              % (r, nm, Jt, base, res.fun, res.fun - base, res2.fun, res2.fun - base))

print('\n[F] 规范不变性自检：地板对 F->sF、F->OF 应当严格不动；率曲线应是同一族的重新参数化')
F = np.atleast_2d(Sun)
th = 0.7
Q = np.array([[np.cos(th), -np.sin(th)], [np.sin(th), np.cos(th)]])
for nm, Fx in [('F', F), ('2F', 2 * F), ('QF', Q @ F), ('D F (D=diag(3,1))', np.diag([3., 1.]) @ F)]:
    p = phi_iter(Fx)[0]
    Js = [rate_at(Fx, t) for t in (45.0, 80.0)]
    print('  %-22s phi_iter=%.6f   I@45=%s  I@80=%s' % (nm, p, *['%7.3f' % x for x in Js]))
print('  => 若前三个 I 相同而第四个不同：率侧设计空间 = 平面 + 增益形状，正定左乘里只有非正交部分起作用')
print('== E62 done ==')
