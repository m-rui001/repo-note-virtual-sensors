r"""E63b = 复刻原文 Fig.1 的红线：原文在锥 K_F 内自由优化信息矩阵，而 C 的 §11.2 只走了一条射线。

E63 量的事实：红线左端 (D=46.99, I=4.924)、右端 (D=89.81, I=1.546)。
按 C 的口径（$V=vI_r$，即 $G=\tau F^\top F$ 这一条**射线**）在 I=1.546 处给出的成本是几百，
完全贴不上红线。原因写在原文 §IV 那句话里："restricts the admissible information matrices to the
cone $K_F$"——原文优化的是
    G = F^T V^{-1} F = Z M Z^T,  M 为该平面上的任意正定阵（r=2 时 3 维，模尺度 2 维形状），
而 C 把 M 钉成 tau*I。地板 Phi_0 只看支撑平面（对 M 盲），率侧完全看 M。
所以：
 [A] 先量射线口径与锥口径在同一条 F_2 上差多少（这就是 §11.2 整张表的系统性偏差方向）；
 [B] 锥口径下按 I 定点最小化 J，看能不能落到数字化红线上（判定我们是否真的复刻了原文算例）；
 [C] 同一套代码跑 blue(r=3)/magenta(r=4) 两条，检查"右端点/左端点"是否都由 I 窗口决定；
 [D] 顺带把 C §11.4 那条"左端点=不可行墙"判掉：墙是 I->inf 的渐近线，图上的左端点是 I=5 的窗口边。
"""
import sys
sys.stdout.reconfigure(encoding='utf-8')
import numpy as np
from scipy.linalg import schur, solve_discrete_are as sdare
from scipy.optimize import minimize
from PIL import Image

np.set_printoptions(precision=4, suppress=True, linewidth=170)
_src = open('p0/exp_c_audit.py', encoding='utf-8').read().split("print(r'== E53")[0]
_ns = {'__name__': 'e53_preamble'}
exec(compile(_src, 'p0/exp_c_audit.py[preamble]', 'exec'), _ns)
A, W, TH, JC, n, sym = _ns['A'], _ns['W'], _ns['TH'], _ns['JC'], _ns['n'], _ns['sym']
phi_iter = _ns['phi_iter']
print('== E63b：锥 K_F 内的率-成本前沿 vs C 的射线口径 ==')

Ts, U = schur(A, output='real', sort=lambda a: abs(a) < 1.0)[:2]
Zs = {r: U[:, n - r:] for r in (2, 3, 4)}          # 原文 span(F) 的正交基（列）


def IJ(F, V):
    """一般 (F,V)：返回 (I, J)。F 为 (r,n)，V 为噪声协方差 (r,r)。"""
    F = np.atleast_2d(F); r = F.shape[0]
    try:
        Pt = sym(sdare(A.T, F.T, sym(W), sym(V)))
    except Exception:
        return None
    if not np.all(np.isfinite(Pt)):
        return None
    FPt = F @ Pt
    P = sym(Pt - FPt.T @ np.linalg.solve(FPt @ F.T + V, FPt))
    sign, ld = np.linalg.slogdet(np.eye(r) + np.linalg.solve(V.T, F @ Pt @ F.T).T)
    if sign <= 0:
        return None
    return 0.5 * ld / np.log(2), JC + float(np.trace(TH @ P))


def J_at_I_ray(Z, I_target):
    """C 的口径：M = tau*I（射线），二分 tau 命中 I_target。"""
    F = Z.T; lo, hi = -6.0, 14.0
    for _ in range(60):
        m = 0.5 * (lo + hi)
        a = IJ(F, np.eye(F.shape[0]) / 10.0 ** m)
        if a is None or a[0] < I_target:
            lo = m
        else:
            hi = m
    return IJ(F, np.eye(F.shape[0]) / 10.0 ** hi)


def J_at_I_cone(Z, I_target, starts=3):
    """原文口径：M = s*(任意正定) —— 用 Cholesky 参数化 r(r+1)/2 维，尺度由 I 的二分吸收。"""
    r = Z.shape[1]
    nz = r * (r + 1) // 2

    def chol_M(x):
        L = np.zeros((r, r)); k = 0
        for i in range(r):
            for j in range(i + 1):
                L[i, j] = np.exp(x[k]) if i == j else x[k]
                k += 1
        return L @ L.T

    def obj(x):
        M = chol_M(x)
        if np.linalg.cond(M) > 1e12:
            return 1e6
        lo, hi = -8.0, 16.0
        for _ in range(34):
            m = 0.5 * (lo + hi)
            a = IJ(Z.T, np.linalg.inv(10.0 ** m * M))
            if a is None or a[0] < I_target:
                lo = m
            else:
                hi = m
        a = IJ(Z.T, np.linalg.inv(10.0 ** hi * M))
        return 1e6 if a is None else a[1]

    best = (np.inf, None)
    rng = np.random.default_rng(7)
    x_iso = np.zeros(nz)
    k = 0
    for i in range(r):
        for j in range(i + 1):
            x_iso[k] = 0.0 if i != j else -2.0
            k += 1
    for x0 in [x_iso] + [rng.normal(0, 1.0, nz) for _ in range(starts - 1)]:
        res = minimize(obj, x0, method='Nelder-Mead',
                       options=dict(xatol=1e-4, fatol=1e-7, maxiter=60 * nz))
        if res.fun < best[0]:
            best = (float(res.fun), res.x)
    return best

# ---------------------------------------------------------------- 数字化红线（精确掩码）
Aimg = np.asarray(Image.open('.work3/fig_p5_img0.png').convert('RGB')).astype(int)
H, Wd, _ = Aimg.shape
x0, x1, y0, y1 = 121.0, 1924.0, 52.0, 1239.0
px2D = lambda px: 40.0 + (px - x0) * 50.0 / (x1 - x0)
px2I = lambda py: 5.0 - (py - y0) * 4.0 / (y1 - y0)
leg = np.zeros(Aimg.shape[:2], bool)
leg[int(y0):int(y0 + 0.34 * (y1 - y0)) + 60, int(x0 + 0.40 * (x1 - x0)) + 110:] = True
box = np.zeros_like(leg); box[int(y0) + 2:int(y1) - 1, int(x0) + 3:int(x1) - 2] = True
COL = {'red': (255, 0, 0), 'blue': (0, 0, 255), 'magenta': (255, 0, 255)}
DIG = {}
print('\n[A] 三条曲线的数字化（精确同色掩码，列中位数）')
for nm, c in COL.items():
    m = np.all(Aimg == np.array(c), axis=2) & box & (~leg)
    xs = np.where(m.any(0))[0]
    pts = []
    for px in xs:
        py = np.where(m[:, px])[0]
        pts.append((px2D(px), px2I(py.mean())))
    pts = np.array(pts)
    DIG[nm] = pts
    o = np.argsort(pts[:, 0])
    print('  %-8s %d 列: D in [%.2f, %.2f]; 左端 (D=%.2f, I=%.3f), 右端 (D=%.2f, I=%.3f)'
          % (nm, len(xs), pts[:, 0].min(), pts[:, 0].max(), pts[o[0], 0], pts[o[0], 1],
             pts[o[-1], 0], pts[o[-1], 1]))
    for It in (1.55, 2.0, 2.5, 3.0, 4.0, 4.9):
        sel = np.abs(pts[:, 1] - It) < 0.03
        print('        I≈%.2f 处 D = %s' % (It, np.round(np.median(pts[sel, 0]), 2) if sel.sum() else '---'))

# ---------------------------------------------------------------- [B] 射线 vs 锥
print('\n[B] 同一条 F_2（原文 Schur 尾 2 列）：射线口径 tau*I vs 锥口径最优')
print('  %-8s %-12s %-12s %-12s' % ('目标 I', '射线 J', '锥 J', '差(射线-锥)'))
for It in (1.546, 2.0, 3.0, 4.0, 4.924, 6.0):
    jr = J_at_I_ray(Zs[2], It)[1]
    jc, x = J_at_I_cone(Zs[2], It)
    L = np.zeros((2, 2)); k = 0
    for i in range(2):
        for j in range(i + 1):
            L[i, j] = np.exp(x[k]) if i == j else x[k]; k += 1
    ev = np.linalg.eigvalsh(L @ L.T)
    print('  %-8.3f %-12.4f %-12.4f %-12.4f   最优形状 M 的特征值比 = %.3f : 1'
          % (It, jr, jc, jr - jc, ev.max() / ev.min()))

print('\n[C] 与数字化红线的对照（锥口径应贴住原文曲线）')
red = DIG['red']
for It in (1.55, 2.0, 2.5, 3.0, 4.0, 4.9):
    sel = np.abs(red[:, 1] - It) < 0.04
    if not sel.sum():
        continue
    Dd = float(np.median(red[sel, 0]))
    jc, _ = J_at_I_cone(Zs[2], It)
    jr = J_at_I_ray(Zs[2], It)[1]
    print('  I=%.2f: 图 D=%.3f | 锥 %.3f (%+.3f) | 射线 %.3f (%+.3f)' % (It, Dd, jc, jc - Dd, jr, jr - Dd))

print('\n[D] blue(r=3)/magenta(r=4)：各自窗口两端')
for nm, r in (('blue', 3), ('magenta', 4)):
    for It in (1.55, 4.9):
        try:
            jr = J_at_I_ray(Zs[r], It)[1]
        except Exception:
            jr = np.nan
        sel = np.abs(DIG[nm][:, 1] - It) < 0.05
        Dd = float(np.median(DIG[nm][sel, 0])) if sel.sum() else np.nan
        print('  %-8s I=%.2f: 图 D=%9.3f | 射线口径 J=%9.3f | 地板=%9.4f'
              % (nm, It, Dd, jr, JC + phi_iter(Zs[r].T)[0]))

print('\n[E] 结论用的三个数')
for r in (2, 3, 4):
    print('  r=%d: 墙(地板)=%.4f | 射线口径 I=4.924 的 D=%.4f | 锥口径同点 D=%.4f'
          % (r, JC + phi_iter(Zs[r].T)[0], J_at_I_ray(Zs[r], 4.924)[1], J_at_I_cone(Zs[r], 4.924)[0]))
print('== E63b done ==')
