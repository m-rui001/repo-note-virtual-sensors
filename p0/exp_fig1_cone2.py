r"""E63c = 把 Fig.1 三条 TRV 曲线全部复刻（锥口径），并给出"曲线从哪条窗口边出去"的预测。

E63b 已经把红线钉住：锥口径 D(I) 与数字化红线的差 |ΔD|<=0.14（I>=2），而 C 的射线口径差 +9~+39。
这里把同一套代码用到 blue(r=3)、magenta(r=4=全状态，原文说它=unrestricted) 上：
 - 若在 I=2.0/2.5/3.0 三格都贴上数字化值，则"我们复刻了原文算例的三条曲线"成立，
   整块板子的植物读数（A,B,W,R,Q=I）就此闭环；
 - 再预测每条曲线的出窗边：D(4.926) 若 >50 则该曲线从顶边出，否则从左边缘 D=40 出。
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
print('== E63c：Fig.1 三条曲线的锥口径复刻 ==')

Ts, U = schur(A, output='real', sort=lambda a: abs(a) < 1.0)[:2]


def IJ(F, V):
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


def cone_J(Z, I_target, starts=5):
    """给定率，最小化成本：M 为该平面上任意正定（Cholesky 参数化），尺度由内层二分命中 I_target。"""
    r = Z.shape[1]
    nz = r * (r + 1) // 2

    def build(x):
        L = np.zeros((r, r)); k = 0
        for i in range(r):
            for j in range(i + 1):
                L[i, j] = np.exp(x[k]) if i == j else x[k]; k += 1
        return L @ L.T

    def obj(x):
        M = build(x)
        try:
            if np.linalg.cond(M) > 1e10:
                return 1e6
        except Exception:
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

    x_iso = np.zeros(nz)
    k = 0
    for i in range(r):
        for j in range(i + 1):
            x_iso[k] = -2.0 if i == j else 0.0
            k += 1
    rng = np.random.default_rng(11)
    best = (np.inf, None)
    for x0 in [x_iso] + [rng.normal(0, 0.8, nz) for _ in range(starts - 1)]:
        res = minimize(obj, x0, method='Powell',
                       options=dict(xtol=1e-4, ftol=1e-7, maxiter=120 * nz, maxfev=120 * nz))
        if res.fun < best[0]:
            best = (float(res.fun), res.x)
    return best


# 数字化三条曲线
import os
_IMG = next(p for p in ('p0/fig/lcss_p5_x327.png', '.work3/fig_p5_img0.png')
               if os.path.exists(p))    # the repro package ships the raster under p0/fig/
Aimg = np.asarray(Image.open(_IMG).convert('RGB')).astype(int)
x0, x1, y0, y1 = 121.0, 1924.0, 52.0, 1238.0   # frame bottom as auto-detected by exp_fig1_digitize.py
px2D = lambda px: 40.0 + (px - x0) * 50.0 / (x1 - x0)
px2I = lambda py: 5.0 - (py - y0) * 4.0 / (y1 - y0)
box = np.zeros(Aimg.shape[:2], bool); box[int(y0) + 2:int(y1) - 1, int(x0) + 3:int(x1) - 2] = True
leg = np.zeros_like(box); leg[int(y0):int(y0 + 0.34 * (y1 - y0)) + 60, int(x0 + 0.40 * (x1 - x0)) + 110:] = True
CURVE = {'red': ((255, 0, 0), 2), 'blue': ((0, 0, 255), 3), 'magenta': ((255, 0, 255), 4)}
DIG = {}
for nm, (c, r) in CURVE.items():
    m = np.all(Aimg == np.array(c), axis=2) & box & (~leg)
    pts = [(px2D(px), px2I(np.where(m[:, px])[0].mean())) for px in np.where(m.any(0))[0]]
    DIG[nm] = np.array(pts)
    print('  %-8s r=%d: D in [%.2f, %.2f], I in [%.3f, %.3f]'
          % (nm, r, DIG[nm][:, 0].min(), DIG[nm][:, 0].max(), DIG[nm][:, 1].min(), DIG[nm][:, 1].max()))

print('\n[A] 锥口径复刻（目标率 -> 最小成本），对照数字化曲线')
ISOP = {}
for nm, (c, r) in CURVE.items():
    Z = U[:, n - r:]
    fl = JC + phi_iter(Z.T)[0]
    print('  %s (r=%d, 地板 %.4f):' % (nm, r, fl))
    for It in (2.0, 2.5, 3.0):
        sel = np.abs(DIG[nm][:, 1] - It) < 0.04
        if not sel.sum():
            print('    I=%.2f: 图上没有这一段' % It); continue
        Dd = float(np.median(DIG[nm][sel, 0]))
        jd, x = cone_J(Z, It)
        # 射线口径（C 的 §11 表所用的各向同性 V=vI）
        lo, hi = -8.0, 16.0
        for _ in range(34):
            m2 = 0.5 * (lo + hi)
            a = IJ(Z.T, np.eye(r) / 10.0 ** m2)
            if a is None or a[0] < It:
                lo = m2
            else:
                hi = m2
        jray = IJ(Z.T, np.eye(r) / 10.0 ** hi)[1]
        ISOP[(nm, It)] = (Dd, jd, jray)
        print('    I=%.2f: 图 D=%8.3f | 锥 %8.3f (%+6.3f, 相对 %+.2f%%) | 射线 %8.3f (%+6.3f)'
              % (It, Dd, jd, jd - Dd, 100 * (jd - Dd) / Dd, jray, jray - Dd))

print('\n[B] 出窗边的预测：图框是 D in [40,90], I in [1,5]；每条曲线在哪个边被截断')
for nm, (c, r) in CURVE.items():
    Z = U[:, n - r:]
    Dtop, _ = cone_J(Z, 4.926)
    print('  %-8s 锥口径 D(I=4.926)=%9.3f；数字化左端 D=%.2f, I=%.3f => %s'
          % (nm, Dtop, DIG[nm][np.argmin(DIG[nm][:, 0]), 0], DIG[nm][np.argmin(DIG[nm][:, 0]), 1],
             '从顶边出窗' if Dtop > 40.5 else '从左边缘 D=40 出窗'))

print('\n[C] 反算数字化曲线的率：给定图上 (D,I) 点，锥口径应给出 J<=D 且 I>=图上 I（不可行区检查）')
for nm, (c, r) in CURVE.items():
    Dd = float(DIG[nm][np.argmin(DIG[nm][:, 0]), 0])
    Id = float(DIG[nm][np.argmin(DIG[nm][:, 0]), 1])
    jn, _ = cone_J(U[:, n - r:], Id)
    print('  %-8s 左端 (D,I)=(%.2f,%.3f): 锥 J(I)=%.3f  差 %+.3f | 该 F 的墙 %.4f'
          % (nm, Dd, Id, jn, jn - Dd, JC + phi_iter(U[:, n - r:].T)[0]))
print('== E63c done ==')
