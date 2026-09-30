r"""E63e = 锥的**闭包**口径：S=Z M Z^T，M=NN^T（N 全自由 -> 自动允许秩亏），G=S^{1/2}，V=I。

为什么要换口径。E63b/E63c 用 Cholesky M=LL^T（正定）+ cond(M)<=1e10 惩罚，结果是
"最优形状特征值比 6.0e7 / 9.6e8"——优化器在顶帽子。而原文的 K_F 定义为
{S 半正定 : range(S) subset range(F^T)}，**边界本来就属于可行集**。
边界上那条零精度方向 = 一条根本不用的通道（V=inf），所以"用了 r 支传感器"
不等于"锥里有 r 个活跃方向"。

数值上用 S=G^T G 的因式（G=S^{1/2} 半正定平方根，秩亏无害）走滤波器 DARE
sdare(A^T, G^T, W, I)，全程不出现 V=S^{-1}，条件数天然好。

自检（必须先过）：blue 的平面 span(U[:,1:]) 真包含 red 的 span(U[:,2:])，
所以同一目标率下必有 J_blue <= J_red。E63c 给的是 62.897 > 62.896，即被帽子截断。
"""
import sys, time
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
print('== E63e：锥闭包口径（秩亏允许）复刻 Fig.1 ==')

Ts, U = schur(A, output='real', sort=lambda a: abs(a) < 1.0)[:2]
In_ = np.eye(n)


def s_half(S):
    ev, EV = np.linalg.eigh(sym(S))
    ev = np.clip(ev, 0.0, None)
    return sym(EV @ np.diag(np.sqrt(ev)) @ EV.T)


def IJ_S(S):
    """S=信息增益矩阵（闭包内任意半正定）；返回 (率 bit/步, 总成本 J)。"""
    G = s_half(S)
    try:
        Pt = sym(sdare(A.T, G.T, sym(W), In_))
    except Exception:
        return None
    if not np.all(np.isfinite(Pt)):
        return None
    P = sym(Pt - Pt @ G @ np.linalg.solve(G @ Pt @ G + In_, G @ Pt))
    sign, ld = np.linalg.slogdet(In_ + sym(Pt @ S))
    if sign <= 0:
        return None
    return 0.5 * ld / np.log(2), JC + float(np.trace(TH @ P))


def closure_J(Z, I_target, starts=3, maxfev=400, bisect=24):
    """min J s.t. I>=I_target，M=NN^T 在 span(Z) 上（含秩亏）。尺度由二分命中。"""
    r = Z.shape[1]

    def obj(N):
        M = sym(N @ N.T)
        lo, hi = -8.0, 16.0
        for _ in range(bisect):
            m = 0.5 * (lo + hi)
            a = IJ_S(Z @ (10.0 ** m * M) @ Z.T)
            if a is None or a[0] < I_target:
                lo = m
            else:
                hi = m
        a = IJ_S(Z @ (10.0 ** hi * M) @ Z.T)
        return 1e6 if a is None else a[1]

    rng = np.random.default_rng(7)
    best = (np.inf, None)
    # 第一个起点 = 各向同性（等价 C 的射线口径），保证闭包值 <= 射线值
    for N0 in [np.eye(r)] + [rng.normal(0, 1, (r, r)) for _ in range(starts - 1)]:
        res = minimize(obj, N0, method='Powell',
                       options=dict(xtol=1e-4, ftol=1e-8, maxfev=maxfev))
        if res.fun < best[0]:
            best = (float(res.fun), res.x)
    N = np.asarray(best[1]).reshape(r, r)
    M = sym(N @ N.T)
    ev = np.clip(np.linalg.eigvalsh(M), 0, None)
    return best[0], ev[::-1] / max(ev.max(), 1e-300)


Aimg = np.asarray(Image.open('.work3/fig_p5_img0.png').convert('RGB')).astype(int)
x0, x1, y0, y1 = 121.0, 1924.0, 52.0, 1239.0
px2D = lambda px: 40.0 + (px - x0) * 50.0 / (x1 - x0)
px2I = lambda py: 5.0 - (py - y0) * 4.0 / (y1 - y0)
box = np.zeros(Aimg.shape[:2], bool)
box[int(y0) + 2:int(y1) - 1, int(x0) + 3:int(x1) - 2] = True
leg = np.zeros_like(box)
leg[int(y0):int(y0 + 0.34 * (y1 - y0)) + 60, int(x0 + 0.40 * (x1 - x0)) + 110:] = True
CURVE = {'red': ((255, 0, 0), 2), 'blue': ((0, 0, 255), 3), 'magenta': ((255, 0, 255), 4)}
DIG = {}
for nm, (c, r) in CURVE.items():
    m = np.all(Aimg == np.array(c), axis=2) & box & (~leg)
    DIG[nm] = np.array([(px2D(px), px2I(np.where(m[:, px])[0].mean()))
                        for px in np.where(m.any(0))[0]])

print('\n[A] 自检：同一天花板下的嵌套性。射线口径 = 起点 N=I 的那一列（Powell 第一步就该 <= 它）')
for nm, (c, r) in CURVE.items():
    Z = U[:, n - r:]
    lo, hi = -8.0, 16.0
    M = np.eye(r)
    for _ in range(24):
        m = 0.5 * (lo + hi)
        a = IJ_S(Z @ (10.0 ** m * M) @ Z.T)
        if a is None or a[0] < 2.0:
            lo = m
        else:
            hi = m
    ray = IJ_S(Z @ (10.0 ** hi * M) @ Z.T)[1]
    print('  %-8s r=%d 射线(闭包写法) I=2.00 -> D=%9.3f' % (nm, r, ray))

print('\n[B] 闭包最优（目标率 -> 最小成本），对照数字化曲线；特征值比=活跃方向数')
RES = {}
for nm, (c, r) in CURVE.items():
    Z = U[:, n - r:]
    fl = JC + phi_iter(Z.T)[0]
    print('  %s (r=%d, 地板 %.4f):' % (nm, r, fl))
    for It in (2.0, 2.5, 3.0, 4.0):
        sel = np.abs(DIG[nm][:, 1] - It) < 0.04
        if not sel.sum():
            print('    I=%.2f: 图上没有这一段' % It)
            continue
        Dd = float(np.median(DIG[nm][sel, 0]))
        t = time.time()
        jd, ratio = closure_J(Z, It)
        RES[(nm, It)] = (Dd, jd)
        print('    I=%.2f: 图 D=%8.3f | 闭包锥 %8.3f (%+6.3f, %+.2f%%) | 形状 %s (%.0fs)'
              % (It, Dd, jd, jd - Dd, 100 * (jd - Dd) / Dd,
                 np.array2string(ratio[:r], precision=2), time.time() - t))

print('\n[C] 嵌套性检查：J_blue(r=3) <= J_red(r=2) <= J_magenta(r=4)？'
      '（span 单调，闭包口径必须成立；违反 = 优化器没收敛）')
for It in (2.0, 2.5, 3.0, 4.0):
    row = ' | '.join('%s %8.3f' % (nm, RES[(nm, It)][1]) for nm in ('red', 'blue', 'magenta')
                     if (nm, It) in RES)
    print('  I=%.2f: %s' % (It, row))

print('\n[D] 出窗边：每条曲线在窗口顶 I=4.926 的闭包成本，与数字化左端对照')
for nm, (c, r) in CURVE.items():
    Z = U[:, n - r:]
    Dtop, ratio = closure_J(Z, 4.926, starts=2)
    i = int(np.argmin(DIG[nm][:, 0]))
    print('  %-8s 闭包 D(4.926)=%9.3f 形状 %s | 数字化左端 (D,I)=(%.2f,%.3f) 墙 %.4f'
          % (nm, Dtop, np.array2string(ratio[:r], precision=2),
             DIG[nm][i, 0], DIG[nm][i, 1], JC + phi_iter(Z.T)[0]))
print('== E63e done ==')
