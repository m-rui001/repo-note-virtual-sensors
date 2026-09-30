r"""E63f = (1) red 平面的锥**穷尽**全局最小（含秩亏边界）；(2) C 的 c33 "S=c*Theta 对偶下界"真伪。

判据先写下：
 [A] 穷尽值 < E63b 的 Powell 值(62.896/54.287/50.602) -> 那条"锥口径复刻数字化红线（差<=0.14）"
     是优化器未收敛的巧合，必须收回；相等 -> 保留。
 [B] 无约束（rank<=4 自由锥）最优 J 若在某个 I 低于 C 的对偶曲线 J_dual(I) -> S=c*Theta 不是
     无约束极小点，他的"下界"身份作废（他的 {P:I(P)>=I_t} 不凸：I 的前项 logdet 凹、后项 -logdet 凸）。
 [C] 无约束最优 S* 的形状与 Theta 的偏离（特征值比、主方向夹角）给出 S=c*Theta 的反例。

路线：S = s*F^T F，F = N^T Z^T（M=NN^T 的秩 k 因式，k 可小于 r => 闭包），V=I_k，
用 sdare(A^T,F^T,W,I) 拿 P~ —— 与 C 的信息型定点 JI 对表到 1e-9；再用闭环谱半径判可检测性，
不可检测就印 inf（这是闭包上唯一会出幻觉的地方）。
"""
import sys, time
sys.stdout.reconfigure(encoding='utf-8')
import numpy as np
from scipy.linalg import schur, solve_discrete_are as sdare
from scipy.optimize import minimize
from PIL import Image

np.set_printoptions(precision=4, suppress=True, linewidth=170)
_src = open('p0/exp_c_audit.py', encoding='utf-8').read().split("print(r'== E53")[0]
_ns = {'__name__': 'p'}
exec(compile(_src, 'p0/exp_c_audit.py[preamble]', 'exec'), _ns)
A, W, TH, JC, n, sym = _ns['A'], _ns['W'], _ns['TH'], _ns['JC'], _ns['n'], _ns['sym']
print('== E63f：锥穷尽（秩亏允许）+ C 的 c33 对偶下界检验 ==')
In_ = np.eye(n)
EPSR = 1e-11


def JI(S, tol=1e-12, maxit=20000):
    """C c33 同构的信息型定点，对任意半正定 S 安全；不收敛=inf。"""
    Pt = W.copy()
    for _ in range(maxit):
        Pm = sym(np.linalg.inv(np.linalg.inv(Pt) + S))
        Pn = sym(A @ Pm @ A.T + W)
        if not np.all(np.isfinite(Pn)) or np.max(np.abs(Pn)) > 1e14:
            return np.inf, np.inf
        if np.max(np.abs(Pn - Pt)) < tol * max(1.0, np.max(np.abs(Pn))):
            Pt = Pn; break
        Pt = Pn
    else:
        return np.inf, np.inf
    return (float(np.trace(TH @ Pm)) + JC,
            float(0.5 * (np.linalg.slogdet(Pt)[1] - np.linalg.slogdet(Pm)[1]) / np.log(2)))


def JIfast(S):
    """S 的秩 k 因式 -> F(k x n), V=I_k -> sdare 拿 P~ -> (J,I)；闭环不稳定则 inf。"""
    ev, EV = np.linalg.eigh(sym(S))
    keep = ev > max(EPSR * ev.max(), 1e-300)
    if not keep.any():
        return np.inf, np.inf
    F = (EV[:, keep] * np.sqrt(ev[keep]))[:, :].T          # k x n
    k = F.shape[0]
    try:
        Pt = sym(sdare(A.T, F.T, sym(W), np.eye(k)))
    except Exception:
        return np.inf, np.inf
    if not np.all(np.isfinite(Pt)):
        return np.inf, np.inf
    FPt = F @ Pt
    try:
        P = sym(Pt - FPt.T @ np.linalg.solve(FPt @ F.T + np.eye(k), FPt))
    except Exception:
        return np.inf, np.inf
    cl = A - A @ Pt @ F.T @ np.linalg.solve(FPt @ F.T + np.eye(k), F)
    if np.max(np.abs(np.linalg.eigvals(cl))) >= 1.0 - 1e-9:
        return np.inf, np.inf        # 不是稳定解 => 该 S 不可检测
    ld = np.linalg.slogdet(Pt)[1] - np.linalg.slogdet(P)[1]
    if ld <= 0:
        return np.inf, np.inf
    return JC + float(np.trace(TH @ P)), 0.5 * ld / np.log(2)


Ts, U = schur(A, output='real', sort=lambda a: abs(a) < 1.0)[:2]
Z2 = U[:, n - 2:]
print('\n[0] 双实现自检（同一 S，信息定点 vs 因式 sdare）+ 我自己的 sym(det) bug 定量')
def Irate_wrong(Pt, S):
    """E62/E63b/c 里我用过的写法：det(I+sym(Pt@S)) —— 对称化改变了行列式。"""
    return 0.5 * np.linalg.slogdet(In_ + sym(Pt @ sym(S)))[1] / np.log(2)


worst = 0.0
t0 = time.time()
for M in [np.diag([3.0, 0.7]), np.diag([1.0, 1e-3]), np.diag([1e-2, 1.0]), np.eye(2)]:
    for sc in (10.0, 1e3):
        S = sym(Z2 @ (sc * M) @ Z2.T)
        a = JI(S); b = JIfast(S)
        d = abs(a[0] - b[0]) if np.isfinite(a[0]) and np.isfinite(b[0]) else np.nan
        worst = max(worst, 0 if not np.isfinite(d) else d)
        ev2, EV2 = np.linalg.eigh(sym(S))
        Fc = (EV2[:, ev2 > 1e-11 * ev2.max()] * np.sqrt(ev2[ev2 > 1e-11 * ev2.max()])).T
        Ptc = sym(sdare(A.T, Fc.T, sym(W), np.eye(Fc.shape[0])))
        iw = Irate_wrong(Ptc, S)
        print('   M=%s x%-6.0e: 定点 J=%10.5f I=%8.5f | 因式 J=%10.5f I=%8.5f | dJ=%.2e | sym-bug 虚高 %+6.3f bit'
              % (np.array2string(np.diag(M), precision=3).replace(' ', ''), sc, a[0], a[1], b[0], b[1], d,
                 iw - b[1]))
print('   最大 |dJ| = %.2e；因式路线均摊 %.1f ms/eval' % (worst, 1000 * (time.time() - t0) / 8))

import os
_IMG = next(p for p in ('p0/fig/lcss_p5_x327.png', '.work3/fig_p5_img0.png')
               if os.path.exists(p))    # the repro package ships the raster under p0/fig/
Aimg = np.asarray(Image.open(_IMG).convert('RGB')).astype(int)
x0, x1, y0, y1 = 121.0, 1924.0, 52.0, 1238.0   # frame bottom as auto-detected by exp_fig1_digitize.py
px2D = lambda px: 40.0 + (px - x0) * 50.0 / (x1 - x0)
px2I = lambda py: 5.0 - (py - y0) * 4.0 / (y1 - y0)
box = np.zeros(Aimg.shape[:2], bool)
box[int(y0) + 2:int(y1) - 1, int(x0) + 3:int(x1) - 2] = True
leg = np.zeros_like(box)
leg[int(y0):int(y0 + 0.34 * (y1 - y0)) + 60, int(x0 + 0.40 * (x1 - x0)) + 110:] = True
DIG = {}
for nm, c in [('red', (255, 0, 0)), ('blue', (0, 0, 255)), ('magenta', (255, 0, 255))]:
    m = np.all(Aimg == np.array(c), axis=2) & box & (~leg)
    DIG[nm] = np.array([(px2D(px), px2I(np.where(m[:, px])[0].mean()))
                        for px in np.where(m.any(0))[0]])


def digD(nm, It, w=0.04):
    sel = np.abs(DIG[nm][:, 1] - It) < w
    return float(np.median(DIG[nm][sel, 0])) if sel.sum() else np.nan


def J_hit(S1, It, lo=-3.0, hi=13.0, steps=26):
    """给定形状 S1，二分尺度命中 I=It；返回最小 J（不可行 -> nan）。"""
    a = JIfast(10.0 ** hi * S1)
    if not np.isfinite(a[1]) or a[1] < It:
        return np.nan
    for _ in range(steps):
        m = 0.5 * (lo + hi)
        b = JIfast(10.0 ** m * S1)
        if np.isfinite(b[1]) and b[1] >= It:
            hi = m
        else:
            lo = m
    return JIfast(10.0 ** hi * S1)[0]


print('\n[A] red 平面 K_{F2} 穷尽：M=R(th)diag(1,rho)R(th)^T，rho=0 就是关掉一条通道')
RHO = [0.0, 1e-4, 1e-2, 1e-1, 0.5, 1.0]
TH_ = np.linspace(0, np.pi, 91)[:-1]
IS = (2.0, 2.5, 3.0, 4.0)
best = {It: (np.inf, None) for It in IS}
t0 = time.time()
for rho in RHO:
    for th in TH_:
        Rv = np.array([[np.cos(th), -np.sin(th)], [np.sin(th), np.cos(th)]])
        S1 = sym(Z2 @ sym(Rv @ np.diag([1.0, rho]) @ Rv.T) @ Z2.T)
        for It in IS:
            j = J_hit(S1, It)
            if np.isfinite(j) and j < best[It][0]:
                best[It] = (j, (rho, th))
POW = {2.0: 62.896, 2.5: 54.287, 3.0: 50.602, 4.0: 47.929}
for It in IS:
    j, arg = best[It]
    Dd = digD('red', It)
    print('  I=%.2f: 穷尽 %9.4f (rho=%.3g, th=%.1f deg) | E63b Powell %8.3f | 穷尽-Powell %+8.4f | 图 %8.3f | 穷尽-图 %+.3f'
          % (It, j, arg[0], np.degrees(arg[1]), POW[It], j - POW[It], Dd, j - Dd))
print('  耗时 %.0f s。穷尽-Powell 若为负 => E63b 的"复刻红线"是未收敛巧合。' % (time.time() - t0))

print('\n[A2] 同一张表里把射线口径（C 的 §11 用的 V=vI，即 M=I）也用**修好的率**重算：'
      '到底哪一种口径贴住数字化红线')
for It in IS:
    jray = J_hit(sym(Z2 @ np.eye(2) @ Z2.T), It)
    jco = best[It][0]
    Dd = digD('red', It)
    print('  I=%.2f: 图 %8.3f | 射线 %9.4f (%+7.3f) | 锥穷尽 %9.4f (%+7.3f)'
          % (It, Dd, jray, jray - Dd, jco, jco - Dd))
for It in (1.55, 4.90):
    sel = np.abs(DIG['red'][:, 1] - It) < 0.02
    if sel.sum():
        Dd = float(np.median(DIG['red'][sel, 0]))
        jr = J_hit(sym(Z2 @ np.eye(2) @ Z2.T), It)
        print('  I=%.2f: 图 %8.3f | 射线 %9.4f (%+7.3f)' % (It, Dd, jr, jr - Dd))

print('\n[B] 无约束（S 任意半正定，rank<=4）vs C 的对偶族 S=c*Theta')
wT, VT = np.linalg.eigh(TH); VT = VT[:, ::-1]
dual_J = {}
for It in (2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0):
    dual_J[It] = J_hit(TH / np.trace(TH), It)
FREE = {}
for It in (2.0, 3.0, 4.0):
    rng = np.random.default_rng(int(It * 37))
    nz = n * (n + 1) // 2

    def obj(x):
        Lm = np.zeros((n, n)); k = 0
        for i in range(n):
            for jj in range(i + 1):
                Lm[i, jj] = np.exp(x[k]) if i == jj else x[k]; k += 1
        S1 = sym(Lm @ Lm.T)
        tr = np.trace(S1)
        if tr <= 0 or not np.isfinite(tr):
            return 1e6
        j = J_hit(S1 / tr, It)
        return 1e6 if not np.isfinite(j) else j

    starts = [np.full(nz, -1.5)]
    d0 = np.log(np.sqrt(np.clip(wT[::-1], 1e-9, None)))
    x0 = np.zeros(nz); k = 0
    for i in range(n):
        for jj in range(i + 1):
            x0[k] = d0[i] if i == jj else 0.0; k += 1
    starts += [x0, 2 * x0, rng.normal(-1.0, 1.2, nz)]
    bb = (np.inf, None)
    for x0 in starts:
        r_ = minimize(obj, x0, method='Powell', options=dict(maxfev=400, xtol=2e-4, ftol=1e-8))
        if r_.fun < bb[0]:
            bb = (float(r_.fun), r_.x)
    Lm = np.zeros((n, n)); k = 0
    for i in range(n):
        for jj in range(i + 1):
            Lm[i, jj] = np.exp(bb[1][k]) if i == jj else bb[1][k]; k += 1
    S1 = sym(Lm @ Lm.T); S1 = S1 / np.trace(S1)
    ev = np.clip(np.linalg.eigvalsh(S1), 0, None)[::-1]
    u1 = np.linalg.eigh(S1)[1][:, -1]
    ang = np.degrees(np.arccos(np.clip(abs(float(np.dot(u1, VT[:, 0]))), 0, 1)))
    FREE[It] = (bb[0], ev / ev[0], ang)
    print('  I=%.2f: 对偶(Theta) %9.4f | 无约束 %9.4f | 无约束-对偶 %+9.4f | rank(>=1e-6) %d | 形状 %s | 顶方向离 Theta_1 %.2f deg'
          % (It, dual_J[It], bb[0], bb[0] - dual_J[It], int((ev > 1e-6).sum()),
             np.array2string(ev[:4], precision=3), ang))
print('  判据：任一行"无约束-对偶"为负 => C 的 S=c*Theta 不是无约束极小点，[C] 表的"下界"身份作废。')

print('\n[B2] 嵌套性门（正确性检查）：span(U[:,1:]) 真包含 span(U[:,2:]) => 同一 I 上 J_blue<=J_red。'
      '\n     另外无约束值必须 <= 任何平面值。三条同时成立才算数。')
Z3 = U[:, n - 3:]
for It in (2.0, 3.0, 4.0):
    rng = np.random.default_rng(int(It * 91))
    nz3 = 6

    def obj3(x):
        Lm = np.zeros((3, 3)); k = 0
        for i in range(3):
            for jj in range(i + 1):
                Lm[i, jj] = np.exp(x[k]) if i == jj else x[k]; k += 1
        S1 = sym(Lm @ Lm.T); tr = np.trace(S1)
        j = J_hit(sym(Z3 @ (S1 / tr) @ Z3.T), It)
        return 1e6 if not np.isfinite(j) else j

    bb = (np.inf, None)
    for x0 in [np.full(nz3, -1.5), np.full(nz3, -0.3)] + [rng.normal(-1.0, 1.2, nz3) for _ in range(3)]:
        r_ = minimize(obj3, x0, method='Powell', options=dict(maxfev=350, xtol=2e-4, ftol=1e-8))
        if r_.fun < bb[0]:
            bb = (float(r_.fun), r_.x)
    ok1 = bb[0] <= best[It][0] + 1e-6
    ok2 = FREE[It][0] <= bb[0] + 1e-6
    print('  I=%.2f: blue 平面 %9.4f | red 穷尽 %9.4f -> 嵌套 %s | 无约束 %9.4f -> 嵌套 %s'
          % (It, bb[0], best[It][0], 'OK' if ok1 else '违反(优化器未收敛)',
             FREE[It][0], 'OK' if ok2 else '违反(优化器未收敛)'))

print('\n[C] 凸性诊断：把 C 的 {P:I(P)>=I_t} 直接当集合测（在 S 空间取两点，看可行集的 J 行为）')
S1 = sym(Z2 @ np.diag([1.0, 1e-3]) @ Z2.T)
S2 = sym(Z2 @ np.diag([1e-3, 1.0]) @ Z2.T)
for al in (0.0, 0.25, 0.5, 0.75, 1.0):
    Smix = sym(al * S1 + (1 - al) * S2)
    lo, hi = -3.0, 13.0
    for _ in range(30):
        m = 0.5 * (lo + hi)
        b = JIfast(10.0 ** m * Smix)
        if np.isfinite(b[1]) and b[1] >= 3.0:
            hi = m
        else:
            lo = m
    j, i_ = JIfast(10.0 ** hi * Smix)
    print('  a=%.2f: 单位形状混合后命中 I=3.0 -> 尺度 %.3e, J=%9.4f, I=%.4f' % (al, 10.0 ** hi, j, i_))
print('== E63f done ==')
