r"""E65 = 只用信息型定点（C 的 c33 同构，无我那个坏掉的判据）复核四件事。

为什么重写：E63f/E63e 的"因式 sdare"路线里我加了一条闭环谱半径判据
`cl = A - A@Pt@F.T@(...)`，在同一株植物上它对有的方向给 rho=21.7、对别的给 rho<1，
方向稍动就翻转 —— 判据不可靠，删掉。定点路线本来就在不可检测时不收敛（返回 inf），
语义正确；而且它和我因式路线在 E63f[0] 的 8 个格子上 J 与 I 都对到 1e-10。

[A] red 平面锥穷尽（去掉判据）—— E63f 的 62.8975/54.2924/50.6023/47.9396 是否站得住
[B] C 的 S=c*Theta 分支：他自己的 10^0.5 网格 + 线性插值 vs 二分命中（他表里 I=2.00 的 70.9922）
[C] 等地板方向配对（打 §12.4"率只通过地板依赖子空间"）：地板差 <0.5%/0.1% 时的率散布 + 交叉对
[D] 三条曲线各自的"射线口径 vs 锥穷尽"：判定原文 Fig.1 画的是哪个口径
"""
import sys, time
sys.stdout.reconfigure(encoding='utf-8')
import numpy as np
from scipy.linalg import schur
from scipy.optimize import minimize
from PIL import Image

np.set_printoptions(precision=4, suppress=True, linewidth=170)
_src = open('p0/exp_c_audit.py', encoding='utf-8').read().split("print(r'== E53")[0]
_ns = {'__name__': 'p'}
exec(compile(_src, 'p0/exp_c_audit.py[preamble]', 'exec'), _ns)
A, W, TH, JC, n, sym = _ns['A'], _ns['W'], _ns['TH'], _ns['JC'], _ns['n'], _ns['sym']
print('== E65：纯定点路线复核 ==')


def JI(S, tol=1e-12, maxit=20000):
    Pt = W.copy()
    try:
        for _ in range(maxit):
            Pm = sym(np.linalg.solve(np.linalg.inv(Pt) + S, np.eye(n)))
            Pn = sym(A @ Pm @ A.T + W)
            if not np.all(np.isfinite(Pn)) or np.max(np.abs(Pn)) > 1e14:
                return np.inf, np.inf
            if np.max(np.abs(Pn - Pt)) < tol * max(1.0, np.max(np.abs(Pn))):
                Pt = Pn; break
            Pt = Pn
        else:
            return np.inf, np.inf
    except np.linalg.LinAlgError:
        return np.inf, np.inf
    return (float(np.trace(TH @ Pm)) + JC,
            float(0.5 * (np.linalg.slogdet(Pt)[1] - np.linalg.slogdet(Pm)[1]) / np.log(2)))


def J_hit(S1, It, lo=-4.0, hi=13.0, steps=26):
    if not np.isfinite(JI(10.0 ** hi * S1)[0]):
        return np.nan
    for _ in range(steps):
        m = 0.5 * (lo + hi)
        j, i_ = JI(10.0 ** m * S1)
        if np.isfinite(i_) and i_ >= It:
            hi = m
        else:
            lo = m
    return JI(10.0 ** hi * S1)[0]


Ts, U = schur(A, output='real', sort=lambda a: abs(a) < 1.0)[:2]
Z = {2: U[:, n - 2:], 3: U[:, n - 3:], 4: U[:, :]}
wT, VT = np.linalg.eigh(TH); VT = VT[:, ::-1]
THETA = {2: VT[:, :2], 3: VT[:, :3], 4: VT[:, :]}
GEO = {}
for r in (1, 2, 3):
    try:
        f = np.load('.work3/c24_bestF_r%d.npy' % r)
        GEO[r] = f if f.shape[0] == r else f.T
    except Exception as e:
        print('  (c24_bestF_r%d 读不到: %s)' % (r, e))

Aimg = np.asarray(Image.open('.work3/fig_p5_img0.png').convert('RGB')).astype(int)
x0, x1, y0, y1 = 121.0, 1924.0, 52.0, 1239.0
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


def digD(nm, It, w=0.03):
    sel = np.abs(DIG[nm][:, 1] - It) < w
    return float(np.median(DIG[nm][sel, 0])) if sel.sum() else np.nan


print('\n[A] red 平面（Schur 尾 2 列 = 原文 F2）锥穷尽，无判据版')
IS = (2.0, 2.5, 3.0, 4.0)
bestA = {It: (np.inf, None) for It in IS}
t0 = time.time()
for rho in [0.0, 1e-4, 1e-2, 1e-1, 0.5, 1.0]:
    for th in np.linspace(0, np.pi, 61)[:-1]:
        Rv = np.array([[np.cos(th), -np.sin(th)], [np.sin(th), np.cos(th)]])
        S1 = sym(Z[2] @ sym(Rv @ np.diag([1.0, rho]) @ Rv.T) @ Z[2].T)
        for It in IS:
            j = J_hit(S1, It)
            if np.isfinite(j) and j < bestA[It][0]:
                bestA[It] = (j, (rho, th))
for It in IS:
    print('  I=%.2f: 定点穷尽 %9.4f (rho=%-6.3g th=%.0f) | E63f 有判据版 %9.4f | 图 %8.3f | 差 %+.3f'
          % (It, bestA[It][0], bestA[It][1][0], np.degrees(bestA[It][1][1]),
             {2.0: 62.8975, 2.5: 54.2924, 3.0: 50.6023, 4.0: 47.9396}[It], digD('red', It),
             bestA[It][0] - digD('red', It)))
print('  %.0f s' % (time.time() - t0))

print('\n[B] C 的 S=c*Theta 分支：他的 10^0.5 网格插值 vs 二分命中')
his = {2.0: 70.9922, 2.5: 54.0421, 3.0: 45.3757, 3.5: 40.9311, 4.0: 37.5685, 4.5: 36.2756, 5.0: 34.9827}
for It, jh in his.items():
    jb = J_hit(TH / np.trace(TH), It)
    print('  I=%.2f: c33 插值 %9.4f | 定点二分 %9.4f | 二分-插值 %+7.4f' % (It, jh, jb, jb - jh))

print('\n[C] 等地板配对：r=1 方向族。命题"率只通过地板依赖 range(S)" <=> 等地板必等率')
rng = np.random.default_rng(5)
Zs = rng.normal(0, 1, (320, n)); Zs /= np.linalg.norm(Zs, axis=1, keepdims=True)
t0 = time.time()
FL, J2s, J5s = [], [], []
for z in Zs:
    f = JI(10.0 ** 13 * np.outer(z, z))[0]
    if not np.isfinite(f):
        continue
    FL.append(f); J2s.append(J_hit(np.outer(z, z), 2.0)); J5s.append(J_hit(np.outer(z, z), 5.0))
FL = np.array(FL); J2s = np.array(J2s); J5s = np.array(J5s)
ok = np.isfinite(J2s) & np.isfinite(J5s)
FL, J2s, J5s = FL[ok], J2s[ok], J5s[ok]
print('  可用 %d 个方向（%.0f s），地板 [%.3f, %.3f]' % (len(FL), time.time() - t0, FL.min(), FL.max()))
for tol in (0.005, 0.001, 0.0002):
    sp, npair = [], 0
    for i in range(len(FL)):
        d = np.abs(FL - FL[i]) / np.maximum(FL, FL[i])
        for k in np.where((d < tol) & (np.arange(len(FL)) > i))[0]:
            npair += 1
            sp += [abs(J2s[i] - J2s[k]) / min(J2s[i], J2s[k]), abs(J5s[i] - J5s[k]) / min(J5s[i], J5s[k])]
    if sp:
        sp = np.array(sp)
        print('  地板相对差 <%.4f: %d 对，率相对极差 中位 %.1f%% / 最大 %.1f%%' % (tol, npair, 100 * np.median(sp), 100 * sp.max()))
    else:
        print('  地板相对差 <%.4f: 无配对' % tol)
for nm, JV in (('I=2', J2s), ('I=5', J5s)):
    cross = 0; ex = None; worst = 0
    for i in range(len(FL)):
        w = np.where((FL > FL[i]) & (JV < JV[i]))[0]
        cross += len(w)
        if len(w):
            k = w[np.argmax(FL[w] - FL[i])]
            if FL[k] - FL[i] > worst:
                worst, ex = FL[k] - FL[i], (i, k)
    if ex is None:
        print('  %s: 无交叉对（该率上地板序=率序）' % nm)
    else:
        print('  %s: 交叉对 %d 个；最极端一例：地板高 %.3f（%.1f%%）却省率 J %.3f -> %.3f（%.1f%%）'
              % (nm, cross, worst, 100 * worst / FL[ex[0]], JV[ex[0]], JV[ex[1]],
                 100 * (JV[ex[1]] - JV[ex[0]]) / JV[ex[0]]))

print('\n[D] 三条曲线：射线口径(V=vI) vs 锥穷尽 vs 数字化，判定原文画的是哪个口径')
for nm, r in (('red', 2), ('blue', 3), ('magenta', 4)):
    Zp = Z[r] if r in Z else U[:, :]
    Its = (2.0, 3.0)
    nzr = r * (r + 1) // 2

    def obj(x, It):
        Lm = np.zeros((r, r)); k = 0
        for i in range(r):
            for jj in range(i + 1):
                Lm[i, jj] = np.exp(x[k]) if i == jj else x[k]; k += 1
        S1 = sym(Lm @ Lm.T); tr = np.trace(S1)
        if not np.isfinite(tr) or tr <= 0:
            return 1e6
        j = J_hit(sym(Zp @ (S1 / tr) @ Zp.T), It)
        return 1e6 if not np.isfinite(j) else j

    print('  %s (r=%d, 平面=Schur 尾 %d 列):' % (nm, r, r))
    for It in Its:
        Dd = digD(nm, It)
        jray = J_hit(sym(Zp @ np.eye(r) @ Zp.T), It)
        bb = np.inf
        for x0 in [np.full(nzr, -1.0), np.full(nzr, -3.0), np.full(nzr, 0.0)]:
            rr = minimize(obj, x0, args=(It,), method='Powell',
                          options=dict(maxfev=300, xtol=2e-4, ftol=1e-8))
            bb = min(bb, float(rr.fun))
        print('    I=%.2f: 图 %8.3f | 射线 %9.4f (%+7.3f) | 锥 Powell %9.4f (%+7.3f)'
              % (It, Dd, jray, jray - Dd, bb, bb - Dd))
print('== E65 done ==')
