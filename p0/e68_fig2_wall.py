r"""E68 = 注记 Fig.2 的数据：一族可行子空间沿各自射线逼近**自己的**墙。

目的只有一个，就是支撑 §20.3 / 注记 Proposition 1(i) 与 Corollary 1：
$D_{\min}(V)=j_c+\Phi_0(V)$ 是 $D(I)$ 的**水平渐近线**（$\tau\to\infty\Rightarrow I\to\infty$，
因为 $\det P_m\sim1/\tau$），所以原文 $I\le5$ 的窗口里墙是不可见的。
不做新的优化、不加小数、不动 0.1832。
"""
import sys, time, csv
sys.stdout.reconfigure(encoding='utf-8')
import numpy as np
from scipy.linalg import schur

np.set_printoptions(precision=4, suppress=True, linewidth=170)
_src = open('p0/exp_c_audit.py', encoding='utf-8').read().split("print(r'== E53")[0]
_ns = {'__name__': 'p'}
exec(compile(_src, 'p0/exp_c_audit.py[preamble]', 'exec'), _ns)
A, W, TH, JC, n, sym = _ns['A'], _ns['W'], _ns['TH'], _ns['JC'], _ns['n'], _ns['sym']
phi_iter = _ns['phi_iter']
print('== E68: 射线沿子空间逼近各自的墙（注记 Fig.2 数据）==')


def JI(S, tol=1e-11, maxit=8000):
    r"""信息型定点：返回 $(D, I)$；不收敛（不可检测）返回 $(\infty,\infty)$。"""
    Pt = W.copy()
    try:
        for _ in range(maxit):
            Pm = sym(np.linalg.solve(np.linalg.inv(Pt) + S, np.eye(n)))
            Pn = sym(A @ Pm @ A.T + W)
            if not np.all(np.isfinite(Pn)) or np.max(np.abs(Pn)) > 1e14:
                return np.inf, np.inf
            if np.max(np.abs(Pn - Pt)) < tol * max(1.0, np.max(np.abs(Pn))):
                return float(np.trace(TH @ Pm)) + JC, \
                    float(0.5 * (np.linalg.slogdet(Pn)[1] - np.linalg.slogdet(Pm)[1]) / np.log(2))
            Pt = Pn
    except np.linalg.LinAlgError:
        return np.inf, np.inf
    return np.inf, np.inf


# 子空间族：原文例子的三种具名选择 + 我缓存的 r=1,2,3 最优支撑
Tu, Uu, _ = schur(A, output='real', sort=lambda a: abs(a) < 1.0)
FAMS = [('schur_tail_1', Uu[:, n - 1:].T),
        ('schur_tail_2', Uu[:, n - 2:].T),
        ('schur_tail_3', Uu[:, n - 3:].T)]
try:
    BEST = np.load('p0/e53_best.npy', allow_pickle=True).item()
    for r in sorted(BEST):
        FAMS.append(('best_support_r%d' % r, np.atleast_2d(BEST[r])))
except Exception as e:
    print('  (缓存最优支撑未载入：%s)' % e)

VS = np.logspace(-3, 12, 151)
rows = []
t0 = time.time()
print('  %-16s %-3s %-12s %-11s %-11s %-11s' %
      ('family', 'r', 'wall D_min', 'I@v=1e-3', 'I@v=1e12', 'gap@I=5'))
for lab, F in FAMS:
    F = np.linalg.qr(F.T)[0].T                      # 行正交化，F^T F 即子空间投影
    r = F.shape[0]
    w, _, _, tag = phi_iter(F)
    if tag != 'conv':
        print('  %-16s %-3d 墙不可用（%s），跳过' % (lab, r, tag))
        continue
    wall = JC + w
    P = F.T @ F
    pts = []
    for v in VS:
        D, I = JI(v * P)
        if not np.isfinite(D):
            continue
        pts.append([lab, r, '%.6f' % wall, '%.6e' % v, '%.6f' % I, '%.6f' % D,
                    '%.6f' % (D - wall), '%.6f' % (100 * (D - wall) / wall),
                    int(I <= 5.0)])
        rows.append(pts[-1])
    win = [x for x in pts if float(x[4]) <= 5.0]
    e5 = 'I=%.2f 处仍差 %.3f%%' % (float(win[-1][4]), float(win[-1][7])) if win else '窗口内无采样点'
    print('  %-16s %-3d %-12.4f %-11.2f %-11.2f  %s' %
          (lab, r, wall, float(pts[0][4]), float(pts[-1][4]), e5))

out = 'p0/note/data/wall_approach.csv'
with open(out, 'w', newline='') as f:
    wcsv = csv.writer(f)
    wcsv.writerow(['family', 'r', 'D_min', 'v', 'I', 'D', 'gap', 'rel_gap_pct', 'in_pub_window'])
    wcsv.writerows(rows)
print('  %d 行 -> %s（%.0f s）' % (len(rows), out, time.time() - t0))

print('\n[A] 墙要可见，I 必须多大（各族相对间隙降到 1% / 0.1% 所需的率）')
for lab, _ in FAMS:
    sub = [x for x in rows if x[0] == lab]
    if not sub:
        continue
    a = b = np.nan
    for x in sub:
        if float(x[7]) <= 1.0 and not np.isfinite(a):
            a = float(x[4])
        if float(x[7]) <= 0.1 and not np.isfinite(b):
            b = float(x[4])
    print('  %-16s 1%%: I>=%s | 0.1%%: I>=%s' %
          (lab, '-' if not np.isfinite(a) else '%.1f' % a,
           '-' if not np.isfinite(b) else '%.1f' % b))
print('== E68 done ==')
