r"""E63 = 原文 Fig.1 四条曲线的"左端点"到底是什么：C 的 §11.4 说红线左端点 46.5 命中
j_c+Phi_0(F_2)=46.1231，并据此写"红线的左端点就是原文自己那支传感器的不可行墙（不是画图裁剪，
也不是数据端点）"。这条我可以直接否证：图框的 I 轴上限若是 5 bit，那么任何曲线的左端点都是
**它自己 I=5 那一行的成本**，与墙无关。判别实验：magenta=全状态(r=4) 的墙是 j_c=31.48，
若它的左端点落在 54 而不是 31.48，就证明左端点由率窗口决定，而不是墙。

顺带把两件原文事实钉住（这是我第一次读到原文 §V 的原文措辞）：
 - 原文的 F 显式写成 F = [0 I_nu] U^T（ordered real Schur，不稳定模放最后），nu=2
   ⇒ C 的"不稳 Schur 尾"确实就是原文那支传感器（他的数没错）；
 - 但原文正文把它叫作 "task representation is selected as the unstable subspace of A" 与
   "unstable modes only"。Schur 尾坐标张成的平面**不是** A-不变平面（我 E61：不变性残差 1.61，
   ‖TT^T−T^T T‖_F=7.92），真正的不稳定不变平面给 Phi=78.6010 而不是 14.6398。
   ⇒ 原文的规则名与其公式不是同一个对象，而且公式反而更好。这条要在注记里写清。
 - 原文 Fig.1 虚线 I_inf=1.169 bit：核对 sum_{|lam|>=1} log2|lam| 是否等于它。

只读 .work3 的 PNG，产物全部落在我自己的 p0/。
跑法：`python -X utf8 p0/exp_fig1_ends.py > p0/e63_out.txt 2>&1`
"""
import sys
sys.stdout.reconfigure(encoding='utf-8')
import numpy as np
from PIL import Image
from scipy.linalg import schur, solve_discrete_are
from scipy import ndimage

np.set_printoptions(precision=4, suppress=True, linewidth=170)
_src = open('p0/exp_c_audit.py', encoding='utf-8').read().split("print(r'== E53")[0]
_ns = {'__name__': 'e53_preamble'}
exec(compile(_src, 'p0/exp_c_audit.py[preamble]', 'exec'), _ns)
A, W, TH, JC, n, sym = _ns['A'], _ns['W'], _ns['TH'], _ns['JC'], _ns['n'], _ns['sym']
phi_iter = _ns['phi_iter']

print('== E63：Fig.1 左端点的成因（率窗口 vs 不可行墙）+ 原文 §V 规则的实体核对 ==')

# ---------------------------------------------------------------- 0) 原文事实的算术核对
lam = np.linalg.eigvals(A)
absu = np.abs(lam[np.abs(lam) >= 1.0])
print('[0] eig(A)=%s' % np.round(lam, 4))
print('    sum_{|lam|>=1} log2|lam| = %.4f bit   (原文 Fig.1 虚线 I_inf = 1.169)' % float(np.sum(np.log2(absu))))
print('    其中 |lam| = %s' % np.round(np.sort(absu)[::-1], 4))

Ts, Us = schur(A, output='real', sort=lambda a: abs(a) < 1.0)[:2]
U = Us                      # A = U T U^T，不稳定块放最后（sort 把 |lam|<1 放前面）
nu = int(np.sum(np.abs(np.linalg.eigvals(A)) >= 1.0))
print('    Schur 排序后 diag(T) 块：%s   nu=%d' % (np.round(np.diag(Ts), 4), nu))
Fs = {}
for k in (2, 3, 4):         # 原文 (i)(ii)(iii)：F = [0 I_k] U^T
    Fs['schur_tail_%d' % k] = U[:, n - k:].T
Fs['eigen_unstable_2'] = None
# 真正的不稳定不变平面（本征向量张成），用于对照
Vl = np.linalg.eig(A)[1]
sel = [i for i in range(n) if abs(np.linalg.eigvals(A)[i]) >= 1.0]
cols = []
for i in sel:
    vv = Vl[:, i]
    if np.abs(vv.imag).max() > 1e-9:
        cols += [vv.real, vv.imag]
    else:
        cols.append(vv.real)
Binv = np.linalg.qr(np.array(cols).T)[0]
Fs['eigen_unstable_2'] = Binv[:, :2].T

from scipy.linalg import solve_discrete_are as sdare


def IJ(F, tau):
    F = np.atleast_2d(F); r = F.shape[0]
    V = np.eye(r) / tau
    try:
        Pt = sym(sdare(A.T, F.T, sym(W), V))
    except Exception:
        return None
    if not np.all(np.isfinite(Pt)):
        return None
    FPt = F @ Pt
    P = sym(Pt - FPt.T @ np.linalg.solve(FPt @ F.T + V, FPt))
    sign, ld = np.linalg.slogdet(np.eye(r) + tau * (F @ Pt @ F.T))
    if sign <= 0:
        return None
    return 0.5 * ld / np.log(2), JC + float(np.trace(TH @ P))


def at_rate(F, Itarget, lo=-4.0, hi=16.0, steps=90):
    for _ in range(steps):
        m = 0.5 * (lo + hi)
        a = IJ(F, 10.0 ** m)
        if a is None or a[0] < Itarget:
            lo = m
        else:
            hi = m
    return IJ(F, 10.0 ** hi)


print('\n[0a] 两个"不稳定子空间"实体的分别核对（原文正文说的是后者，公式给的是前者）')
for tag, S in [('Schur 尾平面 U[:,2:]', U[:, 2:]), ('真不稳定不变平面', Binv[:, :2])]:
    res_iv = np.linalg.norm((np.eye(n) - S @ S.T) @ A @ S)
    p = phi_iter(S.T)[0]
    print('  %-22s  A-不变性残差 ||(I-SS^T)AS|| = %.4f   Phi_0 = %9.4f   j_c+Phi = %9.4f'
          % (tag, res_iv, p, JC + p))
print('  原文 §V 公式 F = [0 I_nu] U^T 展开即第一行；正文措辞 "task representation is selected as'
      '\n  the unstable subspace of A" / "unstable modes only" 指的是第二行。两者在本例差 %.1f 倍地板。'
      % (phi_iter(Binv[:, :2].T)[0] / phi_iter(U[:, 2:].T)[0]))

print('\n[0b] 原文三支架构的地板与"墙上所需率"')
for k, F in Fs.items():
    fl = JC + phi_iter(np.atleast_2d(F))[0]
    print('  %-18s r=%d  地板 j_c+Phi_0 = %10.4f   行范数 %s' %
          (k, F.shape[0], fl, np.round(np.linalg.norm(F, axis=1), 3)))

# ---------------------------------------------------------------- 1) 图的定标（自己做，不用 C 的常数）
Aimg = np.asarray(Image.open('.work3/fig_p5_img0.png').convert('RGB')).astype(int)
H, Wd, _ = Aimg.shape
R, G, B = Aimg[:, :, 0], Aimg[:, :, 1], Aimg[:, :, 2]
dark = (R < 100) & (G < 100) & (B < 100)


def groups(idx, gap=3):
    out = []
    if len(idx) == 0:
        return out
    s = p = idx[0]
    for v in idx[1:]:
        if v - p > gap:
            out.append((s, p)); s = v
        p = v
    out.append((s, p))
    return out


gc = groups(np.where(dark.sum(0) > 0.6 * H)[0]); gr = groups(np.where(dark.sum(1) > 0.6 * Wd)[0])
x0, x1 = gc[0][0], gc[-1][-1]; y0, y1 = gr[0][0], gr[-1][-1]
print('\n[1] 图框像素 x[%g,%g] y[%g,%g]  尺寸 %dx%d' % (x0, x1, y0, y1, Wd, H))

# 浅灰网格线：给刻度位置（不假设轴范围，先量出来）
gray = (np.abs(R - G) < 14) & (np.abs(G - B) < 14) & (R > 140) & (R < 235)
gv = groups(np.where(gray[:, int(x0) + 3:int(x1) - 1].sum(0) > 0.5 * (y1 - y0))[0])
gh = groups(np.where(gray[int(y0) + 3:int(y1) - 1, int(x0) + 3:int(x1) - 1].sum(1) > 0.5 * (x1 - x0))[0])
gv = [int((g[0] + g[1]) / 2) + int(x0) + 3 for g in gv]
gh = [int((g[0] + g[1]) / 2) + int(y0) + 3 for g in gh]
print('    竖向网格线 px=%s  间距=%s' % (gv, np.diff(gv)))
print('    横向网格线 px=%s  间距=%s' % (gh, np.diff(gh)))

# 虚线（原文 I_inf=1.169）：非框、非网格的水平长 run
hline = []
for py in range(int(y0) + 2, int(y1) - 1):
    row = dark[py, int(x0) + 3:int(x1) - 2]
    if row.sum() > 0.45 * (x1 - x0):
        hline.append((py, int(row.sum())))
print('    候选水平虚线行 (py, 黑像素数) = %s' % hline[:8])

# 曲线颜色：箱内非灰主色直方图
inbox = np.zeros_like(dark); inbox[int(y0) + 2:int(y1) - 1, int(x0) + 3:int(x1) - 2] = True
leg = np.zeros_like(dark); leg[int(y0):int(y0 + 0.34 * (y1 - y0)) + 60, int(x0 + 0.40 * (x1 - x0)) + 110:] = True
keep = inbox & (~leg) & (~dark) & (~gray)
pix = Aimg[keep]
colr, cnt = np.unique(pix, axis=0, return_counts=True)
o = np.argsort(-cnt)
print('    箱内前 8 主色 RGB / 计数：')
for i in o[:8]:
    print('      %s : %d' % (colr[i], cnt[i]))

# ---------------------------------------------------------------- 2) 按颜色的曲线端点（宽容掩码）
MAIN = {str(tuple(colr[i])): tuple(colr[i]) for i in o[:6]}
print('\n[2] 每条主色的曲线覆盖（列中位数），左端点 D 与对应 I')
AX = dict(D0=40.0, D1=90.0, I0=5.0, I1=1.0)   # 先用 C 的假设，下面用网格线/虚线自检
px2D = lambda px: AX['D0'] + (px - x0) * (AX['D1'] - AX['D0']) / (x1 - x0)
px2I = lambda py: AX['I0'] + (py - y0) * (AX['I1'] - AX['I0']) / (y1 - y0)
if len(gv) >= 2:
    sp = np.median(np.diff(gv))
    print('    竖网格间距 %.2f px；若 1 格 = 5 单位(D 轴)，则 px2D 斜率核对：%.6f vs %.6f'
          % (sp, 5.0 / sp, (AX['D1'] - AX['D0']) / (x1 - x0)))
if len(gh) >= 2:
    sp2 = np.median(np.diff(gh))
    print('    横网格间距 %.2f px；若 1 格 = 1 单位(I 轴)，则 px2I 斜率核对：%.6f vs %.6f'
          % (sp2, -1.0 / sp2, (AX['I1'] - AX['I0']) / (y1 - y0)))
for py, c in hline[:4]:
    print('    水平虚线 py=%d -> 标定 I=%.4f（原文声称 1.169）' % (py, px2I(py)))

res = {}
for i in o[:6]:
    c = colr[i].astype(float)
    m = np.all(np.abs(Aimg - c) <= 45, axis=2) & keep
    xs = []
    ys = []
    for px in range(int(x0) + 2, int(x1) - 1):
        colsel = np.where(m[:, px])[0]
        if len(colsel):
            xs.append(px); ys.append(np.median(colsel))
    if len(xs) < 3:
        print('    RGB%s: 只有 %d 列，跳过' % (tuple(c.astype(int)), len(xs))); continue
    xs = np.array(xs, float); ys = np.array(ys)
    Ds, Is = px2D(xs), px2I(ys)
    res[tuple(c.astype(int))] = (Ds, Is)
    lef = np.argmin(Ds)
    print('    RGB%s: %d 列, D 覆盖 [%.2f, %.2f], 左端 (D=%.2f, I=%.3f), 右端 (D=%.2f, I=%.3f)'
          % (tuple(c.astype(int)), len(xs), Ds.min(), Ds.max(), Ds.min(), Is[lef], Ds.max(), Is[np.argmax(Ds)]))

# ---------------------------------------------------------------- 3) 理论预测：左端点 = 窗口上沿的率
print('\n[3] 预测：若 I 轴上限为 Imax，则各架构左端点 = J(Imax)，与墙无关')
NAME = {'red_2unstable': 'schur_tail_2', 'blue_2u1s': 'schur_tail_3', 'magenta_2u2s': 'schur_tail_4'}
for Imax in (4.5, 5.0):
    print('  Imax=%.1f bit:' % Imax)
    for tag, key in NAME.items():
        F = np.atleast_2d(Fs[key]); r = F.shape[0]
        a = at_rate(F, Imax)
        fl = JC + phi_iter(F)[0]
        print('    %-16s r=%d  J(Imax)=%9.4f   地板=%9.4f  （若左端点=墙则两者应相等）' % (tag, r, a[1], fl))

print('\n[4] 反证要点：magenta=全状态 r=4，它的地板就是 j_c=%.4f，'
      '而它的左端点在 D 轴 ~54 处（C 的 c1/extract2 数字化）——左端点比地板高 22.7，'
      '所以"左端点=不可行墙"这条推断在同一张图里就被 magenta 否证。' % JC)
print('    若按窗口解释：magenta 在 I=5 的成本 = %.4f，与 54.19 比较。' % at_rate(np.atleast_2d(Fs['schur_tail_4']), 5.0)[1])
print('== E63 done ==')
