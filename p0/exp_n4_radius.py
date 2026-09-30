"""E55：N4 的临界半径 —— 上一条对 D 的预言被自己的抽样否证，换成能算、能写见证的那个。

32.6 我在 E53[7] 押的是"$\\varepsilon$ 一大，球里就出现不可检测的 $\\hat A$，$\\sup$ 恒为 $+\\infty$"。
数跑出来：120 个随机谱范数扰动、$\\varepsilon$ 到 $0.5$，**无有限解的个数 $=0$**，
$\\Phi_{\\max}/$名义只到 $9.95\\times(r=1)$、$18.7\\times(r=2)$。所以那句话**按字面是错的**：
随机抽样看不见危险，因为"$A+E$ 有一条落在 $\\ker F$ 里的本征向量"在 $E$ 空间里是低维曲面，
有限次抽样命中概率为 $0$。

但 $\\sup$ 真的是 $+\\infty$，而且越过哪条线才算数有闭式。要求 $(A+E)v=\\lambda v$、$Fv=0$、$|\\lambda|\\ge1$。
给定实单位 $v$ 与**实** $\\lambda$ 时最小范数解是秩一 $E=(\\lambda I-A)vv^\\top$，范数 $\\|(A-\\lambda I)v\\|$；
对固定实 $\\lambda$ 再对 $v\\in\\ker F$ 取极小就是 $\\sigma_{\\min}\\bigl((A-\\lambda I)Z_{\\ker F}\\bigr)$，
而实矩阵的最小奇异向量是实的，所以整条链留在实数域内，见证可以直接写出来：

$$\\varepsilon^*_{\\rm real}(A,F)=\\min_{\\lambda\\in\\mathbb R,\\ |\\lambda|\\ge1}\\ \\sigma_{\\min}\\bigl((A-\\lambda I)Z_{\\ker F}\\bigr).$$

结论形状：$\\varepsilon<\\varepsilon^*_{\\rm real}$ 时球内处处可检测、$\\Phi$ 连续有界；$\\varepsilon\\ge\\varepsilon^*_{\\rm real}$ 时
$\\sup=+\\infty$，本文件把达到它的 $E$ 印出来并自检 $(A+E)v=\\lambda v$。
复 $\\lambda$ 版本只是下界（对实扰动不必然成立），一并报作对照。

所以 D 的 R3 交付物改成：**报 $(\\bar\\Phi(\\varepsilon),\\ \\varepsilon^*_{\\rm real}(A,F))$ 两个数**，
$\\varepsilon\\ge\\varepsilon^*$ 直接印 $+\\infty$ 并附见证。否则这统计量在抽样时永远"通过"、
在最坏情形下永远是 $+\\infty$，两头都错。

植物、$\\Phi$ 口径、最优支撑沿用 E53（$\\Theta$ 来自控制 DARE 的 $Q=I$；$\\Phi$ 用后验）。
"""
import sys
import numpy as np
sys.stdout.reconfigure(encoding='utf-8')
from scipy.linalg import solve_discrete_are, schur

# 只取 E53 的表头与工具，避免重跑它的扫描段
_src = open('p0/exp_c_audit.py', encoding='utf-8').read().split("print(r'== E53")[0]
_ns = {'__name__': 'e53_preamble'}
exec(compile(_src, 'p0/exp_c_audit.py[preamble]', 'exec'), _ns)
A, B, W, n, R = _ns['A'], _ns['B'], _ns['W'], _ns['n'], _ns['R']
TH, JC = _ns['TH'], _ns['JC']
sym, ctrl, post, phi_iter, pbh_bad = (_ns['sym'], _ns['ctrl'], _ns['post'],
                                      _ns['phi_iter'], _ns['pbh_bad'])
BESTF = np.load('p0/e53_best.npy', allow_pickle=True).item()


def ker_basis(F):
    _, s, Vh = np.linalg.svd(F)
    r = F.shape[0]
    return Vh[r:].T if r < n else np.zeros((n, 0))


def eps_over(F, lams):
    r"""$\\min_{\\lambda\\in\\text{lams}}\\sigma_{\\min}((A-\\lambda I)Z_{\\ker F})$，连同取到它的 $v$。"""
    Z = ker_basis(F)
    if Z.shape[1] == 0:
        return None
    best = None
    for lam in lams:
        u, sv, vh = np.linalg.svd((A - lam * np.eye(n)) @ Z, full_matrices=False)
        v = Z @ vh[-1]
        real = abs(np.imag(lam)) < 1e-14
        if real:
            v = np.real(v)
        c = (float(sv[-1]), v / np.linalg.norm(v), float(np.real(lam)) if real else lam)
        if best is None or c[0] < best[0]:
            best = c
    return best


def lam_real(N=2800):
    g = np.geomspace(1.0, 4.5, N)
    return np.concatenate([-g[::-1], g])


def lam_circle(rad=1.0, M=1440):
    return rad * np.exp(1j * np.linspace(0, 2 * np.pi, M, endpoint=False))


def lam_annulus(M=240, radii=(1.0, 1.05, 1.2, 1.5, 2.0, 3.0, 4.5)):
    out = list(lam_circle(1.0, 6 * M))
    for r in radii[1:]:
        out += list(lam_circle(r, M))
    return out


def witness(F, frac=None):
    r"""实扰动见证：把 $\\ker F$ 里一条射线弯成 $A+E$ 的不稳定本征向量。"""
    e, v, lam = eps_over(F, lam_real())
    E = np.outer(lam * v - A @ v, v)
    if frac is not None:
        E = (frac * e / max(np.linalg.norm(E, 2), 1e-30)) * E
    return E, e, lam, v


def phi_hat(Ah, F):
    bad = pbh_bad(F, Ax=Ah)
    if bad:
        return float('inf'), 'PBH 失败 $|\\lambda|=%.4f$' % bad[0]
    r = F.shape[0]
    try:
        Pt = sym(solve_discrete_are(Ah.T, F.T, W, np.zeros((r, r))))
        _, _, Thh, _ = ctrl(Ah, B, W, np.eye(n))
        res = np.abs(sym(Ah.T @ Pt @ Ah) - Pt -
                     sym(Ah.T @ Pt @ F.T @ np.linalg.solve(F @ Pt @ F.T, F @ Pt @ Ah))).max()
        return float(np.trace(Thh @ post(F, Pt))), '奇异 DARE 残差 %.1e' % res
    except Exception as ex:
        return float('inf'), 'DARE 抛错 %s' % type(ex).__name__


print('== E55：N4 的临界半径与实扰动见证 ==')
print('  植物：C 的 $n=4$ 那株，$\\operatorname{eig}(A)=%s$，$j_c=%.4f$'
      % (np.array2string(np.linalg.eigvals(A), precision=4), JC))
print('  E53[7] 抽样事实：120 个随机扰动到 $\\varepsilon=0.5$ 零失败 —— 随机法看不见这层危险。\n')

print('  [A] 三种半径：实 $\\lambda$（可写见证）、复 $\\lambda$ 单位圆、复 $\\lambda$ 环域（下界对照）')
print('  %-16s %-3s %-15s %-11s %-15s %-15s %-13s' %
      ('支撑 $F$', '$r$', '$\\varepsilon^*_{\\rm real}$', '$\\lambda_{\\rm worst}$',
       '$\\varepsilon^*_{\\rm circ}$', '$\\varepsilon^*_{\\rm ann}$', '$\\Phi$ 名义'))
_, Ut, _ = schur(A, output='real', sort=lambda a: abs(a) < 1.0)
for lbl, F in [('精修最优', BESTF[1]), ('精修最优', BESTF[2]), ('精修最优', BESTF[3])] + \
               [('Schur 尾 %d 列' % k, Ut[:, n - k:].T) for k in (1, 2, 3)]:
    r = F.shape[0]
    cr = eps_over(F, lam_real())
    cc = eps_over(F, lam_circle())
    ca = eps_over(F, lam_annulus())
    v0 = phi_iter(F)[0]
    print('  %-16s %-3d %-15.6f %-11.4f %-15.6f %-15.6f %-13s'
          % (lbl, r, cr[0], cr[2], cc[0], ca[0],
             '%.6f' % v0 if np.isfinite(v0) else 'inf'))

print('\n  [B] 沿实见证走：$\\Phi(A+E,F)$ 随 $\\varepsilon/\\varepsilon^*_{\\rm real}$（$r=1$ 精修最优支撑）')
F1 = BESTF[1]
E1, e1, lam1, v1 = witness(F1)
print('  见证自检：$\\|E_1\\|_2=%.6f$（$\\varepsilon^*=%.6f$）  $\\lambda=%.4f$  $\\max|Fv|=%.1e$  '
      '$\\max|(A+E_1)v-\\lambda v|=%.1e$  $\\max|\\operatorname{Im}v|=%.1e$'
      % (np.linalg.norm(E1, 2), e1, lam1, np.abs(F1 @ v1).max(),
         np.abs((A + E1) @ v1 - lam1 * v1).max(), np.abs(np.imag(v1)).max()))
print('  %-13s %-13s %-16s %-12s %s' %
      ('$\\varepsilon/\\varepsilon^*$', '$\\varepsilon$', '$\\Phi(A+E,F)$', '$\\Phi/j_c$', '复核'))
for f in (0.1, 0.3, 0.5, 0.7, 0.9, 0.95, 0.99, 0.999, 1.0):
    ph, tag = phi_hat(A + f * E1, F1)
    print('  %-13.3f %-13.6f %-16s %-12.3f %s'
          % (f, f * e1, '%.6f' % ph if np.isfinite(ph) else 'inf',
             ph / JC if np.isfinite(ph) else float('inf'), tag))

F2 = BESTF[2]
E2, e2, lam2, v2 = witness(F2)
ph2 = phi_hat(A + E2, F2)[0]
print('  $r=2$ 同口径：名义 $\\Phi=%.6f$，$\\varepsilon^*_{\\rm real}=%.6f$（$\\lambda=%.4f$），'
      '$\\varepsilon=1.0\\times\\varepsilon^*$ 处 $\\Phi=$ %s'
      % (phi_iter(F2)[0], e2, lam2, '%.6f' % ph2 if np.isfinite(ph2) else 'inf'))

print('\n  [C] 抽样 vs 最坏：同一株植物、同一个 $\\varepsilon=0.5$，两条口径给出相反的安全感')
ph, tag = phi_hat(A + 0.5 * E1 / np.linalg.norm(E1, 2), F1)
print('  E53[7]（120 次随机方向，$\\varepsilon=0.5$）：$\\Phi_{\\max}/$名义 $=9.95\\times$，失败 $0/120$ —— 读起来"稳"。')
print('  本模块（实见证方向，$\\varepsilon=0.5=%.1f\\times\\varepsilon^*_{\\rm real}$）：$\\Phi=$ %s，%s'
      % (0.5 / e1, '%.6f' % ph if np.isfinite(ph) else 'inf', tag))
print('  两者不矛盾：抽样命不中那一张低维曲面，$\\sup$ 却是对整个球取的。')

print('\n  [D] 交付形态（给 D 的 N4 的改写）')
print('    $\\bar\\Phi(\\varepsilon)$ 只在 $\\varepsilon<\\varepsilon^*_{\\rm real}(A,F)$ 上有定义；')
print('    $\\varepsilon\\ge\\varepsilon^*$ 直接印 $+\\infty$，并附见证 $E=(\\lambda I-A)vv^\\top$（秩一、实、$v\\in\\ker F$）。')
print('    成本：$\\varepsilon^*$ 是 $O(N_\\lambda)$ 次 $8\\times4$ SVD（这里 $N_\\lambda=5600$，实测一个支撑 $<2$ 秒），')
print('    比 120 次蒙特卡洛便宜，而且是证书不是采样。')

print('\n  [E] 发散的形状：$\\Phi$ 到可检测性边界是不是**单极点**（留数、指数）')
print('      若 $\\Phi\\sim\\kappa/(\\varepsilon^*-\\varepsilon)$，那 $\\Phi\\cdot(\\varepsilon^*-\\varepsilon)\\to\\kappa$ 是常数，')
print('      而 $\\log\\Phi$ 对 $\\log(\\varepsilon^*-\\varepsilon)$ 的回归斜率就是指数（应为 $-1$）。')
print('      这条同时是 Prop C/"$\\Phi<\\infty\\iff$ 可检测"的**定量续集**：边界上怎么发散、发散多快。')
F3 = BESTF[3]
E3, e3, _, _ = witness(F3)
for lbl, F, Ex, ex_ in [('精修最优 $r=1$', F1, E1, e1),
                        ('精修最优 $r=2$', F2, E2, e2),
                        ('精修最优 $r=3$', F3, E3, e3)]:
    gap, vv = [], []
    for f in (0.8, 0.9, 0.95, 0.98, 0.99, 0.995, 0.999, 0.9995):
        ph = phi_hat(A + f * Ex, F)[0]
        if np.isfinite(ph) and ph > 0:
            gap.append(ex_ * (1 - f))
            vv.append(ph)
    gap, vv = np.array(gap), np.array(vv)
    m = np.polyfit(np.log(gap), np.log(vv), 1)
    res = vv * gap
    print('  %-16s $N=$%d  斜率 $=%.3f$（应 $-1$）  $\\Phi\\cdot$间距 $=[%.3f,\\ %.3f]$  末两点比 $=%.4f$'
          % (lbl, len(gap), m[0], res[0], res[-1],
             (np.log(vv[-1]) - np.log(vv[-2])) / (np.log(gap[-1]) - np.log(gap[-2]))))
print('  读法：斜率贴 $-1$ 且留数收敛 $\\Rightarrow$ 单极点，$\\bar\\Phi(\\varepsilon)$ 有显式上界')
print('      $\\kappa/(\\varepsilon^*-\\varepsilon)$；这条比"抽样最大值"强，因为它给的是**随半径的函数形式**。')
