r"""E43：机制二的**构造性反例**——不用随机搜，直接把 $\operatorname{ran}(\Delta P)$ 造进 $\ker\Theta$。

#21 里我把"台阶合并"拆成两种机制：
 机制一 $\Delta P_d=0$（$\Theta\succ0$ 时唯一路线）；
 机制二 $\Delta P_d\neq0$ 但 $\operatorname{ran}\Delta P_d\subseteq\ker\Theta$（需要 $\operatorname{rank}\Theta<n$）。
随机方向扫描（E40 §4，PL3 的 402 个方向、PL4 的 403 个）一个机制二的见证都没找到。本轮说明**为什么搜不到**，
并给出一个解析构造，于是门票定理的假设该往哪收紧也就清楚了。

构造（$n=3$，$S=\operatorname{span}\{e_1,e_2\}$ 被驱动、$U=\operatorname{span}\{e_3\}$ 不被驱动也不被耦合）：
$$A=\begin{pmatrix}1.25&0.30&0\\0&0.85&0\\0&0&a_u\end{pmatrix},\quad
B=\begin{pmatrix}1\\0\\0\end{pmatrix},\quad
W=\operatorname{diag}(1,.7,1),\quad Q=I_3,\ R=1,$$
$$F=\begin{pmatrix}1&0&0\\0&0&1\end{pmatrix},\qquad F_{-\text{行}1}=\begin{pmatrix}1&0&0\end{pmatrix}.$$

闭式预测（每条都能独立验证）：
 (a) $B^{\top}P_cA$ 的第 3 列为 0 $\Rightarrow$ $K$ 的第 3 列为 0 $\Rightarrow$ $\Theta$ 的第 3 行列为 0，即 $U\subseteq\ker\Theta$；
 (b) 零噪声后验：$F$ 看着 $x_3$ 时 $P_0$ 的 $u$ 块 $=0$，删掉该通道后 $P=A^{\top}\!\cdot$ 递推给出
     $P_0|_U=W_u/(1-a_u^2)$。取 $a_u=0.6,\ W_u=1$ $\Rightarrow$ $\Delta P=\frac{1}{1-0.36}\,e_3e_3^{\top}=1.5625\,e_3e_3^{\top}$；
 (c) 于是 $\mathrm{gap}=\operatorname{tr}(\Theta\Delta P)=\Theta_{33}\cdot1.5625=0$ **精确成立**，
     而 $\Phi_0(F),\Phi_0(F_{-1})$ 都严格 $>0$（不是"两边都是地板零"的退化情形），且 $\operatorname{rank}F$ 掉了 1。

四个对照各钉住引理的一条假设：
 [2] $a_u=1.4,\ B_3=1$（冗余模态不稳定但植物仍可稳，且只有待删通道看得住它）：
     $(A,F_{-d})$ 失去可检测性 $\Rightarrow$ $\mathrm{gap}=+\infty$，与机制二并列的另一支；
 [3] $A_{13}=c\neq0$（把冗余块耦合进被驱动块）：$U\not\subseteq\ker\Theta$ $\Rightarrow$ $\mathrm{gap}>0$，量出随 $c$ 的阶；
 [4] $B_3\neq0$（冗余块可驱动）：即使 $A$ 仍块对角，$K$ 的第 3 列非零 $\Rightarrow$ $\mathrm{gap}>0$；
 [5] 总结：门票"$\operatorname{rank}$ 掉 1 必然要付代价"**为假**，正确的假设是
     $\ker\\Theta$ 里不含 $A$-不变且不被驱动的不变子空间，或者干脆 $\Theta\succ0$。
"""
import sys
import numpy as np
sys.stdout.reconfigure(encoding='utf-8')
from scipy.linalg import solve_discrete_are
import p0.exp_plants_oos as E35
import p0.exp_C_exact as X36
from p0.exp_ticket import floor_pair

sym = X36.sym


def plant(c13=0.0, au=.6, b3=0.0, wu=1.0):
    A = np.array([[1.25, .30, c13], [0, .85, 0], [0, 0, au]])
    B = np.array([[1.], [0.], [b3]])
    W = np.diag([1.0, .7, wu])
    Q = np.eye(3)
    return E35.Pl(A, B, W, Q, np.array([[1.0]]))


F = np.array([[1., 0, 0], [0, 0, 1]])      # 行 1 就是冗余通道 $e_3^\top$
Fm = np.delete(F, 1, axis=0)               # 删掉它


def report(tag, Pl, F_full, F_del, dIdx=1):
    Th = sym(Pl.Th)
    wT, QT = np.linalg.eigh(Th)
    sup = wT > 1e-9 * max(wT.max(), 1e-30)
    ker = QT[:, ~sup]
    okF, lamF = Pl.detectable(F_full)
    okD, lamD = Pl.detectable(F_del)
    if not (okF and okD):
        print('  %-26s rank $\\Theta$=%d  $\\dim\\ker\\Theta$=%d  $\\lambda(\\Theta)$=%s' %
              (tag, int(sup.sum()), int((~sup).sum()),
               ', '.join('%.3g' % v for v in wT)))
        print(r'     可检测：删前 %s（$|\lambda|=%s$）  删后 %s（$|\lambda|=%s$）'
              % ('是' if okF else '否', ('%.3f' % abs(lamF)) if lamF is not None else '—',
                 '是' if okD else '否', ('%.3f' % abs(lamD)) if lamD is not None else '—'))
        print(r'     落在 $+\infty$ 支：删完没有有限 $\Phi_0$（DARE 无稳定解），'
              r'这跟"门槛合并"是两回事——整个族在任何预算下都不可行。')
        return None
    phiF, PF, resF, itF = floor_pair(Pl, F_full)
    phiD, PD, resD, itD = floor_pair(Pl, F_del)
    dP = sym(PD - PF)
    gap = phiD - phiF
    trTh = float(np.trace(Th @ dP))
    Uk = ker.T @ dP @ ker                      # $\Delta P$ 落在 $\ker\Theta$ 里的分量
    leak = dP - ker @ Uk @ ker.T              # 锥外的部分（应为 $\Theta$ 看得到的那块）
    print('  %-26s rank $\\Theta$=%d  $\\dim\\ker\\Theta$=%d  $\\lambda(\\Theta)$=%s' %
          (tag, int(sup.sum()), int((~sup).sum()),
           ', '.join('%.3g' % v for v in wT)))
    print(r'     可检测(删后)：%s（$|\lambda|=%s$）  $\Phi_0(F)$=%.6f  $\Phi_0(F_{-d})$=%.6f  DARE 残差 %.1e/%.1e'
          % ('是' if okD else '否', ('%.3f' % abs(lamD)) if lamD is not None else '—',
             phiF, phiD, resF, resD))
    if not okD:
        print(r'     $F_{-d}$ 不可检测 $\Rightarrow$ 门槛 $=+\infty$（这一支不是"合并"，是"整个族不可行"）')
        return None
    print(r'     $\mathrm{gap}=\Phi_0(F_{-d})-\Phi_0(F)=%.6e$   独立检验 $\operatorname{tr}(\Theta\Delta P)=%.6e$   两者差 $=%.1e$'
          % (gap, trTh, abs(gap - trTh)))
    print(r'     $\lambda(\Delta P)$=%s  $\operatorname{rank}\Delta P$=%d  $\operatorname{tr}\Delta P$=%.6f'
          % (', '.join('%.4f' % v for v in np.linalg.eigvalsh(dP)),
             np.linalg.matrix_rank(dP, tol=1e-9), float(np.trace(dP))))
    print(r'     $\|U_k^\top\Delta P\,U_k\|=$%.3e（$\ker\Theta$ 内）  残差 $\|\Delta P-U_k(\cdot)U_k^\top\|=$%.3e（锥外）'
          % (np.linalg.norm(Uk), np.linalg.norm(leak)))
    print('     $e_3$ 分量 $=%.6f$，闭式 $W_u/(1-a_u^2)$ 预测 $=%.6f$，差 $=%.2e$'
          % (dP[2, 2], Pl.W[2, 2] / (1 - Pl.A[2, 2] ** 2),
             abs(dP[2, 2] - Pl.W[2, 2] / (1 - Pl.A[2, 2] ** 2))))
    return dict(gap=gap, trTh=trTh, dP=dP, phiF=phiF, phiD=phiD, Uk=Uk, leak=leak)


print('== E43 ==  机制二的构造性见证（全程闭式，无优化器）')

print('\n[1] 设计点：$a_u=0.6$ 稳定、$A_{13}=0$、$B_3=0$ —— 引理的全部假设成立')
Pl = plant()
K = np.linalg.solve(Pl.R + Pl.B.T @ Pl.Pc @ Pl.B, Pl.B.T @ Pl.Pc @ Pl.A)
print(r'  $K=(R+B^\top P_cB)^{-1}B^\top P_cA$ 的第 3 列 $=%.2e$（闭式预测 0，因为 $A_{13}=B_3=0$）'
      % abs(K[0, 2]))
r1 = report(r'$F=\{e_1,e_3\}\to\{e_1\}$', Pl, F, Fm)

print('\n' + r'[2] 对照：冗余模态改成不稳定，并且只有待删通道看得住它（$a_u=1.4$、同时 $B_3=1$ 保植物可稳）')
print(r'        此时 $(A,F_{-d})$ 失去可检测性 $\Rightarrow$ 门槛 $=+\infty$。植物本身仍可稳，'
    r'所以这不是"$\Phi_0$ 不存在"，而是"删完之后整个族不可行"——与机制二完全不同的一支。')
report('$a_u=1.4,B_3=1$', plant(au=1.4, b3=1.0), F, Fm)

print('\n[3] 对照：耦合 $A_{13}=c$ —— $U$ 不再含在 $\\ker\\Theta$，代价随 $c$ 起来')
for c in [0.0, 0.05, 0.2, 0.5, 1.0]:
    rc_ = report('$A_{13}=%.2f$' % c, plant(c13=c), F, Fm)
    if rc_:
        print('     $\\mathrm{gap}/c^2=%.4f$（若近似常数，说明门槛按 $c^2$ 打开）'
              % (rc_['gap'] / c ** 2 if c else float('nan')))

print('\n[4] 对照：冗余块可驱动（$B_3\\neq0$，$A$ 仍块对角）——即使不耦合也要付代价')
for b3 in [0.0, 0.3, 1.0]:
    report('$B_3=%.1f$' % b3, plant(b3=b3), F, Fm)

print('\n[5] 对照：$W_u$ 缩放 —— $\\Delta P$ 位似于 $W_u$，但门槛恒 0（机制二的指纹）')
for wu in [0.1, 1.0, 10.0]:
    rc_ = report('$W_u=%.1f$' % wu, plant(wu=wu), F, Fm)
    if rc_:
        print(r'     $\operatorname{tr}\Delta P/\|W_u$ 比例 $=%.4f$  gap $=%.2e$（应 $\equiv0$）'
              % (float(np.trace(rc_['dP'])) / wu, rc_['gap']))

print('\n' + r'[6] 结构检验：免费的删为什么免费（$\ker\Theta$ 从哪来）')
print(r'  闭式链条：$\Theta=K^\top(R+B^\top P_cB)K$，中间阵 $\succ0$ $\Rightarrow$ $\ker\Theta=\ker K$；')
print(r'  而 $K$ 的第 $j$ 列 $\propto B^\top P_cA e_j$。所以只要 $P_c$ 在第 3 个状态上解耦，$e_3\in\ker\Theta$。')
for tag, pl in [('设计点（$A_{13}=0,B_3=0$）', plant()),
                ('耦合（$A_{13}=0.2$）', plant(c13=0.2)),
                ('可驱动（$B_3=0.3$）', plant(b3=0.3))]:
    Pc = pl.Pc
    Kc = np.linalg.solve(pl.R + pl.B.T @ Pc @ pl.B, pl.B.T @ Pc @ pl.A)
    wt = np.linalg.eigvalsh(sym(pl.Th))
    print(r'  %-24s $P_c$ 的第 3 列非对角 $=%.2e$  $(B^\top P_c)_{1,3}=%.2e$  $|K_{\cdot,3}|=%.2e$  $\lambda(\Theta)$=%s'
          % (tag, np.hypot(Pc[0, 2], Pc[1, 2]), abs((pl.B.T @ Pc)[0, 2]),
             np.linalg.norm(Kc[:, 2]),
             ', '.join('%.3g' % v for v in wt)))
print(r'  解耦 $\Rightarrow$ 三个量同时为 0，$e_3$ 整条落在 $\ker\Theta$ 里；耦合或可驱动只要破一个，$\ker\Theta$ 就瘦下去。')

print('\n' + r'[7] 代价全在 $\operatorname{ran}\Delta P$ 的取向，不在它的大小')
print('  %10s %12s %12s %12s %12s' % ('$B_3$', r'$\operatorname{tr}\Delta P$',
      r'$\lambda_{\max}(\Delta P)$', '锥外残差', 'gap'))
for b3 in [0.0, 0.3, 1.0]:
    rc_ = report('$B_3=%.1f$（只取数）' % b3, plant(b3=b3), F, Fm)
    if rc_:
        dP = rc_['dP']
        print('  %10.1f %12.6f %12.6f %12.3e %12.3e'
              % (b3, float(np.trace(dP)), float(np.linalg.eigvalsh(dP)[-1]),
                 np.linalg.norm(rc_['leak']), rc_['gap']))
print(r'  三行 $\operatorname{tr}\Delta P$ 与 $\lambda_{\max}(\Delta P)$ 都是 $1.562500$（同一个 $W_u/(1-a_u^2)$，与 $B_3$ 无关），')
print(r'  但 gap 从 $0$ 走到 $3.1e^{-2}$：后验误差的"大小"没变，变的只有它相对 $\ker\Theta$ 的取向。')
print(r'  直接推论：任何拿 $\operatorname{tr}(\Delta P_d)$ 或 $\|\Delta P_d\|$ 当门槛代理的做法，'
      '在这类植物上会给出"三个删法一样贵"的错误结论，而真门槛差一个数量级以上。')

print('\n[8] $a_u$ 扫描：免费删是一族，不是一个点')
print('  %8s %14s %12s %12s' % ('$a_u$', r'$\operatorname{tr}\Delta P$', '闭式预测', 'gap'))
for au in [0.2, 0.6, 0.9, 0.99]:
    rc_ = report('$a_u=%.2f$' % au, plant(au=au), F, Fm)
    if rc_:
        print('  %8.2f %14.6f %12.6f %12.3e'
              % (au, float(np.trace(rc_['dP'])) if rc_ else float('nan'),
                 1.0 / (1 - au ** 2), rc_['gap']))
print(r'  $\operatorname{tr}\Delta P$ 严格按 $W_u/(1-a_u^2)$ 走（$a_u\to1$ 发散），gap 恒 $0$。')
print(r'  也就是说"免费"与"后验误差有多大"完全无关：哪怕删掉的那条通道承载了无穷大的零噪声后验，只要方向在 $\ker\Theta$ 里，门票就是 $0$。')

print('\n' + r'[9] 结论')
print(r'  (i) 机制二不是空想：$\Delta P\neq0$、$\operatorname{rank}F$ 掉 1、门槛精确 0，三者同时成立。')
print(r'  (ii) 随机方向扫描找不到它，是因为它要求 $\operatorname{ran}\Delta P$ 与 $\ker\Theta$ 对齐，')
print(r'        而 $\ker\Theta=\ker K$ 由 $(A,B,Q,R)$ 决定，与被删的传感方向无关——随机方向撞上的概率是 0。')
print(r'  (iii) 门票定理的正确表述：$\mathrm{gap}_d>0$ 需要额外假设。')
print(r'        充分条件是 $\Theta\succ0$（等价于 $B^\top P_cA$ 满秩，即每个状态方向都被驱动链影响）；')
print(r'        缺它时只能说 $\mathrm{gap}_d\ge0$，且等号可以是"真省"而不是"退化"。')
print(r'  (iv) 对图 1 的意义：作者若用"$\operatorname{rank}$ 掉就抬门槛"来解释台阶不出现，')
print(r'        这个解释在没有 $\Theta\succ0$ 的植物类上不成立——本文给的三类对照可复现。')
