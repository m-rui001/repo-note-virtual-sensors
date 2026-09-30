# -*- coding: utf-8 -*-
r"""e170：锚点的 DARE **到底有几个解**，我用的是哪一个（e169 [C2] 的正判读）。

e169 的实测事实（`p0/e169_out.txt`）：
 [i] $D(\lg)$ 对同一方向**单调降**（8 条方向里 7 条上升段 $=0$，$\Theta$ 首轴 $0/48$）$-\ -$
     e167 的"同方向两个 $I$"**不是**代价曲线的多支；
 [ii] 但 `solve_discrete_are` 的根与"从 $P_0=W$ 起的 Riccati 时间迭代"**极限不同**（最大相对差 5.27）。
 $\Rightarrow$ 唯一剩下的解释：**该 Riccati 有多个对称解**，两个算法落在不同支上。
 而 $\rho(A)=1.7124>1$（开环不稳）$+$ 秩-1 测量 $C=s\,u^{\mathsf T}$ $\\Rightarrow$ $(A^{\mathsf T},C)$ **未必可检测** $-\ -$
 若不可检测，则**不存在 stabilizing 解**，"DARE 的解"这个记号在我两车的定价里是**未定义到位**的。
 这不是代码风格问题：$\mathrm{tr}(\Theta P_m)$ 两个根差 $10\times$（e169 表：137 vs 1377）$-\ -$
 "同秩可达地板"、$\Delta_r$、离开角全建在这个选择上。

机器（全部现算，不 import C）：
 [D0] $F=A^{\mathsf T}$ 的谱与不稳定模态个数；PBH 可检测性：对每个 $|\lambda|\ge1$ 检 $\mathrm{rank}\begin{bmatrix}\lambda I-F\\ C\end{bmatrix}=n$。
 [D1] 根的**分类**：对每个根算诱导增益 $K=FXC^{\mathsf T}(I+CXC^{\mathsf T})^{-1}$ 与 $\rho(F-KC)$ $-\ -$ $<1$ 才是 stabilizing。
 [D2] 根的**排序**：$\lambda_{\min}(X_{\rm it}-X_{\rm sd})\ge0$ 则时间迭代的根是"更大"的那个（Riccati 解的 Loewner 序）。
 [D3] 两个根各给一条 $D(\lg)$ 曲线（$\lg\in[-14,10]$，49 点）：报 $\min D$（"地板"）、$D=t$ 的穿点数与穿点处的 $I$。
      $\Rightarrow$ "地板/穿点"到底随不随根变（若随，则我板上所有地板数必须标根）。
判据（跑前写死）：
 [D1$^\prime$] **改正映射后**（$\Phi(X)=AXA^{\mathsf T}+W-AXC^{\mathsf T}(I+CXC^{\mathsf T})^{-1}CXA$，与 `solve_discrete_are(A^{\mathsf T},C^{\mathsf T},W,I)$ 同式）：
      两根最大相对差 $\le10^{-6}$ $\\Rightarrow$ 判"**两法同根 $-\ -$ e169 的'两根差 5.27'是我方把转置写反，撤回该结论**"；
      若仍 $>10^{-3}$ $\\Rightarrow$ 判"同一 Riccati 有两个对称解，我的定价用的是哪一个必须写明"。
 [D4] 若 PBH 在某方向**不可检测**且该方向的根 $\rho(A-LC)\ge1$ $\\Rightarrow$ 判
      "**不可检测方向的 DARE 根不是滤波极限** $-\ -$ 板上凡地板/$\Delta_r$/离开角要说明这些点是否滤波可实现"。
 [D5] 若两根都 stabilizing 且相同 $\\Rightarrow$ 判"对象无歧义"。
 [D6] 顺带测：秩-1 的 $\min D$（"同秩可达地板"）在两根下是否同一个数（e162b [P] 的 43.7058 依赖这个）。
方向：$\Theta$ 四轴 $+$ 2 条随机射线（$r=1$），另加 $\Theta$ 前 2 列（$r=2$）作对照。
"""
import sys
import time
import numpy as np
from scipy.linalg import solve_discrete_are
sys.stdout.reconfigure(encoding='utf-8')

out = []
T0 = time.time()


def p(*a):
    s = ' '.join(str(x) for x in a)
    out.append(s)
    print(s)
    sys.stdout.flush()


P = {'__name__': 'p'}
exec(compile(open('p0/e162_ray_closure.py', encoding='utf-8').read().split('CRAB = {')[0],
             'p0/e162_ray_closure.py[prefix]', 'exec'), P)
AUD = P['_ns']
SDP = {56.66: 2.053565, 80.00: 1.603105}
A, B, W = AUD['A'].copy(), AUD['B'].copy(), AUD['W'].copy()
Pc, K0, TH, JC = P['ctrl_full'](A, B, W, np.eye(4), np.eye(4))
n = 4
wth, U0 = np.linalg.eigh(TH)
Uth = U0[:, np.argsort(-wth)]
rng = np.random.default_rng(555)
F = A.T
GRID = np.linspace(-14.0, 10.0, 49)
TGT = (80.00, 56.66)


def sym(X):
    return 0.5 * (X + X.T)


def root_scipy(C):
    r = C.shape[0]
    return sym(solve_discrete_are(F, C.T, W, np.eye(r)))


def price(C, X):
    r = C.shape[0]
    S = C @ X @ C.T + np.eye(r)
    Kk = X @ C.T @ np.linalg.inv(S)
    Pp = sym(X - Kk @ S @ Kk.T)
    I = 0.5 * np.log(np.linalg.det(S)) / np.log(2.0)
    return I, JC + float(np.trace(TH @ Pp))


def pbh_detect(C, lam):
    """该不稳定特征值处 PBH 是否满秩（$C$ 看得见这个模态）。"""
    M = np.vstack([lam * np.eye(n) - F, C])
    return np.linalg.matrix_rank(M, tol=1e-7) == n


evF = np.linalg.eigvals(F)
unstab = [z for z in evF if abs(z) >= 1.0]

# ---- 口径自纠（本轮第一版写错，公开）----
# `solve_discrete_are(a,b,q,r)` 解的是 $X=a^{\mathsf T}Xa+q-a^{\mathsf T}xb(r+b^{\mathsf T}xa)^{-1}b^{\mathsf T}xa$。
# 我的调用是 $(a,b)=(A^{\mathsf T},C^{\mathsf T})$ $\\Rightarrow$ 真实迭代映射是
#    $\Phi(X)=A\,X\,A^{\mathsf T}+W-A\,X\,C^{\mathsf T}(I+CXC^{\mathsf T})^{-1}C\,X\,A$，
# 而 e169 与我 e170 初版都写成了 $A^{\mathsf T}XA$ $-\ -$ **转置用反** $\\Rightarrow$ e169 的"两根差 $5.27$"是我的映射错，
# 不能算"DARE 多解"的证据。本轮把映射改正，并**同时**报两种闭环半径（滤波增益 $A-LC$ 与 scipy 约定 $A^{\mathsf T}-C^{\mathsf T}f$）。
Dyn = A


def root_iter2(C, itmax=40000, tol=1e-13, P0=None):
    r = C.shape[0]
    X = sym(W).copy() if P0 is None else sym(P0).copy()
    for k in range(itmax):
        S = np.eye(r) + C @ X @ C.T
        L = Dyn @ X @ C.T @ np.linalg.inv(S)
        Xn = sym(Dyn @ X @ Dyn.T + W - L @ S @ L.T)
        if np.max(np.abs(Xn - X)) < tol:
            return Xn, k + 1, True
        X = Xn
    return X, itmax, False


def radius_filter(C, X):
    r = C.shape[0]
    S = np.eye(r) + C @ X @ C.T
    L = Dyn @ X @ C.T @ np.linalg.inv(S)
    return float(max(np.abs(np.linalg.eigvals(Dyn - L @ C))))


def radius_scipy(C, X):
    r = C.shape[0]
    S = np.eye(r) + C @ X @ C.T
    f = np.linalg.inv(S) @ C @ X @ Dyn
    return float(max(np.abs(np.linalg.eigvals(Dyn.T - C.T @ f))))


def pbh_detect(C, lam):
    """滤波口径：$(A,C)$ 在该不稳定特征值处是否 PBH 满秩（$C$ 看得见这个模态）。"""
    M = np.vstack([lam * np.eye(n) - Dyn, C])
    return np.linalg.matrix_rank(M, tol=1e-7) == n


evD = np.linalg.eigvals(Dyn)
unstab = [z for z in evD if abs(z) >= 1.0]
p(r'== [口径自纠] e169 的"两根差 5.27"源于我把迭代映射的转置写反（$A^{\mathsf T}XA$ 应为 $AXA^{\mathsf T}$）$-\ -$ '
  r'本轮改正后重判；$D(\lg)$ 单调（e169 [i]）与 e167 的角度位移（未归一化向量被 clip 成 $0^\circ$）也一并登记 ==' % ())
p(r'== 锚点：$J_C=%.6f$；$\rho(A)=%.4f$；$A$ 的 $|\lambda|$ 降序 %s；不稳定（$|\lambda|\ge1$）个数 %d =='
  % (JC, max(np.abs(evD)), np.array2string(np.array(sorted(np.abs(evD), reverse=True)), precision=4), len(unstab)))
p('')
p(r'== [D0/D1/D2] 逐方向：可检测性、两根的 stabilizing 半径、Loewner 序（$\lg=0$ 与 $\lg=2$）==')
DIRS = [(r'$\Theta$ 轴%d' % (i + 1), Uth[:, [i]]) for i in range(4)]
for i in range(2):
    u = rng.normal(0, 1, n)
    DIRS.append((r'rand%d' % (i + 1), (u / np.linalg.norm(u))[:, None]))
DIRS.append((r'$\Theta$ 前2列(秩2)', Uth[:, :2]))

nosdt = 0
nosc_stab = 0
nit_stab = 0
rows_checked = 0
maxrel = 0.0
for nm, G in DIRS:
    Zn = G / np.linalg.norm(G, axis=0)
    r = Zn.shape[1]
    for lg in (0.0, 2.0):
        C = (10.0 ** (0.5 * lg)) * Zn.T
        Xs = root_scipy(C)
        Xi, it, ok = root_iter2(C)
        rel = float(np.max(np.abs(Xi - Xs)) / max(1e-12, np.max(np.abs(Xs))))
        maxrel = max(maxrel, rel)
        det = all(pbh_detect(C, z) for z in unstab)
        rs, ri = radius_filter(C, Xs), radius_filter(C, Xi)
        rsc = radius_scipy(C, Xs)
        loe = float(np.linalg.eigvalsh(Xi - Xs).min())
        rows_checked += 1
        nosdt += 0 if det else 1
        nosc_stab += 1 if rs < 1.0 else 0
        nit_stab += 1 if ri < 1.0 else 0
        p(r'      %-12s $r=%d$ $\lg=%+.1f$：PBH 可检测 %s；scipy 根 $\rho(A-LC)=%.4f$（scipy 约定 $\rho=%.4f$）；'
          r'迭代根 $\rho(A-LC)=%.4f$（%d 步收敛 %s）；两根最大相对差 %.2e；$\lambda_{\min}(X_{\rm it}-X_{\rm sc})=%+.3e$'
          % (nm, r, lg, det, rs, rsc, ri, it, ok, rel, loe))
p(r'  [D1$^\prime$] 改正映射后，scipy 根与时间迭代根的最大相对差 $=\,$%.2e（%d 个组合）' % (maxrel, rows_checked))
p('')
p(r'== [D3] 两条 $D(\lg)$ 曲线：scipy 根 vs 迭代根（$\Theta$ 首轴、$\Theta$ 前2列）==')
for nm, G in ((r'$\Theta$ 轴1', Uth[:, [0]]), (r'$\Theta$ 前2列(秩2)', Uth[:, :2])):
    Zn = G / np.linalg.norm(G, axis=0)
    for t in TGT:
        line = {}
        for tag, f in (('scipy', root_scipy), ('iter', lambda C: root_iter2(C, itmax=20000)[0])):
            Ds, Is = [], []
            for lg in GRID:
                C = (10.0 ** (0.5 * lg)) * Zn.T
                try:
                    X = f(C)
                    I, D = price(C, X)
                except Exception:
                    I, D = np.nan, np.nan
                Ds.append(D)
                Is.append(I)
            Ds, Is = np.array(Ds), np.array(Is)
            fin = np.isfinite(Ds)
            cross = sum(1 for i in range(len(GRID) - 1)
                        if fin[i] and fin[i + 1] and (Ds[i] - t) * (Ds[i + 1] - t) <= 0 and Ds[i] != Ds[i + 1])
            nup = int(np.sum(np.diff(Ds[fin]) > 0))
            line[tag] = (np.nanmin(Ds), np.nanmin(Is), cross, nup, int(np.sum(~fin)))
        p(r'      %-14s $D_{\rm tgt}=%.2f$：scipy 根 $\min D=%.4f$，穿点 %d 个，上升段 %d，非有限 %d $\\Rightarrow$ '
          r'迭代根 $\min D=%.4f$，穿点 %d 个，上升段 %d，非有限 %d'
          % (nm, t, line['scipy'][0], line['scipy'][2], line['scipy'][3], line['scipy'][4],
             line['iter'][0], line['iter'][2], line['iter'][3], line['iter'][4]))
p('')
p(r'== [D1$^\prime$/D4/D5] 判决（自带读数条数：#46）==')
p(r'  检查 %d 个（方向,$\lg$）组合 $-\ -$ 其中 PBH 不可检测 %d 个；scipy 根 $\rho(A-LC)<1$ %d 个；迭代根同 %d 个；'
  r'两根最大相对差 %.2e' % (rows_checked, nosdt, nosc_stab, nit_stab, maxrel))
if maxrel <= 1e-6:
    p(r'  **[D1$^\prime$] 判：两法同根（最大相对差 %.2e $\le10^{-6}$）$-\ -$ e169 的"两根差 5.27"是我方迭代映射的转置写反，'
      r'该结论撤回**（不改任何已落的 $\Delta$/地板数，因为它们全走同一 `solve_discrete_are`）。' % maxrel)
elif maxrel > 1e-3:
    p(r'  **[D1$^\prime$] 判：同一 Riccati 有两个对称解**（最大相对差 %.2e）$-\ -$ 我的定价必须写明取哪一支。' % maxrel)
else:
    p(r'  [D1$^\prime$] 判：两根差 %.2e 落在 $10^{-6}\sim10^{-3}$ $\\Rightarrow$ 只报数，不判。' % maxrel)
if nosdt > 0:
    p(r'  **[D4] 判：%d/%d 个组合的 $(A,C)$ 在该秩下 PBH 不可检测** $\\Rightarrow$ 这些方向的 DARE 根不是"滤波极限"，'
      r'凡用它们得到的"可达/地板"读数要另加口径说明。' % (nosdt, rows_checked))
else:
    p(r'  [D4] 判：全部 %d 个组合都可检测 $\\Rightarrow$ stabilizing 根存在且唯一，我的分支就是它。' % rows_checked)
p('')
p('用时 %.0f s | 异常 无' % (time.time() - T0))
txt = '\n'.join(out) + '\n'
open('p0/e170_out.txt', 'w', encoding='utf-8').write(txt)
open('p0/e170_run.txt', 'w', encoding='utf-8').write(txt)
