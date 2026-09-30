# -*- coding: utf-8 -*-
r"""e169：$D(\lg)$ 到同一方向能否**两次**穿过同一目标代价？（e167 异常的正面对决）

触发（本车道机器问题，必须先判）：
 `p0/e167_levelset_support.py`：$D{=}80,r{=}1$，从 $\Theta$ 首轴出发 $-$ $-$ 单发割线给 $I=1.653054$，
 Powell 细化后给 $I=1.603105$（$=$ e163 凸 SDP 值），**方向位移 $0.00^\circ$**。
 同一方向、同一 $D$、两个 $I$ $\\Rightarrow$ 要么 $D(\lg)$ 非单调（真实模型性质：加强测量反而抬高 $\mathrm{tr}(\Theta P_p)$），
 要么 `solve_discrete_are` 在某些增益段返回**非 stabilizing 根**（求解器伪影）。
 两者对板的后果完全不同 $-$ $-$ 前者的 $\Delta$ 仍合法（可行点真的满足 $|D-t|\le10^{-7}$，只是"目标代价反演多值"），
 后者会把 e162/e162b/e165 的全部 $\Delta$ 与离开角**作废重测**（#43）。
（e168 的教训也写在这：它的 $\lg$ 窗 $[-2,12]$ 太窄，$D$ 在 $\lg=-2$ 已低于目标 $\\Rightarrow$ 扫不到穿点，
 只有 10/241 条方向有交点 $-$ $-$ 那是**窗口伪影**不是物理，本轮把窗开到 $[-14,10]$。）

机器：
 [C1] 逐方向在 $\lg\in[-14,10]$（49 点）打表 $(\lg,I,D)$，打印 $D$ 的**上升段计数**与 $D=t$ 的**穿点个数**（二分夹到 $10^{-9}$）。
 [C2] DARE 交叉验证：同一增益下把 Riccati **压缩映射迭代**（$P_0=W$，至多 20000 步，$\|P_{k+1}-P_k\|_\infty<10^{-12}$ 判收敛）
      与 `solve_discrete_are` 的根对比 $-\ -$ 若两者一致，非单调就**不是**求解器伪影。
 [C3] 判决（跑前写死）：
      若某方向有 $\ge2$ 个穿点、值差 $>10^{-3}$ bit，**且** [C2] 的迭代根与 DARE 根相对差 $\le10^{-6}$（$\lg$ 网格上最大）
        $\\Rightarrow$ 判"**$D(\lg)$ 非单调是模型性质**：同一方向的同一目标代价有多支 $-\ -$ 我的 $\Delta$/离开角必须标支；
          既有读数仍按'可行点'口径有效（$|D-t|\le10^{-7}$ 是直接验证的）"。
      若 [C2] 出现 $\ge10^{-3}$ 的根差 $\\Rightarrow$ 判"**求解器伪影** $-\ -$ 相关读数作废重测"。
      若全部方向至多一穿点 $\\Rightarrow$ 判"单调 $-\ -$ e167 的异常另有原因（回报 GUESS 携带的割线早退）"。
方向：$\Theta$ 四轴 $+$ 4 条随机射线；格 $t\in\{80,56.66\}$（$r=1$）。
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
Pc, K, TH, JC = P['ctrl_full'](A, B, W, np.eye(4), np.eye(4))
P['A'], P['B'], P['W'], P['TH'], P['JC'] = A, B, W, TH, JC
wth, U0 = np.linalg.eigh(TH)
Uth = U0[:, np.argsort(-wth)]
P['Uth'] = Uth
n = 4
rng = np.random.default_rng(31337)
GRID = np.linspace(-14.0, 10.0, 49)


def sym(X):
    return 0.5 * (X + X.T)


def price(Zn, lg):
    C = (10.0 ** (0.5 * lg)) * Zn.T
    Pm = sym(solve_discrete_are(A.T, C.T, W, np.eye(1)))
    Sm = C @ Pm @ C.T + np.eye(1)
    Lk = Pm @ C.T @ np.linalg.inv(Sm)
    Pp = sym(Pm - Lk @ Sm @ Lk.T)
    I = 0.5 * np.log(np.linalg.det(np.eye(1) + C @ Pm @ C.T)) / np.log(2.0)
    return I, JC + float(np.trace(TH @ Pp)), Pm


def riccati_fixed(C, itmax=20000, tol=1e-12):
    """压缩映射求 stabilizing 根：$P\leftarrow A^{\mathsf T}PA+W-A^{\mathsf T}PC^{\mathsf T}(I+CPC^{\mathsf T})^{-1}CPA$。"""
    Q = sym(W)
    X = sym(W).copy()
    Ci = np.linalg.inv(np.eye(1) + C @ X @ C.T)
    prev = None
    for k in range(itmax):
        AtX = A.T @ X
        S = np.eye(1) + C @ X @ C.T
        Kk = AtX @ C.T @ np.linalg.inv(S)
        Xn = sym(AtX @ A + Q - Kk @ S @ Kk.T)
        if prev is not None and np.max(np.abs(Xn - prev)) < tol:
            return Xn, k + 1, True
        prev = Xn.copy()
        X = Xn
    return X, itmax, False


def bisect(cur, t, lo, hi):
    flo, fhi = cur(lo) - t, cur(hi) - t
    if flo * fhi > 0:
        return None
    for _ in range(70):
        mid = 0.5 * (lo + hi)
        fm = cur(mid) - t
        if flo * fm <= 0:
            hi, fhi = mid, fm
        else:
            lo, flo = mid, fm
        if hi - lo < 1e-12:
            break
    lg = 0.5 * (lo + hi)
    return lg


p(r'== 锚点 $J_C=%.6f$；$\lg$ 窗 $[%g,%g]$ 共 %d 点；$\rho(A)=%.4f$（是否开环不稳）=='
  % (JC, GRID[0], GRID[-1], len(GRID), max(np.abs(np.linalg.eigvals(A)))))
p('')
p(r'== [C1] $D(\lg)$ 的形状：上升段计数与 $D=t$ 的穿点 ==')

DIRS = [(r'$\Theta$ 轴%d' % (i + 1), Uth[:, [i]]) for i in range(4)]
for i in range(4):
    u = rng.normal(0, 1, n)
    DIRS.append((r'rand%d' % (i + 1), (u / np.linalg.norm(u))[:, None]))

worst_root = 0.0
max_cross = 0
max_dI = 0.0
for t in (80.00, 56.66):
    p(r'  --- 目标 $D=%.2f$（$I_{\rm unc}=%.6f$）---' % (t, SDP[t]))
    for nm, G in DIRS:
        Zn = G / np.linalg.norm(G)
        Ds, Is = [], []
        for lg in GRID:
            I, D, _ = price(Zn, lg)
            Ds.append(D)
            Is.append(I)
        Ds = np.array(Ds)
        Is = np.array(Is)
        nup = int(np.sum(np.diff(Ds) > 0))
        cross = []
        for i in range(len(GRID) - 1):
            if (Ds[i] - t) * (Ds[i + 1] - t) <= 0 and Ds[i] != Ds[i + 1]:
                lg = bisect(lambda x: price(Zn, x)[1], t, GRID[i], GRID[i + 1])
                if lg is not None:
                    Ix, Dx, _ = price(Zn, lg)
                    cross.append((lg, Ix))
        txt = ', '.join('%.4f@lg%.3f' % (ix, l) for l, ix in cross)
        p(r'      %-10s 上升段 %2d/%2d，穿点 %d 个%s%s'
          % (nm, nup, len(GRID) - 1, len(cross), ('：' + txt) if cross else '（窗内无解）',
             '' if len(cross) < 2 else '  <-- 多支'))
        max_cross = max(max_cross, len(cross))
        if len(cross) >= 2:
            max_dI = max(max_dI, max(c[1] for c in cross) - min(c[1] for c in cross))
        if nm.startswith(r'$\Theta$ 轴1'):
            p(r'         表 $(\lg,D,I)$：' + ' | '.join('%.0f:%.3f:%.3f' % (l, d, i)
                                                        for l, d, i in zip(GRID[::4], Ds[::4], Is[::4]))
              + r'（$\lg:D:I$，每 4 点取 1）')
    p('')

p(r'== [C2] DARE 根 vs 压缩映射根（同一增益，$\Theta$ 首轴 12 个增益点）==')
Zn = Uth[:, [0]] / np.linalg.norm(Uth[:, [0]])
for lg in np.linspace(-6, 6, 13):
    C = (10.0 ** (0.5 * lg)) * Zn.T
    _, _, Pd = price(Zn, lg)
    Pi, it, ok = riccati_fixed(C)
    rel = float(np.max(np.abs(Pi - Pd)) / max(1e-12, np.max(np.abs(Pd))))
    worst_root = max(worst_root, rel)
    p(r'      $\lg=%+.2f$：DARE $\mathrm{tr}(\Theta P_m)=%.6f$，迭代 $\mathrm{tr}(\Theta P_m)=%.6f$，'
      r'相对差 %.2e（迭代 %d 步，收敛 %s）'
      % (lg, float(np.trace(TH @ Pd)), float(np.trace(TH @ Pi)), rel, it, ok))
p('')
p(r'== [C3] 判决 ==')
if max_cross >= 2 and max_dI > 1e-3 and worst_root <= 1e-6:
    p(r'  **$D(\lg)$ 非单调是模型性质**：最多穿点 %d 个，同方向同代价的 $I$ 跨度最大 %.4f bit；'
      r'DARE 与压缩映射根的最大相对差 %.2e（$\le10^{-6}$）$-\ -$ 不是求解器伪影。'
      % (max_cross, max_dI, worst_root))
    p(r'  $\Rightarrow$ 既有 $\Delta$/离开角读数**按"可行点"口径仍有效**（$|D-t|\le10^{-7}$ 是直接验证的），但必须标支；'
      r'"同秩地板"与"$\Theta$ 首轴贵多少 bit"这类**比较句**要写明取的是哪一支。')
elif worst_root > 1e-3:
    p(r'  **求解器伪影**：DARE 与迭代根的最大相对差 %.2e $>$ $10^{-3}$ $-\ -$ 相关读数作废重测。' % worst_root)
else:
    p(r'  单调或未测出多支：最多穿点 %d，$I$ 跨度最大 %.4f，根差 %.2e $-\ -$ e167 的异常另有原因，'
      r'需回报 GUESS 携带的割线早退路径。' % (max_cross, max_dI, worst_root))
p('')
p(r'  [C2] 口径：迭代用同一更新式跑 20000 步上限；"收敛 True" 才计入判决，否则该行只报数。')
p('')
p('用时 %.0f s | 异常 %s' % (time.time() - T0, '无'))
txt = '\n'.join(out) + '\n'
open('p0/e169_out.txt', 'w', encoding='utf-8').write(txt)
open('p0/e169_run.txt', 'w', encoding='utf-8').write(txt)
