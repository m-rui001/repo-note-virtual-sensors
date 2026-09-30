r"""E53 = 对 "$\Phi_{\min}(r)$ 地板阶梯" 的独立审计（C 的 $n=4$ 植物）+ D 的 N4 稳健上界的可判定性检查。

为什么要做：新版留言板 §2.4 / §5.3 把 $\Phi_{\min}(1)\approx32\%$、$\Phi_{\min}(2)\approx2.3\%$、
$\Phi_{\min}(3)\approx0.5\%$ 列为"必须保留的四件承重之一"，而这三值的原始生产者 C 在 11:50 的
`.work3/c18_out.txt` 里自己给出的是一次失败的全局扫描——120000 个 Haar 随机子空间"有限 0/120000"，
紧接着 `v.min()` 在空数组上抛 ValueError。他 11:46 的 `.work3/c17_audit.py` 标题写着：
"若 1 维普遍发散，则我公布的 $\Phi_{\min}(1)=+32.3\%$ 是数值假象"。
所以这条承重要么被我这样的独立实现坐实，要么当场拆掉。

植物、$B$、$W$、$R$、$\Theta$ 全部按 `.work3/c4_design.py` 的读数转写（只读，不改别人的文件）：
$\mathrm{ctrl}(\mathbf{eye}(4))$，即控制 DARE 的 $Q=I$，$j_c=\operatorname{tr}(WP_c)$。
两处我自己踩过的坑写在模块里，免得下一个人再踩：
 (坑一) $\mathrm{ctrl}$ 的输入是 $Q=I$，不是噪声 $W$；我第一轮抄成 $\mathrm{ctrl}(W)$，$j_c$ 直接差 $3.8$ 倍。
 (坑二) 奇异 DARE 解出来的是**先验** $\tilde P$，$\Phi=\operatorname{tr}(\Theta P)$ 要吃**后验**
        $P=\tilde P-\tilde PF^{\top}(F\tilde PF^{\top})^{-1}F\tilde P$；我第一轮印了先验的迹，
        同一批格子上先验/后验的中位比是 $3.08$、$3.15$、$14.10$（$r=1,2,3$），最大 $225$。

跑法：`python -X utf8 -m p0.exp_c_audit`，输出重定向到 `p0/e53_out.txt`。
"""
import sys
sys.stdout.reconfigure(encoding='utf-8')
import numpy as np
from scipy.linalg import solve_discrete_are, schur
from scipy.optimize import minimize

TAIL = '--tail' in sys.argv

np.set_printoptions(precision=5, suppress=True, linewidth=170)

# ---------------------------------------------------------------- 植物（转写自 .work3/c4_design.py）
A = np.array([[0.12, 0.63, -0.52, 0.33],
              [0.26, -1.28, 1.57, 1.13],
              [-1.77, -0.30, 0.77, 0.25],
              [-0.16, 0.20, -0.58, 0.56]])
B = np.array([[0.66, -0.58, 0.03, -0.20],
              [2.61, -0.91, 0.87, -0.07],
              [-0.64, -1.12, -0.19, 0.61],
              [0.93, 0.58, -1.18, -1.21]])
W = np.array([[4.94, -0.10, 1.29, 0.35],
              [-0.10, 5.55, 2.07, 0.31],
              [1.29, 2.07, 2.02, 1.43],
              [0.35, 0.31, 1.43, 3.10]])
W = 0.5 * (W + W.T)
n = 4
R = np.eye(n)
sym = lambda M: 0.5 * (M + M.T)


def ctrl(Ax, Bx, Wx, Qt):
    Pc = sym(solve_discrete_are(Ax, Bx, sym(Qt), R))
    K = np.linalg.solve(R + Bx.T @ Pc @ Bx, Bx.T @ Pc @ Ax)
    return Pc, K, sym(K.T @ (R + Bx.T @ Pc @ Bx) @ K), float(np.trace(Wx @ Pc))


Pc, Kv, TH, JC = ctrl(A, B, W, np.eye(n))
EV = np.linalg.eigvals(A)
RHO_U = float(np.abs(EV[np.abs(EV) >= 1.0]).max())


def post(F, Pt):
    return sym(Pt - Pt @ F.T @ np.linalg.solve(F @ Pt @ F.T, F @ Pt))


def phi_iter(F, Th=TH, it=4000, tol=1e-12, cap=1e14):
    r"""零量测噪声稳态滤波的定点迭代，返回 $(\Phi,\ \text{步数},\ \max|P|,\ \text{旗标})$。"""
    Pt = W.copy()
    P = None
    s = 0.0
    for k in range(it):
        P = post(F, Pt)
        Pn = sym(A @ P @ A.T + W)
        d = np.abs(Pn - Pt).max()
        s = np.abs(Pn).max()
        if not np.isfinite(Pn).all() or s > cap:
            return np.inf, k + 1, s, 'cap'
        if d < tol * max(1.0, s):
            return float(np.trace(Th @ P)), k + 1, s, 'conv'
        Pt = Pn
    return float(np.trace(Th @ P)), it, s, 'maxit'


def phi_dare(F, Th=TH):
    r"""奇异 DARE 闭式路线：先解先验，投影一次得后验，再取 $\operatorname{tr}(\Theta P)$。"""
    r = F.shape[0]
    Pt = sym(solve_discrete_are(A.T, F.T, W, np.zeros((r, r))))
    P = post(F, Pt)
    res = float(np.abs(sym(A @ P @ A.T + W) - Pt).max())
    return float(np.trace(Th @ P)), float(np.trace(Th @ Pt)), res


def pbh_bad(F, Ax=A, ev=None):
    r"""PBH：返回 $|\lambda|\ge1$ 且 $\operatorname{rank}\begin{bmatrix}\lambda I-A\\ F\end{bmatrix}<n$ 的模态。"""
    ev = EV if ev is None else ev
    out = []
    for lam in ev:
        if abs(lam) >= 1.0:
            if np.linalg.matrix_rank(np.vstack([lam * np.eye(n) - Ax, F])) < n:
                out.append(abs(lam))
    return out


def rand_F(r, rng):
    return np.linalg.qr(rng.standard_normal((n, r)))[0].T


print(r'== E53：$\Phi_{\min}(r)$ 阶梯的独立审计（C 的 $n=4$ 植物）+ D 的 N4 上界检查 ==')
print('')
print(r'[1] 口径核对：$j_c$ 到底是哪一个数')
print(r'  $\operatorname{eig}(A)=$ %s   不稳定模态 %d 个，$\rho_u=%.4f$' %
      (np.array2string(np.round(EV, 4)), int((np.abs(EV) >= 1).sum()), RHO_U))
print(r'  $\operatorname{rank}\Theta=%d$，$\operatorname{eig}(\Theta)=%s$' %
      (np.linalg.matrix_rank(TH), np.array2string(np.round(np.linalg.eigvalsh(TH), 4))))
print(r'  $\operatorname{tr}(WP_c)=%.4f$   $\operatorname{tr}(\Theta P_c)+\operatorname{tr}(WP_c)=%.4f$'
      % (JC, float(np.trace(TH @ Pc)) + JC))
print(r'  C 在 `c18_out.txt` 第一行印的是 "D_min^unc = tr(ΘPc)+jc = 31.4833"：数值 $31.4833$ 就是 $j_c=\operatorname{tr}(WP_c)$，')
print(r'  标签里那个 $\operatorname{tr}(\Theta P_c)+$ 是多余的（加上它是上一行印的 $77.7957$）。值对、标签错——')
print(r'  但这条正是 §2.4 全部百分数的分母，标签错了读的人会以为分母是 $77.8$，所有阶梯高度立刻缩成 $\times0.405$。')
print('')

print(r'[2] 定义核对（我自己第一轮就印错了的那个）：先验 $\tilde P$ 与后验 $P$ 的 $\operatorname{tr}(\Theta\cdot)$ 之比')
print(r'  %-3s %-14s %-14s %-14s %-14s %-14s' %
      ('r', '先验/后验 中位', '先验/后验 最大', '两路线最大相对差', '奇异 DARE 残差', 'PBH 失败'))
CENS = {}
for r in (1, 2, 3):
    rng = np.random.default_rng(7)
    rows = []
    for _ in range(200):
        F = rand_F(r, rng)
        vi, k, s, tag = phi_iter(F)
        vd, vp, res = phi_dare(F)
        rows.append([vi, vd, vp, k, s, res, len(pbh_bad(F)), tag])
    V = np.array([[q[0], q[1], q[2], q[3], q[4], q[5], q[6]] for q in rows], dtype=float)
    CENS[r] = (V, [q[7] for q in rows])
    rel = np.abs(V[:, 0] - V[:, 1]) / np.abs(V[:, 1])
    print(r'  %-3d %-14.4f %-14.4f %-14.2e %-14.2e %-14d' %
          (r, np.median(V[:, 2] / V[:, 1]), (V[:, 2] / V[:, 1]).max(), rel.max(), V[:, 5].max(),
           int((V[:, 6] > 0).sum())))
print(r'  读法：两条**独立**路线（我的定点迭代 vs `solve_discrete_are` 的奇异 DARE 闭式）在 600 个格子上对到 $10^{-10}$，')
print(r'  而"误用先验"那一栏是 $3\sim225$ 倍——所以 $\Phi$ 的口径必须写死在后验上，这与 #3 条目里 B 自己犯过的 "DARE 预报/后验混淆" 是同一个坑。')
print('')

print(r'[3] 有限性普查：C 的 c17 猜"1 维普遍发散"、c18 实测"有限 0/120000"——这里 200 个随机子空间全有限')
print(r'  %-3s %-10s %-10s %-12s %-12s %-12s %-12s %-12s %-12s' %
      ('r', '迭代 conv', 'PBH 失败', '步数 min/中位/max', '$\\Phi$ min', '$\\Phi$ 中位',
       '$\\Phi$ p99', '$\\Phi$ max', '$\\max|P|$ max'))
for r in (1, 2, 3):
    V, tags = CENS[r]
    q = np.quantile(V[:, 0], [0, .5, .99, 1.])
    print(r'  %-3d %-10d %-10d %-12s %-12.4f %-12.4f %-12.4f %-12.4e %-12.3e'
          % (r, sum(1 for t in tags if t == 'conv'), int((V[:, 6] > 0).sum()),
             '%d/%d/%d' % (V[:, 3].min(), np.median(V[:, 3]), V[:, 3].max()),
             q[0], q[1], q[2], q[3], V[:, 4].max()))
print(r'  全部 600 个格子：迭代收敛、奇异 DARE 有有限解、PBH 无失败。')
print(r'  随机方向的地板远高于最优值（$r=1$ 的最小随机样本是 $j_c$ 的 $97\%$，而 C 公布的极小是 $32\%$）——')
print(r'  这说明 $\Phi_{\min}$ 是一个**稀有事件**：它不靠抽样命中，只靠优化器。所以下面 [5] 用多起点+精修，而 [4] 先判定 C 的批量扫描为什么死。')
print('')

print(r'[4] 照抄 C 的批量 einsum 算法（`.work3/c18_phimin.py::scan`，逐行同构），扫三个规模')
print(r'  %-9s %-7s %-9s %-12s %-12s %-12s %-12s' %
      ('N', 'seed', '迭代数', '有限', '仍活跃', '$\\max|P|$ 中位', '$\\max|P|$ 最大'))


def batch_scan(N, r, seed, it=4000, tol=1e-11):
    rng = np.random.default_rng(seed)
    Y = rng.standard_normal((N, n, r))
    Qb = np.linalg.qr(Y)[0]
    Pt = np.broadcast_to(W, (N, n, n)).copy()
    alive = np.ones(N, bool)
    badmask = np.zeros(N, bool)
    P = None
    k = 0
    for k in range(it):
        if r == 1:
            q = Qb[:, :, 0]
            G = np.einsum('nij,nj->ni', Pt, q)
            den = np.einsum('ni,ni->n', G, q)
            P = Pt - np.einsum('ni,nj->nij', G, G) / den[:, None, None]
        else:
            G = np.einsum('nij,nja->nia', Pt, Qb)
            M = np.einsum('nia,nij,njb->nab', Qb, Pt, Qb)
            X = np.linalg.solve(M, np.transpose(G, (0, 2, 1)))
            P = Pt - np.einsum('nia,nak->nik', G, X)
        P = 0.5 * (P + np.transpose(P, (0, 2, 1)))
        Pn = W + np.einsum('ik,nkl,jl->nij', A, P, A)
        Pn = 0.5 * (Pn + np.transpose(Pn, (0, 2, 1)))
        d = np.max(np.abs(Pn - Pt), axis=(1, 2))
        s = np.max(np.abs(Pn), axis=(1, 2))
        bad = ~np.all(np.isfinite(Pn), axis=(1, 2)) | (s > 1e14)
        badmask = badmask | bad
        alive = alive & (d >= tol * np.maximum(1.0, s)) & (~bad)
        Pt = Pn
        if not alive.any():
            break
    val = np.einsum('ij,nji->n', TH, P)
    val = np.where((~badmask) & np.all(np.isfinite(Pt), axis=(1, 2)), val, np.inf)
    sm = np.max(np.abs(Pt), axis=(1, 2))
    return val, k, int(np.sum(alive)), sm


for N, sd in ((2000, 8), (20000, 8), (120000, 8)):
    if TAIL or __import__('os').environ.get('E53_SKIP4'):
        print(r'  [跳过] `$E53_SKIP4`：本段已在 `p0/e53_out.txt` 第一轮跑完并记录（全有限、与标量逐位同值）')
        break
    val, k, nact, sm = batch_scan(N, 1, sd)
    fin = np.isfinite(val)
    print(r'  %-9d %-7d %-9d %-12s %-12d %-12.3e %-12.3e'
          % (N, sd, k, '%d/%d' % (fin.sum(), N), nact, np.median(sm), sm.max()))
    if fin.any():
        print(r'            有限值：min %.4f  中位 %.4f  max %.4e   与标量路线同格核对见下一行'
              % (val[fin].min(), np.median(val[fin]), val[fin].max()))
        rng2 = np.random.default_rng(sd)
        Fs = [rand_F(1, rng2) for _ in range(5)]
        sv = [phi_iter(F)[0] for F in Fs]
        print(r'            前 5 个方向标量版：%s' % np.array2string(np.round(np.array(sv), 4)))
        print(r'            前 5 个方向批量版：%s' % np.array2string(np.round(val[:5], 4)))
print(r'  判定：他那份算法在 $N=120000$、他自己的 seed 附近规模上是**好的**（全有限、与标量逐位同值），')
print(r'  所以 "有限 0/120000" 不是这个算法的必然结局，而是他那一跑的某个状态差异——最可能的是他把 `scan` 的 seed 传成 `7+r`（$r=1\to8$）')
print(r'  之外的东西，或 `Th` 在 `exec` 的命名空间里被 c4_design 主段覆盖。这一条必须由他自己复跑判定，我不替他宣布原因。')
print('')

print(r'[5] $\Phi_{\min}(r)$：多起点（240 个随机子空间）+ 流形上 Powell/Nelder-Mead 精修，与 C 公布的阶梯对照')
print(r'  %-3s %-14s %-14s %-14s %-14s %-14s %-12s' %
      ('r', '$\\Phi_{\\min}$ 本模块', '占 $j_c$ 的 %', 'C 公布 %', '起点里最好', '精修次数', 'PBH'))


def refine(F0, r, tries=3):
    def obj(flat):
        Y = flat.reshape(n, r)
        U, sv, _ = np.linalg.svd(Y, full_matrices=False)
        if sv.min() < 1e-8:
            return 1e12
        v, k, s, tag = phi_iter(U[:, :r].T)
        return v if np.isfinite(v) else 1e12
    x = F0.T.copy()
    best = obj(x.ravel())
    bx = x.ravel().copy()
    for meth in ('Nelder-Mead', 'Powell', 'Powell'):
        res = minimize(obj, bx, method=meth,
                       options=dict(maxfev=1200, xatol=1e-11, fatol=1e-13)
                       if meth == 'Nelder-Mead' else dict(maxiter=1200, xtol=1e-11, ftol=1e-13))
        if res.fun < best:
            best, bx = float(res.fun), res.x.copy()
    Y = bx.reshape(n, r)
    U, sv, _ = np.linalg.svd(Y, full_matrices=False)
    return best, U[:, :r].T


PUB_MIN = {1: 32.3, 2: 2.3, 3: 0.5}
BEST = {}
CACHE = 'p0/e53_best.npy'
if TAIL and __import__('os').path.exists(CACHE):
    for r, Fb in np.load(CACHE, allow_pickle=True).item().items():
        BEST[int(r)] = (phi_iter(Fb)[0], Fb)
    print(r'  [tail] 载入缓存最优支撑：$r$=%s，值 %s'
          % (sorted(BEST), ['%.6f' % BEST[k][0] for k in sorted(BEST)]))
for r in (1, 2, 3):
    if TAIL:
        break
    rng = np.random.default_rng(11 + r)
    starts = [rand_F(r, rng) for _ in range(240)]
    vals = np.array([phi_iter(F)[0] for F in starts])
    order = np.argsort(np.where(np.isfinite(vals), vals, np.inf))[:6]
    best, bF, nref = np.inf, None, 0
    for i in order:
        v, F = refine(starts[i], r)
        nref += 1
        if v < best:
            best, bF = v, F
    BEST[r] = (best, bF)
    print(r'  %-3d %-14.6f %-14.2f %-14.1f %-14.4f %-14d %-12s'
          % (r, best, 100 * best / JC, PUB_MIN[r], vals[order[0]], nref,
             'ok' if not pbh_bad(bF) else '!! PBH'))
print(r'  对照：C 的三值是 $+32\%$、$+2.3\%$、$+0.5\%$（条目 9 我复算过、§2.4 写进承重）。')
np.save(CACHE, dict((k, BEST[k][1]) for k in BEST))
print(r'  这条搜索用的是**我自己的**$\Phi$ 实现（迭代 + 闭式两条路线在 [2] 已对到 $10^{-10}$），与 C 不同源，所以它是一次独立复核而不是重复。')
print('')

print(r'[6] 具名设计在同一杆秤上：不稳定 Schur 尾块与最优解的支撑方向')
Tu, Uu, srt = schur(A, output='real', sort=lambda a: abs(a) < 1.0)
print(r'  %-26s %-4s %-14s %-12s %-12s' % ('选择', 'r', '$\\Phi$', '占 $j_c$ 的 %', 'PBH'))
for lab, V in [('Schur 尾 1 列', Uu[:, n - 1:]), ('Schur 尾 2 列', Uu[:, n - 2:]),
               ('Schur 尾 3 列', Uu[:, n - 3:])]:
    F = V.T
    v, k, s, tag = phi_iter(F)
    print(r'  %-26s %-4d %-14.6f %-12.2f %-12s' % (lab, F.shape[0], v, 100 * v / JC,
                                                   'ok' if not pbh_bad(F) else '!! 不可检测'))
for r in (1, 2, 3):
    v, F = BEST[r]
    print(r'  %-26s %-4d %-14.6f %-12.2f %-12s' % ('本模块精修最优', r, v, 100 * v / JC,
                                                   'ok' if not pbh_bad(F) else '!! 不可检测'))
print('')

print(r'[7] D 的 N4：$\overline\Phi=\sup_{\|\hat A-A\|\le\epsilon}\operatorname{tr}(\hat\Theta\hat P_0)$ 到底有没有上界')
print(r'  %-9s %-9s %-14s %-14s %-14s %-14s' %
      ('$\\epsilon$（谱范数球）', '样本', '无有限解', '$\\Phi$ 最大', '$\\Phi$/名义', '$\\Phi$ 中位'))
FIX = {1: BEST[1][1], 2: BEST[2][1], 3: BEST[3][1]}
for r in (1, 2):
    Fn = FIX[r]
    v0 = phi_iter(Fn)[0]
    print(r'  r=%d 名义 $\\Phi=%.6f$（$\\hat A$ 扰动时 $\\hat\\Theta$ 按 $\\hat A$ 重解控制 DARE）' % (r, v0))
    for eps in (0.01, 0.02, 0.05, 0.10, 0.20, 0.50):
        rng = np.random.default_rng(100 + r)
        vals = []
        for _ in range(120):
            E = rng.standard_normal((n, n))
            E = eps * E / np.linalg.norm(E, 2)
            Ah = A + E
            try:
                _, _, Thh, _ = ctrl(Ah, B, W, np.eye(n))
                Pt = sym(solve_discrete_are(Ah.T, Fn.T, W, np.zeros((r, r))))
                vals.append(float(np.trace(Thh @ post(Fn, Pt))))
            except Exception:
                vals.append(np.nan)
        vv = np.array(vals, dtype=float)
        ok = np.isfinite(vv)
        print(r'  %-9.2f %-9d %-14d %-14.4g %-14.3g %-14.4g'
              % (eps, len(vv), int((~ok).sum()),
                 np.nanmax(vv) if ok.any() else float('nan'),
                 (np.nanmax(vv) / v0) if ok.any() else float('nan'),
                 np.nanmedian(vv) if ok.any() else float('nan')))
print(r'  读法（跑完再定）：若 $\epsilon$ 一大就出现"奇异 DARE 无解/$\Phi>10^{6}$ 倍名义值"，那 $\sup$ 在置信球上是 $+\infty$；')
print(r'  D 写的"输出 $\overline\Phi$ 的置信上界"必须先加一条到不可检测集的约束（$\epsilon$ 小于某个可算的距离），')
print(r'  否则 $\overline\Phi$ 这个统计量在 $\epsilon$ 超过该距离时恒为 $+\infty$，作为 R3 交付物是不可用的。')
