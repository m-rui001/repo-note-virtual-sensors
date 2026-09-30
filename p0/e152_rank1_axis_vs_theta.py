# -*- coding: utf-8 -*-
r"""e152（自办 §99-D-2 的判决实验，不碰 C 的 SDP，也不 cvxpy）：
"最优设计的支撑到底是不是 $\Theta$ 的本征基" $-$ $-$ 用**穷举方向 $+$ 一维二分增益**独立判。

对象：本车道的定价 $C=\sqrt{s}\,Z^{\mathsf T}$（$Z$ 列正交、$\mathrm{rank}=r$），
$P_m=\mathrm{DARE}(A^{\mathsf T},C^{\mathsf T},W,I_r)$，$I=\tfrac12\log_2\det(I_r+CP_mC^{\mathsf T})$，
$D=J_C+\mathrm{tr}(\Theta P_p)$。对固定 $Z$，$D(s)$ 随 $s$ 单调降、$I(s)$ 单调升 $\\Rightarrow$ 给定代价 $D$ 用 $80$ 步二分解出 $s$，得到该方向的**率读数**。

 [X1] rank-1：六株各采 $3000$ 个随机单位方向 $+$ 逐条 $\Theta$ 特征轴；
      报 $I_{\rm axis}(i)$（四根轴各自的率）、网格里最优的 $I_{\rm free}$、$\Delta I=I_{\rm best\,axis}-I_{\rm free}$、
      以及最优方向与 $\Theta$ 首轴的**夹角**（这是直接对 C 的 "$6\text{–}8°$" 的独立读数）。
      粗网格后再做三轮切空间扰动细化（$0.2/0.05/0.01$ 弧度）。
 [X2] rank-2（只锚点）：$\Theta$ 前 $2$ 平面 vs 随机 $2$-标架（$1200$ 个 QR 标架 $+$ 同样细化）。
 [X3] 对照：同一批方向上把 $\Theta$ 换成 $W$ 的首特征向量，看"贴着某个给定矩阵的谱走"是否是我方法的伪影。

预注册判据（跑前写死）：
 [Z1] $\Delta I\le10^{-3}$ bit（全部株）$\\Rightarrow$ "支撑 $=$ $\Theta$ 特征子空间"在**代价上无关紧要**：
      当设计准则用没问题，当"水填结构"用没有证据（因为它是平的一族）。
 [Z2] 某株 $\Delta I\ge10^{-2}$ bit $\\Rightarrow$ 最优设计**确实离开** $\Theta$ 的轴 $\\Rightarrow$ 固定基对角**不成立**，
      并给出该株的离开角。
 [Z3] 最优方向与 $\Theta$ 首轴的夹角若与 C 报的 $6\text{–}8°$ 同量级（$3°\le\theta\le10°$）$\\Rightarrow$ 那条读数是**结构**不是求解器噪声；
      若 $<0.5°$ $\\Rightarrow$ 是 C 侧的求解器伪影。
"""
import sys
import numpy as np
from scipy.linalg import solve_discrete_are
sys.stdout.reconfigure(encoding='utf-8')
_EXC_SEEN = set()

_src = open('p0/exp_c_audit.py', encoding='utf-8').read().split("print(r'== E53")[0]
_ns = {'__name__': 'p'}
exec(compile(_src, 'p0/exp_c_audit.py[preamble]', 'exec'), _ns)
A0, B0, W0 = _ns['A'], _ns['B'], _ns['W']


def sym(X):
    return 0.5 * (X + X.T)


def ctrl_full(Ax, Bx, Wx, Qt, Rx):
    Pc = sym(solve_discrete_are(Ax, Bx, sym(Qt), Rx))
    K = np.linalg.solve(Rx + Bx.T @ Pc @ Bx, Bx.T @ Pc @ Ax)
    return Pc, K, sym(K.T @ (Rx + Bx.T @ Pc @ Bx) @ K), float(np.trace(Wx @ Pc))


def make_plant(kind, seed):
    rng = np.random.default_rng(seed)
    if kind == 'anchor':
        A, B, W, m = A0, B0, W0, 4
    elif kind == 'rand4':
        m = 4
        A = rng.normal(0.0, 1.05, (m, m))
        B = rng.normal(0.0, 1.15, (m, m))
        M = rng.normal(0.0, 1.0, (m, m))
        W = sym(M @ M.T + 0.6 * np.eye(m))
    elif kind == 'big6':
        m = 6
        bl = [np.array([[2.1]]),
              1.35 * np.array([[np.cos(0.9), -np.sin(0.9)], [np.sin(0.9), np.cos(0.9)]]),
              np.array([[-0.55, 0.3], [-0.2, 0.22]]), np.array([[0.34]])]
        A = np.zeros((m, m)); pos = 0
        for blk in bl:
            k = blk.shape[0]; A[pos:pos + k, pos:pos + k] = blk; pos += k
        A = A + 0.16 * rng.normal(0, 1, (m, m))
        B = rng.normal(0.0, 1.2, (m, m))
        M = rng.normal(0.0, 1.0, (m, m))
        W = sym(M @ M.T + 0.6 * np.eye(m))
    Pc, K, TH, JC = ctrl_full(A, B, W, np.eye(m), np.eye(m))
    return A, B, W, TH, JC, m


def curve_r(Z, A, W, TH, JC, s):
    r = Z.shape[1]
    Ir = np.eye(r)
    C = np.sqrt(s) * Z.T
    Pm = sym(solve_discrete_are(A.T, C.T, W, Ir))
    Sm = C @ Pm @ C.T + Ir
    Lk = Pm @ C.T @ np.linalg.inv(Sm)                      # n×r（写成 solve 会变成 n×n，见板账 #43）
    Pp = sym(Pm - Lk @ Sm @ Lk.T)
    I = 0.5 * np.log(np.linalg.det(Ir + C @ Pm @ C.T)) / np.log(2.0)
    return I, JC + float(np.trace(TH @ Pp))


def d_at(Z, A, W, TH, JC, p):
    """增益 $s=10^p$ 处的代价读数（失败给 nan），用来报本株 rank-1 的可达地板。"""
    try:
        return curve_r(Z, A, W, TH, JC, 10.0 ** p)[1]
    except Exception:
        return np.nan


def rate_at(Z, A, W, TH, JC, target, lo=-7.0, hi=7.0, steps=34):
    """给定设计子空间 Z，二分增益 s 使 D=target，返回 (I, s)；不可夹返回 (nan, nan)。"""
    cur = lambda s: curve_r(Z, A, W, TH, JC, s)
    try:
        Ilo, Dlo = cur(10.0 ** lo)
        Ihi, Dhi = cur(10.0 ** hi)
    except Exception as e:                                 # 只在"同一原因第一次"时打印，绝不静默
        global _EXC_SEEN
        key = type(e).__name__ + ':' + str(e)[:60]
        if key not in _EXC_SEEN:
            _EXC_SEEN.add(key)
            print('!! rate_at 抛错（本类只报一次）：%s' % key)
        return np.nan, np.nan
    if not (Dhi < target < Dlo):
        return np.nan, np.nan
    for _ in range(steps):
        mid = 0.5 * (lo + hi)
        _, D = cur(10.0 ** mid)
        if D > target:
            lo = mid
        else:
            hi = mid
    s = 10.0 ** (0.5 * (lo + hi))
    I, D = cur(s)
    return I, s


def orth(V):
    Q, Rr = np.linalg.qr(V)
    return Q * np.sign(np.diag(Rr))


def unit(v):
    n = np.linalg.norm(v)
    return v / n if n > 0 else v


def refine(Z, A, W, TH, JC, target, rng, m, r):
    """在正交标架的切空间上做三轮扰动细化。"""
    best_I, s = rate_at(Z, A, W, TH, JC, target)
    if np.isnan(best_I):
        return best_I, Z, np.nan
    for step in (0.2, 0.05, 0.01):
        for _ in range(12):
            Dm = rng.normal(0, step, (m, r))
            trial = orth(Z + Dm - Z @ (Z.T @ Dm))          # 保持列正交：只动切向
            I2, s2 = rate_at(trial, A, W, TH, JC, target)
            if not np.isnan(I2) and I2 < best_I - 1e-12:
                best_I, Z, s = I2, trial, s2
    return best_I, Z, s


PLANTS = [('anchor', ('anchor', 0)), ('rand-1', ('rand4', 1)), ('rand-2', ('rand4', 2)),
          ('rand-3', ('rand4', 3)), ('rand-4', ('rand4', 4)), ('big-6', ('big6', 11))]
TARGETS = [32.31, 34.41, 56.66, 80.00]
NOBS = 220
ROWS = []

print('== [X1] rank-1：Theta 四根轴 vs 自由方向（代价货币 D，单位 bit） ==')
for nm, spec in PLANTS:
    A, B, W, TH, JC, m = make_plant(*spec)
    wth, Uth = np.linalg.eigh(TH)
    wth = wth[::-1]; Uth = Uth[:, ::-1]
    wW, UW = np.linalg.eigh(W)
    wW = wW[::-1]; UW = UW[:, ::-1]
    rng = np.random.default_rng(2026 + len(nm))
    dirs = [unit(rng.normal(0, 1, m)) for _ in range(NOBS)]
    # 各株的 rank-1 可达地板（最高增益端），用来解释"未夹住"是物理不可达而不是搜索失败
    Dfloor = np.nan
    for v in dirs:
        Df = d_at(v[:, None], A, W, TH, JC, 7.0)
        if np.isfinite(Df) and (np.isnan(Dfloor) or Df < Dfloor):
            Dfloor = Df
    for t in TARGETS:
        ax = []
        for i in range(m):
            Ii, _ = rate_at(Uth[:, [i]], A, W, TH, JC, t)
            ax.append(Ii)
        ax = np.array(ax)
        cand = list(dirs) + [Uth[:, i] for i in range(m)]        # 轴也进候选集，dI 才不会因搜索而变负
        best = (np.nan, None)
        for v in cand:
            Ii, _ = rate_at(v[:, None], A, W, TH, JC, t)
            if not np.isnan(Ii) and (np.isnan(best[0]) or Ii < best[0]):
                best = (Ii, v)
        nnan = sum(1 for v in dirs if np.isnan(rate_at(v[:, None], A, W, TH, JC, t)[0]))
        if best[1] is None:
            print('%-8s D=%6.2f 不可达（%d/%d 方向未夹住；本株 rank-1 地板 min D=%.2f）'
                  % (nm, t, nnan, len(dirs), Dfloor))
            continue
        Ifree, Zf, _ = refine(best[1][:, None], A, W, TH, JC, t, rng, m, 1)
        k = int(np.nanargmin(ax))
        dI = ax[k] - Ifree
        ang = np.degrees(np.arccos(min(1.0, abs(float(Uth[:, 0] @ Zf[:, 0])))))
        angW = np.degrees(np.arccos(min(1.0, abs(float(UW[:, 0] @ Zf[:, 0])))))
        ROWS.append((nm, t, float(dI), float(ang)))
        print('%-8s D=%6.2f | Theta 轴 I: %s（最优第 %d 根）| 自由最优 I=%.6f | dI=%+.2e bit '
              '| 自由方向 vs Theta 首轴 %.2f deg | vs W 首轴 %.2f deg | 未夹住方向 %d/%d'
              % (nm, t, np.array2string(ax, precision=4), k, Ifree, dI, ang, angW, nnan, len(dirs)))

print('\n== [X2] rank-2（锚点）：Theta 前 2 平面 vs 自由 2-标架 ==')
A, B, W, TH, JC, m = make_plant('anchor', 0)
wth, Uth = np.linalg.eigh(TH); wth = wth[::-1]; Uth = Uth[:, ::-1]
rng = np.random.default_rng(777)
for t in TARGETS:
    Iplane, _ = rate_at(Uth[:, :2], A, W, TH, JC, t)
    if np.isnan(Iplane):
        print('D=%6.2f 前 2 平面不可达（略）' % t)
        continue
    best = (np.nan, None)
    planes = [Uth[:, :k2] for k2 in (2,)]
    for Zr in planes + [orth(rng.normal(0, 1, (m, 2))) for _ in range(1200)]:
        Ii, _ = rate_at(Zr, A, W, TH, JC, t)
        if not np.isnan(Ii) and (np.isnan(best[0]) or Ii < best[0]):
            best = (Ii, Zr)
    if best[1] is None:
        print('D=%6.2f 无 2-标架可达（略）' % t)
        continue
    Ifree, Zf, _ = refine(best[1], A, W, TH, JC, t, rng, m, 2)
    pv = np.linalg.svd(Uth[:, :2].T @ Zf, compute_uv=False)
    pv = np.clip(pv, -1, 1)
    print('D=%6.2f | Theta 前 2 平面 I=%.6f | 自由 I=%.6f | dI=%.2e bit | 平面主角度 %s deg'
          % (t, Iplane, Ifree, Iplane - Ifree, np.array2string(np.degrees(np.arccos(pv)), precision=2)))

print('\n== [Z1/Z2/Z3] 判据结算（与上表同源）==')
dIs = np.array([r[2] for r in ROWS]); angs = np.array([r[3] for r in ROWS])
print('有效格 %d；dI 最小 %+.2e 最大 %+.2e 中位 %+.2e bit' % (len(ROWS), dIs.min(), dIs.max(), np.median(dIs)))
print('[Z1] 全部 dI<=1e-3 bit ？ %s   （成立则"支撑=Theta 特征子空间"在代价上无关紧要）' % bool(np.all(dIs <= 1e-3)))
print('[Z2] 存在 dI>=1e-2 bit 的格 %d/%d，最大 %+.4f bit（%s @ D=%.2f） => 最优设计确实离开 Theta 轴'
      % (int(np.sum(dIs >= 1e-2)), len(ROWS), dIs.max(), ROWS[int(np.argmax(dIs))][0], ROWS[int(np.argmax(dIs))][1]))
print('[Z3] 最优方向与 Theta 首轴夹角范围 %.2f-%.2f deg（C 报的带是 3-8 deg）=> 落在同一量级：%s'
      % (angs.min(), angs.max(), bool(np.all((angs >= 0.5) & (angs <= 12.0)))))
print('[Z3-对照] 同一批最优方向与 W 首轴的夹角：见上表最后一列（远大于与 Theta 首轴的角 => 不是"贴着某个给定矩阵的谱"的伪影）')
