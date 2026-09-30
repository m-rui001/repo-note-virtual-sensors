# -*- coding: utf-8 -*-
r"""实验 133（§82-F 预注册）：**只用拟合窗内的信息**能否识别出"$\gamma$ 单值拟合坏了"的设计？

背景（§82-D）：基函数没错（指数赢），唯一失效模式是 `rand-1` 那类植物上 $\gamma$ 沿 $\Delta I$ 不是常数。
若窗内残差就能预警，定价器才敢对没见过的植物说话；若不能，正文只能退回"只对自己拟合过的植物族有效"。

跑前写死的口径：
 [V1] 窗 $[0.30,0.60]$（300 点 $s$ 网格），左半 $[0.30,0.45]$ / 右半 $[0.45,0.60]$ 各 $\ge4$ 点，否则该设计剔除并计数。
      地板 $D_{\rm floor}$ 同 e132（$s=10^6/10^8/10^{12}$，相对漂移 $>1\%$ 剔除）。
 [V2] 不稳定度 $q=\big|\ln(\gamma_{\rm hi}/\gamma_{\rm lo})\big|$，两次拟合都是**带真地板的两参数** $D_{\rm floor}+\tilde K/(e^{\gamma\Delta I}-1)$。
 [V3] "真坏"度量：同一设计用**全窗** $\gamma$（TF2 口径）在 $f\in\{1.83,2.50,3.67\}$（即 $\Delta I\in\{1.1,1.5,2.2\}$）定价，
      取三档的 $\max|e|$；真值 60 步二分。三档的格数必须相等（纪律 26 第 5 款：整列异常先怀疑代码）。
 [V4] 判据：① 池化 Spearman $\rho(q,\max|e|)\ge0.6$；
      ② 存在 $\tau$ 使"标记 $=q>\tau$"的设计中 $\max|e|>8\%$ 的比例 $\ge0.7$，**且**在排除 `rand-1` 的 5 株上误标率 $\le15\%$。
      ①②同时满足 ⇒ 正文可写 in-window flag，并报 $\tau$ 与混淆矩阵；否则明写"窗内不可识别"。
 [V5] 附加（不改变 V4）：逐株 $q$ 的中位/$90\%$ 分位；逐株内部 $\rho(q,\max|e|)$；留一株外推（用其余 5 株定 $\tau$，看 `rand-1` 命中情况）。
"""
import sys
import numpy as np
from scipy.linalg import solve_discrete_are, eigh
from scipy.optimize import curve_fit
from scipy.stats import spearmanr
sys.stdout.reconfigure(encoding='utf-8')

_src = open('p0/exp_c_audit.py', encoding='utf-8').read().split("print(r'== E53")[0]
_ns = {'__name__': 'p'}
exec(compile(_src, 'p0/exp_c_audit.py[preamble]', 'exec'), _ns)
A0, B0, W0 = _ns['A'], _ns['B'], _ns['W']
sym = _ns['sym']
GT = 2.0 * np.log(2.0)
WLO, WMID, WHI = 0.30, 0.45, 0.60
TG_DS = [1.1, 1.5, 2.2]          # f = 1.83 / 2.50 / 3.67
SGRID = 10.0 ** np.linspace(-6.0, 2.5, 300)
FLOORS = [10.0 ** 6, 10.0 ** 8, 10.0 ** 12]
BISECT = 60


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
        A = rng.normal(0.0, 1.05, (m, m)); B = rng.normal(0.0, 1.15, (m, m))
        M = rng.normal(0.0, 1.0, (m, m)); W = sym(M @ M.T + 0.6 * np.eye(m))
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
        M = rng.normal(0.0, 1.0, (m, m)); W = sym(M @ M.T + 0.6 * np.eye(m))
    Pc, K, TH, JC = ctrl_full(A, B, W, np.eye(m), np.eye(m))
    return A, B, W, TH, JC, m


def r_exp(A):
    mod = np.abs(np.linalg.eigvals(A))
    return float(np.log(mod[mod > 1.0]).sum() / np.log(2.0))


def curve(Z, s, A, CW, TH, JC):
    r = Z.shape[1]
    Ir = np.eye(r)
    C = np.sqrt(s) * Z.T
    Pm = sym(solve_discrete_are(A.T, C.T, CW, Ir))
    Sm = C @ Pm @ C.T + Ir
    Lk = Pm @ C.T @ np.linalg.inv(Sm)
    Pp = sym(Pm - Lk @ Sm @ Lk.T)
    I = 0.5 * np.log(np.linalg.det(Ir + C @ Pm @ C.T)) / np.log(2.0)
    return I, JC + np.trace(TH @ Pp)


def true_at(Z, A, W, TH, JC, target, nstep=BISECT):
    lo, hi = -7.0, 5.0
    if not (curve(Z, 10 ** lo, A, W, TH, JC)[0] < target < curve(Z, 10 ** hi, A, W, TH, JC)[0]):
        return None
    for _ in range(nstep):
        mid = 0.5 * (lo + hi)
        if curve(Z, 10 ** mid, A, W, TH, JC)[0] < target:
            lo = mid
        else:
            hi = mid
    return curve(Z, 10 ** (0.5 * (lo + hi)), A, W, TH, JC)[1]


f_exp = lambda di, K, g: K / (np.exp(g * di) - 1.0)


def fit_gam(dif, Df, dfloor):
    p, _ = curve_fit(lambda x, K, g: dfloor + f_exp(x, K, g), dif, Df,
                     p0=[max(Df.mean() * dif.mean(), 1e-6), GT],
                     bounds=([1e-12, 1e-3], [np.inf, 20.0]), maxfev=60000)
    return float(p[0]), float(p[1])


PLANTS = [('anchor', ('anchor', 0)), ('rand-1', ('rand4', 1)), ('rand-2', ('rand4', 2)),
          ('rand-4', ('rand4', 4)), ('rand-3', ('rand4', 3)), ('big-6', ('big6', 11))]
rows, drop = [], {'floor': 0, 'pts': 0, 'fit': 0, 'truth': 0}
print('== e133：窗内左/右半 γ 一致性 q = |ln(γ_hi/γ_lo)|  vs  三档外推的真实误差 ==')
for pi, (lab, spec) in enumerate(PLANTS):
    A, B, W, TH, JC, m = make_plant(*spec)
    REXP = r_exp(A)
    _, Vth = eigh(TH)
    rng = np.random.default_rng(20260930 + 7 * pi)
    pool = []
    for rk in (1, 2, min(3, m - 1)):
        for c in range(6):
            Z, _ = np.linalg.qr(rng.standard_normal((m, rk)))
            pool.append((rk, Z))
        pool.append((rk, Vth[:, :rk]))
        pool.append((rk, Vth[:, m - rk:]))
    for rk, Z in pool:
        fl = [curve(Z, s, A, W, TH, JC)[1] - JC for s in FLOORS]
        if (max(fl) - min(fl)) / max(min(fl), 1e-12) > 0.01:
            drop['floor'] += 1
            continue
        dfloor = fl[1]
        Is, Ds = [], []
        for s in SGRID:
            try:
                I, D = curve(Z, s, A, W, TH, JC)
            except Exception:
                continue
            if np.isfinite(I) and np.isfinite(D):
                Is.append(I); Ds.append(D)
        di = np.array(Is) - REXP
        Df_all = np.array(Ds) - JC
        inw = (di >= WLO) & (di <= WHI)
        lft = (di >= WLO) & (di <= WMID)
        rgt = (di > WMID) & (di <= WHI)
        if min(inw.sum(), lft.sum(), rgt.sum()) < 4:
            drop['pts'] += 1
            continue
        try:
            K_full, g_full = fit_gam(di[inw], Df_all[inw], dfloor)
            _, g_lo = fit_gam(di[lft], Df_all[lft], dfloor)
            _, g_hi = fit_gam(di[rgt], Df_all[rgt], dfloor)
        except Exception:
            drop['fit'] += 1
            continue
        errs = []
        for t in TG_DS:
            Dt = true_at(Z, A, W, TH, JC, REXP + t)
            if Dt is None:
                drop['truth'] += 1
                errs = None
                break
            errs.append(100 * abs((JC + dfloor + f_exp(t, K_full, g_full)) / Dt - 1.0))
        if errs is None:
            continue
        q = abs(np.log(g_hi / g_lo))
        rows.append({'plant': lab, 'rk': rk, 'q': q, 'g_lo': g_lo, 'g_hi': g_hi, 'g_full': g_full,
                     'mx': max(errs), 'e11': errs[0], 'e15': errs[1], 'e22': errs[2]})
    print('  [%s] 累计有效 %d（地板 %d / 点数 %d / 拟合 %d / 真值 %d 剔除）  Rexp=%.4f'
          % (lab, len(rows), drop['floor'], drop['pts'], drop['fit'], drop['truth'], REXP))

Rq = np.array([r['q'] for r in rows]); Rm = np.array([r['mx'] for r in rows])
assert all(len([r['e11'], r['e15'], r['e22']]) == 3 for r in rows), '三档列格数不等 ⇒ 先怀疑代码'
print('\n  有效设计合计 %d；三档误差列的格数按构造相等（同一设计同批算）' % len(rows))

rho = spearmanr(Rq, Rm).statistic
print('\n== [V4-①] 池化 Spearman ρ(q, max|e|) = %.3f ⇒ %s' % (rho, '≥0.6' if rho >= 0.6 else '<0.6 不满足'))

print('\n== [V5] 逐株 q 分布与内部 ρ ==')
for lab, _ in PLANTS:
    sel = [r for r in rows if r['plant'] == lab]
    if len(sel) < 6:
        continue
    q_ = np.array([r['q'] for r in sel]); m_ = np.array([r['mx'] for r in sel])
    print('  %-8s n=%2d  q med=%.3f  q p90=%.3f  max|e| med=%.2f%%  内部 ρ=%.3f'
          % (lab, len(sel), np.median(q_), np.percentile(q_, 90), np.median(m_),
             spearmanr(q_, m_).statistic))

print('\n== [V4-②] 阈值扫描（标记 = q > τ）==')
oth = [r for r in rows if r['plant'] != 'rand-1']
r1 = [r for r in rows if r['plant'] == 'rand-1']
FEAS = []
for tau in np.unique(np.round(np.percentile(Rq, np.linspace(50, 99, 26)), 3)):
    flag = np.array([r['q'] > tau for r in rows])
    if flag.sum() == 0:
        continue
    prec = np.mean([r['mx'] > 8.0 for r, f in zip(rows, flag) if f])
    good_oth = [r for r in oth if r['mx'] <= 8.0]
    fp = (np.mean([r['q'] > tau for r in good_oth]) if good_oth else float('nan'))
    hit_r1 = np.mean([r['q'] > tau for r in r1]) if r1 else float('nan')
    ok = bool(prec >= 0.7 and fp <= 0.15)
    if ok:
        FEAS.append((float(tau), float(prec), float(fp), float(hit_r1), int(flag.sum())))
    print('  τ=%6.3f  标记 %3d/%d  精度(max|e|>8%%)=%.2f  排除rand-1后误标率=%.2f  rand-1命中率=%.2f  %s'
          % (tau, int(flag.sum()), len(rows), prec, fp, hit_r1, '✔满足②' if ok else ''))
if FEAS:
    tb = FEAS[0][0]
    bm = max(FEAS, key=lambda x: x[1])
    print('\n  [V4-②] 结论：存在可行阈值 %d 个。最小可行 τ=%.3f（精度 %.2f、排除 rand-1 后误标率 %.2f、rand-1 命中 %.2f）；'
          '最高精度 τ=%.3f（精度 %.2f、误标率 %.2f、标记 %d/%d）'
          % (len(FEAS), FEAS[0][0], FEAS[0][1], FEAS[0][2], FEAS[0][3], bm[0], bm[1], bm[2], bm[4], len(rows)))
    print('  取 τ=%.3f 的逐株分解（坏 := 三档 max|e|>8%%，误标 := 不差却被标记）：' % FEAS[0][0])
    for lab, _ in PLANTS:
        sel = [r for r in rows if r['plant'] == lab]
        bad = [r for r in sel if r['mx'] > 8.0]
        print('    %-8s 坏 %2d/%2d 中被标记 %2d   误标 %2d'
              % (lab, len(bad), len(sel), sum(1 for r in bad if r['q'] > FEAS[0][0]),
                 sum(1 for r in sel if r['mx'] <= 8.0 and r['q'] > FEAS[0][0])))
else:
    print('\n  [V4-②] 结论：扫描范围内不存在满足 ② 的阈值')

print('\n== [V5] 留一株外推：用其余 5 株的中位 q 作 τ，看被检出的设计 ==')
for lab, _ in PLANTS:
    sel = [r for r in rows if r['plant'] == lab]
    trn = [r for r in rows if r['plant'] != lab]
    if len(sel) < 4 or not trn:
        continue
    tau = float(np.median([r['q'] for r in trn]))
    bad = [r for r in sel if r['mx'] > 8.0]
    print('  留出 %-8s τ(其余5株中位 q)=%.3f  株内 max|e|>8%% 的设计 %d/%d，其中被标记 %d 个'
          % (lab, tau, len(bad), len(sel), sum(1 for r in bad if r['q'] > tau)))

print('\n== 混淆口径（取 τ = 全体 q 的中位）==')
tau = float(np.median(Rq))
tab = {'标记&坏': 0, '标记&不差': 0, '未标记&坏': 0, '未标记&不差': 0}
for r in rows:
    bad = r['mx'] > 8.0
    fg = r['q'] > tau
    tab[('标记' if fg else '未标记') + ('&坏' if bad else '&不差')] += 1
print('  ' + '  '.join('%s=%d' % (k, v) for k, v in tab.items())
      + '   （坏 := 三档 max|e| > 8%）')
print('  判据 [V4]：① %s（ρ=%.3f）  ② %s  ⇒ %s'
      % ('满足' if rho >= 0.6 else '不满足', rho,
         '满足（%d 个可行 τ）' % len(FEAS) if FEAS else '不满足',
         '可写 in-window flag' if (rho >= 0.6 and FEAS) else
         '窗内不可识别 ⇒ 定价器退回"只对自己拟合过的植物族有效"'))
