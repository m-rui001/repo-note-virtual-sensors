# -*- coding: utf-8 -*-
r"""实验 132（§81-D 预注册）：把"定价误差由外推比支配"从**五点单调的观察**变成**可引用的有效域曲线**。

跑之前写死的口径（不许跑完再挑）：
 [W0] 拟合窗固定 $[\Delta I_{\rm lo},\Delta I_{\rm hi}]=[0.30,0.60]$（比 e130/e131 的窗窄一半），窗内点数 $<5$ 的设计剔除并计数。
      真值用二分：先把二分步数从 80 降到 60，**但**要拿锚点一格同时用 60/80 步算真值，相对差须 $<10^{-6}$，否则退回 80 步。
 [W1] 定价点扫描：$\Delta I_{\rm tgt}\in\{0.7,0.9,1.1,1.3,1.5,1.8,2.2,2.6,3.0,3.5,4.0\}$，外推比定义为
      $f=\Delta I_{\rm tgt}/\Delta I_{\rm hi}$（$=1$ 即窗顶点，$f<1$ 是内插）。不可达（$I$ 超出 $s\in[10^{-7},10^{5}]$ 的像）的格子明写 skip。
 [W2] 三列口径同 e130/e131：TF2（真地板 $+$ 自由 $\gamma$）、FIX（真地板 $+$ 钉 $2\ln2$）、HF2（真地板 $+$ 幂律 $\nu$）。
      地板取 $s=10^6/10^8/10^{12}$ 三档、相对漂移 $>1\%$ 剔除。每格三列必须用**同一批设计**（落盘前断言 $n$ 相等，纪律 26 第 5 款）。
 [W3] 交叉点：一株植物上 $f$ 递增扫描中，$\mathrm{med}|e|(\mathrm{HF2})<\mathrm{med}|e|(\mathrm{TF2})$ 第一次成立处的 $f$，
      记 $f_\times$（对数线性内插）；整段扫描都不成立记 $f_\times=+$inf 并明写。
 [W4] **预注册判据**：若所有有效植物的 $f_\times$ 都有限、且 $\max f_\times/\min f_\times\le2.0$、且 $\min f_\times\ge1.2$
      ⇒ 写"指数基在窗外 $f_\times$ 内占优、更远处退化"（给出数字区间）；
      其余任何情形（含"有 $\ge2$ 株永不交叉而其他株交叉"）⇒ 正文只许写
      "两种基都不是长外推的可靠定价器"。
 [W5] 次要结论（不改变 W4）：把全部 $(\text{植物},f)$ 格子的 TF2 误差对 $f$ 做 Spearman，并报每株内部 Spearman；
      若池化 $\rho\ge0.8$ 且各株内部 $\rho$ 中位 $\ge0.8$ ⇒ "误差随外推比单调"可升为跨植物陈述；否则只能逐株陈述。
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
A0, B0, W0, ctrl = _ns['A'], _ns['B'], _ns['W'], _ns['ctrl']
sym = _ns['sym']
GT = 2.0 * np.log(2.0)
WLO, WHI = 0.30, 0.60
TARGETS = [0.7, 0.9, 1.1, 1.3, 1.5, 1.8, 2.2, 2.6, 3.0, 3.5, 4.0]
SGRID = 10.0 ** np.linspace(-6.0, 2.5, 300)  # 90 点时窄窗 [0.3,0.6] 内只有 4 点，会被 [W0] 的 ≥5 门整批拒掉（分辨率问题，不是物理）
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
    ilo = curve(Z, 10 ** lo, A, W, TH, JC)[0]
    ihi = curve(Z, 10 ** hi, A, W, TH, JC)[0]
    if not (ilo < target < ihi):
        return None
    for _ in range(nstep):
        mid = 0.5 * (lo + hi)
        if curve(Z, 10 ** mid, A, W, TH, JC)[0] < target:
            lo = mid
        else:
            hi = mid
    return curve(Z, 10 ** (0.5 * (lo + hi)), A, W, TH, JC)[1]


f_exp = lambda di, K, g: K / (np.exp(g * di) - 1.0)
fpow = lambda di, K, nu: K / di ** nu

PLANTS = [('anchor', ('anchor', 0)), ('rand-1', ('rand4', 1)), ('rand-2', ('rand4', 2)),
          ('rand-4', ('rand4', 4)), ('rand-3', ('rand4', 3)), ('big-6', ('big6', 11))]

# [W0] 二分步数自检：锚点、rank-1 的一个设计，60 步 vs 80 步
Aa, Ba, Wa, THa, JCa, ma = make_plant('anchor', 0)
rea = r_exp(Aa)
Za = eigh(THa)[1][:, ma - 1:]
d60 = true_at(Za, Aa, Wa, THa, JCa, rea + 2.2, 60)
d80 = true_at(Za, Aa, Wa, THa, JCa, rea + 2.2, 80)
rel = abs(d60 - d80) / d80
print(r'== [W0] 二分步数自检：锚点 Θmax、ΔI=2.2 ⇒ 60 步 %.10f / 80 步 %.10f，相对差 %.2e ⇒ %s'
      % (d60, d80, rel, 'OK，用 60 步' if rel < 1e-6 else '不合格，退回 80 步'))
BISECT = 60 if rel < 1e-6 else 80

cells = []
FIT_ERR = []
print('\n== [W1] 每株 × 每个定价点：TF2 / FIX / HF2 的 med|e|（拟合窗 ΔI∈[%.2f,%.2f]）==' % (WLO, WHI))
for lab, spec in PLANTS:
    A, B, W, TH, JC, m = make_plant(*spec)
    REXP = r_exp(A)
    _, Vth = eigh(TH)
    rng = np.random.default_rng(20260930 + 7 * PLANTS.index((lab, spec)))
    pool = []
    for rk in (1, 2, min(3, m - 1)):
        for c in range(6):
            Z, _ = np.linalg.qr(rng.standard_normal((m, rk)))
            pool.append((rk, Z))
        pool.append((rk, Vth[:, :rk]))
        pool.append((rk, Vth[:, m - rk:]))
    fits, dropped = [], {'floor': 0, 'nwin': 0, 'fit': 0}
    for rk, Z in pool:
        fl = [curve(Z, s, A, W, TH, JC)[1] - JC for s in FLOORS]
        if (max(fl) - min(fl)) / max(min(fl), 1e-12) > 0.01:
            dropped['floor'] += 1
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
        ok = (di >= WLO) & (di <= WHI)
        if ok.sum() < 5:
            dropped['nwin'] += 1
            continue
        dif, Df = di[ok], Df_all[ok]
        try:
            pT, _ = curve_fit(lambda x, K, g: dfloor + f_exp(x, K, g), dif, Df,
                              p0=[Df.mean() * dif.mean(), GT],
                              bounds=([1e-12, 1e-3], [np.inf, 20.0]), maxfev=60000)
            base = f_exp(dif, 1.0, GT)
            Kfix = float(np.dot(base, Df - dfloor) / np.dot(base, base))
            pP, _ = curve_fit(lambda x, K, nu: dfloor + fpow(x, K, nu), dif, Df,
                              p0=[Df.mean() * dif.mean(), 1.0],
                              bounds=([1e-12, 0.05], [np.inf, 5.0]), maxfev=60000)
        except Exception:
            if not FIT_ERR:
                import traceback
                FIT_ERR.append(traceback.format_exc())
            dropped['fit'] += 1
            continue
        fits.append((rk, Z, dfloor, pT, Kfix, pP))
    nvalid = len(fits)
    di_max = max([np.nanmax([curve(Z, s, A, W, TH, JC)[0] for s in (SGRID[-1], 10 ** 5.0)]) - REXP
                  for _, Z, *_ in fits] or [float('nan')])
    print('\n  [%s] n=%d/%d（地板不收敛 %d，窗内点<5 %d，拟合失败 %d）  Rexp=%.4f  曲线上 ΔI 上界≈%.2f'
          % (lab, nvalid, len(pool), dropped['floor'], dropped['nwin'], dropped['fit'], REXP, di_max))
    if FIT_ERR:
        print('    [首次拟合失败回溯最后一行] ' + FIT_ERR[0].strip().split(chr(10))[-1])
    print('    f=ΔI/0.60   ' + '  '.join('%6.2f' % (t / WHI) for t in TARGETS))
    row = {}
    for k in ('TF2', 'FIX', 'HF2'):
        row[k] = []
    for t in TARGETS:
        tgt = REXP + t
        etf, efix, ehf = [], [], []
        for rk, Z, dfloor, pT, Kfix, pP in fits:
            Dt = true_at(Z, A, W, TH, JC, tgt)
            if Dt is None:
                continue
            etf.append((JC + dfloor + f_exp(t, *pT)) / Dt - 1.0)
            efix.append((JC + dfloor + f_exp(t, Kfix, GT)) / Dt - 1.0)
            ehf.append((JC + dfloor + fpow(t, *pP)) / Dt - 1.0)
        n = len(etf)
        assert len(efix) == n and len(ehf) == n, 'W2 各列 n 不等'
        if n < 5:
            row['TF2'].append(np.nan); row['FIX'].append(np.nan); row['HF2'].append(np.nan)
            row.setdefault('n', []).append((t, n))
            continue
        row['TF2'].append(100 * np.median(np.abs(etf)))
        row['FIX'].append(100 * np.median(np.abs(efix)))
        row['HF2'].append(100 * np.median(np.abs(ehf)))
        cells.append({'plant': lab, 'di': t, 'f': t / WHI, 'n': n,
                      'TF2': 100 * np.median(np.abs(etf)),
                      'FIX': 100 * np.median(np.abs(efix)),
                      'HF2': 100 * np.median(np.abs(ehf)),
                      'sTF2': 100 * np.median(etf), 'sHF2': 100 * np.median(ehf)})
    for k in ('TF2', 'FIX', 'HF2'):
        print('    %-5s     ' % k + '  '.join(('   skip' if np.isnan(v) else '%6.2f' % v) for v in row[k]))
    fin = [(f, a, b) for f, a, b in zip([t / WHI for t in TARGETS], row['TF2'], row['HF2'])
           if not (np.isnan(a) or np.isnan(b))]
    fx = np.inf
    for i in range(len(fin)):
        if fin[i][2] < fin[i][1]:
            if i == 0:
                fx = fin[0][0]
            else:
                f0, d0a, d0b = fin[i - 1]
                f1, d1a, d1b = fin[i]
                g0, g1 = d0b - d0a, d1b - d1a
                fx = f0 if g1 == g0 else f0 + (f1 - f0) * (0 - g0) / (g1 - g0)
            break
    print('    [W3] 交叉外推比 f× = %.2f%s' % (fx, '（扫描范围内不交叉）' if np.isinf(fx) else ''))
    cells_f = [c for c in cells if c['plant'] == lab]
    if len(cells_f) >= 4:
        rho = spearmanr([c['f'] for c in cells_f], [c['TF2'] for c in cells_f]).statistic
        print('    [W5] 株内 Spearman ρ(f, med|e|TF2) = %.3f  (n=%d 格)' % (rho, len(cells_f)))

print('\n== [W4] 预注册判决 ==')
cross = {}
for lab, _ in PLANTS:
    cc = sorted([c for c in cells if c['plant'] == lab], key=lambda c: c['f'])
    fx = np.inf
    for i in range(len(cc)):
        if cc[i]['HF2'] < cc[i]['TF2']:
            if i == 0:
                fx = cc[0]['f']
            else:
                f0, g0 = cc[i - 1]['f'], cc[i - 1]['HF2'] - cc[i - 1]['TF2']
                f1, g1 = cc[i]['f'], cc[i]['HF2'] - cc[i]['TF2']
                fx = f0 if g1 == g0 else f0 + (f1 - f0) * (0 - g0) / (g1 - g0)
            break
    cross[lab] = fx
print('  各株 f×：' + '  '.join('%s=%.2f%s' % (k, v, '(∞)' if np.isinf(v) else '') for k, v in cross.items()))
fin = [v for v in cross.values() if not np.isinf(v)]
ninfl = len(cross) - len(fin)
consist = (len(fin) >= 3 and ninfl == 0 and max(fin) / max(min(fin), 1e-9) <= 2.0 and min(fin) >= 1.2)
print('  有限交叉 %d 株 / 永不交叉 %d 株；max/min=%.2f，min=%.2f ⇒ 判据 [W4]：%s'
      % (len(fin), ninfl, (max(fin) / min(fin)) if fin else float('nan'),
         (min(fin) if fin else float('nan')),
         '一致 ⇒ 可写"指数基在窗外 f× 内占优、更远处退化"' if consist
         else '不一致 ⇒ 正文只许写"两种基都不是长外推的可靠定价器"'))

print('\n== [W5] 池化：外推比是否比"植物身份"更能解释误差 ==')
if cells:
    f_ = np.array([c['f'] for c in cells])
    e_ = np.array([c['TF2'] for c in cells])
    print('  格数 %d（%d 株）；池化 ρ(log f, log med|e|) = %.3f'
          % (len(cells), len({c['plant'] for c in cells}),
             spearmanr(f_, e_).statistic))
    rel_ = []
    for lab in {c['plant'] for c in cells}:
        cc = [c for c in cells if c['plant'] == lab]
        mu = np.median([c['TF2'] for c in cc])
        rel_ += [c['TF2'] / mu for c in cc]
    print('  逐株归一后再池化 ρ(log f, log(误差/该株中位)) = %.3f  ⇒ 外推比能解释"跨株归一后的"漂移吗'
          % spearmanr(f_, np.array(rel_)).statistic)
    rhos = []
    for lab in {c['plant'] for c in cells}:
        cc = [c for c in cells if c['plant'] == lab]
        if len(cc) >= 4:
            rhos.append(spearmanr([c['f'] for c in cc], [c['TF2'] for c in cc]).statistic)
    print('  株内 ρ：%s ⇒ 中位 %.3f（判据 ≥0.8 才升为跨植物陈述）'
          % (np.round(rhos, 3), np.median(rhos) if rhos else float('nan')))

print('\n== 有效域表（TF2 口径，全部格）==')
print('  plant      ΔI    f   n   med|e|TF2%  med|e|HF2%  med|e|FIX%  有符号TF2%')
for c in sorted(cells, key=lambda c: (c['plant'], c['f'])):
    print('  %-9s %5.2f %5.2f %3d    %7.2f    %7.2f    %7.2f    %+7.2f'
          % (c['plant'], c['di'], c['f'], c['n'], c['TF2'], c['HF2'], c['FIX'], c['sTF2']))
