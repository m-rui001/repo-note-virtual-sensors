# -*- coding: utf-8 -*-
r"""实验 134（§83-F-② 预注册）：in-window flag 到底能不能用在**设计选择**上。

背景：§83 过了两支判据（$\rho=0.777$、18 个可行 $\tau$），但 83-C 把它拆开后剩一句诚实话——
flag 预警的是"$\gamma$ 沿 $\Delta I$ 非常数"这个**机制**，不是"这次定价会错"这个**结果**。
机制本身不值钱，值钱的是能不能拿它挡掉一个坏选择。本实验就把定价器放到"选传感方向"这个任务上，
并且**预先写好它可能没用的分支**。

跑前写死的口径：
 [U1] 植物/设计池/拟合口径与 e133 逐字相同（6 株、每株 24 设计、窗 $[0.30,0.60]$、真地板 $+$ 自由 $\gamma$ 的全窗拟合、
      $q=\lvert\ln(\gamma_{\rm hi}/\gamma_{\rm lo})\rvert$ 取自左/右半窗）。定价点取 $\Delta I\in\{1.1,2.2\}$
      （$f=1.83/3.67$，短/长外推各一档），真值 60 步二分。每株每档都要印有效设计数。
 [U2] 三个选择器（在同一批候选里各挑预测代价最小者）：
      ORACLE $=\arg\min D_{\rm true}$；PRICE $=\arg\min D_{\rm pred}$；
      FILTER $=$ 只在 $q\le\tau$ 的候选里 $\arg\min D_{\rm pred}$，$\tau$ 用**其余 5 株 $q$ 的中位**（留一株外推，同 e133 [V5]）。
 [U3] 读数 = regret $=\big(D_{\rm true}(\text{被选中})-D_{\rm true}(\text{ORACLE})\big)/D_{\rm true}(\text{ORACLE})\times100\%$，
      对每个选择器逐株逐档给。另报：FILTER 剩余候选数、ORACLE 自身被标记的株/档、PRICE 命中 ORACLE 的比例。
 [U4] 判据（**两支都为真才算 flag 有用**）：
      ① PRICE 的 regret 不是本来就可忽略：$\mathrm{med}_{{\rm 株,档}}\lvert e_{\rm regret}\rvert\ge1\%$。
        若 $<1\%$ ⇒ **明写"选择层不存在需要 flag 解决的问题"**，不许靠 flag 抢救一个不存在的需求。
      ② 在 $\ge4/6$ 株上，FILTER 的两档 regret 都不高于 PRICE。
      ①②同时成立 ⇒ 正文可写"窗内一致性可作选择前的验收过滤"；否则降级为"诊断量"，正文只许留在 §83 的口径。
 [U5] 附加（不改变 [U4]）：逐株逐档给 $\rho_{\rm Spearman}(D_{\rm pred},D_{\rm true})$ 跨 24 个设计，
      和池化的"预测排序是否与真排序一致"读数——这决定定价误差是**整体平移**（无害）还是**换序**（有害）。
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
TG_DS = [1.1, 2.2]
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
print('== e134：定价器拿去选设计。窗 [0.30,0.60]，定价点 dI in %s ==' % TG_DS)
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
    n_ok = 0
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
        pv, tv = {}, {}
        bad = False
        for t in TG_DS:
            Dt = true_at(Z, A, W, TH, JC, REXP + t)
            if Dt is None:
                drop['truth'] += 1
                bad = True
                break
            pv[t] = JC + dfloor + f_exp(t, K_full, g_full)
            tv[t] = Dt
        if bad:
            continue
        q = abs(np.log(g_hi / g_lo))
        rows.append({'plant': lab, 'rk': rk, 'q': q, 'g': g_full, 'pv': pv, 'tv': tv})
        n_ok += 1
    print('  [%-8s] 有效 %2d/24（地板 %d / 点数 %d / 拟合 %d / 真值 %d 累计剔除）  Rexp=%.4f'
          % (lab, n_ok, drop['floor'], drop['pts'], drop['fit'], drop['truth'], REXP))

LABS = [l for l, _ in PLANTS]
assert len(rows) >= 130, '有效设计太少 ⇒ 先怀疑代码'
print('\n  合计有效设计 %d' % len(rows))

print('\n== [U2][U3] 逐株逐档：三种选择器的 regret（%，相对本株本档 ORACLE）==')
reg = {'PRICE': [], 'FILTER': []}
hit_price, flagged_oracle = [], []
rho_rows = []
for lab in LABS:
    sel = [r for r in rows if r['plant'] == lab]
    if len(sel) < 8:
        print('  %-8s 有效 %d <8 ⇒ 不参与判决' % (lab, len(sel)))
        continue
    for t in TG_DS:
        tau = float(np.median([r['q'] for r in rows if r['plant'] != lab]))
        T = np.array([r['tv'][t] for r in sel]); P = np.array([r['pv'][t] for r in sel])
        Q = np.array([r['q'] for r in sel])
        i_or = int(np.argmin(T)); i_pr = int(np.argmin(P))
        keep = Q <= tau
        okc = [i for i in range(len(sel)) if keep[i]]
        if len(okc) < 8:
            print('  %-8s dI=%.1f  过滤后仅剩 %d 候选（tau=%.3f）⇒ 过滤过度，该格不计数'
                  % (lab, t, len(okc), tau))
            continue
        i_fl = min(okc, key=lambda i: P[i])
        base = T[i_or]
        rp = 100 * (T[i_pr] - base) / base
        rf = 100 * (T[i_fl] - base) / base
        reg['PRICE'].append(rp); reg['FILTER'].append(rf)
        hit_price.append(i_pr == i_or)
        flagged_oracle.append(not bool(keep[i_or]))
        rho = spearmanr(P, T).statistic
        rho_rows.append((lab, t, float(rho)))
        print('  %-8s dI=%.1f tau=%.3f 候选%d/%d  ORACLE rk=%d rank%.1f | PRICE rk=%d regret=%+.2f%% hit=%s | '
              'FILTER regret=%+.2f%% | rho(pred,true)=%.3f'
              % (lab, t, tau, len(okc), len(sel), sel[i_or]['rk'], sel[i_or]['g'],
                 sel[i_pr]['rk'], rp, i_pr == i_or, rf, rho))

RP = np.array(reg['PRICE']); RF = np.array(reg['FILTER'])
med_abs = float(np.median(np.abs(RP)))
print('\n== [U4-①] PRICE 的 regret：med|regret|=%.3f%%  最差=%.3f%%  命中 ORACLE=%.2f (%d/%d)'
      % (med_abs, float(np.max(RP)), np.mean(hit_price), sum(hit_price), len(hit_price)))
print('   ⇒ 选择层%s需要 flag' % ('不存在待解决的问题 ⇒ ' if med_abs < 1.0 else ''))
better = []
per = {}
for i, (lab, t, rho) in enumerate(rho_rows):
    per.setdefault(lab, []).append(i)      # rho_rows 与 RP/RF 同序 append，索引天然对齐
for lab, idx in per.items():
    ok = all(RF[j] <= RP[j] + 1e-12 for j in idx)
    better.append((lab, ok))
    print('   [%s] 两档 regret：PRICE=%s  FILTER=%s  ⇒ %s'
          % (lab, ['%+.2f' % RP[j] for j in idx], ['%+.2f' % RF[j] for j in idx],
             'FILTER 不劣' if ok else 'FILTER 有劣化'))
nb = sum(1 for _, b in better if b)
print('== [U4-②] FILTER 两档 regret 都不高于 PRICE 的株数 = %d / %d（逐株：%s）'
      % (nb, len(better), '  '.join('%s:%s' % (l, 'OK' if b else 'no') for l, b in better)))
print('   ORACLE 自身被 flag 挡掉的格数 = %d / %d' % (sum(flagged_oracle), len(flagged_oracle)))
print('\n== [U4] 判决：%s' % (
    '①②同时成立 ⇒ 可写"窗内一致性作选择前过滤"' if (med_abs >= 1.0 and nb >= 4) else
    '① 不成立（regret 本就 <1%）⇒ 选择层无问题可救，flag 停在诊断量'
    if med_abs < 1.0 else '② 不成立（%d/6 株）⇒ flag 不降低选择风险，正文只许写"诊断"' % nb))

print('\n== [U5] 预测排序 vs 真排序（跨本株 24 设计的 Spearman）==')
for lab, t, rho in rho_rows:
    print('  %-8s dI=%.1f  rho=%.3f' % (lab, t, rho))
rr = np.array([x[2] for x in rho_rows])
print('  中位=%.3f  最差=%.3f  <0.9 的格数=%d/%d' % (np.median(rr), rr.min(), int((rr < 0.9).sum()), len(rr)))
