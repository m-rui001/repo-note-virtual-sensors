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
 [U6] 新增（e134b）：除 $\rho$ 外，直接报**错序率**（成对不一致比例）与 $D_{\rm pred}/D_{\rm true}$ 的
      中位数与跨距。若跨距 $\ll$ 该列误差本身，说明误差是近似**同乘子**（选择无关），这才配叫"解释"。

本版是 e134 的**修正重跑**，两笔账登记在此（跑前不隐瞒）：
 (a) e134 第 1 版把"过滤后候选 $<8$"的格子**整格跳过**，连带把该格的 PRICE 也一起丢了 ⇒ `rand-1`（唯一
     定价真坏的那株）**24/24 设计全被评估、但两档一格都没进判决**。这属于 §81-E 同源的"悄悄减样本"。
     e134b 改成：ORACLE/PRICE 无条件评估，只有 FILTER 受候选数门槛。
 (b) e134 第 1 版判决行的打印文字自相矛盾（"不存在待解决的问题 ⇒ 需要 flag"），是 §83-D 第 19 号账的同类。
 (c) **[U7] 是跑后追加的读数，不是预注册判据**（第 1 版跑完 12/12 命中后我才想到要问"这是不是白捡的"）。
     它不参与 [U4] 判决，只用于给"命中率"定价；追加登记在此，不许事后改写成跑前就写好。
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
RPc, RFc, ratio_rows, risk_rows = [], [], [], []
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
        base = T[i_or]
        rp = 100 * (T[i_pr] - base) / base
        if len(okc) >= 8:
            i_fl = min(okc, key=lambda i: P[i])
            rf = 100 * (T[i_fl] - base) / base
            over = ''
        else:
            rf = None
            over = ' [过滤过度：剩 %d 候选<8 ⇒ FILTER 不计数，PRICE 照算]' % len(okc)
        RPc.append(rp); RFc.append(rf)
        hit_price.append(i_pr == i_or)
        flagged_oracle.append((not bool(keep[i_or])) if rf is not None else None)
        Rat = P / T
        dP = P[:, None] - P[None, :]; dT = T[:, None] - T[None, :]
        npair = len(T) * (len(T) - 1) // 2
        inv = 100.0 * int(np.triu((dP * dT) < 0, 1).sum()) / npair
        rho = spearmanr(P, T).statistic
        # [U7] 风险读数：预测把多少个"其实更贵"的设计排到了 ORACLE 之下，其中最差的真实 regret 是多少
        below = [i for i in range(len(sel)) if P[i] <= P[i_or] and i != i_or]
        worst_below = max([100 * (T[i] - base) / base for i in below] or [0.0])
        gap2 = 100 * (np.sort(T)[1] - base) / base
        ovl3 = len(set(np.argsort(P)[:3]).intersection(set(np.argsort(T)[:3]))) / 3.0
        risk_rows.append((lab, t, len(below), worst_below, gap2, ovl3))
        rho_rows.append((lab, t, float(rho)))
        ratio_rows.append((lab, t, float(np.median(Rat)), float(np.percentile(Rat, 90) - np.percentile(Rat, 10)),
                           float(Rat.max() - Rat.min()), inv))
        reg['PRICE'].append(rp)
        if rf is not None:
            reg['FILTER'].append(rf)
        print('  %-8s dI=%.1f tau=%.3f 候选%2d/%d  ORACLE rk=%d g=%.2f | PRICE regret=%+.2f%% hit=%-5s | '
              'FILTER regret=%8s | rho(pred,true)=%.4f 错序率=%5.2f%% | pred/true med=%.4f 跨距=%.4f%s'
              % (lab, t, tau, len(okc), len(sel), sel[i_or]['rk'], sel[i_or]['g'],
                 rp, str(i_pr == i_or), ('%+.2f%%' % rf) if rf is not None else 'n/a',
                 rho, inv, np.median(Rat), Rat.max() - Rat.min(), over))

RP = np.array(RPc)
med_abs = float(np.median(np.abs(RP)))
print('\n== [U4-①] PRICE 的 regret（全 %d 格，含 rand-1）：med|regret|=%.3f%%  最差=%+.3f%%  命中 ORACLE=%.2f (%d/%d)'
      % (len(RPc), med_abs, float(np.max(RP)), np.mean(hit_price), int(sum(hit_price)), len(hit_price)))
print('   ⇒ %s' % ('regret 中位 <1% ⇒ **选择层不存在需要 flag 来解决的问题**（预注册的 ① 否证分支）'
                   if med_abs < 1.0 else 'PRICE 的选择确实被定价误差伤到 ⇒ 再看 ②'))

per = {}
for i, (lab, t, rho) in enumerate(rho_rows):
    per.setdefault(lab, []).append(i)      # rho_rows 与 RPc/RFc 同序 append ⇒ 索引天然对齐
better = []
for lab, idx in per.items():
    use = [j for j in idx if RFc[j] is not None]
    if not use:
        better.append((lab, None))
        print('   [%s] PRICE regret=%s  FILTER=无可评价格（过滤过度）' % (lab, ['%+.2f' % RPc[j] for j in idx]))
        continue
    ok = all(RFc[j] <= RPc[j] + 1e-12 for j in use)
    better.append((lab, ok))
    print('   [%s] PRICE regret=%s  FILTER regret=%s ⇒ %s'
          % (lab, ['%+.2f' % RPc[j] for j in idx], ['%+.2f' % RFc[j] for j in idx],
             'FILTER 不劣' if ok else 'FILTER 有劣化'))
nb = sum(1 for _, b in better if b is True)
n_eval = sum(1 for _, b in better if b is not None)
print('== [U4-②] FILTER 各档 regret 都不高于 PRICE 的株数 = %d / 可评 %d（全 %d 株）'
      % (nb, n_eval, len(better)))
fo = [x for x in flagged_oracle if x is not None]
print('   ORACLE 自身被 flag 挡掉的格数 = %d / %d（被挡 ⇒ FILTER 必然付出代价）' % (sum(fo), len(fo)))
print('\n== [U4] 判决：%s' % (
    '①②同时成立 ⇒ 可写"窗内一致性作选择前过滤"' if (med_abs >= 1.0 and nb >= 4) else
    ('① 不成立 ⇒ flag 停在诊断量，不许写成选择器' if med_abs < 1.0
     else '② 不成立（%d 株）⇒ flag 不降低选择风险' % nb)))

print('\n== [U5] 预测排序 vs 真排序（跨本株 24 设计）==')
for (lab, t, rho), rr in zip(rho_rows, ratio_rows):
    print('  %-8s dI=%.1f  rho=%.4f  错序率=%5.2f%%  pred/true med=%.4f  p90-p10=%.4f  跨距=%.4f'
          % (lab, t, rho, rr[5], rr[2], rr[3], rr[4]))
rr_ = np.array([x[2] for x in rho_rows]); iv = np.array([x[5] for x in ratio_rows])
sp = np.array([x[4] for x in ratio_rows])
print('  rho 中位=%.4f 最差=%.4f | 错序率 中位=%.2f%% 最差=%.2f%% | pred/true 跨距 中位=%.4f 最差=%.4f'
      % (np.median(rr_), rr_.min(), np.median(iv), iv.max(), np.median(sp), sp.max()))

print('\n== [U7] 12/12 是不是白捡的：预测把 ORACLE 排到第几名之下、以及那些设计真贵多少 ==')
print('   below := 预测价 <= ORACLE 预测价 的其它设计数；worst := 这些里最差的真实 regret；')
print('   gap2 := 真值第 2 便宜与第 1 的间距；ovl3 := 预测 top-3 与真 top-3 的重合率')
for lab, t, nb, wb, g2, o3 in risk_rows:
    print('  %-8s dI=%.1f  below=%2d  worst_below=%+7.2f%%  gap2=%6.2f%%  top3重合=%.2f'
          % (lab, t, nb, wb, g2, o3))
wb_ = np.array([x[3] for x in risk_rows]); nb_ = np.array([x[2] for x in risk_rows])
g2_ = np.array([x[4] for x in risk_rows]); o3_ = np.array([x[5] for x in risk_rows])
print('  ⇒ worst_below 最大=%.2f%%（中位 %.2f%%）；below 合计 %d，均值 %.1f/格；'
      'gap2 中位=%.2f%% 最小=%.2f%%；top3 重合 中位=%.2f 最小=%.2f'
      % (wb_.max(), np.median(wb_), int(nb_.sum()), nb_.mean(), np.median(g2_), g2_.min(),
         np.median(o3_), o3_.min()))
print('  [U7] 口径：若 below 普遍 >0 且 worst_below 明显大于 med|regret| ⇒ "12/12 命中"含运气成分，'
      '正文只许写 regret 的**分布上界**而不是命中率。')
