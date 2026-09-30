# -*- coding: utf-8 -*-
r"""实验 143（§83-F-1 预注册，跑前写死，跑后不改）：**留一株外推的阈值标定 + ROC**。

动机：§83-B 的 $\tau=0.083$ 是在**全部 6 株**上扫出来的，正文若报它就是"用全集选阈值"。
本节把它换成干净口径：每一株的 $\tau$ 只允许用**其余五株**标定，然后打到留出株上算召回/误标。

复用：e133 的第 1 段（数据构造）用**纯增字段**的方式扩展，$q$、三档误差 $\max|e|$ 的算式一字未动；
[M0] 把重算出的池化 $\rho$ 与逐株 $q$ 分布按同一格式打印，与 `p0/e133_out.txt` 用 diff 对账，
不一致就是 bug 不是结果。

口径（写死）：
 [S1] 设计级样本与 e133 相同；"坏" := 三档 $\max|e|_{\rm TF2}>8\%$（沿用 §82-F 的定义，不重新调线）。
 [S2] 两个基线预测量（同为窗内信息，零额外 DARE 解）：
        gw   $=\lvert\ln(\gamma_{\rm full}/2\ln2)\rvert$（单条全窗拟合偏离 $2\ln2$）；
        rmed $=$ 全窗拟合的相对残差中位 $\mathrm{med}\,\lvert D_{\rm model}-D_{\rm true}\rvert/D_{\rm true}$（窗内点）。
      被检验对象仍是 $q=\lvert\ln(\gamma_{\rm hi}/\gamma_{\rm lo})\rvert$；基线用来问"旗标是否只是换了个说法"。
 [M1] **留出 AUC**（无阈值）：逐株只用留出株自身的 (分数, bad) 算 Mann–Whitney AUC；三个预测量各一条。
 [M2] **留出工作点**：对留出株 $p$，只用其余五株的样本池，按 §82-F 的 [V4-②] 同一条规则取**最小可行** $\tau_p$
      （标记精度 $\ge0.7$ 且训练池误标率 $\le0.15$，候选格 $=$ 训练池 $q$ 的 50–99 分位 26 档）；
      再把 $\tau_p$ 打到 $p$ 上，报召回与误标率。
 [M3] **池化留出 ROC**：每设计的安全裕度 $m=q/\tau_p$（$\tau_p$ 来自它自己的留出标定）合并成一条 ROC，
      报 AUC（梯形）与曲线上的实际工作点。
 [M4] **泄漏对照**：同一批留出株改用全集 $\tau=0.083$（§83-B）再算一遍召回/误标，量化"用全集选阈值"虚高了多少。

判决分支（穷举；先写后跑）：
 (a) 传递成立：召回有定义的留出株中 $\ge4$ 株召回 $\ge0.7$，**且**池化留出误标率 $\le0.15$
     ⇒ 正文可写"$\tau$ 由其余植物标定"，并报留出召回与 ROC。
 (b) 只认极端株：召回 $\ge0.7$ 的株 $\le3$ ⇒ 旗标降级为"只识别 $\gamma$ 剧烈漂移这一类机制"，不得写"阈值可迁移"。
 (c) 标定本身失败：$\ge2$ 株在训练池上**不存在**可行 $\tau_p$ ⇒ 以 (c) 为准写正文（与 (a)(b) 可同命中）。
 (d) 留出株无坏设计（bad$=0$）⇒ 该株召回**无定义**，从 (a)(b) 的分母剔除并单独打印。
 (e) 留出株 AUC$(q)<0.6$ ⇒ 单独列为"无排序力"，不改 (a) 的字面判定，但必须同写。
特异性质检（不改变 a–e，单独一句）：逐株报三者 AUC；若 $q$ 的中位 AUC $\le \max($基线中位 AUC$)+0.05$，
就明写"旗标不比自带的全窗拟合残差/偏离 $2\ln2$ 提供更多窗内信息"。
"""
import sys
import numpy as np
from scipy.stats import spearmanr, mannwhitneyu
sys.stdout.reconfigure(encoding='utf-8')

GT = 2.0 * np.log(2.0)
TAU_FULL = 0.083          # §83-B 的全集阈值，仅用于 [M4] 泄漏对照
SRC = 'p0/e133_inwindow_flag.py'
src = open(SRC, encoding='utf-8').read()

anchor = "        q = abs(np.log(g_hi / g_lo))\n"
assert src.count(anchor) == 1, 'rows.append 前锚点不唯一'
add = ("        Dfit = dfloor + f_exp(di[inw], K_full, g_full)\n"
       "        relres = np.abs(Dfit - Df_all[inw]) / np.maximum(Df_all[inw], 1e-12)\n"
       "        gw = abs(np.log(g_full / GT))\n"
       "        rmed = float(np.median(relres))\n" + anchor)
dictline = "'mx': max(errs), 'e11'"
assert src.count(dictline) == 1, 'rows.append 字典锚点不唯一'
dictnew = "'mx': max(errs), 'gw': gw, 'rmed': rmed, 'e11'"
cut = src.index("Rq = np.array([")
body = src[:cut].replace(anchor, add).replace(dictline, dictnew)
ns = {'__name__': 'e133part', 'GT': GT}
exec(compile(body, SRC + '[part1+diag]', 'exec'), ns)
rows, PLANTS = ns['rows'], ns['PLANTS']
assert all('gw' in r and 'rmed' in r for r in rows), '诊断字段未落到每一行'

print('\n== [M0] 一致性对账：以下各行应与 p0/e133_out.txt 逐字符相同 ==')
Rq = np.array([r['q'] for r in rows]); Rm = np.array([r['mx'] for r in rows])
rho = spearmanr(Rq, Rm).statistic
print('  %-8s n=%2d  q med=%.3f  q p90=%.3f  max|e| med=%.2f%%  内部 ρ=%.3f'
      % ('合计', len(rows), np.median(Rq), np.percentile(Rq, 90), np.median(Rm), rho))
for lab, _ in PLANTS:
    sel = [r for r in rows if r['plant'] == lab]
    if len(sel) < 6:
        continue
    q_ = np.array([r['q'] for r in sel]); m_ = np.array([r['mx'] for r in sel])
    print('  %-8s n=%2d  q med=%.3f  q p90=%.3f  max|e| med=%.2f%%  内部 ρ=%.3f'
          % (lab, len(sel), np.median(q_), np.percentile(q_, 90), np.median(m_),
             spearmanr(q_, m_).statistic))


def badf(rs):
    return np.array([r['mx'] > 8.0 for r in rs])


def auc(score, y):
    s = np.asarray(score, float)
    if y.all() or not y.any():
        return float('nan')
    n1, n0 = int(y.sum()), int((~y).sum())
    return float(mannwhitneyu(s[y], s[~y], alternative='two-sided').statistic) / (n1 * n0)


def tau_on_train(tr):
    qs = np.array([r['q'] for r in tr]); yb = badf(tr)
    feas = []
    for tau in np.unique(np.round(np.percentile(qs, np.linspace(50, 99, 26)), 3)):
        fg = qs > tau
        if fg.sum() == 0:
            continue
        prec = float(yb[fg].mean())
        good = ~yb
        fp = float((qs[good] > tau).mean()) if good.any() else float('nan')
        if prec >= 0.7 and (not np.isfinite(fp) or fp <= 0.15):
            feas.append(float(tau))
    return (min(feas) if feas else None), len(feas)


print('\n== [M1]+[M2] 留一株外推：τ 只由其余五株标定，打到留出株 ==')
per, pool_m, pool_y = [], [], []
pool_y_full, pool_fg_full = [], []
n_c = n_d = 0
for lab, _ in PLANTS:
    te = [r for r in rows if r['plant'] == lab]
    tr = [r for r in rows if r['plant'] != lab]
    if len(te) < 6 or not tr:
        continue
    tp, nfeas = tau_on_train(tr)
    yb = badf(te)
    qq = np.array([r['q'] for r in te])
    a_q = auc(qq, yb); a_gw = auc([r['gw'] for r in te], yb)
    a_rm = auc([r['rmed'] for r in te], yb)
    fl_full = qq > TAU_FULL
    for r, f in zip(te, fl_full):
        pool_y_full.append(r['mx'] > 8.0); pool_fg_full.append(bool(f))
    if tp is None:
        n_c += 1
        per.append((lab, len(te), int(yb.sum()), None, None, None, a_q, a_gw, a_rm))
        print('  %-8s n=%2d 坏=%2d  训练池可行 τ：%d 个 ⇒ **无留出工作点**（分支 c）'
              '   AUC q=%.3f gw=%.3f rmed=%.3f' % (lab, len(te), int(yb.sum()), nfeas, a_q, a_gw, a_rm))
        continue
    fg = qq > tp
    rec = float(yb[fg].sum() / yb.sum()) if yb.any() else None
    fp = float((fg & ~yb).sum() / (~yb).sum()) if (~yb).any() else None
    if rec is None:
        n_d += 1
    # [M4] 泄漏对照：全集 τ 打同一株
    fl_full = qq > TAU_FULL
    rec_full = float(yb[fl_full].sum() / yb.sum()) if yb.any() else None
    fp_full = float((fl_full & ~yb).sum() / (~yb).sum()) if (~yb).any() else None
    per.append((lab, len(te), int(yb.sum()), tp, rec, fp, a_q, a_gw, a_rm))
    for r, f in zip(te, fg):
        pool_m.append(r['q'] / tp); pool_y.append(r['mx'] > 8.0)
    print('  %-8s n=%2d 坏=%2d  τ_留=%.3f  召回=%s  误标率=%s   ‖ 全集τ=%.3f 时召回=%s 误标=%s'
          '   AUC q=%.3f gw=%.3f rmed=%.3f'
          % (lab, len(te), int(yb.sum()), tp,
             '—(无坏设计)' if rec is None else '%.2f(%d/%d)' % (rec, int(yb[fg].sum()), int(yb.sum())),
             '—' if fp is None else '%.2f' % fp, TAU_FULL,
             '—' if rec_full is None else '%.2f' % rec_full,
             '—' if fp_full is None else '%.2f' % fp_full, a_q, a_gw, a_rm))

defined = [p for p in per if p[4] is not None]
passes = sum(1 for p in defined if p[4] >= 0.7)
pool_m = np.array(pool_m); pool_y = np.array(pool_y, bool)
fp_pool = float((pool_m > 1.0)[~pool_y].mean()) if (~pool_y).any() else float('nan')
rec_pool = float((pool_m > 1.0)[pool_y].mean()) if pool_y.any() else float('nan')
auc_pool = auc(pool_m, pool_y)
aqs = np.array([p[6] for p in per]); ags = np.array([p[7] for p in per]); ars = np.array([p[8] for p in per])
print('\n  [M3] 池化留出（裕度 m=q/τ_留，n=%d，坏 %d）：召回=%.2f  误标率=%.2f  AUC(m)=%.3f'
      % (len(pool_m), int(pool_y.sum()), rec_pool, fp_pool, auc_pool))
pyf, pff = np.array(pool_y_full, bool), np.array(pool_fg_full, bool)
print('  [M4] 泄漏对照汇总：留出 τ 召回合计 %.2f（%d/%d 坏设计）；全集 τ=%.3f 在同一批设计上召回合计 %.2f（%d/%d）、'
      '误标率 %.2f ⇒ 阈值泄漏只虚高 %d 个坏设计（%.3f）；留出误标率 %.2f'
      % (rec_pool, int((pool_m > 1.0)[pool_y].sum()), int(pool_y.sum()), TAU_FULL,
         float(pff[pyf].mean()) if pyf.any() else float('nan'),
         int((pff & pyf).sum()), int(pyf.sum()),
         float((pff & ~pyf).sum() / (~pyf).sum()) if (~pyf).any() else float('nan'),
         int((pff & pyf).sum() - (pool_m > 1.0)[pool_y].sum()),
         float(pff[pyf].mean() - rec_pool) if pyf.any() else float('nan'), fp_pool))
print('  AUC 中位：q=%.3f  gw=%.3f  rmed=%.3f   胜场（q 为最大者）：%d/6'
      % (float(np.nanmedian(aqs)), float(np.nanmedian(ags)), float(np.nanmedian(ars)),
         int(sum(1 for i in range(len(aqs)) if aqs[i] >= max(ags[i], ars[i]) - 1e-12))))

x = np.unique(np.round(pool_m, 4))
pts = []
for t in x:
    f = pool_m > t
    if (~pool_y).any() and pool_y.any():
        pts.append((float((f & ~pool_y).sum() / (~pool_y).sum()),
                    float((f & pool_y).sum() / pool_y.sum())))
pts.sort()
auc_roc = float(np.trapezoid([p[1] for p in pts], [p[0] for p in pts])) if len(pts) > 1 else float('nan')
print('  ROC 点数 %d（阈值扫 %s）⇒ 梯形 AUC=%.3f' % (len(pts), 'm=q/τ_留', auc_roc))
print('  无排序力株（AUC(q)<0.6，分支 e）：%s'
      % ([per[i][0] for i in range(len(per)) if aqs[i] < 0.6] or '无'))

print('\n== 判决（预注册分支，穷举）==')
print('  召回有定义的留出株 %d/6，其中召回≥0.7 的 %d 株；无坏设计株 %d；训练池无可行 τ 的株 %d'
      % (len(defined), passes, n_d, n_c))
print('  [V-a] 两联子句逐条打印：passes≥4 ? %s（实际 %d）   池化留出误标率≤0.15 ? %s（实际 %.2f）'
      % (passes >= 4, passes, fp_pool <= 0.15, fp_pool))
if n_c >= 2:
    print('  ⇒ **(c) 标定本身失败**（≥2 株在其余五株上找不到可行 τ）：留一株口径下连阈值都给不出')
elif passes >= 4 and fp_pool <= 0.15:
    print('  ⇒ **(a) 传递成立**：阈值可由其余植物标定，池化留出误标率 %.2f≤0.15' % fp_pool)
else:
    print('  ⇒ **(b) 只认极端株**：召回≥0.7 的仅 %d 株（首条子句不成立）⇒ 旗标不得写成"阈值可迁移"'
          % (passes,))
if float(np.nanmedian(aqs)) <= max(float(np.nanmedian(ags)), float(np.nanmedian(ars))) + 0.05:
    print('  特异性质检：q 的中位 AUC %.3f ≤ 基线中位 AUC %.3f/%.3f + 0.05 ⇒ **旗标不比自带的全窗拟合'
          '残差/偏离 2ln2 提供更多窗内信息**' % (float(np.nanmedian(aqs)), float(np.nanmedian(ags)), float(np.nanmedian(ars))))
else:
    print('  特异性质检：q 的中位 AUC %.3f 高于两基线（%.3f / %.3f）超过 0.05 ⇒ 旗标有独立信息'
          % (float(np.nanmedian(aqs)), float(np.nanmedian(ags)), float(np.nanmedian(ars))))
print('EXIT=0')
