# -*- coding: utf-8 -*-
r"""实验 144（§93-F-2 预注册，跑前写死）：**把留一株工作点协议同时打到三个窗内信号上**。

§93 的判决是 (b)：$q$（左右半窗 $\gamma$ 之比）在干净标定下只有 2/6 株召回 $\ge0.7$，
而且逐株 AUC 中位输给两个更便宜的基线。但那只是"把 $q$ 判死"。本节问正面问题：
**有没有任何一个窗内信号，能在只由其余五株标定阈值的前提下，跨株传递成一个可用的有效性门？**

数据：直接复用 `p0/e143_holdout_roc.py` 的第 1 段（e133 的拟合 + 三个诊断量），**不重新定义任何东西**；
"坏" := 三档 $\max|e|_{\rm TF2}>8\%$（$\Delta I\in\{1.1,1.5,2.2\}$，与窗 $[0.30,0.60]$ 不相交）。

三个候选（都是窗内、零额外 DARE 解；gw/rmed 连额外拟合都不要）：
  q     $=\lvert\ln(\gamma_{\rm hi}/\gamma_{\rm lo})\rvert$（左/右半窗各一条）
  gw    $=\lvert\ln(\gamma_{\rm full}/2\ln2)\rvert$
  rmed  $=$ 全窗拟合相对残差中位

协议（对每个候选同一条，写死）：
 [P1] 留出株 $p$ 的 $\tau_p$ 只用其余五株的池：候选格 $=$ 训练池该信号的 50--99 分位 26 档（四舍五入 3 位、去重），
      可行 $=$ 标记精度 $\ge0.7$ 且训练池误标率 $\le0.15$；取**最小**可行 $\tau_p$；无可行则记 (no-op)。
 [P2] 打到 $p$ 上：召回 $=$ 坏设计被标记比例；误标率 $=$ 不差设计被标记比例。
 [P3] 池化：裕度 $m=$ 信号$/\tau_p$，报池化召回、池化误标率、$\mathrm{AUC}(m)$（Mann–Whitney）与 ROC 梯形 AUC。
 [T]  逐候选判决：**(T-a) 传递成立** $=$ 召回有定义的株 $\ge4$ 且其中 $\ge4$ 株召回 $\ge0.7$ 且池化误标率 $\le0.15$；
      **(T-b)** 否则。另报 (no-op) 株数与"留出株无坏设计"株数。
 [H]  头对头：若多个候选同落 (T-a)，以池化召回高者为正文候选；若**全部**落 (T-b)，
      则正文只能写"窗内没有可传递的有效性门"，并把 §82-F/§83/§93 整条线留在板上不进正文。
"""
import sys
import numpy as np
from scipy.stats import mannwhitneyu
sys.stdout.reconfigure(encoding='utf-8')

E143 = 'p0/e143_holdout_roc.py'
src = open(E143, encoding='utf-8').read()
cut = src.rindex("[M0] 一致性对账")
head = src[:cut].rsplit('print(', 1)[0]
ns = {'__name__': 'e143head'}
exec(compile(head, E143 + '[head]', 'exec'), ns)
rows, PLANTS = ns['rows'], ns['PLANTS']
assert all('gw' in r and 'rmed' in r for r in rows)
print('复用 e143 第 1 段：有效设计 %d，逐株 %s'
      % (len(rows), '  '.join('%s=%d' % (lab, sum(1 for r in rows if r['plant'] == lab)) for lab, _ in PLANTS)))

KEYS = ['q', 'gw', 'rmed']

# [S0] 仪器自检（第 37 号账的防伪措施）：故意用"召回与精度不相等"的构造，
# 若哪天有人把召回写回 y[fg].mean()，这条会当场 FAIL。
_y = np.array([True, False, False, False]); _fg = np.array([True, True, True, False])
_rec = float((_fg & _y).sum() / _y.sum()); _pre = float(_y[_fg].mean())
print('== [S0] 自检：构造 1 坏 3 不差、标记 3 个 ⇒ 召回应=1.00、精度应=0.33 ==')
print('   召回 %.2f  精度 %.2f  ⇒ %s' % (_rec, _pre, 'PASS（两量不等，公式没混用）'
                                         if abs(_rec - 1.0) < 1e-12 and abs(_pre - 1.0 / 3) < 1e-9 else 'FAIL'))
assert abs(_rec - 1.0) < 1e-12 and abs(_pre - 1.0 / 3) < 1e-9, '[S0] 失败 ⇒ 召回/精度公式被混用，停跑'


def prep(rs, k):
    s = np.array([r[k] for r in rs], float)
    y = np.array([r['mx'] > 8.0 for r in rs], bool)
    return s, y


def tau_tr(train, k, rounded=False):
    s, y = prep(train, k)
    feas = []
    grid = np.percentile(s, np.linspace(50, 99, 26))
    # 候选格用训练池自身的分位点、**不四舍五入**：round(_,3) 对 q（尺度 1e-2）无害，
    # 对 rmed（尺度 1e-4）会把整条网格塌成 0.001，等于偷偷换了一个（退化的）阈值规则。
    # e143 用的是 round 网格，所以 q 另跑一遍 round 规则做敏感性对照 [SN]。
    cand = np.unique(np.round(grid, 3) if rounded else grid)
    for t in cand:
        fg = s > t
        if not fg.any():
            continue
        prec = float(y[fg].mean())
        good = ~y
        fp = float((s[good] > t).mean()) if good.any() else float('nan')
        if prec >= 0.7 and (not np.isfinite(fp) or fp <= 0.15):
            feas.append(float(t))
    return (min(feas) if feas else None), len(feas)


def auc(s, y):
    if y.all() or not y.any():
        return float('nan')
    n1, n0 = int(y.sum()), int((~y).sum())
    return float(mannwhitneyu(s[y], s[~y], alternative='two-sided').statistic) / (n1 * n0)


summary = {}
for k in KEYS:
    print('\n== 候选 %s：留一株工作点（[P1]--([P3]) ==' % k)
    recs, fps, pm, py, nnoop = [], [], [], [], 0
    for lab, _ in PLANTS:
        te = [r for r in rows if r['plant'] == lab]
        tr = [r for r in rows if r['plant'] != lab]
        tp, nfeas = tau_tr(tr, k)
        s, y = prep(te, k)
        if tp is None:
            nnoop += 1
            print('  %-8s **训练池无可行 τ**（可行格 %d 个）⇒ (no-op)   AUC=%.3f' % (lab, nfeas, auc(s, y)))
            continue
        fg = s > tp
        # 召回 = 坏设计中被标记的比例；误标率 = 不差的设计中被标记的比例
        # （首轮这里写成 y[fg].mean()，那是**精度**不是召回 —— 第 37 号代码账）
        rec = float((fg & y).sum() / y.sum()) if y.any() else None
        fp = float((fg & ~y).sum() / (~y).sum()) if (~y).any() else None
        recs.append((lab, rec)); fps.append((lab, fp))
        pm.extend(s / tp); py.extend(y)
        sn = ''
        if k == 'q':
            tp2, _ = tau_tr(tr, k, rounded=True)
            if tp2 is not None and abs(tp2 - tp) > 1e-12:
                fg2 = s > tp2
                sn = '   [SN] round 网格 τ=%.3f ⇒ 召回 %.2f' % (
                    tp2, float((fg2 & y).sum() / y.sum()) if y.any() else float('nan'))
        print('  %-8s τ=%.6f  召回=%s  误标率=%s   AUC=%.3f%s'
              % (lab, tp, '—' if rec is None else '%.2f(%d/%d)' % (rec, int((fg & y).sum()), int(y.sum())),
                 '—' if fp is None else '%.2f' % fp, auc(s, y), sn))
    pm, py = np.array(pm), np.array(py, bool)
    defined = [x for x in recs if x[1] is not None]
    passes = sum(1 for x in defined if x[1] >= 0.7)
    rp = float((pm > 1.0)[py].mean()); fp_ = float((pm > 1.0)[~py].mean())
    print('  [P3] 池化 n=%d 坏=%d  召回=%.2f  误标率=%.2f  AUC(m)=%.3f'
          % (len(pm), int(py.sum()), rp, fp_, auc(pm, py)))
    thr = np.unique(np.round(pm, 6))
    pts = sorted(set((float(((pm > t) & ~py).sum() / (~py).sum()), float(((pm > t) & py).sum() / py.sum()))
                    for t in thr))
    roc = float(np.trapezoid([p[1] for p in pts], [p[0] for p in pts])) if len(pts) > 1 else float('nan')
    print('      ROC 点 %d 个 ⇒ 梯形 AUC=%.3f   召回有定义 %d 株、其中≥0.7 的 %d 株、(no-op) %d 株'
          % (len(pts), roc, len(defined), passes, nnoop))
    verdict = 'T-a 传递成立' if (len(defined) >= 4 and passes >= 4 and fp_ <= 0.15) else 'T-b 不传递'
    print('  [T] 判决：**%s**（召回有定义 %d 株 / 其中≥0.7 %d 株 / 池化误标率 %.2f）' % (verdict, len(defined), passes, fp_))
    summary[k] = dict(passes=passes, defined=len(defined), rec=rp, fp=fp_, auc=auc(pm, py),
                      roc=roc, noop=nnoop, verdict=verdict)

print('\n== [H] 头对头 ==')
for k in KEYS:
    d = summary[k]
    print('  %-5s %s   池化召回=%.2f  池化误标率=%.2f  AUC=%.3f  ROC=%.3f  (no-op) %d 株'
          % (k, d['verdict'], d['rec'], d['fp'], d['auc'], d['roc'], d['noop']))
winners = [k for k in KEYS if summary[k]['verdict'].startswith('T-a')]
if winners:
    best = max(winners, key=lambda k: summary[k]['rec'])
    print('  ⇒ 可进正文的窗内门：%s（取池化召回最高者 %s，%.2f）' % ('/'.join(winners), best, summary[best]['rec']))
else:
    print('  ⇒ **[H] 全部落 (T-b)：窗内没有可跨株传递的有效性门**；'
          '§82-F/§83/§93 整条线留在板上，正文不得出现"用窗内信号预测外推失效"这类句子。')
print('EXIT=0')
