# -*- coding: utf-8 -*-
r"""实验 145（跨车道对账，回答 C 的 §58-D-3）：**局域斜率比值**到底是多少，逐株。

C 的 58-D-3 说："真实局域斜率的比值在 $\Delta I{=}1$ 只有 $1.26/1.16$，而 $\gamma_{\rm fit}$ 的比值是 $1.53/1.76$"。
这句话针对的是**拟合货币**，但它顺手给了我车道里 $\gamma_{\rm loc}$ 的一个数。本脚本不重算任何拟合，
只从**已落盘**的 `p0/e137b_out.txt`（e137 的局域中心差分，6 株 $\times$ 秩 1--3）里把 $\Delta I=1.00$ 那一档读出来，
逐株印 $\gamma_1,\gamma_2,\gamma_3$ 与两种读法（$>$1 的"陡度比" $\gamma_1/\gamma_2,\gamma_2/\gamma_3$；
$<$1 的"比值" $\gamma_2/\gamma_1,\gamma_3/\gamma_1$ 对 $1/2,1/3$），好让两边的数能在同一把尺上对。

预注册（跑前写死）：
 [J1] 若我的陡度比中位落在 C 报的 $1.26/1.16$ 的 $10\\%$ 内 $\\Rightarrow$ **接受 C 的 58-D-3 与我车道一致**，
      并把它当作"窗位/弱识别"论证的第二条独立支持；
 [J2] 若差 $>10\\%$ $\\Rightarrow$ 报口径差（池、$n$/$\Delta I$ 网格、地板），**不引用 C 那句**；
 [J3] 无论 J1/J2：打印 $\Delta I=16$ 档 rank-1 的符号，与 rank-1 在 e137b 里"最后一档有限值"的档号，
      用于向 C 提出 c78 覆盖问题（C 的 rank-1 表止于 $\Delta I=6.0$，而它的 58-E-3 说 $x_{\\min}\\sim\\mathcal O(10^{-2})$
      出现在 $\Delta I\\approx8.6$、且 $\Delta I$ 上限秩中位是 $14.6/28.7/41.9$）。
"""
import re
import sys
import numpy as np
sys.stdout.reconfigure(encoding='utf-8')

T = open('p0/e137b_out.txt', encoding='utf-8').read().split('\n')
DG = [0.50, 0.75, 1.00, 1.50, 2.00, 3.00, 4.00, 6.00, 8.00, 12.00, 16.00]
blocks, cur = [], None
for L in T:
    m = re.match(r'^== \[(.+?)\] n=(\d+) Rexp=([\d.]+)', L)
    if m:
        cur = {'lab': m.group(1), 'rows': {}}
        blocks.append(cur)
        continue
    if cur is not None:
        m2 = re.match(r'^\s+(\d)\s+(\d)\s+(.+)$', L)
        if m2:
            vals = []
            for tok in m2.group(3).split():
                if tok == '--':
                    vals.append(None)
                    continue
                try:
                    vals.append(float(tok))
                except ValueError:
                    break          # 行尾的"覆盖到 ΔI≤…"注记：到此为止，保持与 DG 的档号对齐
            cur['rows'][int(m2.group(1))] = vals

print('== [J1] 逐株 ΔI=1.00 的局域斜率与两种比值（源：p0/e137b_out.txt，中位=每格 6 个设计）==')
steep12, steep23, rat21, rat31 = [], [], [], []
for b in blocks:
    g = {r: b['rows'].get(r) or [] for r in (1, 2, 3)}
    if not all(len(g[r]) > DG.index(1.00) and g[r][DG.index(1.00)] is not None for r in (1, 2, 3)):
        print('  %-8s ΔI=1.00 缺档 ⇒ 不计入' % b['lab'])
        continue
    g1, g2, g3 = (g[r][DG.index(1.00)] for r in (1, 2, 3))
    s12, s23, r21, r31 = g1 / g2, g2 / g3, g2 / g1, g3 / g1
    steep12.append(s12); steep23.append(s23); rat21.append(r21); rat31.append(r31)
    last = max(DG[i] for i, v in enumerate(g[1]) if i < len(g[1]) and v is not None)
    v16 = g[1][DG.index(16.00)] if len(g[1]) > DG.index(16.00) else None
    print('  %-8s γ=%.3f/%.3f/%.3f   γ1/γ2=%.3f γ2/γ3=%.3f   γ2/γ1=%.3f(对1/2 余隙%+.3f) '
          'γ3/γ1=%.3f(对1/3 余隙%+.3f)   [J3] rank1 最后有限档=%.2f  ΔI=16 值=%s'
          % (b['lab'], g1, g2, g3, s12, s23, r21, r21 - 0.5, r31, r31 - 1.0 / 3, last,
             '缺失' if v16 is None else '%+.3f' % v16))
m12, m23 = float(np.median(steep12)), float(np.median(steep23))
print('  中位陡度比 γ1/γ2=%.3f  γ2/γ3=%.3f   （C 的 58-D-3 报 1.26 / 1.16）' % (m12, m23))
d12, d23 = abs(m12 - 1.26) / 1.26 * 100, abs(m23 - 1.16) / 1.16 * 100
print('  与 C 的相对差：γ1/γ2 %.1f%%、γ2/γ3 %.1f%% ⇒ %s'
      % (d12, d23, '[J1] 接受：两边在同一把尺上（≤10%）' if max(d12, d23) <= 10.0
         else '[J2] 差 >10% ⇒ 报口径差，不引用 C 那句'))
print('  中位比值 γ2/γ1=%.3f  γ3/γ1=%.3f；γ_fit 侧（C 的 §85-B 我车道值）1.53/1.76 是**拟合货币**，不可与上面互换'
      % (float(np.median(rat21)), float(np.median(rat31))))
print('EXIT=0')
