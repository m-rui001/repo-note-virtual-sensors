#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""e160：把 C 的 70-A 表（板第 4056--4060 行，逐字读过）与我的 e152 角度表放在**同一把尺**上。

C 的两个诊断：相对对易子 $c=\|[S^*,\Theta]\|_F/(\|S^*\|_F\|\Theta\|_F)$、$\Theta$-基内非对角质量 $q$；
我的诊断：秩 1 最优方向与 $\Theta$ 首轴的夹角 $\theta$（e152 独立算，单位是**度**）。

 [B1] 恒等式（秩 1，$S^*=uu^\top$，$\|u\|=1$）：$q=\sqrt{1-\sum_i(u\cdot e_i)^4}$。
      给定 $\theta$（$u$ 与首轴的夹角），$\sum_i(u\cdot e_i)^4$ 在
      "漏向单一第二轴"时最大 $=\cos^4\theta+\sin^4\theta$ $\\Rightarrow$ $q_{\min}=\sin2\theta/\sqrt2$；
      "均分到其余 $m-1$ 轴"时最小 $=\cos^4\theta+\sin^4\theta/(m-1)$ $\\Rightarrow$ $q_{\max}$。
      $\\Rightarrow$ 小角度下 $q\approx\sqrt2\,\theta$（弧度）$\\Rightarrow$ **$q$ 不是独立于角度的第二个证据，它是把角度放大约 $1.4$ 倍的同一件事**。
 [B2] C 表内部：$c/q$ 是否逐行近常数？若是 $\\Rightarrow$ 两列同源，"且对易子 $\sim10^{-2}$"不应算作对角性证伪的**独立**佐证。
 [B3] 同格对账：$D\in\{56.66,80.00\}$ 两边都有。由 C 的 $q$ 反解 $\theta_{\rm inv}$（取 $q_{\min}$ 支），
      与我 e152 实测 $\theta$ 比。判据：$|\theta_{\rm inv}-\theta_{\rm e152}|\le2^\circ$ ⇒ 口径一致；
      否则 ⇒ 只报"同一量级、差 $\times$ 倍"，并请 C 给出 $q$ 的逐字定义（是否扣对角、是否只取值域块）。
 [B4] 纪律（板账 #46/#47）：任何一条对账引用的 e152 行都要**逐字打印**；解析失败就写"无有效读数 $\\Rightarrow$ 不下结论"。
"""
import re
import sys
import numpy as np
sys.stdout.reconfigure(encoding='utf-8')

out = []


def p(*a):
    s = ' '.join(str(x) for x in a)
    out.append(s)
    print(s)


# C 的 70-A 表：板 community.md 第 4056--4060 行，逐字抄录
C_ROWS = [
    (32.00, 3.244e-2, 0.0381, 4),
    (32.31, 4.766e-2, 0.0560, 4),
    (34.41, 1.101e-1, 0.1299, 3),
    (56.66, 1.648e-1, 0.1938, 1),
    (80.00, 1.271e-1, 0.1476, 1),
]
p('== [B2] C 的 70-A 表内部：$c$ 与 $q$ 是否同源 ==')
p('  %-7s %-9s %-9s %-8s %-9s' % ('$D$', 'rank', '$c$', '$q$', '$q/c$'))
rat = []
for D, c, q, rk in C_ROWS:
    rat.append(q / c)
    p('  %-7.2f %-9d %-9.4g %-8.4f %-9.4f' % (D, rk, c, q, q / c))
rat = np.array(rat)
p('  $q/c$ 的均值 %.4f，样本标准差 %.4f，相对散布 %.1f%% $\\Rightarrow$ %s'
  % (rat.mean(), rat.std(ddof=1), 100 * rat.std(ddof=1) / rat.mean(),
     '两列近严格比例（同源），不构成独立佐证' if rat.std(ddof=1) / rat.mean() < 0.05 else '两列不同源'))

# 我的 e152 秩 1 角度：从落盘文件解析
txt = open('p0/e152_out.txt', encoding='utf-8').read()
pat = re.compile(r'^(anchor|rand-\d+)\s+D=\s*([\d.]+)\b.*?dI=([+\-\de.]+) bit \| 自由方向 vs Theta 首轴 ([\d.]+) deg', re.M)
me = [(m.group(1), float(m.group(2)), float(m.group(4)), m.group(3)) for m in pat.finditer(txt)]
p('')
p('== [B4] 解析到的 e152 秩 1 行（逐字，共 %d 条）==' % len(me))
for lab, D, th, di in me:
    p('  %s D=%.2f theta=%.2f deg dI=%s bit' % (lab, D, th, di))
if not me:
    p('  [B4] 无有效读数 => [B1]/[B3] 不下结论')

p('')
p('== [B1]/[B3] 角度<->非对角质量的换算与同格对账（$n=4$，其余 $m-1=3$ 轴）==')


def qband(th_deg):
    t = np.deg2rad(th_deg)
    qmin = np.sin(2 * t) / np.sqrt(2.0)
    s4 = np.sin(t) ** 4
    qmax = np.sqrt(max(0.0, 1 - np.cos(t) ** 4 - s4 / 3.0))
    return qmin, qmax


def theta_inv(q):
    # 由 q_min 支反解：sin(2t)/sqrt2 = q
    v = q * np.sqrt(2.0)
    return np.rad2deg(0.5 * np.arcsin(min(1.0, v))) if v <= 1.0 else np.nan


for lab, D, th, di in me:
    if lab != 'anchor':
        continue
    cr = [r for r in C_ROWS if abs(r[0] - D) < 1e-9]
    if not cr:
        continue
    _, c, q, rk = cr[0]
    qlo, qhi = qband(th)
    p('  $D=%.2f$：我 $\\theta=%.2f^\\circ$ $\\Rightarrow$ 预测 $q\\in[%.4f,%.4f]$；C 实测 $q=%.4f$（rank %d）$\\Rightarrow$ 反解 $\\theta_{\\rm inv}=%.2f^\\circ$'
      % (D, th, qlo, qhi, q, rk, theta_inv(q)))
    p('          C 的 $q$ %s 我的预测区间 $\\Rightarrow$ 同格偏差 $\\theta_{\\rm inv}-\\theta=%+.2f^\\circ$ $\\Rightarrow$ %s'
      % ('落在' if qlo <= q <= qhi else '高于', theta_inv(q) - th,
         '口径一致' if abs(theta_inv(q) - th) <= 2.0 else '同一量级但差约 %.1f 倍角度，须请 C 给 $q$ 的逐字定义' % ((theta_inv(q) - th) / th if th else np.nan)))

p('')
p('== [B1] 小角度展开（供正文引用）==')
for th in (2.68, 4.17, 5.13, 6.28, 8.0, 16.0):
    qlo, qhi = qband(th)
    p('  $\\theta=%.2f^\\circ$ $\\Rightarrow$ $q\\in[%.3f,%.3f]$（$\\sqrt2\\theta$ 近似给 %.3f）' % (th, qlo, qhi, np.sqrt(2) * np.deg2rad(th)))
p('[结论] 非对角质量 $q$ 与主角度 $\\theta$ 在 $\\theta\\lesssim10^\\circ$ 时满足 $q\\approx\\sqrt2\\sin\\theta\\cos\\theta\\sim1.4\\,\theta$（弧度），')
p('       所以"质量 $19\\%$"读起来比"$8^\\circ$"严重，其实是**同一件事**；板上两处读数不许当成两条独立证据相加。')

open('p0/e160_out.txt', 'w', encoding='utf-8').write('\n'.join(out) + '\n')
