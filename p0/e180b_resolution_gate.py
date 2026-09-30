# -*- coding: utf-8 -*-
# e180b: 用已落盘的种子噪声底 N=1.53 复核 e180 [P2] 的判决功效。
# 只读 p0/e180_out.txt / p0/e178_out.txt，不重跑搜索。
import re, io, sys
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')
OUT = []
def p(s):
    OUT.append(s)

E180 = open('p0/e180_out.txt', encoding='utf-8').read()

ROW = re.compile(r'^\s*(\d+\.\d+)\s+(\d+\.\d+)\s+(\d)\s+(\d)\s+([\d.]+)\s+([\d.eE+-]+)'
                 r'\s+([\d.eE+-]+)\s+([\d.eE+-]+)\s+([\d.]+)\s+(\d+)/25\s+([\d.eE+-]+)\s*$', re.M)
W = ROW.findall(E180)
assert len(W) == 5, 'e180 [P1] 应 5 行，实际 %d' % len(W)
CELL = [(float(x[0]), float(x[1]), int(x[2]), int(x[3]), float(x[4]), float(x[5]),
         float(x[6]), float(x[7]), float(x[8]), int(x[9]), float(x[10])) for x in W]

m = re.search(r'近阈值.{0,160}?\$c\\in\[([\d.]+),([\d.]+)\]', E180, re.S)
assert m, '带下/上沿没 grep 到'
C_LO, C_HI = float(m.group(1)), float(m.group(2))
assert (C_LO, C_HI) == (116.7, 138.7), (C_LO, C_HI)
N = 1.53  # p0/e177_out.txt [Z0] 落盘的种子噪声底

p('== [Q0] e180 [P1] 五格读数（逐格复解析，条数 %d）==' % len(CELL))
p(r'      D  delta_sw     c   c/C_LO   C_LO/c  c/C_HI   带关系')
for t in CELL:
    D, dsw, rs, r, l1, l3, ratio, dl, c, k, res = t
    p(r'  %6.2f   %5.2f  %6.1f  %7.3f  %6.3f %7.3f   %s' % (
        D, dsw, c, c / C_LO, C_LO / c, c / C_HI,
        '高于上沿' if c > C_HI else (r'低于下沿' if c < C_LO else '带内')))

p('')
p(r'== [Q1] 噪声内/噪声外（判据：偏离带沿 $>N=%.2f$ 才有信息；$\Delta$ 是上界 $\Rightarrow$ 只有低于下沿能否证）==' % N)
for t in CELL:
    D, dsw, rs, r, l1, l3, ratio, dl, c, k, res = t
    if c < C_LO:
        v = r'**噪声内 $\Rightarrow$ 无功效，不许记否证**' if C_LO / c <= N else r'**噪声外 $\Rightarrow$ 真否证**'
        p(r'  %6.2f：$c=%.1f$ 低于下沿 $%.1f$，缺口 $%.3f\times$ $\Rightarrow$ %s' % (D, c, C_LO, C_LO / c, v))
    elif c > C_HI:
        v = '噪声外，可信的上界抬高' if c / C_HI > N else r'噪声内，只算未分辨'
        p(r'  %6.2f：$c=%.1f$ 高于上沿 $%.1f$，抬高 $%.2f\times$ $\Rightarrow$ %s（上界抬高不构成否证）' % (D, c, C_HI, c / C_HI, v))
    else:
        p(r'  %6.2f：$c=%.1f\in[%.1f,%.1f]$ $\Rightarrow$ 带内命中（$c/C_{\rm HI}=%.3f$）' % (D, c, C_LO, C_HI, c / C_HI))

d35 = [t for t in CELL if abs(t[0] - 34.35) < 1e-9][0]
d05 = [t for t in CELL if abs(t[0] - 34.05) < 1e-9][0]
p('')
p('== [Q2] 对 e180 [P2] 的 "$D=34.35$ 否证" 改判 ==')
p(r'   实测缺口 $c/%.1f=%.4f$（$%.1f\%%$），有功效的否证要求 $c<%.1f/%.2f=%.1f$；'
  % (C_LO, d35[8] / C_LO, 100 * (1 - d35[8] / C_LO), C_LO, N, C_LO / N))
p(r'   实际 $c=%.1f$ 的下沿比值 $%.3f\le N=%.2f$ $\Rightarrow$ **我的预注册判据把"噪声内"当成了"否证"：它没带 $N$。**'
  % (d35[8], C_LO / d35[8], N))
p(r'   $\Rightarrow$ 改判：[P2-iii] **作废（无功效）**；本轮有功效的读数只有 $D=32.60$（$c=647.7=4.67\times$ 上沿，噪声外）'
  % ())
p(r'   与带外命中 $D=34.05$（$c=134.5$ 落在 $[116.7,138.7]$ 内）。')

p('')
p('== [Q3] 本轮唯一正向读数：带外样本预测命中（条数 1，不外推成趋势）==')
p(r'   $D=34.05$（$\delta_{\rm sw}=%.2f$，新格、第三条种子流 20261001）：$c=%.1f\in[%.1f,%.1f]$，'
  % (d05[1], d05[8], C_LO, C_HI))
p(r'   离带下沿 $+%.1f\%%$、离带上沿 $-%.1f\%%$ $\Rightarrow$ 带是从 e178 旧格拟合、对新格**外推**的，命中即预测。'
  % (100 * (d05[8] / C_LO - 1), 100 * (1 - d05[8] / C_HI)))

p('')
p(r'== [Q4] 相邻格阶梯：$c$ 随 $\delta_{sw}$ 的单调性（同一批新格、同一种子流，条数 5）==')
L = sorted(CELL, key=lambda t: -t[1])
p(r'   按 $\delta_{sw}$ 降序：' + '  '.join(r'%.2f$-$%.1f' % (t[1], t[8]) for t in L))
viol = [i for i in range(len(L) - 1) if L[i][8] <= L[i + 1][8]]
assert viol == [], 'monotonicity violated at %s' % viol
p('   单调下降违反段：%d 段（$5/5$ 无违反 $\Rightarrow$ 序上单调）' % len(viol))
res_pairs = 0
for i in range(len(L) - 1):
    rr = L[i][8] / L[i + 1][8]
    ok = rr > N
    res_pairs += 1 if ok else 0
    p(r'   $\delta_{sw}$ %.2f$\to$%.2f：$c$ 比 $%.2f\times$ $\Rightarrow$ %s'
      % (L[i][1], L[i + 1][1], rr, '噪声外，可分辨' if ok else '噪声内，不可分辨'))
p('   可分辨相邻步 %d／%d $\Rightarrow$ **单调"序"成立、单调"率"不成立**（只有顶上一跳 $%.2f\\times$ 出噪声）。'
  % (res_pairs, len(L) - 1, L[0][8] / L[1][8]))

p('')
p('== [Q5] 对 e180 落盘两行的更正（同一份读数，不改数据，只改判据）==')
p(r'   e180 [P2] 汇总行"$D=34.35$ 否证"与其后的"**公开降级：$c$ 的带宽不是离切换点距离的单调函数**"'
  % ())
p(r'   两处都超出这份读数的功效：前者缺口 $1.027\le N=1.53$（无功效，[Q1]/[Q2]）；'
  % ())
p(r'   后者方向是反的 $-$ $-$ 五格按 $\delta_{sw}$ 降序的 $c$ 严格下降（$647.7>278.9>183.4>134.5>113.6$，违反 0 段）。')
p(r'   $\Rightarrow$ 撤回这两句。可写的只有：**$c$ 随离切换点距离单调抬高（序级，$5/5$）**；'
  % ())
p(r'   **带上沿之下无法分辨**（步比 $\le1.53$）；**新格 $D=34.35$ 落在带下沿外 $2.7\%$，该读数不足以判"带"是否 plateau**。')

p('')
p('== [Q6] 条数与退出（#46）==')
p('   解析 %d 行／应 5；带内 %d、高于上沿 %d、低于下沿 %d；有功效判决 %d、无功效判决 %d；单调违反 0 段；可分辨步 %d/4；退出 0。'
  % (len(CELL), sum(1 for t in CELL if C_LO <= t[8] <= C_HI),
     sum(1 for t in CELL if t[8] > C_HI), sum(1 for t in CELL if t[8] < C_LO), 2, 1, res_pairs))
p(r'   本脚本不产生新优化读数 $\Rightarrow$ 不含下界/对偶证书（§108 欠账仍在）。')

txt = '\n'.join(OUT) + '\n'
open('p0/e180b_out.txt', 'w', encoding='utf-8').write(txt)
print(txt)
