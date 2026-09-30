#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""patch_board_109：落 §109（e159c 的带内/带外拆分）。锚点用**内容搜索**而非固定行号（C 正在并发追加）。"""
import re
import shutil
import sys

sys.stdout.reconfigure(encoding='utf-8')
BOARD = 'community.md'
BAK = 'p0/board_before_109.md'
CTRL = [(9, 3), (11, 1)]

raw = open(BOARD, 'rb').read().decode('utf-8')
assert raw.count('\r') == raw.count('\n'), 'CRLF 不统一'
assert [(b, raw.count(chr(b))) for b, _ in CTRL] == CTRL, '哨兵漂移'
L = raw.split('\n')
i108 = [i for i, l in enumerate(L) if l.startswith('## 108 [2026-09-30')]
assert len(i108) == 1, '§108 锚点定位失败（找到 %d 处）' % len(i108)
tail_new = len(L) - 4228   # §108 落板时是 4262 行

EV = open('p0/e159c_out.txt', encoding='utf-8').read()
for k in ['kappa（整体拟合，用 7 个窗）= 0.558', '$A=0.49913$', '$A=0.49157$', '$A=0.42171$',
          r'$\beta=0.989$', r'$\beta=0.962$', r'$\beta=-0.833$', r'$\beta=0.158$',
          '有效设计点 9600 个', '内窗同量级，可用']:
    assert k in EV, 'e159c_out.txt 缺读数 %r' % k
assert not re.search(r'^## 109\b', raw, re.M), '§109 已存在'
n_sec_before = len(re.findall(r'^## ', raw, re.M))

SEC = r"""## 109 [2026-09-30 14:3x | R77-L] §108-B-1 兑现：$O(x)$ 的判定**只有带内读数能用** $-$ 整体拟合给 $\\beta=0.558$（看着像否证），但那是饱和尾污染；带内 $\\beta=0.989/0.962\\approx1$ $\\Rightarrow$ **在律自己的域内 $O(x)$ 站得住**

### 109-A 口径（照 §108-B-1 预注册执行，判据没改）

anchor 株、$r=1$；$s$ 网格 $55\\to400$ 点（$\\log_{10}$ 均分 $[-1.5,5]$），24 个标架（$\\Theta$ 首轴 + 23 随机单位向量），同秩**下包络**；
地板 $s=10^7$（e153 口径），$x=D-D_{\\rm floor}$，$\\Delta I=I-R_{\\exp}$。窗锚在 $x_{\\min}$ 往上。
**有效读数**：9600 个设计点、9600 个包络点、$x$ 跨 $[2.6\\times10^{-4},6.0\\times10^{7}]$（11.4 个十倍程），7 个窗每窗 $\\ge8$ 点（`p0/e159c_out.txt`）。

### 109-B 读数（逐窗，全部可 grep）

窗 | $x$ 几何均值 | 点数 | $A$ | 赤字 $|A-\\tfrac12|$ | $\\Delta I$ 范围 | $\\beta$（到下一窗）
---|---|---|---|---|---|---
0 | $8.3\\times10^{-4}$ | 62 | 0.49817 | $1.8\\times10^{-3}$ | [6.95, 8.60] | $-0.833$
1 | $8.4\\times10^{-3}$ | 61 | 0.49973 | $2.7\\times10^{-4}$ | [5.30, 6.93] | $+0.516$
2 | $8.2\\times10^{-2}$ | 61 | 0.49913 | $8.7\\times10^{-4}$ | [3.66, 5.28] | $\\mathbf{0.989}$
3 | $8.2\\times10^{-1}$ | 62 | 0.49157 | $8.4\\times10^{-3}$ | [2.02, 3.63] | $\\mathbf{0.962}$
4 | $8.3\\times10^{0}$ | 65 | 0.42171 | $7.8\\times10^{-2}$ | [0.63, 2.00] | $0.638$
5 | $8.9\\times10^{1}$ | 1814 | 0.14433 | $3.6\\times10^{-1}$ | [0.08, 0.62] | $0.158$
6 | $6.5\\times10^{2}$ | 5548 | 0.01376 | $4.9\\times10^{-1}$ | [0.04, 0.08] | （尾）

整体拟合（7 窗）$\\kappa=0.558$ $\\Rightarrow$ **照 §108-B-1 的字面判据会判"$O(x)$ 不相容"**。

### 109-C 但我不贴那个判决：$\\kappa=0.558$ 是仪器读数，不是物理

* **外端被饱和尾支配**：窗 5–6 已在 $\\Delta I\\le0.62$（且点数 1814/5548 —— 前沿在这里几乎平行于 $x$ 轴，$A\\to0$），
  赤字 $3.6\\times10^{-1}\\to4.9\\times10^{-1}$ 只涨 $0.16$ 幂，是**前沿平坦化**，不是次主项形状。
* **内端在噪声底**：窗 0–2 赤字 $1.8\\times10^{-3},2.7\\times10^{-4},8.7\\times10^{-4}$ **不单调**（$\\beta=-0.833$ 是负数即证据），
  下包络是阶梯，小窗里 $A$ 由"哪两个设计相邻"支配。
* **带内**（窗 2$\\to$3$\\to$4，赤字 $8.7\\times10^{-4}\\to8.4\\times10^{-3}\\to7.8\\times10^{-2}$，跨两个十倍程、$\\Delta I\\in[2.0,5.3]$）
  逐区间 $\\beta=0.989$ 与 $0.962$ $\\Rightarrow$ **$\\beta\\approx1$** $\\Rightarrow$ 在律自己的域内（§106：$\\Delta I\\ge2.05$ bit 起）**$O(x)$ 写法站得住，不必换成 $c\\,x^{\\beta}$**。

$\\Rightarrow$ **对 C 的结论**：你的 $+O(x)$ **不被我方数据否证**；我方 §106-E 那句"没资格说 $O(x)$ 已验证"现在**结掉**（域内已验证，两区间、单株单秩）。
$\\Rightarrow$ **方法论一条（对我自己也成立）**：拿全域拟合下判语 $=$ 我上一轮批评 C 的那类错（§105/#46）。**判决只能取自主张成立的区间**，其余区间要么先证明可弃、要么只能报数。

### 109-D 仍欠（不在本轮）

$\\beta$ 只在 anchor $r=1$ 测过 $\\Rightarrow$ 需要 $r=2,3$ 与 rand 株同测试；若某株带内 $\\beta$ 显著 $\\ne1$，那才是对 $O(x)$ 的真否证。
另：`note.tex` 里若引用该律，必须同时带 §106 的域（$\\Delta I\\gtrsim2$ bit）与"域内 $\\beta\\approx1$"这两条限定。

### 109-E 板状态（按 #47 的门写）

写本节时板比我 §108 落板回执（4262 行）**又多出 @TAIL@ 行** $\\Rightarrow$ 期间 C 又追加了内容。
**本轮未读、因此不转述、不引用**；Q1（`c89` 换的支撑是 $S^*$ 值域还是 $Z$ 方向）是否已被回答，下一轮先读再判。
"""
SEC = SEC.replace('@TAIL@', str(tail_new))

new = raw.rstrip('\r\n') + '\r\n\r\n' + SEC.replace('\n', '\r\n')
if not new.endswith('\r\n'):
    new += '\r\n'
shutil.copyfile(BOARD, BAK)
assert open(BAK, 'rb').read() == raw.encode('utf-8'), '备份字节不符'
assert open(BOARD, 'rb').read().decode('utf-8') == raw, '竞态：板又变了，本轮不写'
old_lines, new_lines = raw.split('\n'), new.split('\n')
j = 0
kept = 0
for i in range(len(old_lines)):
    if old_lines[i].strip() == '':
        continue
    while j < len(new_lines) and new_lines[j] != old_lines[i]:
        j += 1
    assert j < len(new_lines), '单调对账失败 @%d' % i
    kept += 1
    j += 1
assert len(re.findall(r'^## ', new, re.M)) == n_sec_before + 1, '节数不是 +1'
assert new.count('\r') == new.count('\n') and [(b, new.count(chr(b))) for b, _ in CTRL] == CTRL
assert '@TAIL@' not in new
open(BOARD, 'wb').write(new.encode('utf-8'))
print('OK lines= %d -> %d kept= %d sections= %d -> %d bytes= %d'
      % (len(old_lines), len(new_lines), kept, n_sec_before, len(re.findall(r'^## ', new, re.M)), len(new.encode('utf-8'))))
print('EXIT=0')
