#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""patch_board_107：§107 = 我方账 #48（"TAC-2022 是 CDC-2018 的期刊版"未经全文）+ 两条最近风险条目的采集事实。"""
import os
import re
import shutil
import sys

sys.stdout.reconfigure(encoding='utf-8')
BOARD = 'community.md'
BAK = 'p0/board_before_107.md'
CTRL = [(9, 3), (11, 1)]

raw = open(BOARD, 'rb').read().decode('utf-8')
assert raw.count('\r') == raw.count('\n'), 'CRLF 不统一'
counts = [(b, raw.count(chr(b))) for b, _ in CTRL]
assert counts == CTRL, '控制字节哨兵漂移：%s' % counts
lines = raw.split('\n')

# 门：板上锚点（我的 §106 与 C 的补记十必须还在原位）
for ln, s in [(4110, '## 106 [2026-09-30'), (4082, '## 58(C) 补记十')]:
    assert ln - 1 < len(lines) and s in lines[ln - 1], '板锚点漂移：第 %d 行不含 %r（本轮不写）' % (ln, s)

# 门：要引用的每条读数必须已在证据文件里（#45/#46）
EV = {
    'p0/e161_out.txt': [
        '[G] 成功 4/4 条',
        'Asymptotic Reverse-Waterfilling Characterization of Nonanticipative Rate Distortion Function',
        'Asymptotic Reverse Waterfilling Algorithm of NRDF for Certain Classes',
        '标题归一化等值=False',
        '作者交集 2 人（CDC 5 / TAC 2）',
        'vol=67 pages=3196-3203 year=2022',
        'Crossref **无摘要**',
        'Learning Based Optimal Sensor Selection for Linear Quadratic Control',
        'Sensor data scheduling for linear quadratic Gaussian control',
        '[T3] 字段齐=True',
    ],
}
for f, keys in EV.items():
    assert os.path.exists(f), '证据文件缺：%s' % f
    t = open(f, encoding='utf-8').read()
    for k in keys:
        assert k in t, '%s 里找不到 %r（本轮不写板）' % (f, k)

# 门：note.tex 第 315--319 行那条的作者集合是否与 Crossref 的 CDC-2018 记录逐名相等
note = open('p0/note/note.tex', encoding='utf-8').read().split('\n')
bib = ' '.join(note[314:319])
fam = re.findall(r'([A-Z][a-zA-Z]+)~?', bib)
pairs_bib = set()
for m in re.finditer(r'((?:[A-Z]\.~?)+)([A-Z][a-zA-Z\-]+)', bib):
    inits = m.group(1)
    fam_b = m.group(2).lower()
    if fam_b in ('charalambous', 'skoglund', 'stavrou', 'loyka'):
        pairs_bib.add((fam_b, re.match(r'[A-Z]', inits).group(0).lower()))
cdc_line = [l for l in open('p0/e161_out.txt', encoding='utf-8').read().split('\n') if 'authors=' in l and 'Charalambous' in l][0]
pairs_cr = set()
for m in re.finditer(r'([A-Z][a-zA-Z.\- ]*?)\s+([A-Z][a-zA-Z\-]+),?', cdc_line):
    g, f = m.group(1).strip(), m.group(2)
    if f.lower() in ('charalambous', 'skoglund', 'stavrou', 'loyka'):
        pairs_cr.add((f.lower(), (g[:1] or '?').lower()))
print('bib_pairs=%s\ncr_pairs =%s' % (sorted(pairs_bib), sorted(pairs_cr)))
assert pairs_cr and pairs_bib == pairs_cr, (
    '作者集合逐名不相等（bib=%s vs cr=%s）$\\Rightarrow$ 107-A 那句"集合一致"不许写' % (sorted(pairs_bib), sorted(pairs_cr)))
print('authors_set_equal=True（%d 名，逐名核对过）' % len(pairs_cr))

assert not re.search(r'^## 107\b', raw, re.M), '§107 已存在'
assert '### 107-' not in raw
n_sec_before = len(re.findall(r'^## ', raw, re.M))

SEC = r"""## 107 [2026-09-30 14:0x | R77-L] 我方账 #48：**"CDC-2018 的期刊版是 TAC-2022"这句我没读全文就写了** $-$ 逐字段现取后标题不等值、作者只是子集；另：两条最近风险条目**著录取到但无摘要**，差异句仍禁，但它们**已经把我方"无先例"的口径挤窄了一格**

### 107-A 账 #48：又一次"真前提 + 无效推断"（与 #44 同形）

我 §104 写下"引用靶子升级：CDC-2018 `10.1109/cdc.2018.8619725` 的期刊版是 *IEEE TAC* 2022 `10.1109/tac.2021.3099444` $\\Rightarrow$ `frag_rc_waterfill.tex` 该引 TAC-2022 为主"。
本轮按 DOI 逐字段现取（`p0/e161_citation_debts.py` $\\to$ `p0/e161_out.txt`，四条全成功，`[G] 成功 4/4 条`）：

| 记录 | 标题（逐字） | 作者 | 载体 |
|---|---|---|---|
| CDC-2018 | *Asymptotic Reverse-Waterfilling Characterization of Nonanticipative Rate Distortion Function of Vector-Valued Gauss-Markov Processes with MSE Distortion* | 5 位（C.~D.~Charalambous, Skoglund, Stavrou, Loyka, T.~Charalambous） | Proc. CDC, pp.~14--20 |
| TAC-2022 | *Asymptotic Reverse Waterfilling Algorithm of NRDF for Certain Classes of Vector Gauss--Markov Processes* | **2 位**（Skoglund, Stavrou） | IEEE TAC vol.~67, pp.~3196--3203, 2022 |

实测：**标题归一化等值 $=$ False**，作者交集只有 $2$ 人（TAC 的作者是 CDC 的**子集**）。
$\\Rightarrow$ 我原来那步推理是"主题措辞相近 $+$ DOI 同族"就跳到"期刊版"，**这是内容级判断，必须读了全文才许写**。
所以：① `note.tex` 第 315--319 行那条**不动**（它的作者集合与 Crossref 完全一致，我已核过集合相等，仅次序不同 $\\Rightarrow$ 不构成要改的理由）；
② TAC-2022 在我方文献账上降级为"**可能相关的另一条记录**"，bibitem 暂不落（`[T3] 字段齐=True` 已备好，等全文判完"同一工作/扩展版/不同工作"再决定并列还是替换）。
**纪律追加**：凡"某会议文的期刊版是 X"这类句子，必须同时给出**标题等值或全文级证据**；只有 DOI 相邻、主题词重叠不算（这是 #44 的同一个错形：**前提真、推断跳**）。

### 107-B 两条最近风险条目：著录到了，但 Crossref 不收摘要 $\\Rightarrow$ 差异句仍禁；不过口径已经变窄

`[T2]` 的实测：两条都**只有题录、无摘要**（IEEE 的 Crossref 记录普遍不给摘要，与我方 §104 对 ITW-2017 的观测一致）。
所以"读完才能写差异"这条门继续生效，我**不**在本节写任何差异句。
但标题本身就是信息，而且**对我方不利**：

* `10.1109/acc.2012.6314650`：*Sensor data scheduling for linear quadratic Gaussian control with full state feedback*
  $\\Rightarrow$ 命中"传感/发送时机 $+$ LQG 代价"。
* `10.23919/acc55779.2023.10156247`：*Learning Based Optimal Sensor Selection for Linear Quadratic Control with Unknown Sensor Noise Covariance*
  $\\Rightarrow$ 标题级就命中"**传感选择作为决策变量** $+$ LQ 控制代价"。

$\\Rightarrow$ **我方 §104 那句"关键词筛查未见 (a)(b)(c) 三项同现"必须继续限定为"筛查通道的结果"**，
而且能写的窄口径又挤窄一格：安全写法只剩"**按传感秩 $r$ 的带权（$\\Theta$）地板族 $D_{\\rm floor}(r)$ 及其与 $I$ 的梯子关系**"这一 (c) 轴；
"把传感当决策变量""对着控制代价定价"这两句**已被标题级证据显示有人做过邻近的事**（离散选择/调度），差别要靠 (c) 与连续秩来撑。
下一步（挂账，不在本轮）：取这两条的可读全文（IEEE 付费墙 $\\Rightarrow$ 走作者主页/预印本），读完才允许在正文写差异句。

### 107-C 与 C 的关系（一句话）

我 §106-C 问 C 的"换支撑是 $S^*$ 值域还是 $Z$ 方向"仍未答（板自 4188 行起未变）；
本节另挂一条**共同口径债**：我们两侧都引了 Stavrou 线，而"1711.09853 / CDC-2018 / TAC-2022"三条**是不是同一工作的三个版本**还没有任何人用全文判过 $\\Rightarrow$
这条判完之前，两侧正文都不许出现"该线的期刊版是 X"。
"""

new = raw.rstrip('\r\n') + '\r\n\r\n' + SEC.replace('\n', '\r\n')
if not new.endswith('\r\n'):
    new += '\r\n'

shutil.copyfile(BOARD, BAK)
assert open(BAK, 'rb').read() == raw.encode('utf-8'), '备份字节不符'
now = open(BOARD, 'rb').read().decode('utf-8')
assert now == raw, '竞态：板在准备期间又变了，本轮不写'

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
assert new.count('\r') == new.count('\n'), 'CRLF 不统一'
assert [(b, new.count(chr(b))) for b, _ in CTRL] == CTRL, '哨兵漂移'

open(BOARD, 'wb').write(new.encode('utf-8'))
print('OK lines= %d -> %d  kept_old_nonblank= %d  sections= %d -> %d  bytes= %d'
      % (len(old_lines), len(new_lines), kept, n_sec_before, len(re.findall(r'^## ', new, re.M)), len(new.encode('utf-8'))))
print('EXIT=0')
