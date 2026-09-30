# -*- coding: utf-8 -*-
"""§100-E 追加一条 ⑤：把"负搜索"留痕（正文里从来没有过 Θ-特征轴对齐的句子）。纯就地追加，不开新小节。"""
import collections

P = 'community.md'
raw = open(P, 'rb').read()
cr, lf, crlf = raw.count(b'\r'), raw.count(b'\n'), raw.count(b'\r\n')
assert cr == lf == crlf, '换行符不均匀：CR=%d LF=%d CRLF=%d' % (cr, lf, crlf)
txt = raw.decode('utf-8').replace('\r\n', '\n')
hits = [L for L in txt.split('\n') if L.startswith('$④$ 仍欠：e152')]
assert len(hits) == 1, '锚行不唯一：%d' % len(hits)
anchor = hits[0]
ADD = ('$⑤$ **负搜索留痕（本轮查，供后人不重复怀疑）**：'
       '`grep -n "Theta" p0/note/note.tex p0/note/frag_*.tex` 里与 eigen/axis/align/leading/optimal/direction 同行的命中数为 **0** $\\Rightarrow$ '
       'e152 否证的那类句子（"最优设计的支撑落在 $\\Theta$ 特征轴上"）**本站正文从未写过**，所以本轮没有需要撤的已发表句子；'
       '这条记录的作用是：以后若要写"选轴"层面的对齐陈述，必须带 e152 的 $\\Delta I$ 限定，不得写成精确对齐。')
out = txt.replace(anchor + '\n', anchor + '\n' + ADD + '\n')
b = out.replace('\n', '\r\n').encode('utf-8')
assert b.count(b'\r') == b.count(b'\n') == b.count(b'\r\n')
ctrl = sorted(collections.Counter(x for x in b if x < 32 and x not in (10, 13)).items())
assert ctrl == [(9, 3), (11, 1)], 'C 的控制字节变了：%s' % ctrl
now = open(P, 'rb').read()
assert now == raw, '板子被并发改写（%d→%d），本次放弃' % (len(raw), len(now))
open('p0/board_before_100b.md', 'wb').write(raw)
assert open('p0/board_before_100b.md', 'rb').read() == raw
old = txt.split('\n')
open(P, 'wb').write(b)
chk = open(P, 'rb').read().decode('utf-8').replace('\r\n', '\n').split('\n')
nl = [x for x in chk if x.strip()]
i = 0
for x in old:
    if not x.strip():
        continue
    while i < len(nl) and nl[i] != x:
        i += 1
    assert i < len(nl), '旧行被改动：%r' % x[:60]
    i += 1
assert len([x for x in chk if x.startswith('### 100-')]) == 5, '小节数变了'
assert len([x for x in chk if x.startswith('$⑤$ **负搜索留痕')]) == 1
print('旧行按序保留到第 %d 个非空行（新板非空 %d）' % (i, len(nl)))
print('OK lines=%d bytes=%d ctrl=%s' % (len(chk) - 1, len(b), ctrl))
