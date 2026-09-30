# -*- coding: utf-8 -*-
"""就地修 §92 索引行与 92-A 的两处抄录不准：bibitem 条数（三条，不是两条）、
`frag_rc_waterfill.tex:233$` 多一个美元符。规矩照 §91 的更正记录走。"""
P = 'community.md'
raw = open(P, 'rb').read()
cr, lf, crlf = raw.count(b'\r'), raw.count(b'\n'), raw.count(b'\r\n')
assert cr == lf == crlf
txt = raw.decode('utf-8').replace('\r\n', '\n')
FIX = [('两条 bibitem $+$ 三处限定之后', '三条 bibitem $+$ 两处限定之后'),
       ('`frag_rc_waterfill.tex:233$ 记的数', '`frag_rc_waterfill.tex:233` 记的数')]
for a, b in FIX:
    assert txt.count(a) == 1, '锚点命中 %d 次：%r' % (txt.count(a), a[:36])
    txt = txt.replace(a, b)
bts = txt.replace('\n', '\r\n').encode('utf-8')
assert bts.count(b'\r') == bts.count(b'\n') == bts.count(b'\r\n')
import collections
ctrl = sorted(collections.Counter(x for x in bts if x < 32 and x not in (10, 13)).items())
assert ctrl == [(9, 3), (11, 1)], ctrl
assert open(P, 'rb').read() == raw, '板子又被改写，本次放弃'
open(P, 'wb').write(bts)
chk = bts.decode('utf-8').split('\r\n')
assert len([x for x in chk if x.startswith('### 92-')]) == 4
assert len([x for x in chk if x.startswith('| §92 | R77-L |')]) == 1
print('OK lines=%d bytes=%d' % (len(chk) - 1, len(bts)))
