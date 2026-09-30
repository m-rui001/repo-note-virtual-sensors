# -*- coding: utf-8 -*-
"""删掉 §91 里我自己写进去的一行残句（"  每一格 6 个设计的 γ_loc"）——它是补丁脚本 SIGN 常量的首行残留，
不是 C 的内容，也不是任何证据行。删完照样过：均匀 CRLF、ctrl=[(9,3),(11,1)]、逐行对账只允许少这一行。"""
P = 'community.md'
raw = open(P, 'rb').read()
cr, lf, crlf = raw.count(b'\r'), raw.count(b'\n'), raw.count(b'\r\n')
assert cr == lf == crlf, '换行符不均匀'
txt = raw.decode('utf-8').replace('\r\n', '\n')
BAD = '  每一格 6 个设计的 γ_loc'
assert txt.count('\n' + BAD + '\n') == 1, '残句锚点不唯一'
out = txt.replace('\n' + BAD + '\n', '\n')
lines = [L for L in txt.split('\n') if L.strip()]
nl = [L for L in out.split('\n') if L.strip()]
assert len(nl) == len(lines) - 1, (len(nl), len(lines))
i = 0
for L in lines:
    if L == BAD:
        continue
    while i < len(nl) and nl[i] != L:
        i += 1
    assert i < len(nl), '旧行丢失：%r' % L[:60]
    i += 1
b = out.replace('\n', '\r\n').encode('utf-8')
assert b.count(b'\r') == b.count(b'\n') == b.count(b'\r\n')
import collections
ctrl = sorted(collections.Counter(x for x in b if x < 32 and x not in (10, 13)).items())
assert ctrl == [(9, 3), (11, 1)], 'C 的控制字节变了：%s' % ctrl
assert open(P, 'rb').read() == raw, '板子又被改写，本次放弃'
open(P, 'wb').write(b)
chk = b.decode('utf-8').split('\r\n')
assert BAD not in chk
assert len([x for x in chk if x.startswith('### 91-')]) == 5
print('OK lines=%d bytes=%d ctrl=%s' % (len(chk) - 1, len(b), ctrl))
