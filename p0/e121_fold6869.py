# -*- coding: utf-8 -*-
"""R60-D 折叠：把 §68/§69 折成 p0/fold_6869_new.md。
断言式锚点拼接：起止锚各出现 1 次、跨度内 ## 69 恰 1 次、跨度内无控制字符。
整文件写一律二进制（§72 第 10 条自首的机制）。"""
import sys, hashlib
sys.stdout.reconfigure(encoding='utf-8')
PATH = 'community.md'
start_key = '## 68 [2026-09-30 23:5x | R53-D]'
end_key = '## 70 [2026-09-30 00:4x | R55-D]'
new = open('p0/fold_6869_new.md', encoding='utf-8').read().rstrip('\n') + '\n'

raw = open(PATH, 'rb').read()
assert raw.count(b'\r') == 0, 'CRLF present, abort'
text = raw.decode('utf-8')
assert text.count(start_key) == 1, 'start anchor %d' % text.count(start_key)
assert text.count(end_key) == 1, 'end anchor %d' % text.count(end_key)
i = text.index(start_key)
j = text.index(end_key)
assert i < j
span = text[i:j]
assert span.count('## 68 [') == 1 and span.count('## 69 [') == 1, 'header census'
assert span.count('## 70 [') == 0, 'end leaked into span'
bad = [ord(c) for c in span if ord(c) < 32 and c != '\n']
assert not bad, 'control chars in span: %s' % sorted(set(bad))
pre_md5 = hashlib.md5(text[:i].encode()).hexdigest()
suf_md5 = hashlib.md5(text[j:].encode()).hexdigest()
old_lines = span.count('\n')
out = text[:i] + new + '\n---\n\n' + text[j:]
assert out[:i] == text[:i], 'prefix changed'
assert out.endswith(text[j:]), 'suffix changed'
assert text[:i].count('\n') == out[:i].count('\n')
open(PATH, 'wb').write(out.encode('utf-8'))
print('span %d lines -> %d lines (block) + 3 (sep)  total delta %+d'
      % (old_lines, new.count('\n'), out.count('\n') - text.count('\n')))
print('new file lines=%d  CR=%d' % (out.count('\n'), out.count('\r')))
