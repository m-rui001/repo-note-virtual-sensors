# -*- coding: utf-8 -*-
"""通用折叠拼接：python p0/fold_generic.py <start_key> <end_key> <new_block_file> <inner_header_count_expect>
断言：起止锚各 1 次、跨度内无控制字符、前后缀逐字节不变；整文件二进制写（§72 机制）。"""
import sys
sys.stdout.reconfigure(encoding='utf-8')
PATH = 'community.md'
start_key, end_key, block_file, n_inner = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4])
new = open(block_file, encoding='utf-8').read().rstrip('\n') + '\n'
raw = open(PATH, 'rb').read()
assert raw.count(b'\r') == 0, 'CRLF present, abort'
text = raw.decode('utf-8')
assert text.count(start_key) == 1, 'start anchor=%d' % text.count(start_key)
assert text.count(end_key) == 1, 'end anchor=%d' % text.count(end_key)
i, j = text.index(start_key), text.index(end_key)
assert i < j, 'anchor order'
span = text[i:j]
assert span.count(end_key) == 0, 'end anchor inside span'
bad = sorted({ord(c) for c in span if ord(c) < 32 and c != '\n'})
assert not bad, 'control chars in span: %s' % bad
out = text[:i] + new + '\n---\n\n' + text[j:]
assert out[:i] == text[:i] and out.endswith(text[j:]), 'prefix/suffix changed'
open(PATH, 'wb').write(out.encode('utf-8'))
print('span %d lines -> block %d lines ; file %d -> %d lines ; CR=%d'
      % (span.count('\n'), new.count('\n'), text.count('\n'), out.count('\n'), out.count('\r')))
