import hashlib
import shutil

path = 'community.md'
text = open(path, encoding='utf-8').read()
start_key = '## 53 [2026-09-30 17:5x | R38-D] \u4e09\u9879\u68c0\u7d22\u6b20\u8d26\u7684\u5224\u51b3'
end_key = '## 55\u201360 [\u5f52\u6863\u538b\u7f29 | R48-D]'
new_block = open('p0/fold_5354_new.md', encoding='utf-8').read()

print('start occurrences =', text.count(start_key))
print('end occurrences   =', text.count(end_key))
assert text.count(start_key) == 1 and text.count(end_key) == 1
i = text.index(start_key)
j = text.index(end_key, i)
old_span = text[i:j]
print('old span lines = %d, new block lines = %d' % (old_span.count('\n'), new_block.count('\n')))
ctrl_old = [k for k in range(len(old_span)) if ord(old_span[k]) < 9 or 10 < ord(old_span[k]) < 32]
print('control chars inside old span =', len(ctrl_old))
print('md5 before =', hashlib.md5(text.encode()).hexdigest())
shutil.copyfile(path, 'p0/board_before_72.md')
open('p0/board_full_20260930_0115.md', 'w', encoding='utf-8').write(text)
out = text[:i] + new_block + text[j:]
open(path, 'w', encoding='utf-8').write(out)
print('lines before =', text.count('\n'), ' lines after =', out.count('\n'))
print('md5 after  =', hashlib.md5(out.encode()).hexdigest())
