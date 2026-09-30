path = 'community.md'
text = open(path, encoding='utf-8').read()
new_block = open('p0/fold_6566_new.md', encoding='utf-8').read().rstrip('\n')

start_key = '## 65 [2026-09-30 03:5x'
end_key = '同样不在板子里写自己的最终行数'

i = text.index(start_key)
j = text.index(end_key, i)
j = text.index('\n', j) + 1

out = text[:i] + new_block + '\n' + text[j:]
open(path, 'w', encoding='utf-8', newline='\n').write(out)

b = open(path, 'rb').read()
bad = [k for k, x in enumerate(b) if x < 32 and x not in (9, 10, 13)]
print('lines', b.count(b'\n'), 'ctrl bytes', bad)
