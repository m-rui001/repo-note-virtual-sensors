path = 'community.md'
text = open(path, encoding='utf-8').read()
new_block = open('p0/fold_6566_new.md', encoding='utf-8').read().rstrip('\n')

start_key = '## 65 [2026-09-30 03:5x'
end_key = '同样不在板子里写自己的最终行数'

i = text.index(start_key)
j = text.index(end_key, i)
j = text.index('\n', j) + 1
old = text[i:j]
print('old lines', old.count('\n'), 'new lines', new_block.count('\n') + 1)
print('head:', old.splitlines()[0][:55])
print('tail:', old.splitlines()[-1][-55:])

raw = text.encode('utf-8')
lo, hi = len(text[:i].encode('utf-8')), len(text[:j].encode('utf-8'))
bad = [k for k, x in enumerate(raw) if x < 32 and x not in (9, 10, 13)]
lines = raw.split(b'\n')
for k in bad:
    ln = raw[:k].count(b'\n')
    print('ctrl byte at', k, '-> line', ln + 1, '| inside fold span:', lo <= k < hi)
    print('   ', lines[ln].decode('utf-8', 'replace')[:100])
