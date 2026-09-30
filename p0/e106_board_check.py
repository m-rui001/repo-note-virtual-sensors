import re, sys
b = open('community.md', 'rb').read()
bad = [i for i, x in enumerate(b) if x < 32 and x not in (9, 10, 13)]
print('ctrl bytes count', len(bad), bad[:12])
print('bytes', len(b), 'lines', b.count(b'\n'))
for nm, path in [('before67', 'p0/board_before_67.md'), ('now', 'community.md')]:
    t = open(path, encoding='utf-8').read()
    heads = [(m.start(), m.group(0).splitlines()[0][:40]) for m in re.finditer(r'^## \d+.*$', t, re.M)]
    print('---', nm)
    for k, (pos, h) in enumerate(heads):
        end = heads[k + 1][0] if k + 1 < len(heads) else len(t)
        print(f'  {h:<42} lines={t[pos:end].count(chr(10))}')
