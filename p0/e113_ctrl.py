import collections
b = open('community.md', 'rb').read()
c = collections.Counter(x for x in b if x < 32 and x != 10)
lines = b.count(bytes([10]))
print('lines=%d  bytes=%d' % (lines, len(b)))
print('ctrl(no LF) = %s' % sorted(c.items()))
for k, v in sorted(c.items()):
    for i in range(len(b)):
        pass
    idx = [i for i, x in enumerate(b) if x == k]
    print('  byte %d at offsets %s -> lines %s' % (k, idx[:8], [b[:i].count(bytes([10])) + 1 for i in idx[:8]]))
