import hashlib
import os

for path in ('community.md', 'p0/board_full_20260930_0115.md'):
    b = open(path, 'rb').read()
    n_crlf = b.count(bytes([13, 10]))
    n_cr = b.count(bytes([13]))
    fixed = b.replace(bytes([13, 10]), bytes([10]))
    rest_cr = fixed.count(bytes([13]))
    import collections
    ctrl = sorted((k, v) for k, v in collections.Counter(
        x for x in fixed if x < 32 and x != 10).items())
    if fixed != b:
        open(path, 'wb').write(fixed)
    print('%s: CR=%d CRLF=%d -> 剩余 CR=%d ; lines=%d ; ctrl=%s ; md5=%s'
          % (path, n_cr, n_crlf, rest_cr, fixed.count(bytes([10])), ctrl,
             hashlib.md5(fixed).hexdigest()))
    print('   size=%d' % os.path.getsize(path))
