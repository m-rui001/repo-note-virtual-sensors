# -*- coding: utf-8 -*-
"""第二次压板（用户 2026-09-30 要求）：只替换**我自己那三处已被后续实测取代的表格**，
结论句、小节标题、C 的段落一律不动。每处替换都留下"冻结备份文件 + 原行号"，随时可查回。
三处：
  A. §88-D 的 12 行 Pearson 表 $\to$ `e137c` 扩展窗的重算结果（数字可在 p0/e137c_out.txt grep）；
  B. §90-A 的"三种性质"表 $\to$ §92-D 的合并文献账；
  C. §90-D 的 $\Delta I=8$ 对表 $\to$ §91-B 的 $\Delta I=12/16$ 对表。
"""
import collections

P = 'community.md'
BK = 'p0/board_before_squeeze2.md'
raw = open(P, 'rb').read()
cr, lf, crlf = raw.count(b'\r'), raw.count(b'\n'), raw.count(b'\r\n')
assert cr == lf == crlf, '换行符不均匀'
txt = raw.decode('utf-8').replace('\r\n', '\n')
L = txt.split('\n')


def L_(n):
    return L[n - 1]


def run(start_prefix, min_rows=5):
    hits = [n for n in range(1, len(L) + 1) if L_(n).startswith(start_prefix)]
    assert len(hits) == 1, '表头 %r 命中 %d 次' % (start_prefix, len(hits))
    a = hits[0]
    b = a
    while b + 1 <= len(L) and L_(b + 1).startswith('|'):
        b += 1
    assert b - a + 1 >= min_rows, '表只有 %d 行，不像被压缩对象' % (b - a + 1)
    assert L_(b + 1).strip() == '' or not L_(b + 1).startswith('|'), b + 1
    return a, b


A = run('| 格 | Pearson | 判决 |', 12)
B = run('| 记录 | 状态 | 判 |', 5)
C = run('| 秩 | 我 $\\gamma_{\\rm loc}(8.0)$', 3)
print('表位：A=%s B=%s C=%s（行号）' % (A, B, C))
for (a, b) in (A, B, C):
    print('  表前一行：%s' % L_(a - 1)[:60])
    print('  表后一行：%r' % L_(b + 1)[:60])

PA = ('（此处原有 12 行 Pearson 表已被 `p0/e137c_gamma_dI16_fixed.py` 的扩展窗重算**取代**：'
      '$\\Delta I$ 网格加到 $12,16$ 之后，$11/12$ 格的 Pearson 落在 $-0.888\\!\\sim\\!-0.931$，'
      '全部不过我预注册的 $|r|\\ge0.95$ 门 $\\Rightarrow$ 不报斜率；唯一过门的是 rand-4 $r2$，'
      '斜率 $-2.617$、Pearson $-0.9502$，脚本自己标注"拟合窗 10 点，另有 1 个有限点因 $\\le 1/r$ 被排除"'
      '（即该斜率只在 $\\Delta I\\le12$ 内有效）。原 12 行见 `%s` 第 %d--%d 行。）' % (BK, A[0], A[1]))
PB = ('（此处原有"三种性质"5 行表已被 §92-D 的合并文献账**取代**：真欠 $=2$ 条 IEEE 记录 $+$ 3 条只读到摘要的 '
      'Stavrou 线预印本，$a\\wedge b$ 撞车仍为 $0$。原表见 `%s` 第 %d--%d 行。）' % (BK, B[0], B[1]))
PC = ('（此处原有的 $\\Delta I=8$ 三行对表已被 §91-B 的 $\\Delta I=12$（18 格全在 $+0.9\\%%$ 内）与 '
      '$\\Delta I=16$（rank 2/3 各 6/6 株 $\\le0.2\\%%$）**取代**；"差 $\\le1.7\\%%$"那句话的强度由 §91-B 顶上。'
      '原表见 `%s` 第 %d--%d 行。）' % (BK, C[0], C[1]))

# 自下而上替换，行号才不会串
for (a, b), new in sorted([(C, PC), (B, PB), (A, PA)], key=lambda x: -x[0][0]):
    L[a - 1:b] = [new]
out = '\n'.join(L)

nl = out.split('\n')
old_nonblank = [x for x in txt.split('\n') if x.strip()]
new_nb = [x for x in nl if x.strip()]
assert len(new_nb) == len(old_nonblank) - (sum(b - a for a, b in (A, B, C))) + 0
# 逐段对账：三个表之外的每一行都必须原样按序存在
skip = set()
for a, b in (A, B, C):
    skip.update(range(a, b + 1))
i = 0
orig = txt.split('\n')
for n, x in enumerate(orig, 1):
    if n in skip:
        continue
    while i < len(nl) and nl[i] != x:
        i += 1
    assert i < len(nl), '表外某行被改动：%r' % x[:60]
    i += 1
print('表外逐行对账通过（%d 行原样保留），板子 %d → %d 行'
      % (len(orig) - len(skip), len(orig), len(nl)))

bts = out.replace('\n', '\r\n').encode('utf-8')
assert bts.count(b'\r') == bts.count(b'\n') == bts.count(b'\r\n')
ctrl = sorted(collections.Counter(x for x in bts if x < 32 and x not in (10, 13)).items())
assert ctrl == [(9, 3), (11, 1)], 'C 的控制字节变了：%s' % ctrl
assert open(P, 'rb').read() == raw, '板子又被改写，本次放弃'
open(BK, 'wb').write(raw)
assert open(BK, 'rb').read() == raw, '备份未通过 cmp 等价校验'
open(P, 'wb').write(bts)
chk = bts.decode('utf-8').split('\r\n')
assert len([x for x in chk if x.startswith('### 88-')]) == 5
assert len([x for x in chk if x.startswith('### 90-')]) == 5
assert len([x for x in chk if x.startswith('### 91-')]) == 5
assert len([x for x in chk if x.startswith('### 92-')]) == 4
print('OK lines=%d bytes=%d ctrl=%s' % (len(chk) - 1, len(bts), ctrl))
