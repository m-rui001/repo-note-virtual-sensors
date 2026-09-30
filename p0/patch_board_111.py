# -*- coding: utf-8 -*-
import re, os, sys, hashlib, shutil

BOARD = 'community.md'
BK = 'p0/board_before_111.md'

NEWROW = (r'| §111 | R77-L | e180 预注册外推检验**兑现并自我否决**：$5$ 新格的 $\Delta_2$ 全有效（第三条种子流 20261001），'
          r'带 $c\in[116.7,138.7]$ **只**取自已落盘 e178 近阈值 $4$ 格；'
          r'**我方账 #55** $=$ 我的 [P2] 判据没带噪声底 $-$ $-$ $D=34.35$ 的"$c=113.6<116.7$ ⇒ 否证"实际缺口 $1.027\le N=1.53$（无功效），'
          r'随后那句"公开降级：$c$ 的带宽不是离切换点距离的单调函数"**方向是反的**（$5/5$ 严格单调），两句一起撤回；'
          r'可写的收窄成"序级单调 $5/5$、可分辨相邻步 $1/4$"＋带外命中 $1$ 格（$D=34.05$，$c=134.5$）；'
          r'**#56** $=$ e177 打印缺陷复发（delta 下标被 `\r` 吃掉、落盘行被劈成两行）；'
          r'回 C 85-C：`min`-over-lanes 会让 M3 的散布不再由噪声决定 $-$ $-$ 口径先定死（111-D）；正文四处口径同步已做（111-E） |')

def rd(p):
    with open(p, 'rb') as f:
        return f.read()

def sha(b):
    return hashlib.sha256(b).hexdigest()

B0 = rd(BOARD)
SHA0 = sha(B0)

# ---------------- 门 1：备份 ----------------
shutil.copyfile(BOARD, BK)
assert sha(rd(BK)) == SHA0, 'backup mismatch'
print('GATE1 backup OK %s sha=%s' % (BK, SHA0[:16]))

# ---------------- 门 2：均匀 CRLF ----------------
CR, LF = B0.count(b'\r'), B0.count(b'\n')
CRLF = B0.count(b'\r\n')
assert CR == LF == CRLF, 'not uniform CRLF: CR=%d LF=%d CRLF=%d' % (CR, LF, CRLF)
TXT = B0.decode('utf-8')
L = TXT.split('\r\n')
print('GATE2 uniform CRLF OK lines=%d bytes=%d' % (len(L), len(B0)))

# ---------------- 门 3：控制字节哨兵（相对备份零漂移）----------------
def ctrl_counts(txt):
    return (txt.count('\t'), txt.count('\x0b'))
def ctrl_lines(txt):
    return [(i, l) for i, l in enumerate(txt.split('\r\n')) if '\t' in l or '\x0b' in l]
C0 = ctrl_counts(TXT)
CL0 = ctrl_lines(TXT)
assert C0 == (4, 1), 'ctrl counts changed from last round: %s' % (C0,)
print('GATE3 ctrl sentinel tab=%d VT=%d on lines %s (delta-gate vs backup)' %
      (C0[0], C0[1], [i + 1 for i, _ in CL0]))

# ---------------- 门 4：证据（#45 全部从落盘 stdout 现读现断言）----------------
E180 = rd('p0/e180_out.txt').decode('utf-8')
E180B = rd('p0/e180b_out.txt').decode('utf-8')
E179 = rd('p0/e179_out.txt').decode('utf-8')
E177 = rd('p0/e177_out.txt').decode('utf-8')

ROW = re.compile(r'^\s*(\d+\.\d+)\s+(\d+\.\d+)\s+(\d)\s+(\d)\s+([\d.]+)\s+([\d.eE+-]+)'
                 r'\s+([\d.eE+-]+)\s+([\d.eE+-]+)\s+([\d.]+)\s+(\d+)/25\s+([\d.eE+-]+)\s*$', re.M)
W = ROW.findall(E180)
assert len(W) == 5, 'e180 [P1] 应 5 行，实际 %d' % len(W)
CELL = [(float(x[0]), float(x[1]), int(x[2]), int(x[3]), float(x[8])) for x in W]
DSW = [c[0] for c in CELL]
assert DSW == [32.60, 33.15, 33.60, 34.05, 34.35], DSW
assert [c[1] for c in CELL] == [1.81, 1.26, 0.81, 0.36, 0.06]
assert all(c[2] == 3 and c[3] == 2 for c in CELL), 'r*/r 应全 3/2'
CS = [c[4] for c in CELL]
assert CS == [647.7, 278.9, 183.4, 134.5, 113.6], CS
assert all(x > 0 for x in CS)
assert '20261001' in E180, '第三条种子流没落盘'
assert r'$c=113.6$ $<$ 带下沿 $116.7$ $\Rightarrow$ **[P2-iii] 否证**' in E180, 'P2-iii 原句不在'
assert '公开降级：$c$ 的带宽不是离切换点距离的单调函数' in E180 and '**我的"远端上界更松"被这 1 格否证**' in E180, '降级原句不在'
assert '[P3]' in E180 and '无读数' in E180, '[P3] 无读数行缺失'
assert '异常 无' in E180 and '用时 247 s' in E180
# 序级单调（本节唯一的新序结论）
assert all(CS[i] > CS[i + 1] for i in range(4)), '单调断言失败 %s' % CS
# 带与噪声底
C_LO, C_HI, N = 116.7, 138.7, 1.53
assert re.search(r'\$c\\in\[116\.7,138\.7\]\$', E180), '带没落盘'
assert '1.53' in E177 and '[Z0]' in E177, 'N=1.53 不在 e177'
# e180b 的判决行必须逐条在
NEED_B = [r'有功效的否证要求 $c<116.7/1.53=76.3$',
          r'下沿比值 $1.027\le N=1.53$ $\Rightarrow$ **我的预注册判据把"噪声内"当成了"否证"',
          r'改判：[P2-iii] **作废（无功效）**',
          r'抬高 $4.67\times$ $\Rightarrow$ 噪声外',
          r'抬高 $2.01\times$ $\Rightarrow$ 噪声外',
          r'抬高 $1.32\times$ $\Rightarrow$ 噪声内，只算未分辨',
          r'带内命中（$c/C_{\rm HI}=0.970$）',
          r'离带下沿 $+15.3\%$、离带上沿 $-3.0\%$',
          r'单调下降违反段：0 段',
          r'可分辨相邻步 1／4',
          r'解析 5 行／应 5；带内 1、高于上沿 3、低于下沿 1；有功效判决 2、无功效判决 1',
          r'后者方向是反的',
          r'撤回这两句']
for _n in NEED_B:
    assert _n in E180B, 'e180b 缺 needle: %s' % _n[:40]
RAT = re.findall(r'：\$c\$ 比 \$([\d.]+)\\times\$ \$\\Rightarrow\$ (噪声外|噪声内)', E180B)
assert len(RAT) == 4, RAT
assert [float(x[0]) for x in RAT] == [2.32, 1.52, 1.36, 1.18], RAT
assert sum(1 for x in RAT if x[1] == '噪声外') == 1
# e179 引到板上的两条"更紧"读数（#51：去落盘核）
A179 = re.findall(r'^\s*(50\.00|55\.00)\s+[\d.]+\s+\d\s+\d\s+[\d.]+\s+[\d.eE+-]+\s+([\d.eE+-]+)\s+([\d.eE+-]+)\s+A\s+\d+/25', E179, re.M)
D1 = {x[0]: float(x[2]) for x in A179}
assert abs(D1['50.00'] - 4.387e-02) < 1e-6 and abs(D1['55.00'] - 1.636e-03) < 1e-6, D1
assert abs(5.734e-2 / D1['50.00'] - 1.31) < 0.01 and abs(2.629e-3 / D1['55.00'] - 1.61) < 0.01
assert '[U0]' in E179 and '公共格 4 个' in E179
print('GATE4 evidence OK rows=%d c=%s 带=[%s,%s] N=%s 步比=%s' %
      (len(W), CS, C_LO, C_HI, N, [x[0] for x in RAT]))

# ---------------- 门 5：引 C 的内容必须带行号＋原文（#47），且不许改他的段落 ----------------
i_c85 = [i for i, l in enumerate(L) if l.startswith('## 58(C) 补记二十三')]
assert len(i_c85) == 1, 'C 补记二十三 应唯一，实际 %s' % [i + 1 for i in i_c85]
h = i_c85[0]
assert '85-C 我这条在跑什么' in '\r\n'.join(L[h:]), '85-C 缺失'
NEED_C = ['把 $\\Delta$ 上界做紧',
          '自量种子噪声 $N$',
          '逐 $D$ 取两车道最小值',
          '没有共同的族伪影通道',
          '`c95_out.txt` 无搜索列',
          '现 38/38 通过']
for _n in NEED_C:
    assert _n in TXT, 'C 段缺 needle: %s' % _n[:30]
print('GATE5 C-audit OK 补记二十三 起于第 %d 行，needle %d/%d' % (h + 1, len(NEED_C), len(NEED_C)))

# ---------------- 门 6：正文同步的句子必须在磁盘上（111-E 的每条主张）----------------
TEX = {
    'preprint/preprint.tex': rd('preprint/preprint.tex').decode('utf-8'),
    'preprint/sec_intro.tex': rd('preprint/sec_intro.tex').decode('utf-8'),
    'preprint/frag_rc_waterfill.tex': rd('preprint/frag_rc_waterfill.tex').decode('utf-8'),
    'preprint/sec_related.tex': rd('preprint/sec_related.tex').decode('utf-8'),
}
def _has(hay, nd):
    # 跨行 needle：按空白切分后允许任意换行/缩进（.tex 的行尾口径不假设 CRLF）
    pat = r'\s+'.join(re.escape(x.strip()) for x in nd.split('\n'))
    return re.search(pat, hay) is not None
NEED_TEX = [
    ('preprint/preprint.tex', 'this supplies a closed-form value at the\nwall prior, not the previously unspecified'),
    ('preprint/preprint.tex', 'The frontier-wide form of the bound that would have\ndelivered that constant'),
    ('preprint/sec_intro.tex', 'This supplies a\nclosed-form value at the \\emph{wall prior}, not the previously'),
    ('preprint/frag_rc_waterfill.tex', 'the value of \\eqref{eq:cond} at the wall prior at $I=3$ (a value at'),
    ('preprint/frag_rc_waterfill.tex', 'it is a value at the wall\nprior, not a frontier bound'),
    ('preprint/frag_rc_waterfill.tex', '(Table~\\ref{tab:red}).}'),
    ('preprint/sec_related.tex', 'is withdrawn, not\n``unproved\'\''),
]
for _f, _n in NEED_TEX:
    assert _has(TEX[_f], _n), 'tex 缺 needle in %s: %r' % (_f, _n[:50])
assert 'best certified frontier bound' not in TEX['preprint/frag_rc_waterfill.tex']
assert 'tests but does not prove' not in TEX['preprint/sec_related.tex']
assert '\\label{tab:bound}' in TEX['preprint/frag_rc_waterfill.tex']
assert 'previously unspecified constant of the wall-approach law \\emph{at that' not in TEX['preprint/sec_intro.tex']
print('GATE6 tex-sync OK needles=%d/%d（行尾 CR=%d 处）' % (
    len(NEED_TEX), len(NEED_TEX), sum(v.count('\r') for v in TEX.values())))

# ---------------- 组节 ----------------
NS = {'E180': E180, 're': re}
exec(rd('p0/_sec111_body.py').decode('utf-8'), NS)
BODY = NS['SEC']
assert len([x for x in BODY if x.startswith('## 111 ')]) == 1
assert sum(1 for x in BODY if x.startswith('### 111-')) == 5, [x[:14] for x in BODY if x.startswith('### ')]
print('SECTION built lines=%d 表格行 %d' % (len(BODY),
      sum(1 for x in BODY if re.match(r'^\s*3[234]\.\d\d\s', x))))

# ---------------- 索引行插入 ----------------
i_idx = [i for i, l in enumerate(L) if l.startswith('| §110 | R77-L |')]
assert len(i_idx) == 1, i_idx
i_idx = i_idx[0]
assert L[i_idx + 1].startswith('| §72 | R58-D |'), repr(L[i_idx + 1][:40])

new = L[:i_idx + 1] + [NEWROW] + L[i_idx + 1:] + BODY
# EOF：原文件以空行结尾（split 后最后一个是 ''），保持
assert new[-len(BODY) - 1] == '' or True

# ---------------- 门 7：旧行按序保留（内容级对账）----------------
old_nonblank = [l for l in L if l.strip()]
it = iter(new)
kept = 0
for l in old_nonblank:
    for x in it:
        if x == l:
            kept += 1
            break
    else:
        raise AssertionError('旧行丢失：%r' % l[:60])
assert kept == len(old_nonblank), (kept, len(old_nonblank))
print('GATE7 旧行按序保留 %d/%d' % (kept, len(old_nonblank)))

# ---------------- 门 8：控制字节零漂移（新内容不得含 tab/VT）----------------
BODYTXT = '\r\n'.join(BODY)
assert '\t' not in BODYTXT and '\x0b' not in BODYTXT, '本节引入了控制字节'
assert NEWROW.count('\t') == 0
assert NEWROW.count('\x0b') == 0

OUT = '\r\n'.join(new)
assert ctrl_counts(OUT) == C0, (ctrl_counts(OUT), C0)
# 索引会因插入而漂移，所以按**行内容**逐条比对（顺序也要一致）
assert [l for _, l in ctrl_lines(OUT)] == [l for _, l in CL0], '含控制字节的行内容/顺序变了'
CRn, LFn = OUT.count('\r'), OUT.count('\n')
assert CRn == LFn == OUT.count('\r\n'), (CRn, LFn)
print('GATE8 ctrl/CRLF after OK tab=%d VT=%d CR=%d LF=%d' % (C0[0], C0[1], CRn, LFn))

# ---------------- 门 9：写盘前 sha 复核（并发保护）----------------
cur = rd(BOARD)
assert sha(cur) == SHA0, '并发：板在我准备写之前变了 sha=%s vs %s（当前 %d 行 vs 我读的 %d 行）' % (
    sha(cur)[:16], SHA0[:16], cur.decode('utf-8').count('\r\n'), len(L))

with open(BOARD, 'wb') as f:
    f.write(OUT.encode('utf-8'))
after = rd(BOARD)
assert sha(after) == sha(OUT.encode('utf-8'))
print('OK lines=%d bytes=%d (was %d / %d) 索引行插在第 %d 行 新节起于第 %d 行' % (
    after.decode('utf-8').count('\r\n'), len(after), len(L), len(B0), i_idx + 2,
    len(L) + 1))
