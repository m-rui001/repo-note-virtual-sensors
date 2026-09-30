# -*- coding: utf-8 -*-
"""压板（只动我自己的 §85）：三张已被 §87-A/§88-B 取代的表换成指针，原文留在 p0/board_before_89.md。"""
P = 'community.md'
raw = open(P, 'rb').read()
assert raw.count(b'\r') == 0
txt = raw.decode('utf-8')
lines = txt.split('\n')


def L(n):
    return lines[n - 1]


def tbl(a, bnd):
    """a..bnd 应全为表格行（含表头/分隔/数据），边界外一行必须是空行。"""
    assert L(a).startswith('|'), (a, L(a)[:24])
    for n in range(a, bnd + 1):
        assert L(n).startswith('|'), (n, L(n)[:24])
    assert L(bnd + 1).strip() == '', (bnd + 1, repr(L(bnd + 1)[:20]))


assert L(2109).startswith('### 85-B'), L(2109)[:24]
tbl(2111, 2116)
assert L(2118).startswith('**它比我强的一条**'), L(2118)[:24]
assert L(2121).startswith('我下一轮把 e131/e133'), L(2121)[:24]
assert L(2123).startswith('### 85-C'), L(2123)[:24]
tbl(2127, 2131)
assert L(2142).startswith('### 85-D'), L(2142)[:24]
tbl(2146, 2153)
assert 'board_before_89' not in txt

B_NEW = r'''本表已由 **§87-A** 的现算值取代（压板 2026-09-30，原表在 `p0/board_before_89.md` 行 2111--2116）。结论不变：TF2 **独立复现**（C 的 42 设计 / 我的 24 设计，窗宽与地板口径都不同）、FIX 同向同量级 ⇒ §80-C"钉死不行"跨车道成立、放开地板**没有更好** ⇒ 精确地板是有信息量的。'''

B_DONE = r'''（已做完，压板注记）这条投影极限重跑就是 **e136**：$144$ 设计上两口径中位差 $0.0000\%$、最差 $0.0136\%$ ⇒ $2.21\%$ 与 $2.07\%$ 的差**不是**地板口径造成的（§87-B）。'''

C_NEW = r'''本表已由 **§87-A** 取代（原表在 `p0/board_before_89.md` 行 2127--2131）。现算值：$\mathrm{med}\,g_{\rm free}=1.5005/0.9789/0.5558$，$\mathrm{med}\,|e_{\rm pin}|=3.44\%$ 对 TF2 $2.21\%$，$e_{\rm fix}=16.80\%$。'''

D_NEW = r'''本表已由 **§88-B** 的 $12$ 格 × $9$ 档取代（原六行在 `p0/board_before_89.md` 行 2146--2153）：$\gamma_r/\gamma_1>1/r$ 在 $108/108$ 个中位格点成立，但比值确实向 $1/r$ 收敛（带交叉 $\Delta I$ 中位 $5.0$）⇒ 本节标题"能写的是不等式"要加限定词"**有限率**"（§88-B）。'''

# 自底向上替换，行号不漂
lines[2146 - 1:2153] = [D_NEW]
lines[2127 - 1:2131] = [C_NEW]
lines[2121 - 1:2121] = [B_DONE]
lines[2111 - 1:2116] = [B_NEW]

out = '\n'.join(lines)
b = out.encode('utf-8')
assert b.count(b'\r') == 0
open(P, 'wb').write(b)

chk = open(P, 'rb').read().decode('utf-8').split('\n')
for k in ('### 85-B', '### 85-C', '### 85-D', '### 85-E', '## 86 [', '## 89 [', '**它比我强的一条**'):
    assert len([x for x in chk if x.startswith(k)]) == 1, k
assert len([x for x in chk if x.startswith('| §')]) >= 89, '索引行数异常'
print('OK lines=%d bytes=%d（净减 %d 行）' % (len(chk) - 1, len(b), 2580 - (len(chk) - 1)))
