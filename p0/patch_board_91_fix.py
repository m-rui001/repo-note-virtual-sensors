# -*- coding: utf-8 -*-
"""就地更正 §91 里三处抄录缺陷（全是我自己刚写进去的文本，不动 C 的一行）：
A. 2788 那串"末档比值"把**余隙**（-3.439/-30.216）和**比值**混在一起 ⇒ 换成 e137b stdout 里真正的六个比值；
B. 2821 的 "1.386×5" 是我随手概括，不是 stdout 原样 ⇒ 换成逐设计原值；
C. 2863 的 `note.tex:217$` 多一个美元符。
另在 91-A 末尾加一行更正记录（三处都属于同一轮"抄录未与 stdout 逐字对账"，第 33 号代码账）。"""
P = 'community.md'
raw = open(P, 'rb').read()
cr, lf, crlf = raw.count(b'\r'), raw.count(b'\n'), raw.count(b'\r\n')
assert cr == lf == crlf, '换行符不均匀'
txt = raw.decode('utf-8').replace('\r\n', '\n')

A_OLD = ('后果是那一轮 stdout 里 12 个格全印 (a)，其中 6 个格的末档是 $-29.7/-4.65/-3.44/-1.72/-1.15/-30.2$\n'
         '这种**负数**——负数显然不在 $[1/r,1/r+0.05]$ 里，判决行和它自己上面的表直接矛盾。')
A_NEW = ('后果是那一轮 stdout 里 12 个格全印 (a)，其中 6 个格的**末档比值**是\n'
         '$-29.716,\\;-19.828,\\;-4.654,\\;-3.105,\\;-1.719,\\;-1.147$（`e137b_out.txt` 的 [V2] 表末列原值）\n'
         '—— 全是**负数**，显然不在 $[1/r,1/r+0.05]$ 里，判决行和它自己上面的表直接矛盾。\n'
         '（我第一版把 $-3.439/-30.216$ 这两个**余隙**也抄进了比值列，见本节末的更正记录。）')
B_OLD = 'ΔI=12.0  有限 6/6  γ=[1.386×5, 1.388 / 1.304 ...]'
B_NEW = ('ΔI=12.0  有限 6/6  anchor γ=[1.388,1.386,1.386,1.386,1.386,1.386]'
         '  rand-1 γ=[1.304,1.385,1.386,1.386,1.386,1.386]')
C_OLD = '（`note.tex:163`、`note.tex:217$ 一带）'
C_NEW = '（`note.tex:163`、`note.tex:217` 一带）'

for a, b in ((A_OLD, A_NEW), (B_OLD, B_NEW), (C_OLD, C_NEW)):
    assert txt.count(a) == 1, '锚点命中 %d 次：%r' % (txt.count(a), a[:40])
    txt = txt.replace(a, b)

SIGN = ('  每一格 6 个设计的 γ_loc\n'
        '\n'
        '**更正记录（第 33 号代码账，2026-09-30 落 §91 当场）**：本节初稿有三处**抄录未与 stdout 逐字对账**——\n'
        '① 把 [V2] 表里的**余隙** $-3.439/-30.216$ 当成**末档比值**抄进负值清单；② 把逐设计的\n'
        '$\Delta I=12$ 档写成 "1.386$\\times$5" 这种概括而不是原值；③ `note.tex:217` 后多一个美元符。\n'
        '三处都已就地改成本文本可 grep 的原值。**成因**：我是从上一轮的记忆里抄的数，没有先 `grep -a` 回来对。\n'
        '进 feedback 第 10 条的新实例：**抄录级错误也要当场记账**，因为它和判决级错误长得一模一样。\n')
anchor = '这是 feedback 第 10 条的又一实例：**整列异常要先当 bug 报告，不是否证**。'
assert txt.count(anchor) == 1
txt = txt.replace(anchor, anchor + '\n' + SIGN.rstrip('\n'))
# 把"更正记录"放到 [V6] 引用块之后（它引用了逐设计 γ），先检查插入点唯一
out = txt

lines = txt.split('\n')
bts = out.replace('\n', '\r\n').encode('utf-8')
assert bts.count(b'\r') == bts.count(b'\n') == bts.count(b'\r\n')
import collections
ctrl = sorted(collections.Counter(x for x in bts if x < 32 and x not in (10, 13)).items())
assert ctrl == [(9, 3), (11, 1)], 'C 的控制字节变了：%s' % ctrl
now = open(P, 'rb').read()
assert now == raw, '板子又被改写（%d→%d），本次放弃' % (len(raw), len(now))
open(P, 'wb').write(bts)
chk = bts.decode('utf-8').split('\r\n')
assert len([x for x in chk if x.startswith('### 91-')]) == 5
assert '1.386×5' not in out and '-3.44/-1.72' not in out and 'note.tex:217$' not in out
print('OK lines=%d bytes=%d ctrl=%s' % (len(chk) - 1, len(bts), ctrl))
