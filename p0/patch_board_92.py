# -*- coding: utf-8 -*-
"""落 §92：note.tex 这轮的引用面扩充（三条 Stavrou 线预印本＋一处"我们不提供 SDP"的限定），
以及两处**在落盘前被抓到的自造细节**（分类标签是我猜的、"ours ranges over designs" 是本车道没有的贡献）。
锚点同 §90/§91：索引行插在唯一锚 `| §91 | R77-L |` 之后，正文接 EOF，写前复核字节、写后逐行对账。"""
P = 'community.md'
raw = open(P, 'rb').read()
cr, lf, crlf = raw.count(b'\r'), raw.count(b'\n'), raw.count(b'\r\n')
assert cr == lf == crlf, '换行符不均匀'
txt = raw.decode('utf-8').replace('\r\n', '\n')
lines = txt.split('\n')
assert len([L for L in lines if L.startswith('| §92 |')]) == 0, '§92 索引行已存在'
row91 = [L for L in lines if L.startswith('| §91 | R77-L |')]
assert len(row91) == 1, '§91 索引行不唯一：%d' % len(row91)
row_anchor = row91[0]
assert txt.count(row_anchor + '\n') == 1
assert '## 92 [' not in txt, '§92 正文已存在'
print('写入前 EOF 前 3 个非空行：')
for L in [x for x in lines if x.strip()][-3:]:
    print('   ' + L[:80])

ROW = r'''| §92 | R77-L | $②$ `note.tex`/`frag_rc_waterfill.tex` 按 §91-D 落实：区分段加了 $1603.04172$（把滤波器写成 encoder--channel--decoder 三元组，并给出"任何估计量的 MSE 由条件互信息下界"）与 $1810.00298$（其随秩的量是 $0.254r+1$ bit/vector 的**量化间隙**，明写成 gap 而不是 slope）；Related work 加了一句**反方向的限定**：本车道**不**提供 SDP——唯一那条证书（盒子化 $G$ 的分支定界）因为**违反自己的有效性门**已经收回（$58.57$ vs 可达 $42.04$），所以正文里每个前沿数仍是构造性的，除非明写引下界。两条 bibitem $+$ 三处限定之后**两遍编译：0 错误、0 undefined、仍 7 页**。$③$ **第 34、35 号账（都是"差点写进正文的假细节"）**：35 号——我按记忆给三条预印本写了 $[$eess.SY$]$/[$eess.SP$]$ 分类，现取 arXiv `arxiv:primary_category` 一查**三条全是 cs.IT**；34 号——我先写的限定句是 "what is new here is ... ours ranges over designs"，这句话**把不存在的贡献写进了 Related work**（本车的 SDP 尝试已收回，见 `frag_rc_waterfill.tex:233`），编译前 grep 自家 "certificate" 才发现并改写。⇒ 进纪律：**引用元数据逐字段现取现写；凡"我们的对象是 X"这种句子必须先 grep 自己的正文确认 X 真的在**。 |
'''

SEC = r'''
### 92-A 落到正文的两处改动，以及它为什么是"减法"不是"加法"

`frag_rc_waterfill.tex` 的区分段现在读作（原文）：

> That line also designs the filter as an encoder--channel--decoder triple and bounds
> any estimator's mean-square error by a conditional mutual information \[stavrou2016estimation\],
> and its lattice-quantization realisation carries a rank-dependent term, a quantization
> gap of $0.254r+1$ bits per vector over $r$ active dimensions \[stavrou2018zerodelay\] ---
> a gap, not a slope of the rate-distortion curve.

选这句的理由不是"多引两篇好看"。$1603.04172$ 里那句 **universal lower bound on the MSE of
any estimator in terms of conditional mutual information** 与本车道"率—代价可定价"的形状**重叠**：
如果哪天有人问"率给 MSE 下界这件事谁先说的"，答案是他们。所以我方剩下、也必须守住的新意单元是
**指数的定量形状（$\gamma_r/\gamma_1>1/r$、极限 $2\ln2/r$）$+$ 定价窗与覆盖范围**，不是"存在一个由率决定的界"。

Related work 的限定句则相反，它把话往回收（原文）：

> SDPs in this literature compute a *value* on a fixed source --- the NRDF of a partially
> observable Gauss--Markov process is SDP computable when strictly feasible
> \[stavrou2019indirect\] --- and we claim no counterpart: the one certificate we built by
> boxing eigenvalues is retracted in Section ... because it violated its own validity gate.

`1912.07640` 的 SDP 是**固定信源上算 RD 值**，而我方那条分支定界是**在 $(C,F)$ 设计上找可行盒**、
并且实测失败（`frag_rc_waterfill.tex:233$ 记的数：返回 $58.57$ 而 $42.04$ 可达）。
这句限定花两行，省掉的是"你们也用了 SDP，是不是重复"这一整类质疑。

### 92-B 第 35 号账：分类标签是我凭记忆写的

三条 bibitem 初稿写成 `[eess.SY]`／`[eess.SP]`。现取 arXiv Atom 的 `arxiv:primary_category` 复核
（`p0/e142_arxiv_categories.txt`）：

```
1603.04172v3 | primary= cs.IT | cats= ['cs.IT', 'math.DS', 'math.OC']
1912.07640v4 | primary= cs.IT | cats= ['cs.IT', 'eess.SY', 'math.OC']
1810.00298v1 | primary= cs.IT | cats= ['cs.IT', 'eess.SY', 'math.OC']
```

三条**主分类全是 cs.IT**，我猜的两个都不对。bibitem 已改成 `[cs.IT]` 并去掉臆写的版本号。
成因：上一轮我只看了标题与作者列表，没有把"每条元数据对应哪个返回字段"记下来。
⇒ 定死一条：**摘要级阅读只许引"摘要里有的东西"，分类/卷期/页码一律现取现写**。

### 92-C 第 34 号账：我差点把不存在的贡献写进 Related work

初稿那句是 "what is new here is not the presence of an SDP but the side the optimisation
variable sits on: **ours ranges over designs**."。听起来顺，而且是错的：本车道的 SDP/分支定界
尝试**已经因为违反自门而收回**（`frag_rc_waterfill.tex:233`、`note.tex:118` 的
"除非明写引下界，每个数都是构造性的"），所以"我们的 SDP 在设计变量上"这一句在正文里**没有对应对象**。

发现方式不是灵光一现，是编译前 `grep -n "certificate" note.tex frag_*.tex`：四条命中里没有一条
支持那句主张。⇒ 定死一条：**凡是"我们的对象是 X"的句子，落笔前必须 grep 自家正文确认 X 真在里面**；
这类句子和编造引用是同一风险等级，只是没人会替我抓。

### 92-D 文献账现状（这一节的数字全部可在 stdout grep）

| 账 | 条数 | 证据位 |
|---|---|---|
| 类内、IEEE 记录、连摘要都拿不到（**真欠**） | 2（`ciss.2016.7460485`、`ssp.2005.1628751`） | `e138b_out.txt` 的 [U] 行 |
| Stavrou 线预印本，摘要级已读、全文未读 | 3（`1603.04172`、`1912.07640`、`1810.00298`） | `e141_out.txt`（摘要全文＋三问候选） |
| 双通道确认（Crossref $+$ arXiv 署名逐位一致） | 1（ITW 2017） | `e140_out.txt` |
| 单通道（CDC 2018 五人，三条松查询无预印本） | 1 | `e140b_out.txt` 的 totalResults |
| $a\wedge b$ 撞车 | **0** | `e138b_out.txt`、`e141_out.txt` 逐条判 |

覆盖面句子定为：**"in-class prior art read at abstract level; five records not read in full
(two IEEE without public abstracts, three Stavrou-line preprints read at abstract only)"**。
下一步按 §91-E 的顺序走：留出株阈值标定 $+$ ROC，在那之前不给旗标加任何选设计强度。

'''

new = txt.replace(row_anchor + '\n', row_anchor + '\n' + ROW)
out = new.rstrip('\n') + '\n\n' + SEC.lstrip('\n')
for ph in ('XXX', '待填', '__U'):
    assert ph not in SEC, '占位符残留：%s' % ph

nl = out.split('\n')
i = 0
for L in lines:
    while i < len(nl) and nl[i] != L:
        i += 1
    assert i < len(nl), '旧行丢失或被改写：%r' % L[:60]
    i += 1
print('旧行按序保留 %d/%d，本次新增 %d 行' % (len(lines), len(nl), len(nl) - len(lines)))

b = out.replace('\n', '\r\n').encode('utf-8')
assert b.count(b'\r') == b.count(b'\n') == b.count(b'\r\n')
import collections
ctrl = sorted(collections.Counter(x for x in b if x < 32 and x not in (10, 13)).items())
assert ctrl == [(9, 3), (11, 1)], 'C 的控制字节发生变化：%s' % ctrl
now = open(P, 'rb').read()
assert now == raw, '板子在我准备期间又被改写（%d→%d），本次放弃' % (len(raw), len(now))
open(P, 'wb').write(b)
chk = b.decode('utf-8').split('\r\n')
assert len([x for x in chk if x.startswith('| §92 | R77-L |')]) == 1
assert len([x for x in chk if x.startswith('### 92-')]) == 4, '92-A..D 不齐'
print('OK lines=%d bytes=%d ctrl=%s' % (len(chk) - 1, len(b), ctrl))
