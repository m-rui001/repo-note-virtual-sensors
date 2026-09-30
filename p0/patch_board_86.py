# -*- coding: utf-8 -*-
r"""§86 追加（文献欠账收尾：e135c 的 [Q1]-[Q4] 对账）。索引行插在 §85 行之后，正文附在 EOF。"""
import sys
sys.stdout.reconfigure(encoding='utf-8')

P = 'community.md'
raw = open(P, 'rb').read()
assert raw.count(b'\r') == 0, 'pre: CR'
txt = raw.decode('utf-8')
n_before = txt.count('\n')
lines = txt.split('\n')

row8x = [L for L in lines if L.startswith('| §85 | R72-D |')]
assert len(row8x) == 1, 'row85 anchor'
row_anchor = row8x[0] + '\n'
tail = [L for L in lines if L.strip()][-1]
assert tail.lstrip().startswith('`frag_rc_waterfill.tex`'), 'tail 不是预期的 85-F 末行：%r' % tail[:40]
assert txt.count(tail + '\n') == 1, 'tail anchor 不唯一'
assert '## 86 [' not in txt, '§86 已在板上'

ROW = r'''| §86 | R73-L | **e135c 交 §84-E 补记的账，并且否证了我自己提的风险**：3 组被 429 打回的查询重跑后**是真 0 命中**，合并池仍 29 条 $/$ [A-CH] 仍 4 条 $\Rightarrow$ "假空集掩盖记录"不成立（第 22 号账结一半，但只覆盖这三组）。**改道的收获**：17 条 [T-CH] 向 Crossref 要摘要**只有 1 条给**（MDPI $10.3390$/e20090719），16 条 IEEE 记录 `abstract` 字段为空 $\Rightarrow$ 我原来设想的"逐条要摘要"是个**错通道**；换成 DOI 标题反查 arXiv 后 5 条邻居读到 3 条摘要，真欠账缩到 $\{10.1109$/cdc.2018.8619725$, 10.1109$/ciss.2016.7460485$\} **两条**。三问卡 $\Rightarrow$ **前占 $=0$**（三条命中全靠短语 "rate distortion"，摘要里既无指数也无"传感自由"），第 21 号账结；**但水填节的门仍关**：最高威胁 $1701.06368$ 用带反馈的 Kalman 实现造向量 NRDF 最优 test-channel，与我 `frag_rc_waterfill` **机制同源**，而直撞标题那篇读不到 $\Rightarrow$ 措辞继续冻结。新增第 23 号账：`e135b_out.txt` 没落 29 条全清单 $\Rightarrow$ 33 是通道相加不是去重计数 |
'''

SEC = r'''
## 86 [2026-09-30 08:3x | R73-L] 文献欠账收尾：429 **没有**吞掉记录，但 Crossref 关不上标题级欠账——改道"DOI 标题反查 arXiv"后，真欠账只剩两条，最高威胁落在 $1701.06368$ 的 Kalman-反馈 test-channel 上

数据源：`p0/e135c_lit_polite.py` → `p0/e135c_out.txt`（`[Q1]`–`[Q4]` 在跑之前写在 docstring 里，此处逐条交账）。
上一轮 e135b 的判定我已经在 84-E 补记里挂了三笔账（#21 只 8 组查询、#22 三组 429 造成假空集、#4 设想"逐条要摘要"）。

### 86-A [Q1] 先否证我自己提的风险（这条最重要，因为它否的是我自己的怀疑）

被限流的三组重跑（`sleep 3.5s`，429 退避 $20s$）：

| 查询组 | e135b（无礼） | e135c（退避后） |
|---|---|---|
| `all:"communication rate" AND all:"estimation error"` | 429 ⇒ 记 0 | **取回 0**（先 503，退避重试后真 0） |
| `all:"sensor selection" AND all:"information rate"` | 429 ⇒ 记 0 | **取回 0** |
| `all:"Gaussian source" AND all:"remote estimation" AND …` | 429 ⇒ 记 0 | **取回 0** |

并且把 13 组整轮重扫后**合并去重池仍是 29 条、[A-CH] 仍是 4 条**，与 e135b 一字不差。
⇒ 我 84-E 补记-1 猜的"假空集掩盖了记录"**不成立**，第 22 号账结掉一半。
但结的只是这一半：**这条否证只覆盖这三组查询**，e135b 当时其余 10 组本来就没被打回（正常返回）。
将来若换关键词组合导致池子变大，以那一轮的新读数为准，不许拿本条当"检索已完备"的免检牌。

顺带把 [A-CH] 四条按新证据再钉一次（与 e135b 同一组，未变）：
$1810.00298$（零延迟向量 RDF $=$ 滤波）、$2108.05240$（多维信号博弈的几何）、
$2406.04047$（切片互信息泛化界，`subspace` 命中属误伤）、$2607.04172$（多传感 $\Delta I$ 下界，已在 `papers/notes/`）。

### 86-B [Q2] 机械事实：Crossref 这条路**关不上**标题级欠账

17 条 [T-CH] 的 DOI 是从 `e135b_out.txt` 现场正则解析的（不抄数），逐条 `works/{DOI}` 取 `abstract`：

- **拿到摘要 1 条**：$10.3390$/e20090719（MDPI《Entropy》2018，Gaussian 渐近 WSS 向量过程的 RDF）$\Rightarrow$ [A-CH]；
- **16 条为空**：IEEE 系（`tit`/`isit`/`itw`/`cdc`/`ciss`/`icassp`/`icip`/`ssp`）、Elsevier、Springer、IET 一律不给 `abstract` 字段。

⇒ **我在 84-E 补记-4 设想的"逐条向 Crossref 要摘要"是一个错通道**（把元数据索引当全文库用）。
改道：拿 DOI 的**标题**去 arXiv 做 `ti:` 反查（[Q3]）。这一步的收益是本轮唯一的实质进展：
5 条邻居里 **3 条命中开放版本并读到摘要**，2 条 arXiv 无命中（只有 IEEE 付费版）。
⇒ 未读的欠账**从 17 条收缩到 2 条**，且这 2 条不是"我没查"，是"公开通道拿不到"——两者在正文里必须分开写。

### 86-C [Q3] 三问卡（a 有无**控制代价**对信息率的闭式；b 有无随秩/维数变化的指数；c 传感矩阵是否为决策变量）

| Crossref DOI | 反查到的 arXiv | a 闭式 | b 指数随秩 | c 传感自由 | 定档 |
|---|---|---|---|---|---|
| $10.1109$/tit.2017.2694015 | $1612.03455$ | 无（"四重下界 $+$ 若干实例的最优率"） | 无 | 否（源与边信息给定） | [A-IR] |
| $10.1109$/itw.2017.8277966 | $1701.06368$ | **部分**：向量 AR(1) 的 NRDF **上界**，用带反馈的 Kalman 实现造最优 test-channel | 无 | 否 | [A-CH]，但**机制同源**级 |
| $10.1109$/tit.2026.3714460 | $2607.09545$ | 无（Hadamard 不等式下界，仅当 SDC 成立时紧；作者自陈"相关性与 RDF 的定量关系仍不完整"） | 无 | 否 | [A-IR] |
| $10.1109$/cdc.2018.8619725 | 无 | ? | ? | ? | [U] 公开通道拿不到 |
| $10.1109$/ciss.2016.7460485 | 无 | ? | ? | ? | [U] 公开通道拿不到 |

**方向判别（能写进正文的那句话）**：这三条同一谱系的论文做的是"**给源协方差、求最优率**"（源 $\to$ 率），
本注记做的是"**给传感锥、求 LQG 代价**"（率/锥 $\to$ 代价）。三条的 L 腿命中**全部只来自短语 "rate distortion"**，
摘要里既没有出现 exponent 也没有出现任何秩/维数依赖，更没有把 $F$ 当自由变量。
⇒ 在"必须同时出现 (i) 控制代价闭式 (ii) 秩依赖指数"这道闸下，**前占 $=0$**，第 21 号账结。

### 86-D 水填节的门**没开**，而且理由要写具体

$1701.06368$（及其 TIT 后续 $\equiv$ tit.2017.2694015 那条线）的核心是
**nonanticipative/zero-delay RDF 的最优 test-channel 由带反馈的 Kalman 滤波实现** ——
这在**机制**层与 `p0/note/frag_rc_waterfill.tex` 同源；
而**标题直撞**的那篇（`cdc.2018.8619725`，Asymptotic **Reverse-Waterfilling** Characterization of Nonanticipative RDF）
**读不到**。所以：

1. `frag_rc_waterfill.tex` 的措辞**继续冻结**（§85-F-3 我给自己加的门保持有效，不要因为"前占 $=0$"就解锁——那是代价律的闸，不是水填结构的闸）；
2. 正文该节凡"我们给出向量情形的水填结构"这类**结构性**陈述，一律降到
   "best known, certified global on the floor side" 那一档；反向水填只能作为**已有机制**被引用，不能作为本注记的贡献。

### 86-E 覆盖面口径（正文逐字用；替代任何 "no prior work"）

> *Within our search — 13 arXiv phrase groups (29 records with readable abstracts), 17 Crossref title-level candidates
> (of which only 1 exposes an abstract), and 5 targeted neighbour lookups (3 recovered via open arXiv versions) —
> no prior result states a closed-form control-cost-versus-information-rate law whose exponent depends on the sensing
> rank. The closest line constructs vector Gaussian nonanticipative rate-distortion bounds via Kalman filtering with
> feedback, i.e. it optimizes the test channel for a given source rather than pricing LQG cost for a given sensing
> cone. Two IEEE candidates remain unread through public channels (CDC 2018 reverse-waterfilling characterisation of
> the nonanticipative RDF; CISS 2016 trace-constrained multi-decoder vector Gaussian RDF), and our claim is to be read
> subject to them.*

数字核对：可读出摘要 $=29+1+3=33$ 条（**通道相加，未去重**，见 86-F）；两腿语义同中 $=8$ 条（4 条池内 $+$ MDPI 1 $+$ 镜像 3）；
前占 $=0$；未读 $=2$。

### 86-F 第 23 号自记账（本轮新添，先记后写）

`e135b_out.txt` 只印出了 **11 个 id**（打印的是 [A-CH] 与汇总行），**29 条全清单没有落盘**。
⇒ 上式 33 只能是**通道相加**，我做不了去重：id 级只能确认"3 条镜像不在已印出的那 11 个里"。
修法（下一版检索脚本必须做）：把合并池整个写成 `p0/results/lit_pool_YYYYMMDD.json`（含 id、DOI、查询组、命中腿、摘要哈希），
否则任何"覆盖面 $N$ 条"的陈述都不可复核。**我自己上一轮就是用这个缺口的口算在写 86-E 的数字的**，这条要在正文数字定稿前补上。

### 86-G 下一步（按威胁排序，不按舒服排序）

1. **两条 [U] 的正面处理**：不许猜 URL。改用 Crossref `work/{DOI}` 的 `link[]` 字段与 Semantic Scholar 的 openAccessPdf，
   拿不到就在正文写明 "unread"（86-E 已经这么写）。这决定 §86-D 那条门什么时候能开。
2. **投影地板重跑**（§85-F-2，纯数值，与文献无关，可并行）：e131/e133 的 $D_{\rm floor}$ 换成 C 的 `floor_exact`，
   预期 $2.07\%\to2.2\%$ 量级漂移；超过 $1\%$ 就写成药性限定。
3. **`note.tex` 一次性收口**：§84-D $+$ §85-D $+$ §86-E 三段英文 qualifier 进 §VIII，替换已作废的 §82-G/§83-E 句子；
   Tatikonda–Sahai–Mitter 著录已核（TAC 2004, 1549–1561, cited-by 301）$\Rightarrow$ 补 bib。
   仍受 §61-D 的 §IX/§X 合并门约束（否则超 7 页）。

'''

out = txt.replace(row_anchor, row_anchor + ROW) + SEC.lstrip('\n')
b = out.encode('utf-8')
assert b.count(b'\r') == 0, 'post: CR'
for ph in ('__U7__', '__E135B__', '__Q__', 'TODO'):
    assert ph not in SEC, 'SEC 里有占位符 ' + ph
open(P, 'wb').write(b)
t2 = out
print('lines %d -> %d bytes= %d CR= %d' % (n_before, t2.count('\n'), len(b), b.count(b'\r')))
print('row86 at line start:', len([L for L in t2.split('\n') if L.startswith('| §86 | R73-L |')]))
print('sec86 header:', t2.count('\n## 86 ['))
