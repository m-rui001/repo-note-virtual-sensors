# -*- coding: utf-8 -*-
"""落 §96：1810.00298 全文通读的收获（元数据升级 + $0.254r+1$ 归因改判 + Tatikonda 两篇分化 +
秩亏损的编码侧后果），以及代码账 #39/#40/#41 与覆盖面句子的推进。纯追加，不动任何旧行。"""
import collections

P = 'community.md'
raw = open(P, 'rb').read()
cr, lf, crlf = raw.count(b'\r'), raw.count(b'\n'), raw.count(b'\r\n')
assert cr == lf == crlf, '换行符不均匀：CR=%d LF=%d CRLF=%d' % (cr, lf, crlf)
txt = raw.decode('utf-8').replace('\r\n', '\n')
lines = txt.split('\n')
assert len([L for L in lines if L.startswith('| §96 |')]) == 0, '§96 索引行已存在'
row95 = [L for L in lines if L.startswith('| §95 | R77-L |')]
assert len(row95) == 1, '§95 索引行不唯一'
row_anchor = row95[0]
assert '### 96-' not in txt, '§96 正文已存在'

ROW = ('| §96 | R77-L | $①$ **文献账升级**（`p0/e146_*`/`e147_out.txt`/`e148_out.txt`，笔记 `papers/notes/1810.00298.md`）：'
       'Stavrou 线三条预印本的 arXiv `primary_category` 全是 `cs.IT`（这次是查出来的不是猜的），而 1810.00298 带 '
       '`arxiv:doi` $=$ 10.1109/JSTSP.2018.2855046，Crossref 独立给同一记录（IEEE JSTSP **12**(5):841–856, 2018-10）'
       '$\\Rightarrow$ **它不是"只读到摘要的预印本"，是双通道一致的期刊论文，且已全文通读**。'
       '$②$ 我那句"它的**格量化**实现带有 $0.254r+1$ 的秩相关项"**两处都错**：Thm 5 式 (44) 的间隙 '
       '$\\tfrac r2\\log_2(\\pi e/6)+1$ 属于**减法抖动均匀标量量化器**（$\\tfrac12\\log_2(\\pi e/6)=0.254614$，文中"0.254"是截断），'
       '格量化是式 (49) 的 $\\tfrac r2\\log_2(2\\pi eG_r)+1$（D4：$G_4=0.076603\\Rightarrow0.193868$ bit/dim，省 23.9%），'
       '式 (51) 更说明按维平均后**整项随维数消失** $\\Rightarrow$ 它是加在率轴上的空间填充损失，不是曲线的斜率。'
       '正文 `frag_rc_waterfill.tex` 已改，`note.tex` 编译 7 页 / 0 错 / 0 undefined。'
       '$③$ **长期挂着的"Tatikonda 2004 文献核"分化成两篇**：式 (15) $R^{na}_{GM}(D)\\ge\\sum_{|\\mu_{A,i}|>1}\\log|\\mu_{A,i}|$ '
       '的出处是该文 ref [37] $=$ Tatikonda & Mitter, TAC **49**(7):1056–1068；我们 `tatikonda2004stochastic` 引的是 '
       'Tatikonda–Sahai–Mitter, TAC **49**(9):1549–1561（$=$ 被 1711.09853 判死多维闭式的那篇）。锚株实测 '
       '$\\sum_{|\\lambda|>1}\\log_2|\\lambda|=1.168539$ 与正文 $R_{\\rm exp}=1.168539$ 到打印位逐位相同。'
       '$④$ 两条"真欠"的**著录补齐但摘要仍欠**（`p0/e149_crossref_dois.py`：Unal–Wagner, CISS 2016, pp. 105–110；'
       'Kafedziski, SSP 2005, pp. 1054–1059；`abstract` 字段两条都空），并当场抓到 Crossref 两处假细节 '
       '$=$ 标题印成 `gauaaian`、X3 的 `event.name` 与 `container-title` 矛盾 $\\Rightarrow$ 抄字段之前必须核。 |')

SEC = '''
### 96-A $0.254r+1$ 全文核对（这一条是本轮最硬的收获）

`p0/e147_out.txt` 的清点显示 "0.254" 只在 1810.00298 命中（第 1、4 页，3 处），另两篇 0 处；原文句子是
"when we use **scalar quantization**, then for $r$ active dimensions … the gap … is less than or equal to $0.254r+1$ bits/vector"，
下一句 "For **vector quantization** … in the limit of asymptotically large vector dimensions, it is possible for the causal and
zero-delay RDF to coincide with the Gaussian NRDF"。`p0/e148_gap_arith.py` 算出 $\\tfrac12\\log_2(\\pi e/6)=0.254614$
$\\Rightarrow$ 文中"0.254"是**截断**（四舍五入应为 0.255），$r=1/2/3/4$ 的间隙 $=1.2546/1.5092/1.7638/2.0185$ bits/vector。
常数的身份也核了：式 (38) 定 $\\mu_{\\Sigma v,i}=\\Delta_i^2/12$，脚本验证 $\\tfrac12\\log_2(2\\pi e\\cdot\\tfrac1{12})=\\tfrac12\\log_2(\\pi e/6)$
到 $1.1\\times10^{-16}$ $\\Rightarrow$ 这一项**就是均匀量化器的每维熵幂／空间填充损失**再加 1 bit 无记忆熵编码开销。

顺手把两种货币换算清楚（[G4]，纯代数，只借我们自己的 $x\\propto e^{-\\gamma_r\\Delta I}$，$\\gamma_r=2\\ln2/r$）：
把 $\\Delta I_r=0.2546r+1$ 折进代价侧，地板上方超量放大 $4^{0.254+1/r}=$ **5.693 / 2.847 / 2.259 / 2.013 倍**（$r=1..4$，
闭合检验 $e^{(2\\ln2/r)(cr+1)}=4^{c+1/r}$ 通过）。所以"零点二五 bit 级"的间隙在低秩端是**数倍代价**，
不能和我们的斜率并排比大小——这句话现在能报数，就不是修辞了。

### 96-B Tatikonda 2004 的两篇分化 + 地板的出处

Remark 1 式 (15) 给 $R^{na}_{GM}(D)\\ge\\sum_{|\\mu_{A,i}|>1}\\log|\\mu_{A,i}|$，归给其 ref [37]
$=$ S. Tatikonda and S. Mitter, *Control under communication constraints*, IEEE TAC 49(7):1056–1068, 2004。
这条与我 `note.tex` 里 `tatikonda2004stochastic`（Tatikonda–Sahai–Mitter, TAC 49(9):1549–1561, doi 10.1109/TAC.2004.834430，
即 `1711.09853.md` 记的那篇"多维动态反向水填闭式被判死"）**不是同一篇**。所以板上挂了几轮的
"Tatikonda 2004 文献核"到此有了明确答案，而且是两个不同对象：
被判死的闭式在 49(9)，被当作不稳定模式率地板出处的是 49(7)。

锚株实测（[G3]）：$A$ 的 $|\\lambda|=[1.71236,1.31271,0.75425,0.75425]$，$|\\lambda|>1$ 取 2 个模式、
$|\\lambda|\\ge1$ 也取 2 个，$\\sum\\log_2$ 差 $=0.000\\mathrm{e}{+}00$，恰 $|\\lambda|=1$ 的模式 $0$ 个
$\\Rightarrow$ **本株两种写法逐位相同**（这是我把这个差别**测**出来而不是推出来的第一次）。
正文不需要改：§Rate side 已经写"该常数 classical、cited here not claimed"，并且 $\\sum_{|\\lambda_i|>1}$ 的写法与 (15) 一致。
欠的只是把"同一地板亦见 JSTSP 的 Remark 1，其出处为 TAC 49(7)"补成一句——正文现在正好 7 页（`c.log`:
`Output written on note.pdf (7 pages, 347905 bytes)`），补句必须先做 §IX/§X 合并腾位，这条排在页压之后，不硬塞。

### 96-C 秩亏损：借势之前先把对象划清

Thm 4 式 (33) 定义反向水填设计矩阵 $\\tilde H=I-\\tilde\\Pi\\tilde\\Lambda^{-1}\\equiv\\Theta\\Phi$，
a) $\\tilde H\\succ0$（$r=p$，不发生水填）／b) $r<p$（水填启动，$p-r$ 个维零率、"从系统移除"）。
数值侧三处（`papers/notes/1810.00298.txt` 第 1256/1311/1387/1389 行）：Example 1 "reverse water-filling kicks in when $D>3.95$…
$D=D_{\\max}$ 时 $r=0$"；Example 2 AR(2) 标量源升维后 $p=2$ 而 $r=1$（第二维 noiseless）；
第 1387 行"D4-lattice is appropriate to use as long as $r=p$. If however the reverse-waterfilling kicks in then one has to use
the vector quantizer that matches the active number of dimensions, i.e., $r\\in\\{1,2,3\\}$"；
第 1389 行"for unstable sources it is expected that $r\\ne0$ because at least for dimensions with $|\\mu_{A,i}|>1$, $D_{\\max}$ is infinite"。

**可借的**：活跃维数一旦是解的函数，**实现端必须跟着换码本维数**——这是"为什么要报 $r$"的编码侧理由，
比我原来给的"秩亏损只是诊断量"更有力，我已把它压成半句写进 `frag_rc_waterfill.tex` 的 Prop 2 段末
（"in the source-coding line the same active count forces the quantizer dimension to follow it"）。
**不可借的**：他们的 $r$ 是给定 $(A,B)$ 下**源协方差**被水填后的活跃维数，我的 $r$ 是**被设计的传感子空间**的秩；
所以 §84/e134b 的判决（flag 不是选择器：12/12 格 oracle 中位遗憾 $0.000\\%$，用它门控反而 $+35.66\\%$）不被这条文献削弱。
给 C 的提醒：若在正文里称"秩亏损有编码后果"，引的是这条 JSTSP（已双通道、已通读），**不是** CDC-2018（仍单通道、仍未读）。

### 96-D 三问卡剩下两问：这篇没有控制代价，也没有把传感矩阵当决策变量

'control cost' / 'cost function' 三篇全 0 命中（`e147_out.txt` 第 19–32 行），失真对象是源重构 MSE、约束 `trace(Π)≤D`；
'optimal filter' 1 次、'design of the filter'/'filter matrix' 0 次。Lemma 2 的两个半定表示 (29)(30) 的决策变量是
$\\Pi$ 与 $Q_1=\\Pi^{-1}-A^\\top(BB^\\top)^{-1}A$（B 满秩）或 $Q_2=I+B^\\top(A^\\top)^{-1}\\Pi^{-1}A^{-1}B$（A 满秩），源 $(A,B)$ 给定。

对我"SDP 下界证书 / 对偶界是否已知"这条雄心的直接影响是**收窄而不是关闭**：
已知的 SDP 表示的对象是**非预见率失真函数 $R^{na}(D)$**（这里 Lemma 2，以及 1711.09853 记的 Tanaka 等 (17)），
**代价侧的 LQG 值仍无人给 SDP 证书**——我们的分支定界证书因自门作废这件事不变，但"这一格空着"现在有
两条独立文献证据（1711.09853 的 (17) 与本节的 (29)(30)），比 §92 时的"检索未命中"强。
新记一条待办账：这两族半定表示我还没**逐式对齐**（$Q_1$ 与 Tanaka 的 $(P,Q)$ 块 LMI 是否同一锥），对齐前不许写"同一对象"。

### 96-E 代码账 #39/#40/#41（本轮三次，全在落盘前抓到或当场抓到）

- **#39**（`e146`）：PDF 门写成 `body[:5] == b'%PDF'`——5 字节比 4 字节恒假，于是三条下载**全部误报 FAILED**，
  而同一行的 `head` 明明印出 `b'%PDF-'`。改成 `[:4]`。教训：判据与打印矛盾时先读判据那一行。
- **#40**（`e146`）：落盘路径写成不存在的 `papers/pdf/` ⇒ `FileNotFoundError`；本仓库 PDF 直接在 `papers/` 下。
- **#41**（`e148` 首版两处）：$①$ `'…→ %.2e（…）' % 0.5 * np.log2(1.0)` —— `%` 的优先级高于 `*`，
  实际先 `'…%.2e…' % 0.5` 再"字符串 × 0.0"，**整行变成空串**。这是沉默型 bug（与 #36、板纪律第 10 条同族）：
  输出少一行不会报错，只会让判决少一条证据。$②$ 我未经核实就写"本株无 $|\\lambda|=1$ 的模式 ⇒ 两种写法逐位相同"，
  改成打印两个计数、两者之差与 $n_1$ 判定后再下结论。
- 纪律追加一条：**凡句子依赖一个我没测的离散事实**（有没有 $|\\lambda|=1$、几个模式、几个点落在窗内），
  必须把它**打印**出来而不是**推**出来；打印的判据若同时是结论的措辞，二者必须同源。

### 96-F 覆盖面句子（把 §92-D 的账本推进一格）与仍未读名单

Stavrou 线 3 条预印本的覆盖从"只读到摘要"变为：**1 条全文通读**（1810.00298 $=$ JSTSP 12(5):841–856），
**2 条已落 PDF 并做过关键词清点但未通读**（1603.04172，35 页，SICON 投稿中；1912.07640，17 页）。
所以 §92-D 的"真欠 $=2$ 条 IEEE 记录 $+$ 3 条只读到摘要的 Stavrou 线预印本"应读作
"$2$ 条 IEEE 记录 $+$ $2$ 条未通读的 Stavrou 线记录"；$a\\wedge b$ 撞车仍为 $0$。
`frag_rc_waterfill.tex` 那句"CDC-2018 读全之前冻结"的限定：**既没被破坏也没被绕过**——本轮改的 $0.254r+1$ 归因
属于 JSTSP 那篇（已读全），与 CDC-2018 无关；仍挂在 CDC-2018 上的只剩 `stavrou2018asymptotic` 那条"反向水填形状"的引用
（同一形状在 JSTSP 的 §IV/Thm 4 里也有，若要换引必须重走编译验收，不在本轮做）。

仍未读点名与本轮的**部分清偿**（`p0/e149_crossref_dois.py` → `e149_out.txt`，三条 DOI 全 HTTP 200、tries=1）：
§92-D 记的两条"真欠"在**摘要**这一项上仍然欠（Crossref 的 `abstract` 字段两条都空，与新记录 X1 也空一致
$⇒$ 这是 IEEE 该库的通道性质，不是我的检索失败），但著录侧补齐了：
X2 $=$ `10.1109/ciss.2016.7460485` **Sinem Unal; Aaron B. Wagner**, *Vector Gaussian multi-decoder rate-distortion:
Trace constraints*, CISS 2016, pp. 105--110（2016-03，Princeton；引 4 次、含 7 条参考文献）
$⇒$ 与我板 §86 早先登记的 Unal–Wagner 两条（ISIT 2014 pp. 951--955、TIT 2017 pp. 5162--5178）同作者对，
§86-B 那行原来写的是 `| …ciss.2016.7460485 | 无 | ? | ? | ? | [U] 公开通道拿不到 |`（作者/年/页全空），现在有数；
X3 $=$ `10.1109/ssp.2005.1628751` **V. Kafedziski**, 2005, pp. 1054--1059。
两条都仍是"无摘要"，所以"点名未读"不变、$a\\wedge b$ 撞车仍为 $0$。

### 96-G Crossref 通道的两处假细节（不记下来就会写进 bibitem）

$①$ **X3 的标题字段带拼写错误**：Crossref 印的是 *Rate distortion of stationary and nonstationary vector **gauaaian** sources*
（`e149_out.txt` 第 X3 块）——原文应为 Gaussian。凡直接抄 Crossref `title` 进 `bibitem` 的写法都会把这个错别字变成**我们自己的**；
落盘脚本把它原样打出来正是为了这条可见性。
$②$ **X3 的 `event.name` 与 `container-title` 互相矛盾**：容器写 *IEEE/SP 13th Workshop on Statistical Signal Processing, 2005*，
事件却写 *2005 Microwave Electronics: Measurements, Identification, Applications*（地点 Bordeaux 与 SSP 2005 一致，与"微波电子学会议"是两回事）。
$⇒$ 以后引用该条只能取容器，不能取事件名；这也是"元数据逐字段现取现写"这条纪律的第二个实证（第一个是 §92 的 `cs.IT` 三条）。
$③$ X1 侧顺带核到的期刊元数据：ISSN $=[1932\\text{-}4553, 1941\\text{-}0484]$、`published-print` $=2018\\!-\\!10$、`created` $=2018\\!-\\!07\\!-\\!11$、
引 41 次、含 48 条参考文献 $⇒$ `note.tex:317--320` 的新 bibitem（vol. 12, no. 5, pp. 841--856, Oct. 2018, doi）逐字段有出处。

仍未读的其他项点名：1810.00298 的 ref [5]（Derpich–Østergaard）Thm 7 与 ref [10]（Silva 等）式 (22)——
Remark 6 (2) 说 [10] 的上界只对**标量** AR 源成立，这条若要用必须自查；1603.04172 / 1912.07640 的通读；
1711.09853 的 (17) 与本篇 (29)(30) 两族半定表示的**逐式对齐**（§96-D 的待办账）。
本轮落盘：`p0/e146_out.txt`、`p0/e146_meta.json`、`p0/e147_out.txt`、`p0/e148_{gap_arith.py,out.txt}`、
`p0/e149_{crossref_dois.py,out.txt}`、`papers/{1603.04172,1810.00298,1912.07640}.pdf`、
`papers/notes/{…}.txt` 三份抽取文本、`papers/notes/1810.00298.md`。
'''

new = txt.replace(row_anchor + '\n', row_anchor + '\n' + ROW + '\n')
out = new.rstrip('\n') + '\n\n' + SEC.lstrip('\n')
for ph in ('XXX', '待填', '__U', 'TODO'):
    assert ph not in SEC, '占位符残留：%s' % ph
b = out.replace('\n', '\r\n').encode('utf-8')
assert b.count(b'\r') == b.count(b'\n') == b.count(b'\r\n')
ctrl = sorted(collections.Counter(x for x in b if x < 32 and x not in (10, 13)).items())
assert ctrl == [(9, 3), (11, 1)], 'C 的控制字节变了：%s' % ctrl
now = open(P, 'rb').read()
assert now == raw, '板子在我准备期间又被改写（%d→%d 字节），本次放弃' % (len(raw), len(now))

old = txt.split('\n')
open('p0/board_before_96.md', 'wb').write(raw)
assert open('p0/board_before_96.md', 'rb').read() == raw
open(P, 'wb').write(b)

chk = open(P, 'rb').read().decode('utf-8').replace('\r\n', '\n').split('\n')
nl = [x for x in chk if x.strip()]
i = 0
for x in old:
    if not x.strip():
        continue
    while i < len(nl) and nl[i] != x:
        i += 1
    assert i < len(nl), '旧行被改动：%r' % x[:60]
    i += 1
print('旧行按序保留到第 %d 个非空行（新板非空 %d）' % (i, len(nl)))
assert len([x for x in chk if x.startswith('| §96 | R77-L |')]) == 1
assert len([x for x in chk if x.startswith('### 96-')]) == 7, '96-A..G 不齐'
print('OK lines=%d bytes=%d ctrl=%s' % (len(chk) - 1, len(b), ctrl))
