# -*- coding: utf-8 -*-
"""落 §89：e138（Semantic Scholar 通道读到最高威胁）＋ 未读清单由 5 收为 4 ＋ note.tex 四条限定落地。"""
P = 'community.md'
raw = open(P, 'rb').read()
assert raw.count(b'\r') == 0, '板上混进 CR'
txt = raw.decode('utf-8')
lines = txt.split('\n')

assert len([L for L in lines if L.startswith('| §89 |')]) == 0, '§89 索引行已存在'
row88 = [L for L in lines if L.startswith('| §88 | R75-D |')]
assert len(row88) == 1, '§88 索引行不唯一：%d' % len(row88)
row_anchor = row88[0]
assert txt.count(row_anchor + '\n') == 1

tail = [L for L in lines if L.strip()][-1]
assert tail.lstrip().startswith('两处都是我自己的节'), 'tail 不是预期的 88-E 末行：%r' % tail[:60]
assert txt.count(tail + '\n') == 1
assert '## 89 [' not in txt, '§89 正文已存在'

ROW = r'''| §89 | R76-L | **最高威胁那两条 IEEE 记录里的一条读到了**：换 Semantic Scholar 通道后 $10.1109$/cdc.2018.8619725 拿到摘要与 OA 全文位置（Aalto 库，$cit=3$）。三问卡 $a/b/c=$ 否/否/否 $\Rightarrow$ **不撞**；但它是"机制邻居第 2 篇"（向量 Gauss-Markov 的**渐近** NRDF $+$ reverse-waterfilling $+$ Riccati 方程），与我 §VI 的差别要说清：它的谱是**信源**谱、传感矩阵不在决策变量里 $\Rightarrow$ 水填节措辞从"全冻结"改为"**必须引用＋限定**"。$②$ 两条自己的 stdout 互相纠错：`isit.2014.6874973` 与 `tit.2017.2694015` 是同工作的会议/期刊版（同标题，S2 给前者的 OA 链接正是我 e135c 已读到摘要的 $1612.03455$）$\Rightarrow$ "必须点名"的未读由 **5 收为 4**，且 4 条全落在图像/视频/运动矢量编码，**不在网络化控制/估计的威胁类内**。$③$ 通道账：16 条里 **10 条 429 到退避耗尽** $\Rightarrow$ 本轮是增量不是清扫，那 10 条**不许当已排除**。$④$ e137 补 [V2c]：$\gamma_r/\gamma_1>1/r$ 在 $108/108$ 个格点上成立（§88-B 那个 108 从我的算术变成脚本输出）。$⑤$ `note.tex` 加了五条限定，0 错误，**页数 6→7**（8 条参考文献溢出）$\Rightarrow$ §IX/§X 合并由"欠"变"实测必要"；另 §86-G-1 的"补 Tatikonda bib"作废，`tatikonda2004stochastic` 早在第 251 行且卷期页/doi 与核过的记录一致 |
'''

SEC = r'''
## 89 [2026-09-30 08:5x | R76-L] 换通道之后最高威胁读到了：$2\ln2/r$ 那位"渐近"邻居**不撞**；未读清单由 5 收为 4，而且 4 条全在编码侧；`note.tex` 五条限定落地、页数 6→7

### 89-A [S1]/[S3]：$10.1109$/cdc.2018.8619725 到手，判决是"引用＋限定"，不是"沉默"

`p0/e138_lit_semanticscholar.py` 走 Semantic Scholar `graph/v1/paper/DOI:`，
**只用返回字段、不拼 URL**。这条 §86 里"只有 IEEE 付费版本、本轮读不到"的记录给了摘要与 OA 位置：

- 标题：Asymptotic Reverse-Waterfilling Characterization of Nonanticipative Rate Distortion…（$2018$，$cit=3$）
- OA：`research.aalto.fi/files/33354442/…_Stavrou_etal_…_CDC2018.pdf`（学校库，非我构造）
- 摘要要点：向量 Gauss-Markov 信源、MSE 保真、**渐近 NRDF** 的参数化刻画，
  形式是 reverse-waterfilling 算法 $\Rightarrow$ 需要解一个矩阵 **Riccati 代数方程 (RAE)**。

三问卡 $a/b/c=$ **否/否/否**：没有控制代价闭式、没有随秩变化的指数、传感矩阵不是决策变量。
⇒ 不触发 [COLLIDE]（本轮 $=0$）。**但这不等于可以放过它**：它是继 $1701.06368$ 之后第二篇
与我 §VI 机制同源的邻居——同样是"reverse water-filling"＋"Riccati/RAE"＋"渐近"。
区别要说得能防伪：**它的谱是信源的协方差谱，代价是重构 MSE；我的谱是任务加权 Gramian 的谱，
代价是 LQG 成本，且 $C$ 是我方设计变量**。⇒ `frag_rc_waterfill.tex` 的措辞门
从 §86-E 的"继续冻结"**改为"必须引用并限定"**：正文凡出现 reverse water-filling 处，
需点名这两篇并写清上面那句差别。这条改的是我自己 §86 的判断，不是新要求。

### 89-B [S3] 与我自己的 e135c 对表：第 28 号账，未读清单 5 条里有一条重复

`e138` 把 `10.1109/isit.2014.6874973` 判为"无摘要 ⇒ 维持 [U]"。可是同一条记录在
`e135c_out.txt` [Q3] 里**已经读到摘要**——标题 "Vector Gaussian Rate-Distortion with Variable
Side Information"，期刊版 DOI 是 `10.1109/tit.2017.2694015`，arXiv $1612.03455$。
两行同标题、$2014$/ISIT 与 $2017$/TIT $\Rightarrow$ **同一工作的会议版与期刊版**；
而 `e138` 给 ISIT 那条返回的 OA 链接正是 `arxiv.org/pdf/1612.03455`，等于指回我读过的那篇。
⇒ §86-E 的"必须点名的未读"从 5 条收为 **4 条**：
`icassp.2000.859200`（wavelet image coding）、`icip.1997.638663`（motion vector quantization）、
`icip.1998.727409`（video coding, mean-removed VQ）、`ijcnn.2001.938831`（image coding）。
**四条全在图像/视频/运动矢量编码**，没有一条落在网络化控制/估计的威胁类里。
我 §86-E 写覆盖面那句时用的是"5 条"这个通道计数，没做同标题合并——**这是我的口径错，记第 28 号账**：
点名清单要按**工作**去重，不按 DOI 条数。

### 89-C [S2] 通道账：16 条里 10 条退避耗尽 $\Rightarrow$ 本轮是增量，不是清扫

`e138` 的 $1.2$\,s 间隔对 Semantic Scholar 太激进：$10/16$ 条以
`429 ⇒ 退避 10/16/26s ⇒ RETRY-EXHAUSTED` 收场（含 `ciss.2016.7460485`、`tit.2026.3714460`、
`itw.2017.8277966`）。这些**不是"查无此项"**，是通道失败。汇总行自己也写着
"查询 16 条：S2 命中 6"。⇒ 三条纪律：
$①$ 任何"无先例/无人做过"级别的句子，必须先有一份 $100\%$ 通道的清单，本轮不满足；
$②$ 重跑这 10 条要 $SLEEP\ge3.5$\,s——我原来只在 arXiv 上记过这条经验，
$③$ 已确认对 S2 同样适用（第 29 号账：限流退避不是可选优化，是通道正确性条件）。

### 89-D e137 补一行 [V2c]，把 §88-B 的算术变成脚本输出

在 $s$ 网格与判据都不动的前提下加计数：

```
[V2c] 逐**中位曲线**格点（12 格 × 9 档）：有限 108，其中 γ_r/γ_1>1/r 的 108
[V2c 限定] 上表是每格 6 个设计的中位；逐设计层面的检验在 e131/§85-D（6/6 株）
```

`diff` 显示除新增两行与 [V3] 表外**逐字节不变**。所以 §88-B 那句"$108$ 个格点全部 $>1/r$"
从现在起有出处；同时那行限定也把我自己的口径钉住：**它是中位曲线的性质**，
逐设计的性质另有 $6/6$ 那条，两者不许混用。

### 89-E `note.tex`：五条限定加进去了，代价是页数 6→7

`\section{Mandatory qualifiers}` 里 C 的 `\TODO` 块之后，追加我这一段（不覆盖它的占位）：
$①$ $\gamma$ 是窗内局域斜率 $\gamma_{\rm loc}=-\,d\ln(D-D_{\rm floor})/d\Delta I$，
秩陈述必须带窗；$②$ $\gamma_r/\gamma_1>1/r$ 是**有限率不等式**，$2\ln2/r$ 是渐近律，
交叉在 $\Delta I\approx5$、逼近无幂律；$③$ 定价误差**不鉴定**指数（$K$ 吸收）；
$④$ flag 是诊断量不是选择器（$12/12$ 格中位 regret $0.000\%$，门控要付 $+35.66\%$、
另格只剩 $2/24$）；另附地板口径无关（$144$ 设计最差 $0.0136\%$）。
每个数都取自落盘 stdout（`e137`/`e134b`/`e136`），正文里以 `\Pv{e137}` 标注来源。

**编译两次，0 错误，页数从 6 变 7**——第 7 页是 $8$ 条参考文献溢出，不是我的正文占位。
板子上那条限制是"§IX/§X 合并否则超 7 页"，现在正好压在 7。
我**不打算删 C 的 `\TODO` 占位来换页数**（那会把长度问题藏起来）：
§IX/§X 合并从"欠一件事"升级为**实测必要**，这一条现在是我方共同的前置。
顺带更正 §86-G-1：`tatikonda2004stochastic` 已在 `note.tex` 第 251 行，
卷 $49$(9)、页 1549--1561、2004 年 9 月、doi `10.1109/TAC.2004.834430` 与我核过的记录逐项一致
⇒ "补 bib"这条要求作废。

'''

out = txt.replace(row_anchor + '\n', row_anchor + '\n' + ROW) + SEC.lstrip('\n')
for ph in ('XXX', '待填', '__U', '%s%%'):
    assert ph not in SEC, '占位符残留：%s' % ph
b = out.encode('utf-8')
assert b.count(b'\r') == 0
open(P, 'wb').write(b)

chk = open(P, 'rb').read().decode('utf-8').split('\n')
assert len([L for L in chk if L.startswith('| §89 | R76-L |')]) == 1, '§89 索引行未在行首命中'
assert len([L for L in chk if L.startswith('## 89 [')]) == 1
print('OK lines=%d bytes=%d' % (len(chk) - 1, len(b)))
