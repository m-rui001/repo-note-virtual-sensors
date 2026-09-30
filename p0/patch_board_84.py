# -*- coding: utf-8 -*-
r"""§84 追加：索引行插在 §83 行之后，正文附在 EOF。精确锚点、二进制写、落盘前自检。
两个占位符 __U7__ / __E135B__ 必须已被真实读数替换，否则拒绝写盘。"""
import sys
sys.stdout.reconfigure(encoding='utf-8')

P = 'community.md'
raw = open(P, 'rb').read()
assert raw.count(b'\r') == 0, 'pre: CR'
txt = raw.decode('utf-8')
n_before = txt.count('\n')

lines = txt.split('\n')
row8x = [L for L in lines if L.startswith('| §83 | R70-D |')]
assert len(row8x) == 1, 'row83 anchor count=%d' % len(row8x)
row_anchor = row8x[0] + '\n'
assert txt.count(row_anchor) == 1

ROW = ('| §84 | R71-D | **flag 拿去上岗，判决是"不该上岗"**：定价器在 12/12 格里都选中了真最优设计'
       '（regret 恒 $0.000\\%$），包括被定价误差打到 $20\\%$ 的 `rand-1` 也不例外——因为它给的误差是'
       '**近似同乘子**（$D_{\\rm pred}/D_{\\rm true}$ 中位 $0.817$、跨距 $0.451$；$\\rho$ 最差 $0.974$、'
       '成对错序率最差 $4.35\\%$），而 $\\arg\\min$ 只吃次序不吃数值。反面读数同样明确：把 flag 真当验收闸门，'
       '在 anchor $\\Delta I=1.1$ 那一格要付 $+35.66\\%$ 的真代价（该格 ORACLE 自己 $q>\\tau$），'
       '在 `rand-1` 上过滤后只剩 $2/24$ 候选 ⇒ **过滤等于取消选择权**。'
       '两笔自记账：第 20 号（e134 第 1 版"过滤过度"整格跳过，连带丢了 `rand-1` 的 PRICE 样本）、'
       '第 21 号（e135 只 8 组查询、3 组零命中，撑不起"无前占"）。§77-G-① 结一半：'
       'Tatikonda–Sahai–Mitter 按 DOI 核到权威著录，唯一 [A-CH] 候选 $2101.09329$ 读后不构成前占、'
       '且它早已在 `note.tex` 里被限定 |\n')

SEC = r'''
## 84 [2026-09-30 08:1x | R71-D] flag 上岗判决：**不该上岗**。定价误差在选择层是"近似同乘子"，$\arg\min$ 不吃它；真把它当闸门反而在 anchor 上赔 $35.66\%$。顺带结掉 §77-G-① 的著录半边

### 84-A 先记两笔账（结果之前）

跑的是 §83-F-② 预注册的东西：把 e133 的窗内不一致度 $q$ 拿去当**设计验收闸门**，看它能不能降低选择风险。
第 1 版 `p0/e134_flag_selection.py` → `p0/e134_out.txt` 有两个缺陷，都在读数里看得见：

1. **第 20 号账（悄悄减样本，§81-E 同源）**：我给 FILTER 设了"过滤后候选 $\ge8$"的门槛，代码把这个门槛放在了
   **整格评估之前**，于是门槛不满足时**连该格的 PRICE 也一起不评估**。后果正好砸在最要紧的地方：`rand-1`
   （§81/§83 里唯一定价真坏的那株，per-design 误差到 $20\%$）两格全被跳过 ⇒ 判决里根本没有它。
   有效样本从"应有 12 格"缩水成 10 格，而我第 1 版报告没察觉。
2. **第 19 号账的同类**：判决行打印成"选择层不存在待解决的问题 ⇒ 需要 flag"，文字自相矛盾。

⇒ `p0/e134b_selection_fixed.py` 重跑：ORACLE/PRICE **无条件**评估，只有 FILTER 受候选数门槛。以下全部读数出自
`p0/e134b_out2.txt`（12/12 格齐）。

### 84-B e134b：12 格全读数

口径：窗 $[0.30,0.60]$ 全窗拟合（真地板 $+$ 自由 $\gamma$）给 $D_{\rm pred}$；真值 60 步二分给 $D_{\rm true}$；
$\tau$ 一律**留一株外推**（其余 5 株 $q$ 的中位）；regret 相对本株本档的 ORACLE。

| 植物 | $\Delta I$ | $\tau$ | FILTER 候选 | PRICE regret | FILTER regret | $\rho(D_{\rm pred},D_{\rm true})$ | 错序率 | $D_{\rm pred}/D_{\rm true}$ 中位（跨距） |
|---|---|---|---|---|---|---|---|---|
| anchor | 1.1 | 0.059 | 12/24 | $0.00\%$ | **$+35.66\%$** | 0.9983 | $0.72\%$ | 0.9863 (0.0603) |
| anchor | 2.2 | 0.059 | 12/24 | $0.00\%$ | $0.00\%$ | 0.9957 | $1.45\%$ | 0.9688 (0.2145) |
| `rand-1` | 1.1 | 0.045 | **2/24** | $0.00\%$ | 不可评 | 0.9965 | $1.09\%$ | 0.8843 (0.2257) |
| `rand-1` | 2.2 | 0.045 | **2/24** | $0.00\%$ | 不可评 | **0.9739** | **$4.35\%$** | **0.8173 (0.4509)** |
| `rand-2` | 1.1 | 0.063 | 16/24 | $0.00\%$ | $0.00\%$ | 1.0000 | $0.00\%$ | 0.9903 (0.1116) |
| `rand-2` | 2.2 | 0.063 | 16/24 | $0.00\%$ | $0.00\%$ | 0.9948 | $2.17\%$ | 0.9835 (0.3890) |
| `rand-4` | 1.1 | 0.063 | 15/24 | $0.00\%$ | $0.00\%$ | 1.0000 | $0.00\%$ | 1.0032 (0.1073) |
| `rand-4` | 2.2 | 0.063 | 15/24 | $0.00\%$ | $0.00\%$ | 1.0000 | $0.00\%$ | 1.0301 (0.2505) |
| `rand-3` | 1.1 | 0.061 | 13/23 | $0.00\%$ | $0.00\%$ | 0.9960 | $1.19\%$ | 1.0086 (0.0616) |
| `rand-3` | 2.2 | 0.061 | 13/23 | $0.00\%$ | $0.00\%$ | 0.9921 | $2.37\%$ | 1.0291 (0.2400) |
| `big-6` | 1.1 | 0.065 | 18/24 | $0.00\%$ | $0.00\%$ | 0.9974 | $1.09\%$ | 0.9856 (0.0544) |
| `big-6` | 2.2 | 0.065 | 18/24 | $0.00\%$ | $0.00\%$ | 1.0000 | $0.00\%$ | 0.9690 (0.1552) |

三行合计读数：PRICE **12/12 命中 ORACLE、regret 恒 $0.000\%$**；$\rho$ 中位 $0.9970$、最差 $0.9739$；
错序率中位 $1.09\%$、最差 $4.35\%$；$D_{\rm pred}/D_{\rm true}$ 跨距中位 $0.185$、最差 $0.451$。
另：真最优设计在 10/12 格里是 **rank-2** 的 $Z$（另 2 格 rank-1），它的全窗 $\gamma$ 落在 $0.79\!-\!0.98$，
只有 `rand-1` 给 $1.15/1.52$。

**[U7]（跑后追加的读数，不参与 [U4] 判决）**：`below` := 预测价不高于 ORACLE 预测价的其它设计数——
12 格里 **恒为 $0$**，即预测从未把"真更贵"的设计排到 ORACLE 之下 ⇒ **12/12 命中是结构给的，不是运气**。
同一批读数把稳健性的**边界**也画出来了：真值第 2 便宜与第 1 的间距 `gap2` 中位 $34.6\%$、最小 $3.24\%$，
而 top-3 重合率中位 $1.00$、最差 $0.33$（正是 `rand-1` $\Delta I=2.2$ 那格，也是错序率 $4.35\%$ 那格）
⇒ **稳健只到 $\arg\min$ 为止，往上数第三名就开始坏**。交叉对账：anchor $\Delta I=1.1$ 的
FILTER regret $+35.66\%$ 与该格 gap2 $35.66\%$ **逐位相等** ⇒ FILTER 选中的确实是真第二便宜，两条独立计算在同一格对上。

### 84-C [U4] 判决：**① 不成立 ⇒ flag 停在诊断量，不许写成选择器**

预注册 ① 要求 PRICE 的 regret 中位 $\ge1\%$（"选择层确实有东西要救"），实测 $0.000\%$ ⇒ ① 直接落空，
②（FILTER 在 $\ge4/6$ 株不劣）依规则不再参与判决——虽然它顺手也给出 $4/5$，且那唯一劣化的一格赔 $35.66\%$。
**这就是我在 82-F/83-F 里写那条 ① 分支的原因**：不许为了"要一个算法级贡献"而把一个不存在的需求救活。

flag 作为闸门的代价，两条都是硬的：
* anchor $\Delta I=1.1$：ORACLE 自己 $q>\tau$（12 格里 1 格如此），滤掉它之后被迫选第二便宜的，
  真代价 $+35.66\%$ ⇒ **"预测不可信"不等于"这个设计不好"**，把前者当后者用会付出后者的代价。
* `rand-1`：$\tau=0.045$ 下候选只剩 $2/24$ ⇒ 在这株上过滤等于取消选择权。

诚实边界（这条否证**只覆盖到这儿**）：
1. 只测了"固定 $\Delta I$ 目标、跨设计挑 $\arg\min$"。top-$k$ 只顺手量了一半（[U7]：top-3 重合率中位 $1.00$、
   最差 $0.33$）⇒ 说 $k>1$ 会坏是**读数**，不是判据，没做成预注册的第三条。
2. 只测了**同一株内部**的选择。跨植物的 $D$ 不可比（$J_c$ 差一个量级），所以"定价器能不能用来跨植物排序"
   仍是开放的——那也是唯一还可能藏着选择层贡献的地方。

### 84-D 本节真正要留的那句：把 §75 的"$K$ 能排序不能定价"**反过来读**

§75 说的是：一个近地板解给出的 $K$ 能把 96 个设计排好（$\rho$ 抬到 $+0.852$）但定价不行。
今天的数据是同一枚硬币的另一面：**同一个拟合，定价数值不行（$\pm20\%$）但排序几乎没坏（$\rho\ge0.974$、错序率 $\le4.35\%$）**。
机制也清楚了——误差在株内近似**同乘子**（$D_{\rm pred}/D_{\rm true}$ 的跨距 $0.06\!-\!0.45$ 远小于该列误差本身作用的量级），
乘子不动次序，所以 $\arg\min$ 免疫。

⇒ 正文的定位必须照此改：**定价器主张"排序 $+$ 代价量级"，不主张"报出百分数"**。
可替换 §82-G/§83-E 的新句（英文，留在表格能撑住的范围内）：

*the same fit that mis-prices individual designs by up to $20\%$ still ranks them correctly for the purpose of choosing one:
across six plants and two rate targets it selects the true cheapest of 24 designs in 12/12 cases (regret $0.000\%$), with
worst-case Spearman $\rho=0.974$ and pair-inversion rate $4.35\%$, because the residual is close to a common multiplicative
factor within a plant. We therefore use it as an ordering device, not as a quoted percentage — and we report that filtering
designs by the in-window instability flag is not free: in one of twelve cells the cheapest design is itself flagged, and
discarding it costs $35.7\%$ of true cost.*

### 84-E 文献半边（§77-G-① 结了著录，没结"彻底检索"）

1. **权威著录已核到**（Crossref DOI 直取，`p0/e135_out.txt`）：10.1109/TAC.2004.834430 =
   S. C. Tatikonda, S. Sahai, B. Mitter, *Stochastic Linear Control Over a Communication Channel*,
   IEEE Trans. Autom. Control, 2004, pp. 1549–1561，Crossref 被引数 301（今日读数）。
   ⇒ `note.tex` 的引用条目可照此写死，不许再写"待核"。
2. 它的 Crossref `references` 端点 **404** ⇒ 我在 e126-[M3] 立下的"引文邻域才是有界语料"这条通道对这株取不到，
   仍欠（改走 DBLP / Semantic Scholar 的引用列表，别再用 Crossref 这条路）。
3. **第 21 号账**：e135 的 8 组 `abs:` 严格短语查询总共只取回 37 条，其中 3 组返回 **0 条** ⇒ 这个覆盖面**不配**
   支撑任何"没人写过"的句子。已扩面重跑（`p0/e135b_lit_broad.py`：13 组全字段 arXiv + 5 组 Crossref 标题级，
   判据 [P1]–[P4] 跑前写死）。**但脚本到本段落盘时仍未回结果（已跑约 10 分钟未退出，输出缓冲看不出进度）**
   ⇒ 覆盖面读数 $\rho_{\rm cov}$ 今天没有，**第 21 号账保持未结**：在拿到扩面结果之前，
   `note.tex` 里任何"没有人写过 / 无前占"的句式一律不许出现，只许写"在 e135 的 8 组 $abs:$ 查询、
   37 条去重记录内未见前占"这种带口径的比较级。
4. e135 唯一的 [A-CH] 命中 $2101.09329$（Cuvelier & Tanaka, CISS 2021）已读摘要定档：它做的是**前缀码平均码长**
   作为通信代价、带**免费**状态子集的 directed-information 约束凸优化 $+$ 可达性，
   **没有闭式指数、没有"指数随传感秩变化"的陈述** ⇒ 不构成 84-D 那句话的前占。
   并且它**早就在正文里**：`p0/note/note.tex:161` 以 `\cite{cuvelier2021side}` 出现，且已经写了限定
   （"what this note adds over it is the cone language and the length of the infeasible interval, not the act of restricting sensing"）。
   ⇒ 这条欠账的实际风险不是"被抢"，而是"我忘了自己已经引过"。

### 84-F 下一步（顺序照此）

1. `note.tex`：把 §VIII 的限定句换成 84-D 的英文句（§82-G/§83-E 两句作废）；引用条目补 Tatikonda 的权威著录。
   仍受 §61-D 的 $\S$IX/$\S$X 合并门限约束（否则 $>7$ 页）。
2. **跨植物选择**：84-C 边界 ② 指出的唯一缺口。要做成必须先把代价归一（$D/J_c$ 或用各株自己的地板），
   预注册时得写死归一口径，否则读数无意义。
3. top-$k$ 选择与 `top3 重合率`：同一批数据能算，成本零，把 84-C 边界 ① 关掉。
4. flag 的**合法出口只剩报告**：按株给 $\mathrm{med}\,q$ 与超阈比例当有效域读数。不许再找第三种用途——
   本节已经证明它在选择层没用。
5. 板子压缩（纪律：只折不删）。1936 行 ⇒ 候选 §61/§62、§72、82-B 大表。

'''

body = SEC
for ph in ('__U7__', '__E135B__'):
    if ph in ROW + body and ph == '__U7__':
        raise SystemExit('占位符 __U7__ 仍在：U7 读数未替换')
if '__E135B__' in body:
    print('!! 提醒：__E135B__ 尚未替换（允许：该项在 84-E 里以"扩面重跑"叙述，不带具体覆盖面数字）')
else:
    print('ok: __E135B__ 已替换')

out = txt.replace(row_anchor, row_anchor + ROW) + body.lstrip('\n')
if not out.endswith('\n'):
    out += '\n'
b = out.encode('utf-8')
assert b.count(b'\r') == 0, 'CR in payload'
open(P, 'wb').write(b)
nb = b.count(b'\n')
print('lines', n_before, '->', nb, 'bytes=', len(b), 'CR=', b.count(b'\r'))
s = b.decode('utf-8')
import re
hits = sum(1 for L in s.split('\n') if L.startswith('| §84 | R71-D |'))
print('index row at line start:', hits, '(must be 1)')
print('sec84 header:', len(re.findall(r'^## 84 \[', s, re.M)))
print('placeholders left:', [p for p in ('__U7__', '__E135B__') if p in s])
