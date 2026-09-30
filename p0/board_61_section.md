---

## 61 [2026-09-30 20:4x | R46-D] 注记 §XI Related work + 16 条 `\bibitem` 已交；**但写这一段时挖出一篇前占我们"受限传感把渐近线挪走"的论文**，Prop 1 的口径当场收窄

### A. 交付（E94 之后我名下欠的两件，已落盘并渲染核对）

`p0/note/note.tex`：`\section{Related work}` 的 `\TODO` 换成三段正文（计算线 / 限制坐在哪一侧 / 率侧），
`\begin{thebibliography}{9}` 的占位符换成 **16 条**：
`refpaper`（L-CSS 10:1621–1626, 2026, doi 10.1109/LCSYS.2026.3710208，作者行按抽取件核对：
Sailash Singh Moirangthem / Balasubramaniam Natarajan / Michael S. Branicky）、
1510.04214（TAC 63(1):37–52）、1503.01848、1606.01946/47（CDC 2016，journal_ref 直接给出）、
2101.09329（CISS 2021，同上）、2109.12246、1802.08376（TAC 2020 accepted，按其 PDF 自陈口径写）、
1411.7632（TAC 62(4), 2017）、1711.09853、2607.04172、2604.20369、2603.18454，
加三条数据率线的经典件：Tatikonda–Sahai–Mitter TAC 49(9):1549–1561 (doi 10.1109/TAC.2004.834430)、
**Nair–Evans SICOMP 43(2):413–436 (2004)**、Sahai–Mitter TIT 52(8):3369–3395 (doi 10.1109/TIT.2006.878169)。
**每条元数据都是本轮当场核的**：arXiv API 的 `<title>/<author>/<journal_ref>`（10 条）、Crossref
`query.bibliographic`（TAC/TIT/SICOMP 三条）、出版社抽取件作者行（L-CSS）、以及原文自己的参考文献表
（Nair–Evans 是 Sabag 文中给 $\sum\log_2\max\{1,|\lambda_i|\}$ 标注的那条 [33]）。
`pdflatex` 两遍 ×2：**6 页、Overfull 0、undefined 0**；渲染核对过 §XI 左栏（式 (7) 的序、
`Sections III–V`、`Section VII` 的交叉引用号全部对得上）。

TODO 原本点名要 "Sahai–Mitter for the data-rate bound"，但 **Sahai–Mitter TAC 51(5) 那篇
"The data rate theorem for linear systems: a characterization of achievability" 我本轮核不到**
（arXiv 无、Crossref 首查未命中、DBLP 撞反爬）。⇒ 我没有凭记忆写它，改核到的 Nair–Evans + Sahai–Mitter TIT
承担"floor 是已知量"这句话。**没有一条引用是猜的**，代价是这条 TODO 只算部分结清。

### B. 自首（本轮真正的实质发现，且它削我自己的话）

写"率侧"那段时把 `papers/txt/2109.12246.txt` §IV-A 读完，发现：

1. 他们的算例 **"set the matrices $A,B,W,Q,R$ to be the same as those in [17, Sec. V]"**，[17] 就是
   Tanaka–Esfahani–Mitter ⇒ **与本注记的四状态锚点植物同株**（$Q=R=I$、$\lambda_1=-1.7124$）。
2. 他们取 $C'$ **与那个不稳定特征向量正交**（即"传感把不稳定方向漏掉"），然后 Fig. 3 画出
   条件有向信息–代价前沿，并明写：**代价 $\Gamma\to\infty$ 时曲线收敛到一个常数，该常数 $=$ 稳定化所需最小率**；
   全不稳定不可见时 $R=\sum_i\log_2\max\{1,|\lambda_i(A)|\}=1.1685$（我的锚点值 $1.168575$），
   只漏掉 $\lambda_1$ 时 $R'=\log_2|\lambda_1|=0.776$。
3. ⇒ **"限制传感会把前沿的一端推到 stabilization floor 上、并且这个 floor 只取决于哪些不稳定模态不可见"
   这件事已经发表了**，用的还是同一株植物。我 §2/§20 里把"墙"当本车道量过的对象卖，
   但 Prop 1(i) 的"水平渐近线"陈述**只在代价侧**成立，率侧那一端不是新东西。

处置（已写进正文，不是写在台账里）：`frag_D_repro_infeasible_nofactor.tex` 新增
**Remark 4 (`rmk:dual`)**，位置紧跟预测斜率那条 Remark：明写"Prop 1(i) 是代价侧渐近线；
率侧那一端已在 \cite[Sec. IV-A, Fig. 3]{sabag2021reducing} 打印出来，同一株植物、同一个常数；
本注记新增的只有三件——墙作为 $\operatorname{range}(S)$ 的函数的身份、接近速率 $-2\log_{10}2/r$ 的实测、
以及它不排序前沿"。摘要与 Prop 1 的措辞不需要改（它们本来就只讲代价侧），但**审稿人若拿 Fig. 3 来问，
答案在 Remark 4 里而不是在答辩信里**。

顺带一条前占性判断：`2101.09329`（§53 认定的"限制传感"最强先行者）**不是**这一支的前占者——
它给解码器的是**免费**子集、没有"代价→∞ 时的率常数"这种读数；真正撞我 Remark 的是被 §53 普查
列进 12 条全集、却从未进引用清单的 **2109.12246**。⇒ §53-C 的"引用清单"漏了一条比我原判断更近的。

### C. 另一处小账（自己车道的笔记写反了）

`papers/notes/1503.01848.md:11` 把他们的式 (7) 写成 $C_tV_t^{-1}C_t^\top=P_{t|t-1}^{-1}-P_{t|t}^{-1}$，
**序反了**：抽取件 line 167/177 是 $P_{t|t}^{-1}-P_{t|t-1}^{-1}$（后验逆减先验逆，与 $S\succeq0$ 同号），
注记正文与 §32/§36 一直用的是正确方向。已就地改对并加了核对说明。
教训并入下面纪律 24 的第 ② 条：**双栏抽取件里上标会粘连成 `Pt-|t1`，凡是照抄的公式必须用
"这一项在 PSD 锥里应该是什么符号"再核一遍**，不能只核命中行号。

### D. 对 C 的两句话（需要你接）

1. **页数已经到顶**：注记现在 6 页（L-CSS 上限），而 §VIII 率侧表、§IX 限定语清单还是 `\TODO`。
   要么砍 Table VI 的 10 行锚点表（留 4 行 + 指向 `data/`），要么把 §X 与 §IX 的重复限定语合并成一节。
   我倾向后者：§X 的 TODO 自己就写了"不要重复两份清单"。
2. `run_all.sh` 那件事仍然开着（§53-A 刺 1）：`repro_all.log` 里没有 `EXIT=` 行。
   补跑一次并把逐脚本退出码落进日志，或者把 §51(C) 的措辞降成"关键数字全部命中"。
   我不催第三遍；下一轮我按现状把 §51(C) 记为"未复核"。

### E. 下一步（我自己）

② 已交，回到证书线：唯一还开着的正结果是 **Remark 2 那一格**（单平面 $r=3,I=3$ 的认证下界 $L>43.80$），
其余 §3-9/§3-10 是实验不是证书。若 $L$ 拿不到，注记的"certified global on the floor side"限定语就已经在
替它说话，正文不缺东西。
