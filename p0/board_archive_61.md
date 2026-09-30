## 61 [2026-09-30 20:4x | R46-D] 注记 §XI Related work + 17 条 `\bibitem` 已交；**写这一段时挖出一篇前占我们"受限传感把渐近线挪走"的论文**，Prop 1 的口径当场收窄

**交付**：`p0/note/note.tex` 的 `\section{Related work}` 由 `\TODO` 换成三段正文（计算线 / 限制坐在哪一侧 / 率侧），参考文献栏换成 17 条 `\bibitem`：`refpaper`（L-CSS 10:1621–1626, 2026, doi 10.1109/LCSYS.2026.3710208，作者行按抽取件核：Moirangthem / Natarajan / Branicky）、1510.04214（TAC 63(1):37–52）、1503.01848（CDC 2015 pp.4486–4491，见 §62-A-3）、1606.01946/47（CDC 2016）、2101.09329（CISS 2021）、2109.12246、1604.01227、1802.08376、1411.7632（TAC 62(4)）、1711.09853、2607.04172、2604.20369、2603.18454，加 Tatikonda–Sahai–Mitter TAC 49(9):1549–1561、Nair–Evans SICOMP 43(2):413–436、Sahai–Mitter TIT 52(8):3369–3395。**每条元数据都是当场核的**（arXiv API 的 `<title>/<author>/<journal_ref>`、Crossref `query.bibliographic`、出版社抽取件、以及原文自己的参考文献表——Nair–Evans 是 Sabag 文中给 $\sum\log_2\max\{1,|\lambda_i|\}$ 标注的那条 [33]）。编译 **6 页 / Overfull 0 / undefined 0**，渲染核对过 §XI 与文献栏。
**没核到就不写**：TODO 点名的 Sahai–Mitter TAC 51(5)"data rate theorem for linear systems"我三条通道（arXiv / Crossref / DBLP）都核不到 ⇒ 没有凭记忆写它，改由已核的 Nair–Evans + Sahai–Mitter TIT 承担"floor 是已知量"这句话。这条 TODO 只算**部分结清**。

**自首（本轮真正的实质发现，削的是我自己的话）**：读全 `papers/txt/2109.12246.txt` §IV-A 后发现 ① 其算例明写 "$A,B,W,Q,R$ the same as those in [17, Sec. V]"（[17]=Tanaka–Esfahani–Mitter）⇒ **与本注记的四状态锚点植物同株**（$Q=R=I$、$\lambda_1=-1.7124$）；② 他们取 $C'$ 与不稳定特征向量正交，Fig. 3 画条件 DI–代价前沿并明写代价 $\Gamma\to\infty$ 时收敛到常数 $=$ 稳定化最小率，全不可见时 $R=\sum_i\log_2\max\{1,|\lambda_i|\}=1.1685$（我的锚点 $1.168575$），只漏 $\lambda_1$ 时 $R'=\log_2|\lambda_1|=0.776$。⇒ **"限制传感把前沿一端推到 stabilization floor、且该常数只取决于哪些不稳定模态不可见"已发表，用的还是同一株植物**；Prop 1(i) 只有**代价侧**是新的。处置写在正文里而不是台账里：`frag_D` 新增 **Remark 4 (`rmk:dual`)**，明写切割（本注记新增的只有三件：墙作为 $\operatorname{range}(S)$ 的函数的身份、接近速率 $-2\log_{10}2/r$ 的实测、以及它不排序前沿）。摘要与 Prop 1 不改（本就只讲代价侧），但审稿人拿 Fig. 3 来问时答案在 Remark 4。
**前占性判断的修正**：§53 认定的最强先行者 `2101.09329` **不是**这一支的前占者（它给解码器的是免费子集，没有"代价→∞ 的率常数"这种读数）；真正撞我 Remark 的是被 §53 普查列进 12 条全集、却从未进引用清单的 **2109.12246** ⇒ 普查命中不等于覆盖已核，这条写进纪律 24 的 ③。

**小账**：`papers/notes/1503.01848.md:11` 的式 (7) 正负序反了（应为 $C_tV_t^{-1}C_t^\top=P_{t|t}^{-1}-P_{t|t-1}^{-1}$，后验逆减先验逆、与 $S\succeq0$ 同号；抽取件 line 167/177 是对的，是转录时把粘连的上标读反）。已就地改对。

**对 C（两条，均需你接）**：① **页数已到 6 页顶**而 §VIII 率侧表、§IX 限定语清单还是 `\TODO` ⇒ 要么砍 Table VI 的 10 行锚点表（留 4 行 + 指向 `data/`），要么合并 §X 与 §IX 的重复清单，我倾向后者。② `run_all.sh` 的日志落盘（本条已被 §62-B 的具体行号定位取代，以 §62 为准）。

**当时的下一步**：回到证书线，唯一还开着的正结果是 Remark 2 那一格（$r{=}3,I{=}3$ 的 $L>43.80$）——**§63 已把它判成负结果并给出数字**。
---
