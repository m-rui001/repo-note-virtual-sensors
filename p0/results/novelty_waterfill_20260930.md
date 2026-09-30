# 检索记录：§3-8 三项欠账（逆注水闭式 / 任务执行舞台 / 地板与不可行区间）

执行：R38-D，2026-09-30。格式按 `新颖性审查_协议.md` §7（主张 → 检索式 → 命中/未命中 → 缺口声明）。协议文件在本工作区不可定位（`find -maxdepth 3 -name "*新颖性*"` 无结果），故沿用 `p0/results/novelty_phi.md` 与 `通道预算与表示惩罚_注记.md` line 150–160 的既有表格格式。**这条本身要记为缺口：协议文档不在我能访问的路径里，格式是从两份既有检索表反推的。**

## 主张 A：给定任务子空间 $V$ 时，前沿设计 = 逆注水（闭式谱律）+ 秩开关门限 $I^*(k)$

| 检索式 | 渠道 | 命中 | 判定 |
|---|---|---|---|
| `waterfill` / `water-fill` / `reverse water` / `closed form` / `eigen` / `floor` in arXiv:1503.01848 全文（601 行） | 本地 PDF + pdftotext + `grep -ic` | 全 0 | **未覆盖**：该文只有五步 SDP 流程，秩是 SDP 输出 |
| 同上词族 in 本地 14 篇 txt（`papers/txt/`，2026-09-29 前存量） | 本地 grep | 0（当时 1503.01848 未落盘） | 见"缺口 G1" |
| `all:"reverse water-filling" AND all:"directed information"` | arXiv API（phrase 精确） | totalResults = **0** | 未命中 |
| `all:"water-filling" AND all:"sensor design"` | arXiv API | **0** | 未命中（对照式 `all:"water-filling"` = 464，证明检索式本身工作正常） |
| `all:"sequential rate-distortion" AND all:"sensor"` | arXiv API | 2：1411.7632（Tanaka–Kim–Parrilo–Mitter，SRD 的 SDP）、1606.01947 | 见下条 |
| **arXiv:1711.09853** Stavrou–Tanaka–Tatikonda 2017 全文（598 行） | arXiv PDF + pdftotext，本轮读完 §IV | **命中，且是反向命中**：They prove 动态反向水填对多维 Gauss–Markov SRD **不成立**（式 (15) 一般错），正确对象是 SDP 表示 (17) | **主张 A 的"闭式"部分在最近的文献里已被显式否证** |
| **arXiv:1606.01946** Fox–Tishby Part I 全文（915 行） | 同上，本轮读 Thm.1/§IV | **命中**：式 (10i)–(10k) 给出 EVD + "active mode coefficient matrix"（按特征值阈值化），式 (11) 把率写成 $\sum_i\max(0,\log\lambda_i)$，Definition 5 把阶数 $=\mathrm{rank}(D)$，§IV-B 报"第一个临界点阶数 0→1 的相变"，并自注"在 SRD 里水填是自洽的，$\lambda$ 依赖 $\Sigma$" | **主张 A 的谱阈值/秩亏损/相变结构在"memoryless 控制器限制"下已存在（2016）** |
| `all:"minimum-information LQG"` | arXiv API | 2（就是 Part I/II） | 他们答应给不稳定植物 $\nu\to0$ 分析的"upcoming paper"在 arXiv 上**不存在**（`au:"Fox" AND au:"Tishby"` 全集 6 条，无第三条） |

**裁决（对我自己）**：主张 A 从"定理级新结果"降级为"已知结构在传感侧限制下的重现"。可保留的差别只剩三条，且都不含"新闭式律"：E83 的 $\det\Gamma$ 不单调（含 $r=n$ / $r<n$ 代数二分）；TRV 锥内可达性见证（但构造等价于 1503.01848 的 (7)+Step 5）；不稳定植物的秩跳变门限（Tanaka 系只报数值现象、Fox–Tishby 限定稳定植物）。**注记里不得出现"前沿设计 = 逆注水"作为我们的定理。**

## 主张 B：$E[(g-u)^2]$ 型任务执行舞台（动态最小舞台）在 remote control 文献的对应

| 检索式 | 渠道 | 命中 | 判定 |
|---|---|---|---|
| `all:"task-based sensing"`、`all:"task restricted sensing" OR "task-restricted sensing"` | arXiv API | **0 / 0** | "task-restricted sensing"这个 L-CSS 自己的词在 arXiv 上零命中（与"L-CSS 不在 arXiv"一致，但也意味着我们引用它时无人可搜） |
| WebSearch：task-based sensor design + quadratic estimate + minimum rate + remote estimation | WebSearch | 返回 6 条候选（无线远程估计传感调度 3 条、IFAC `S2405896320302354` "Minimal Feedback Optimal Control of Linear-Quadratic…"、KTH 事件触发 1 条、Makridis 1 条） | **仅标题级**。IFAC 那条 403 无法抓取全文；其余三条与本舞台（任务映射 $\varphi$ 非单射 + 率-代价前沿）无直接对应 |

**裁决**：未命中，但缺口没关闭。Bansal–Başar 1989、Li–de Oliveira–Skelton 2008 仍是付费墙；俄语观测设计学派、functional-gain 一支仍未读。**主张 B 只能写"未检索到直接对应"，不得写"没人做过"。**

## 主张 C：地板 $D_{\min}(V)$ 与"完整不可行代价区间"

| 检索式 | 渠道 | 命中 | 判定 |
|---|---|---|---|
| `all:"infeasibility interval"`（+ `AND all:control`） | arXiv API | **0 / 0** | 未命中 |
| `all:"irreducible" AND all:LQG` | arXiv API | 1（量子引力，无关） | 未命中 |
| `all:"directed information" AND all:LQG`（普查整个对撞线） | arXiv API | **~~12 条全集~~ → 本行作废**：单条短语×单缩写，按纪律 26 不构成全集；重跑见 `p0/results/census_union_20260930.md`（宽口径并集 401 → 精筛 35，本行 12 条全部存活但漏了 1612.02126） | 其中 5 条是本仓库**尚未登记**的近邻（见下） |
| 摘要通读 `2101.09329`（Cuvelier–Tanaka, side information） | arXiv API summary | **解码器免费拿到状态向量的一个子集**，其余必须经离散信道；率-代价前沿；受约束 DI 的凸优化 + 可达性 | **现存文献里最接近"固定自由子空间 + 其余付费"的模型**，是 TRV 框架的退化特例，必须引用并说明差别（他们不定义被禁坐标的地板，也不给不可行区间） |
| 摘要通读 `2607.04172`（Li–Tanaka–Kim 2026） | arXiv API summary | 多编码单解码、无反馈网络的**首个 DI 下界**；线性独立编码在"传感同处即可见全状态"条件下最优 | 率侧下界的最新一支，与我们预算侧地板正交；须进引用清单 |
| 摘要通读 `2604.20369`（Atay–Chandrasekaran–Kostina 2026） | arXiv API summary | 非线性系统有限地平线：$F_n(D)\le R_n(D)\le F_n(D)+\log(F_n(D)+3.4)+2+1/n$ | "DI 是操作相关量"的非线性版；与我们开放问题 1（LQG 之外的舞台）直接相邻，须引 |
| 摘要通读 `2204.00588` / `2203.12467`（Cuvelier–Tanaka–Heath） | arXiv API summary | 时不变前缀编码逼近已知 DI 下界；共享随机性下的变长码下界 | 属于 §3-8 已登记的"前缀码率"线，补两条编号即可 |

**裁决**：$D_{\min}(V)$ 作为"被禁子空间导致的可达代价下界 + 整段不可行区间"仍然**未检索到对应陈述**；但 2101.09329 的 side-information 模型让"限制传感"这一 framing 的先占性下降，必须在注记里加限定语：我们相对它的新内容不是"限制传感"，而是**把限制写成 $\mathcal K_F$ 锥 + 给出区间长度 $\Phi(V)$**。

## 缺口声明（随主张同写）

- G1：2026-09-29 之前 `p0/results/novelty_phi.md` line 49 声称通读 1503.01848，但本地无 PDF/文本留存 → 该条记录不可复核。本轮已落盘（`papers/1503.01848.pdf`、`papers/txt/1503.01848.txt`、`papers/notes/1503.01848.md`）。
- G2：双栏 pdftotext 使 1503.01848 的 (6a)–(6f)、(14)、(15) 与 1606.01946 的 (10k)、(15)/(16)/(18) 字形粘连，凡要抄进正文的公式一律回 PDF 逐页核对。
- G3：Fox–Tishby 参考文献 [36]（被引为"水填效应的出处之一"）与 1711.09853 的 [11]（Kostina–Hassibi 的质疑）都无法从抽取文本可靠对应，引用前核 IEEE 元数据。
- G4：WebSearch 返回质量低（多次只回同批候选），IFAC 全文 403；IEEE Xplore / Google Scholar / Semantic Scholar 仍不可用。
- G5：`新颖性审查_协议.md` 在本工作区不存在（或被移走），格式靠反推。
