# 普查重跑：DI×LQG 对撞线的"全集"应当怎么数

执行：R39，2026-09-30。脚本 `p0/e97_census_union.py`（EXIT=0），输出 `p0/e97_out.txt`，原始数据 `p0/e97_census_raw.json`。
格式沿 `p0/results/novelty_waterfill_20260930.md`（主张 → 检索式 → 命中/未命中 → 缺口声明）；`新颖性审查_协议.md` 在本工作区仍不可定位（`find . -maxdepth 3 -iname "*新颖*"` 与 `/e/pdf/topics` 下 `-iname "*协议*"` 均无结果），这条缺口延续未关闭。

**本文件作废 `novelty_waterfill_20260930.md` line 35 的"12 条全集"表述。** 那条用的检索式 `all:"directed information" AND all:LQG` 是单条短语 × 单缩写，按我自己的纪律 26 不构成"全集"。

## 主张：旧清单漏了旗舰论文 1612.02126，且漏检机制可以实测

先把归因做成数字，而不是事后解释。三条对照式（同一 API、同一 `all:` 字段，只差词形）：

| 放开哪一轴 | 检索式 | totalResults | 1612.02126 是否命中 |
|---|---|---|---|
| 都不放（旧式） | `all:"directed information" AND all:LQG` | 12 | 否 |
| 只放信息轴 | `all:"directed mutual information" AND all:LQG` | **0** | 否 |
| 只放控制轴 | `all:"directed information" AND (all:LQR OR all:"linear quadratic")` | 12 | **否** |
| 两轴同放 | 见下 G1 | 14 | **是** |

裁决：漏检的真凶是**词形必须成对放开**。1612.02126 的标题是 "Rate-cost tradeoffs in control"，作者写 "directed mutual information"，控制侧只写 "control"/"quadratic cost" 而不写 LQG；单放信息轴时控制轴仍写死 LQG 所以归零，单放控制轴时信息短语仍写死所以那 12 条里不含它。**单独任何一条同义扩展都救不回旗舰论文**——这条以前我只是断言，现在是数字。

顺带修正 §64-A 的归因：当时是从抽出的**全文** grep 反推 `all:` 覆盖范围。实际 arXiv API 的 `all:` 只覆盖元数据（标题/摘要/comments/作者），不覆盖全文。结论方向没变，机制写错了。

## 宽口径并集（纪律 26：≥3 组同义扩展取并集，逐条落盘）

信息轴固定为 `directed information` ∪ `directed mutual information` ∪ `causal rate-distortion`，四个正交对象轴分别配对：

| 组 | 第二轴 | totalResults | 实抓 | [N2] 分页完整 |
|---|---|---|---|---|
| baseline | `LQG` | 12 | 12 | 是 |
| G1 | 控制/代价：LQR, linear quadratic, quadratic cost, control cost | 14 | 14 | 是 |
| G2 | 传感/观测：sensor, actuator, sensor selection, partial observation, task-based sensing, measurement | 334 | 334（200+134 两页） | 是 |
| G3 | 对偶/下界：converse, lower bound, dual, Lagrangian, semidefinite | 71 | 71 | 是 |
| G4 | 率失真：rate-distortion, rate distortion, data rate, communication rate | 32 | 32 | 是 |

宽口径并集 **401** 条（`e97_out.txt` 第 41–446 行逐条列举，含 G 归属标记）。
精筛 [N6]（标题+摘要须同时含信息词族与控制/代价词族）：**401 → 35**。
[N6] 双向自检：baseline 的 12 条**全部**存活（精筛没有误杀旧清单），1612.02126 存活（精筛接住了旧漏检）。

## 35 条精筛清单（本地覆盖状态 + 摘要关键词证据）

★ = 摘要含"对偶/下界 ∩ 秩/凸"，即与我方论断直接对撞。

| arXiv | 日期 | 本地 | 关键词证据 | 标题 |
|---|---|---|---|---|
| 1201.2334 | 2012-01 | 未入库 | exact | Universal Estimation of Directed Information |
| 1203.5572 | 2012-03 | 未入库 | – | Causal conditioning and instantaneous coupling in causality graphs |
| 1210.3889 | 2012-10 | 未入库 | – | Spatio-temporal Granger causality: a new framework |
| 1510.04214 | 2015-10 | PDF+note | sdp,exact | LQG Control with Minimum Directed Information: SDP Approach |
| 1604.01056 | 2016-04 | 未入库 | dual,exact | Capacity Achieving Distributions … LQG Theory of DI-Part II |
| ★1604.01227 | 2016-04 | 未入库 | lowerbound,sdp,exact,rank | Rate of Prefix-free Codes in LQG Control Systems |
| 1612.02126 | 2016-12 | PDF+note | converse,lowerbound | Rate-cost tradeoffs in control |
| 1702.06445 | 2017-02 | 未入库 | lowerbound,exact | Interplay Between Transmission Delay, Average Data Rate, and Performance |
| 1704.05370 | 2017-04 | 未入库 | – | On Information Transfer in Control Dynamical Systems |
| 1705.02802 | 2017-05 | 未入库 | exact | Directed Information as Privacy Measure in Cloud-based Control |
| 1711.08516 | 2017-11 | 未入库 | – | $k$-NN Estimation of Directed Information |
| 1803.08558 | 2018-03 | 未入库 | – | On Data-Driven Computation of Information Transfer for Causal Inference |
| 2001.03813 | 2020-01 | 未入库 | lowerbound | Fundamental Limits of Prediction, Generalization, and Recursion |
| 2003.04179 | 2020-03 | 未入库 | exact | Capacity of Continuous Channels with Memory via DI Neural Estimators |
| 2004.02356 | 2020-04 | 未入库 | sdp,exact | Scalable Synthesis of Minimum-Information Linear-Gaussian Control by Distributed Optimization |
| 2006.04683 | 2020-06 | 未入库 | – | Causal Structure Identification from Corrupt Data-Streams |
| 2101.09329 | 2021-01 | 未入库 | sdp,rank | Rate of Prefix-free Codes in LQG Control Systems with Side Information |
| 2109.12246 | 2021-09 | PDF+note | sdp,exact,rank | Reducing the LQG Cost with Minimal Communication |
| 2109.13146 | 2021-09 | 未入库 | sdp,rank | Dynamic Allocation of Visual Attention for Vision-based Autonomous Navigation |
| 2111.07401 | 2021-11 | 未入库 | – | Neural Capacity Estimators: How Reliable Are They? |
| 2203.11793 | 2022-03 | 未入库 | – | A Perspective on Neural Capacity Estimation |
| 2203.12467 | 2022-03 | 未入库 | lowerbound | A Lower-bound for Variable-length Source Coding in LQG Control |
| 2204.00588 | 2022-04 | 未入库 | lowerbound | Time-invariant Prefix Coding for LQG Control |
| 2301.00621 | 2023-01 | 未入库 | – | Data-Driven Optimization of Directed Information over Discrete Alphabets |
| 2403.09243 | 2024-03 | 未入库 | – | A simple reconstruction method to infer nonreciprocal interactions |
| 2510.02696 | 2025-10 | 未入库 | rank | Mutual Information-Driven Visualization and Clustering for Core KPI Selection |
| 2510.14096 | 2025-10 | 未入库 | rank | TENDE: Transfer Entropy Neural Diffusion Estimation |
| 2601.12782 | 2026-01 | 未入库 | – | Sensing-Limited Control of Noiseless Linear Systems Under Nonlinear Observations |
| 2602.09711 | 2026-02 | 未入库 | exact | Directed Information: Estimation, Optimization and Applications |
| 2604.20369 | 2026-04 | 未入库 | exact | Rate-Cost Tradeoffs in Nonlinear Control |
| 2606.21754 | 2026-06 | 未入库 | rank | EPSTE: Embedded Polygon Symbolic Transfer Entropy |
| 2606.31396 | 2026-06 | 未入库→本轮读 | lowerbound | How Much Sensing Information Is Needed to Control an Unstable Linear System? |
| ★2607.04172 | 2026-07 | 未入库 | lowerbound,sdp,exact,rank | Lower Bound of Networked Control with Multiple Sensors and One Controller |
| ★2608.26917 | 2026-08 | 未入库 | lowerbound,sdp | Minimum Rate For Partially Observable Linear System with Side Information |
| 2609.27580 | 2026-09 | 未入库 | exact,rank | Action-Directed Information for Distributed Control and Agentic Interaction |

本地覆盖：**3/35 已入库**（1510.04214、1612.02126、2109.12246 均为 PDF+note），2606.31396 本轮新增读数（`papers/notes/2606.31396.md`），其余 31 条仅摘要级。旧 baseline 的 12 条里 9 条未入库（`e97_out.txt` §5）。

## 本轮普查立刻付账的四个新条目

1. **2606.31396 / 2601.12782**（同组，作者表本轮 API 核实：2601 = Li/Liu/Xiong/Xu/Liu 五人；2606 = 同五人 + Yuan/Jin。同组假设成立，2601 是 2606 的 noiseless 前身）。Theorem 1 给出对**任意**观测律 $p(y_t|x_t)$（非线性、非高斯、秩亏）的必要性 $I(z\to y)\ge R_{\exp}=\sum_{|\lambda_i|\ge1}\log_2|\lambda_i|$。这条直接打掉我 §64-D 写的"文献中唯一发表过的显式部分观测反证界"。
2. **2004.02356**（G1×G3 交集，旧清单无）：摘要原文 "optimal control and sensing policies can be synthesized jointly by solving a SDP"。**传感与控制器联合合成的 SDP** 已有发表，任何"我们把传感设计并进 SDP"式表述都是二手。
3. **2608.26917**（2026-08-27，最新）：部分可观测 + side information 的最小率，条件 DI 下界 + 线性策略充分性 + 标例凸性。
4. **2607.04172**（2026-07-05）：多编码器无反馈网络的**首个** DI 下界 + SDP 公式化。

## [N7] 第二通道（Crossref）与跨通道一致性

三条 bibliographic 式各返回 100 条（API 单页上限），标题相关 98/16/58：
- CR1 `directed information LQG`：命中目标论文本体 `10.1109/lcsys.2026.3710208`（Moirangthem–Natarajan–Branicky，L-CSL 2026, 1621–1626，与 C 评的论文同题同号——通道可用性的正证）。
- CR1 命中 `10.1109/cdc42340.2020.9304490` = "The Minimal Directed Information Needed to Improve the LQG Cost"，作者 Sabag/Tian/Kostina/Hassibi，CDC 2020, 1842–1847。已核实这就是清单里 **2109.12246 的会议前身**（`papers/notes/2109.12246.md` line 3 早已登记）。所以第二通道这次**没有**捞出 arXiv 漏项，而是给 arXiv 清单做了一次独立复核：核心线（1510.04214 CDC16/TAC17、2109.12246 CDC20、L-CSL26）在两个索引里互指。
- 补充式 `au:"Sabag" AND all:"directed information"` total=4，全部已在清单内。

## 缺口声明

- G1：401 条宽口径里没有精筛掉的 366 条**未逐条人工判读**，只走了 [N6] 词形规则。规则可能误杀"标题用 information-theoretic 而不用 DI 词族"的论文；本轮未估计误杀率。
- G2：Crossref 三条查询都撞 100 行上限，`items-page` 未翻页 → 第二通道的全集**没有拿到**。已声明，不能声称 Crossref 侧穷尽。
- G3：31/35 条只有摘要级证据；1811.11792（G3 相关，PDF+txt 已落盘）**仍未通读**。
- G4：arXiv 侧漏检无法用 arXiv 自身完全审计（元数据盲区只能靠外部索引发现）。目前只有 Crossref 一个外部索引；DBLP / IEEE Xplore 未接。
- G5：`新颖性审查_协议.md` 不在可访问路径，格式仍是从两份既有检索表反推。

## 落盘件

`p0/e97_census_union.py`、`p0/e97_out.txt`（581 行，含 401 条逐条列举）、`p0/e97_census_raw.json`（per_query / union / precise 三段）。
