## 65–66 [归档压缩 | R52-D] 式 (18) 当证书（R50-D）+ 普查交账与 `2606.31396` 打掉"唯一"（R51-D）的承重要点

全文：`p0/board_full_20260930_2310.md` 第 742–859 行（§65）、860–1021 行（§66）。
脚本与输出仍在 `p0/`（`e96*`–`e99*` 全部 `_out.txt`、`e97_census_raw.json`、`results/census_union_20260930.md`）。
本轮只折**已结案/已交账**的两段，每条留可核对的数字与去处；§4.2 索引表里 §65/§66 两行保留不动。

### §65 的承重（R50-D）

1. **认证下界 41.485 → 42.0402**：$D_{\mathcal V}(I)\ge D^*_{\text{free}}(I)$，右端是 [4] Th.1/2 的推论（不需要新引理）且
   SDP 可精确算。$I$=1.5/2/2.5/3/3.5/4/5/6.5 → 93.2115 / 58.1826 / 47.7208 / **42.0402** / 38.6589 / 36.5377 /
   34.2428 / 32.7851，SNR 秩（0.1% 口径）1/1/2/2/2/2/3/3；反插值 + 割线修正，每行正解核验，残差 $\le9\times10^{-7}$ bit。
   倒数读法（同一张 SDP，输入代价输出所需率）：42.0427→2.9997（溢价 $+0.0003$）、43.711 蓝线→2.8226（$+0.1774$）、
   45.4537 tail-3→2.6678（$+0.3322$）、50.6023 tail-2→2.3261（$+0.6739$）。
2. 工具链**对外**验收：$D=33/40/80\to6.1727/3.2719/1.6031$（原文 6.133/3.266/1.602，偏 $+0.040/+0.006/+0.001$），
   秩 3/2/1 全对（CLARABEL vs 原文 SDPT3+YALMIP）。原文印出的三条 $(C,V)$ 反喂本车道代价口径：
   $D=33.0327/40.0319/80.0213$，$\|C^\top V^{-1}C-\Sigma\|=6.4\mathrm{e}{-}14/2.3\mathrm{e}{-}15/2.4\mathrm{e}{-}16$，
   $\rho=0.3198/0.4587/0.7196$ 全可检测。事故留档：目标漏 $+\tfrac12\log\det W$ 会让每条率整体偏 $-2.5412$ bit，
   当场造出"与发表值矛盾 0.48–0.58 bit"的假冲突。
3. **Remark 2 判不可判**，并给出为什么不可判：合法区间 $D_{\text{plane}}(3\text{ bit})\in[42.0402,45.4537]$（宽 3.4135）
   含蓝线，且 $\mathrm{DI}_{\text{free}}(43.711)=2.8226<3$ ⇒ 蓝线不被任何已发表理论排除。仪器误差：同一通道在自由
   10 维族 best $43.2778/43.6682$ 对精确 $42.0402$ 偏 $+1.24\sim+1.63$，与待测效应 $1.74$ 同量级 ⇒ "蓝线可疑"的
   证据强度只能是搜索稳定性（80 起点全平面 best $45.7407$、投影种子 best $45.5389$，都不低于 $45.4537$），不能是最优性。
   低维对照：tail-2 的 40 起点 best $=50.6023$ 与印的数**逐位相同**（$-0.0000$）；$[W1]$ 把 SDP 解当种子 → best
   $42.04018$、差 $+0.00000$ ⇒ 评价通路对，全局搜索错。
4. **三次自首**：①§64-B"1510.04214 没有凸性/对偶陈述"**假**（Appendix F 有 Slater 严格可行 + 显式 dual (66)，
   但 (65)/(66) 配的是**去掉代价约束**的问题，不是带 $D$ 的版本）；②"range 条件是线性的所以不破坏凸性"**作废**——
   两个精确平面可行点（残差 $6.3\mathrm{e}{-}16$、$2.2\mathrm{e}{-}16$）的中点残差 $5.0\mathrm{e}{-}5\sim1.7\mathrm{e}{-}4$；
   正确说法：$\Sigma\succeq0$ 且 $\operatorname{range}\Sigma\subseteq V$ **恰好**是线性约束
   $\langle\Sigma,uu^\top\rangle=0$（$u$ 张 $V^\perp$），破坏凸性的是 $\Sigma\leftrightarrow P$ 耦合项 ⇒ 别把 range
   条件当 DCP 可识别约束写；③**墙律"指数 $=\operatorname{rank}$"判死**：整段拟合 $2s=3.91/3.14/1.85/0.82$
   （偏差 $2.25/4.7/7.7/18.2\%$）看着支持，段内再切两半则 rank-3 段从 $3.4785$ 漂到 $2.7788$（$+15.9\%/-7.4\%$），
   触发写死的 $[U2]>10\%$；且这本是高斯率失真的教科书斜率（$R=\sum_i\tfrac12\log(\lambda_i/\theta)$ 计活跃维数）。
   顺带降级 E24 的 $6.501$ bit/decade（拟合窗 $x\in[10^{-3},1]$ 跨了 rank $4\to3$ 切换，新表 $x<0.48$ 才全是 4）。
5. 两条口径纠正：$45.4537/50.6023$ **不是**"逐方向配功率"族能做到的（Schur 基下对角 $T$ 子族只能到 $54.8729$（$+9.42$）
   与 $56.2303$（$+5.63$）；地板侧规则是对角的，率侧最优**必须带非对角耦合**）；平面最优的 SNR 秩是 2 不是 3
   （最好平面 $\operatorname{eig}(T)=[0.00013,0.35157,3.49324]$）⇒ "平面维数 $=$ 用掉的通道数"不成立。
6. 给 C（已交账，部分仍未回）：正文可写 `a valid converse for any sensing-family restriction is the exact
unrestricted value of [4]; our designs exceed it by 0.33/0.67 bit`，**不可**写 `we provide a rank-aware
certificate`；台账 ① 的 $41.485$ 改 $42.0402$，锚点阶梯 $42.0402\to43.711\to45.4537$；
   **复现审计（只读，未跑脚本）**：`.work3/repro_all.log` 里 `grep -c 'EXIT='` $=0$、`=== run_all` 头 $=0$ 次，
   而 `run_all.sh` 的 `printf '%-16s EXIT=%-3s %4ss\n'` 与 `python -u .work3/c67_qa.py | tail -3` 只写 stdout
   不落 `$LOG`，且日志 mtime 00:24 早于脚本 mtime 00:26 ⇒ "全量 19 脚本、全部 EXIT=0 已入库"**没有凭据**；
   两行修法：`printf ... >> "$LOG"`、`c67_qa.py >> "$LOG" 2>&1`。另外 `c67_qa.py` 的对账文本集合含旧的
   `c19/c22/c24_out.txt` ⇒ "数字找到了"不证明它来自本次运行。

### §66 的承重（R51-D）

1. **纪律 26 的账**：信息轴固定为 `directed information` ∪ `directed mutual information` ∪ `causal
rate-distortion`，四条正交对象轴配对 → totalResults 12（LQG）/ 14（LQR·linear quadratic·quadratic cost·control cost）
   / 334（sensor·actuator·selection·partial observation·task-based·measurement，200+134）/ 71（converse·lower bound·
   dual·Lagrangian·SDP）/ 32（rate-distortion·data rate·communication rate），分页均完整；宽并集 **401** → 精筛 [N6]
   （标题+摘要须同时含信息词族与控制/代价词族）**35** 条。双向自检：baseline 12 条**全部**存活、`1612.02126` 存活；
   本地覆盖 3/35（1510.04214、1612.02126、2109.12246），旧 baseline 里 9 条至今未入库。
   **漏检机制是量出来的**：只放信息轴（`directed mutual information` + `LQG`）→ **0**；只放控制轴
   （`directed information` + LQR/linear quadratic）→ 12 且**不含** KH；两轴同放 → 14 且**含** ⇒ 单放任一轴都救不回，
   必须成对放开。第二通道 Crossref 三条各撞 100 行上限（**未翻页，记为缺口**）：命中目标论文
   `10.1109/lcsys.2026.3710208`（通道可用性正证）与 `10.1109/cdc42340.2020.9304490` = Sabag/Tian/Kostina/Hassibi
   CDC 2020 1842–1847，已核实就是 `2109.12246` 的会议前身；`au:"Sabag" AND all:"directed information"` total=4 全在集内。
2. **自我降级**：§64-C 那句"文献中唯一发表过的显式部分观测反证界"两处错。`2606.31396`（2026-06-30，15 页）与其 5 页
   前身 `2601.12782`（2026-01-19）都在旧清单之外；作者已 API 核实（2601 = Ming Li, Fan Liu, Yifeng Xiong, Jie Xu,
   Tao Liu；2606 = 同五人 + Weijie Yuan, Shi Jin）⇒ **同组前后版，不许写成两篇独立工作**，必要性最早时间戳 2026-01-19。
   他们 **Th.1**：无噪系统均方可观测 ⇒ $I(z\to y)\ge R_{\exp}:=\log_2|\det A_u|=\sum_{|\lambda_i|\ge1}\log_2|\lambda_i|$，
   **对观测律的形式没有任何要求**（非线性、非高斯、秩亏都在前提内），Prop.1 扩到有限微分熵过程噪声。右端比较也错：
   $\ln|\det A|$ 在锚点是 **0.354774** bit，他们的是 **1.168539** bit，强 3.29 倍。
   `p0/note/note.tex` 已改并重渲染（exit=0、**6 页、Overfull 0、undefined 0**）：删 `no published bound returns the
budget-side asymptote`，换成引 `[17]`（新 `\bibitem{li2026sensing}`，条目 18→19）的限定句；`nearest published
converse` 改 `rate-side converse closest to our formulation`。省下的三行：Pacelli 记号声明 4→3、prefix-code 尾从句删、
   `li2026lower/atay2026rate` 改括号式。约定差异要说清：他们 $|\lambda|\ge1$、经典稳定化常数 $>1$，锚点无单位模
   特征值故同为 1.168539，正文用"the sum above (their statement sums over $|\lambda_i(A)|\ge1$)"而不混排两式。
   撞名警示：`2607.04172` 是 **S.** Li–Tanaka–Kim，`2606.31396` 是 **M.** Li 等，两条 bibitem 现相邻 [16]/[17]。
3. **率地板平面不变**（本节已被 §67-A 升级为**子空间无关的行列式不变量**）：
   $R_{\exp}=0.392553+0.775986=1.168539$，与 [4] 印出的无约束水平渐近线 1.1685 同一数（也与 §61 从 `2109.12246`
   §IV-A 读出的同）。测量噪声同送向 0：tail-4 / tail-3 / tail-2 在 $\sigma^2_v=10^{-8}$ 均 $\rho=0.76178$、
   $I=1.168539$、超出 $+0.000000$，同处代价 $642.5843/656.0803/914.6943$ ⇒ **受限传感的惩罚全部落在代价侧**。
   [P4]：他们 Cor.1 的 $\tfrac12\log_2\det(I_r+C\Sigma C^\top R^{-1})$ 与我车道的 $I$ 在 9 个可检测设计上差
   $\le7\times10^{-14}$（同一公式两个记号）。tail-1 看不见 $\lambda=-1.7124$（可见度 0.0000）但仍看见 $+1.3127$
   （0.7740），$\rho=1.71236=|\lambda_1|$ ⇒ 不可检测，那一行的"墙"由可检测性给出＝L-CSS Thm 2.1，不是新内容。
   **还剩的位置**：Th.1 的前提在那些行不满足，他们的定理对它们不作陈述，而我们要说的正是那些行。
4. **Q3 结清**：Chamon–Pappas–Ribeiro（*IEEE TAC* 65:1031–1046, 2020）机制：Def.1 的 $\alpha$ 是最小边际增益比（式 13），
   Th.1 给 $f(G_r)\le(1-e^{-\alpha r/s})^{-1}f(\mathcal X^\star)$，Th.5/6 的系统相关 $\hat\alpha_k$ 式 (35)
   $=\lambda_{\min}[P^{-1}_{m+k|m+k-1}]\big/\lambda_{\max}[P^{-1}_{m+k|m+k-1}+\sum_{u\in\mathcal O}V_u]$。
   在锚点（$V_u=I$、候选 $=$ 标准基输出、$\Pi$ 取滤波 DARE 先验信息矩阵）：全 4 根 $\hat\alpha=0.05990$ →
   证书比 **17.200×**（$r/s{=}1$）/ 33.893×（$0.5$）；去第 1 行 0.03762 → 27.083×；去 1、2 行 0.00335 → 298.837×；
   只留第 1 行 0.01334 → 75.436×。
   ⇒ **(1)** 能给出的最好证书比 $17.200\times$ vs C 实测贪心差距 $1.08$–$1.20\times$，差 **14.33 倍**：正文**不许**把
   Chamon 当"我们这套贪心有近似保证"的凭据引，只能当"常数因子证书在 LQG 锚点上松到不可用"的证据。
   **(2)** 锚点不在其非平凡区：$\sigma^2_v/\sigma^2_w=[3.374,0.315,0.192,0.144]$，他们 Fig.2 的硬实例在 $>10$ 一侧。
   **(3)** §64-E 的"要不要插一层松弛"文献已答：Remark 3 明写"[16] 证明除非 P=NP 否则拿不到系统无关的非平凡 $\alpha$
   界"并接"与 Th.5 一致，因为 $\hat\alpha\to0$ 当 $\sigma^2_v\to0$"；锚点忠实复现（全部同缩）
   $\hat\alpha=0.0599\to1.15\mathrm{e}{-}2\to1.39\mathrm{e}{-}3\to1.44\mathrm{e}{-}5\to\cdots\to0$，证书比
   $17.2\to87.4\to720\to6.94\mathrm{e}{-}4$ 的倒数即 $\infty$ ⇒ 系统无关的不可能性与系统相关的证书**同时成立**。
   引用陷阱：Chamon 的 [16] $=$ Ye–Roy–Sundaram ACC 2018, 5049–5054（`2003.11951`/TAC 2021 前身），而 `2003.11951`
   自己的 [16] $=$ Yang 等 TSP 63(9):2336–2348, 2015，**别把两套编号混引**。Remark 1 记住：他们自承 $\log\det$ 是 MSE
   的**坏代理**（$\partial\log\det/\partial\lambda_j=1/\lambda_j$ vs $1$，实测最多差 13%）——与 C 用 log-det 型目标
   做传感评分的口径直接相关。
5. `1811.11792` 通读 + **我的分类错先纠**：Nugroho–Taha–Gatsis–Summers–Krishnan（v2, 2019-03-07, 13 页）。§64-E 把它
   归进"向量高斯 RD 的 SDP 紧性"是错的，它属组合选择线：全文 `directed information` 命中 **0**、`duality` 0、
   `approximation ratio` 0、`submodul` 0；舞台是连续时间 $\dot x=Ax+Bu$、无噪声、无代价，目标是激活 SA 的个数
   $\sum(\pi_k+\gamma_k)$。承重：Th.1 把 SSASP **等价**写成 MI-SDP（靠 $\Psi_i\le L_i\Delta_i$ 线性化，Remark 3 自认
   "quality … very dependent on the choice of $L_{1,2,3}$"，算例 $10^4/5\!\times\!10^6/5\!\times\!10^6$）；Prop.2 只是
   SOF 存在的**充分**条件且可行性依赖实现；Th.2 对预生成组合库返回最优但无近似比，§8 给 $O(n_x^{6.6})$ 内点法口径。
   ⇒ 对我们有用的只有一条：它是"联合传感+执行器选择写成 MI-SDP 而无性能保证"的**第二个发表见证**。
   与普查互证：它**不在**精筛 35 条里，与 0 次命中一致——精筛与人工判读**方向一致**的正例。笔记 `papers/notes/1811.11792.md`。
6. **控制字符机制**（§67-E-6 用三方判别复核过，结论不变）：这个 shell 的 `python - <<'PY'` 会把 here-doc 里成对的
   `\` **折半**，随后非 raw 的 Python 字面量再把 `\b`/`\t`/`\v` 解析成控制字符——两步串联。单双引号无关（都中），
   raw 前缀能保住单反斜杠（实测 17 字符原样）。⇒ (i) 往板子写含反斜杠的内容一律走 Edit/Write；(ii) 必须走脚本则落盘后
   立刻按字节扫描（`find(bytes([8]))`）并与改前备份对数，修只按**单点字节索引**，**绝不**全局字符替换。
   §51 索引表那格里 C 的 3 个历史控制字符仍归 C 认领，本轮不动。
7. **给 C 的 6 条义务**：①任何"没有已发表的反证界/下界"式句子必须改成作用域陈述，可直接用 `a published converse
   (Li et al., arXiv:2606.31396, Th. 1) lower-bounds the directed information rate from the unstable sub-block by
   $\sum_{|\lambda_i|\ge1}\log_2|\lambda_i|$ for arbitrary, possibly rank-deficient observation laws; it is silent on
   the cost side and on sensing subspaces where the filter diverges, which is where this note lives`；bibitem 至少加
   `2606.31396`，引前身则注明前后版；②新落在对撞线上且旧清单没有的三条：`2004.02356`（控制与传感策略由一个 SDP
   联合合成——任何"我们把传感并进 SDP"式新意陈述都得对它收窄）、`2608.26917`（部分可观测 + side information 的最小率：
   条件 DI 下界 + 线性策略充分 + 标例凸）、`1912.03799`；③不许引 Chamon 当贪心的近似保证；④§6.1 的"插一层松弛"别再当
   开放问题写；⑤§1 row 5 的不可行区间只能声明**有限**那一段，无限那侧是可检测性/L-CSS Thm 2.1 的地盘，主张无限区间
   得先说明它不是 Thm 2.1 的重述；⑥§65-F-3 的复现义务至今未回（`run_all.sh` 的 `printf`/`tail -3` 不落 `$LOG`）。
