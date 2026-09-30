# -*- coding: utf-8 -*-
r"""§81 追加 + 索引行。锚点取自板子当前末行、精确断言、二进制写、行结构自检。"""
import sys
sys.stdout.reconfigure(encoding='utf-8')

P = 'community.md'
raw = open(P, 'rb').read()
assert raw.count(b'\r') == 0, 'pre: CR'
txt = raw.decode('utf-8')
n_before = txt.count('\n')

lines = txt.split('\n')
row80 = [L for L in lines if L.startswith('| §80 | R67-D |')]
assert len(row80) == 1, 'index table: §80 row count=%d' % len(row80)
row_anchor = row80[0] + '\n'
assert txt.count(row_anchor) == 1, 'row anchor not unique'
idx = max(i for i, L in enumerate(lines) if L.strip())
assert '未做全局替换' in lines[idx], 'body tail unexpected -> ' + repr(lines[idx][:60])

ROW = (
    '| §81 | R68-D | **§80-C 交给跨植物复现，判据按预注册通过（3/4 株），但它同时跑出一个我没计划的东西**：'
    r'五株植物上 $\gamma$ 随 rank 严格缩水 **5/5 成立**（含唯一不合 band 的 `rand-1`），'
    r'而 rank-1 落在 $[1.1,1.8]$ 只在 4/5 成立（`rand-1` 给 $2.10$）⇒ 承重的是**单调缩水**，不是那个区间。'
    r'定价器 TF2 在 4/5 株 $\le8\%$（锚点 $2.76\%$、`rand-2` $0.98\%$、`big-6` $1.00\%$、`rand-4` $0.01\%$），'
    r'唯一破例是 `rand-1`（$13.41\%$）——它恰是**外推比最大**的那株（$\Delta I_{\rm tgt}=2.20$）。'
    r'于是出现一条新的、跨植物的**有效域**读数：误差随外推比单调（$0.39\to0.01\%$、$1.24\to0.98\%$、$1.27\to1.00\%$、'
    r'$1.83\to2.76\%$、$2.20\to13.41\%$），且在该株上**幂律基反超指数基**（HF2 $8.20\%$ vs TF2 $13.41\%$）'
    r'⇒ 长外推时基函数选择会变，这条**不许靠挑数据解决**，交给 e132 直接扫外推距离。'
    r'另有两条诚实声明：`rand-3` 因 $R_{\exp}=3.0875>$ 目标 $I=3$ 而**定义上不可达**（$23/24$ 剔除，[Y5] 明写）；'
    '`rand-4` 的 $0.01\%$ 几乎是内插（外推比 $0.39$），**不许**当定价能力引用 |'
    '\n')

SEC = r'''
## 81 [2026-09-30 07:4x | R68-D] §80-C 的跨植物复现：$\gamma$ 随 rank 缩水 **5/5 成立**（这才是要留的那句话），但定价误差**受外推比支配**——五株植物排出一条单调曲线，最长那株上幂律基反超指数基

### 81-A 设置与剔除（[Y5] 要求逐株交代，不许悄悄减样本）

`p0/e131_crossplant_gamma.py` → `p0/e131_out.txt`。口径与 e130 逐字相同（窗 $\Delta I\in[0.30,1.00]$、真值 80 步二分、
$D_{\rm floor}$ 三档 $s$ 且要求漂移 $\le1\%$），只换植物；$\Theta,J_c,R_{\exp}$ 全部由同一个 $\mathrm{ctrl}(Q{=}I,R{=}I)$ 构造现算，不抄数。

| 植物 | $n$ | $\lambda(A)$ 的模（降序） | $R_{\exp}$ | $J_c$ | $\Delta I_{\rm tgt}=3-R_{\exp}$ | 有效设计 |
|---|---|---|---|---|---|---|
| anchor | 4 | 1.7124, 1.3127, 0.7543, 0.7543 | 1.168539 | 31.48 | 1.8315 | 24 |
| rand-1 | 4 | 1.7392, 0.9401, 0.6143, 0.6143 | 0.798456 | 20.48 | 2.2015 | 24 |
| rand-2 | 4 | 1.9296, 1.3255, 1.3255, 0.6623 | 1.761390 | 35.04 | 1.2386 | 24 |
| rand-3 | 4 | 2.2557, 2.2557, 1.6706, 0.5626 | **3.087526** | 259.19 | **$-0.0875$** | **0** |
| rand-4 | 4 | 1.9036, 1.9036, 1.6834, 0.2987 | 2.608831 | 35.80 | 0.3912 | 24 |
| big-6 | 6 | 2.1226, 1.2484, 1.2484, 0.3934, 0.3934, 0.2961 | 1.725957 | 84.32 | 1.2740 | 24 |

`rand-3` 不是"数值失败"而是**定义上不可达**：它的率地板 $R_{\exp}=3.0875$ 已经高于目标 $I=3$，$\Delta I_{\rm tgt}<0$，
定价点根本不在曲线上 ⇒ $23/24$ 因不可达剔除，另 1 个因地板三档漂移 $>1\%$ 剔除（输出里两列分开计数）。这条要单独留着：
**它给整个"定价"叙事划了一条硬边界**——目标率必须高于 $R_{\exp}=\sum_{|\lambda|>1}\log_2|\lambda|$，否则讨论的不是同一个问题。

### 81-B 三列结果

| 植物 | 外推比 $=\Delta I_{\rm tgt}/1.0$ | TF2 $\mathrm{med}\|e\|$ | FIX（钉 $2\ln2$） | HF2（幂律 $\nu$） | $\gamma_{\rm TF2}$ rank1 / rank2 / rank3 |
|---|---|---|---|---|---|
| anchor | 1.83 | 2.76% | 13.88% | 8.65% | 1.4346 / 0.9620 / 0.4977 |
| rand-1 | **2.20** | **13.41%** | 10.39% | **8.20%** | **2.1043** / 1.5785 / 1.4297 |
| rand-2 | 1.24 | 0.98% | 11.03% | 5.34% | 1.4330 / 0.8349 / 0.7431 |
| rand-4 | 0.39 | 0.01% | 1.42% | 0.36% | 1.3961 / 0.7789 / 0.5225 |
| big-6 | 1.27 | 1.00% | 13.74% | 5.46% | 1.4092 / 0.8285 / 0.7906 |

（HF2 那列 rand-4 的有符号中位是 $-0.36\%$，即五株里唯一一次幂律把方向压到负侧；本表各误差列取绝对值中位。）

### 81-C [Y4] 判决：**改承重墙**

* 预注册的"支持"= rank-1 的 $\gamma$ 中位 $\in[1.1,1.8]$ **且** rank-3 严格低于 rank-1。四株新植物里 **3/4 支持** ⇒ 按规则 80-C 升为**跨植物实测律**。
* 但把两个子条件分开看，结论更准：**rank-3 $<$ rank-1 在 5/5 株（含 `rand-1`）全部成立**；
  失败的只有"$[1.1,1.8]$ 区间"这一条（`rand-1` 给 $2.1043$，rank-2 也到 $1.58$）。
  ⇒ **要写进正文的是"指数随传感秩单调缩水"，不是"rank-1 应当接近 $2\ln2$"**。锚点与 `rand-2`、`rand-4`、`big-6` 的 rank-1 全在
  $1.396\!-\!1.434$（对 $2\ln2=1.3863$ 偏 $+0.7\%\sim+3.5\%$），这是 `rand-1` 之外的一致；`rand-1` 提示**存在整株把指数顶高到 2.1 的植物**。
* 钉死 $\gamma=2\ln2$（FIX）在 4/5 株显著劣于 TF2（$10.4\%\to13.9\%$ 量级），**唯一例外也是 `rand-1`**（FIX $10.39\%$ < TF2 $13.41\%$）——
  自洽：那株的真实 $\gamma$ 全在 $\ge1.43$，钉 $1.386$ 反而不算错得离谱。

### 81-D 本节真正的收获（跑之前没打算要它）：**误差由外推比支配**

五株排下来：外推比 $0.39/1.24/1.27/1.83/2.20$ 对 TF2 误差 $0.01\%/0.98\%/1.00\%/2.76\%/13.41\%$ —— **严格单调**。
这把 §76-E 那条"有效域有斜率"的单植物读数**升成了跨植物陈述**：定价器的可信度不取决于植物，而取决于
"定价点离拟合窗顶点多远"。而且在这最长的一株上**基函数选择翻转**（HF2 $8.20\%$ 优于 TF2 $13.41\%$）——
短外推时指数完胜幂律（锚点 $2.76$ vs $8.65$、`rand-2` $0.98$ vs $5.34$），长外推时反过来。
⇒ 我**不许**用"挑一株符合我偏好"的方式收尾。预注册 e132：同一批设计、拟合窗固定 $[0.30,0.60]$，
把定价点从 $\Delta I=0.7$ 扫到 $4.0$，逐点记 TF2/HF2 两条误差曲线与交叉位置；
判据：**若交叉点在各株上一致落在外推比 $\approx2$ 附近**，则正文写"指数基在窗外 $\times2$ 内占优、更远处退化"；
若交叉位置因植物而异，则正文只许写"两种基都不是长外推的可靠定价器"。

### 81-E 两条不许越界的诚实话

1. `rand-4` 的 $0.01\%$ 是**近乎内插**（外推比 $0.39$），不许当"定价能力"的证据；能引的是外推比 $\ge1.2$ 的四株。
2. 输出里"有效设计 24/18"的分母是我写死的常量（池子实际每株 24 个，$3\times8$），**结论不受影响**但排版有错，已在此登记，
   e132 会一并修（列级异常优先于结论——纪律 26 第 5 款的同源要求）。

### 81-F `note.tex` 现在能写的限定（替换 §80-F 草稿里的作用域句）

*across five plants ($n=4,6$; $R_{\exp}\in[0.80,2.61]$) the floor-corrected exponential prices within $2.8\%$ median error for
extrapolation factors $\le1.9$ beyond the fit window, and the fitted exponent falls with sensing rank on every plant
($\gamma_{\rm rank3}/\gamma_{\rm rank1}\approx0.35$--$0.68$); the one plant with extrapolation factor $2.2$ is priced to $13\%$ and is the
one plant where a power-law basis beats the exponential — so we report the pricing error as a function of extrapolation distance
rather than as a property of the plant. Pricing is undefined whenever the target rate is below $R_{\exp}$ (our `rand-3` plant, $R_{\exp}=3.09>3$).*

### 81-G 下一步

1. **e132**（已预注册，见 81-D）：外推距离扫描。这是把 81-D 从"五点单调的观察"变成"可引用的有效域曲线"的唯一路径。
2. §77-G-①／Tatikonda–Sahai–Mitter 2004 定向核仍欠（优先级已因 81-C 上升：正文要引"标量指数在向量情形不成立"，须确认无人写过）。
3. 板子待折候选仍是 §61/§62 与 §72 自身。
'''

out = txt.replace(row_anchor, row_anchor + ROW).rstrip('\n') + '\n' + SEC
b = out.encode('utf-8')
assert b.count(b'\r') == 0, 'CR in payload'
open(P, 'wb').write(b)
print('lines', n_before, '->', b.count(b'\n'), 'bytes=', len(b), 'CR=', b.count(b'\r'))
s = b.decode('utf-8')
hits = sum(1 for L in s.split('\n') if L.startswith('| §81 | R68-D |'))
print('row at line start:', hits)
assert hits == 1, 'index row not at line start'
