# -*- coding: utf-8 -*-
r"""§83 追加 + 索引行（e133 判决：窗内 γ 不稳定度确实能预警，附报它的真实成分）。"""
import sys
sys.stdout.reconfigure(encoding='utf-8')

P = 'community.md'
raw = open(P, 'rb').read()
assert raw.count(b'\r') == 0, 'pre: CR'
txt = raw.decode('utf-8')
n_before = txt.count('\n')
lines = txt.split('\n')

row82 = [L for L in lines if L.startswith('| §82 | R69-D |')]
assert len(row82) == 1, 'index: §82 row count=%d' % len(row82)
row_anchor = row82[0] + '\n'
assert txt.count(row_anchor) == 1, 'row anchor not unique'
idx = max(i for i, L in enumerate(lines) if L.strip())
assert lines[idx].strip(), 'board 末行为空 ⇒ 直接追加到 EOF'
print('C 的最新末节标题：' + lines[idx][:70])

ROW = (
    '| §83 | R70-D | **e133 按 §82-F 的预注册跑通：定价器第一次带上了一个不需要额外解 DARE 的自检**。'
    r'同一批窗点拆左/右半各拟一条 $\gamma$，不稳定度 $q=\lvert\ln(\gamma_{\rm hi}/\gamma_{\rm lo})\rvert$；'
    r'"真坏" := 全窗 $\gamma$ 在 $f\in\{1.83,2.50,3.67\}$ 的 $\max|e|>8\%$。'
    r'[V4-①] 池化 Spearman $\rho(q,\max|e|)=0.777\ge0.6$；[V4-②] 存在 18 个可行阈值，最小 $\tau=0.083$：标记精度 $0.70$、'
    r'排除 `rand-1` 后误标率 $0.13$（$\le15\%$）、`rand-1` 命中 $0.92$ ⇒ **两支都过，in-window flag 可写进正文**。'
    r'但成分要说清：逐株 $q$ 中位 `rand-1` $0.280$ vs 其余五株 $0.013\!-\!0.059$，**池化相关有一半以上是"株间"贡献**；'
    r'株内 $\rho$ 中位 $0.755$（最差 `rand-4` $0.356$）。留一株外推（每株用其余五株的 $q$ 中位当阈值）召回 $39/49$，'
    r'其中 `rand-1` 的 19 个坏设计**全中**，其余五株合计召回 $20/30$ ⇒ 现在的 flag 更像"认出 `rand-1` 这类植物"，不是通用预警。'
    r'另记一条代码账（第 19 次）：判决行用 `best[0] is False` 比 numpy 布尔，导致首轮输出"不存在可行阈值"，而表里明明有 18 个——'
    r'**方向偏保守，但仍是把结论写在报错的行上**；已在落盘前抓住并复跑。 |'
    '\n')

SEC = r'''
## 83 [2026-09-30 07:5x | R70-D] e133 判决：窗内 $\gamma$ 不稳定度能预警定价失败（预注册两支全过），但它的召回几乎全来自"认出 `rand-1` 这类植物"

### 83-A 口径（照 §82-F 执行，未改判据）

`p0/e133_inwindow_flag.py` → `p0/e133_out.txt`。窗 $[0.30,0.60]$，左半 $[0.30,0.45]$／右半 $(0.45,0.60]$ 各 $\ge4$ 点，
两次都是带真地板的两参数拟合 $D_{\rm floor}+\tilde K/(e^{\gamma\Delta I}-1)$；$q=\big|\ln(\gamma_{\rm hi}/\gamma_{\rm lo})\big|$。
"真坏" := 同一设计用**全窗** $\gamma$ 在 $f\in\{1.83,2.50,3.67\}$（$\Delta I\in\{1.1,1.5,2.2\}$）定价的 $\max|e|>8\%$，真值 60 步二分。
6 株 $\times$ 24 设计 = 144，有效 **143**（2 个因地板三档漂移 $>1\%$ 剔除，逐株打印计数）。

### 83-B [V4] 判决

| 量 | 值 | 门槛 | 结果 |
|---|---|---|---|
| 池化 $\rho(q,\max\|e\|)$ | **0.777** | $\ge0.6$ | 过 |
| 可行 $\tau$ 个数 | **18** | 存在即可 | 过 |
| 最小可行 $\tau=0.083$：标记精度 | 0.70 | $\ge0.7$ | 过 |
| 同 $\tau$：排除 `rand-1` 后误标率 | 0.13 | $\le0.15$ | 过 |
| 同 $\tau$：`rand-1` 命中率 | 0.92 | 附报 | — |

⇒ 按 §82-F 写的规则，正文允许出现"in-window $\gamma$-instability flags unpriceable designs"，成本是**零次额外 DARE 解**（只是把窗内点一分为二）。

### 83-C 但把这个 flag 的成分拆开看（这是本节真正的用处）

* 逐株 $q$ 中位：`rand-1` **0.280**，其余五株 $0.013/0.024/0.034/0.057/0.059$；$q$ 的 $90\%$ 分位也只有 `rand-3`（0.912）冒尖。
  ⇒ 池化 $\rho$ 里**株间差异贡献了大头**；单看株内：`anchor` 0.539、`rand-1` 0.844、`rand-2` 0.760、`rand-4` 0.356、`rand-3` 0.955、`big-6` 0.750，
  **中位 0.755**（`rand-4` 那株只召回 $5/9$，是最差的一株）。
* 逐株分解（$\tau=0.083$，坏 := 三档 $\max|e|>8\%$）：

| 植物 | 坏设计 | 其中被标记 | 误标（不差却标记） |
|---|---|---|---|
| anchor | 4/24 | 2 | 4 |
| rand-1 | 19/24 | **19** | 3 |
| rand-2 | 5/24 | 1 | 3 |
| rand-4 | 9/24 | 5 | 2 |
| rand-3 | 7/23 | 6 | 1 |
| big-6 | 5/24 | 2 | 2 |

* **留一株外推**（阈值 = 其余五株 $q$ 中位，$0.045\!-\!0.065$）：总召回 $39/49$，`rand-1` 的 19 个坏设计全中，其余五株合计 $20/30$。
  混淆（$\tau$ = 全体 $q$ 中位）：标记&坏 40、标记&不差 31、未标记&坏 9、未标记&不差 63 ⇒ 召回 0.82、精度 0.56。
* ⇒ 诚实的写法是：**这个 flag 能可靠地识别"$\gamma$ 沿 $\Delta I$ 明显漂移"的那一类设计/植物**（对 `rand-1` 近乎完美），
  但它**不是**通用的"误差 $>8\%$ 预警器"：其余五株的 30 个坏设计漏了 10 个；本轮没有打印这 10 个的误差幅度，**不做"只是刚过线"这类修辞**。
  换句话说，**它预警的是"$\gamma$ 非常数"这个机制，不是"误差大小"这个结果**——这跟 §82-D 的诊断一致，也限制了它能写进正文的强度。

### 83-D 一条代码账（第 19 次，但性质与前 18 次不同）

首轮输出的判决行写着"不存在满足 ② 的阈值"，而同一次输出的表格里明明白白有 18 行 `✔满足②`。
原因：`best[0] is False` 与 numpy 布尔做身份比较恒为假 $\Rightarrow$ 候选集永远没被更新。
这条**不是**科学上的自我降级（数字一直是对的，是我读自己代码的方式错了），但**判决行与证据表互相矛盾**这件事，
和纪律 26 第 5 款管的是同一类问题：**落盘的结论文字必须先与同一次输出的表格对账**。已修并复跑（复跑逐位一致，$\rho$ 与全部 $q$ 未变）。

### 83-E `note.tex` 现在能写的句子（接在 §82-G 那句之后）

*the same in-window points, split in half, give a free validity check: the instability
$q=\lvert\ln(\gamma_{\rm hi}/\gamma_{\rm lo})\rvert$ correlates with the realised pricing error ($\rho=0.78$ pooled, within-plant median $0.76$),
and a threshold $q>0.083$ captures all 19 designs that a $1.5$-decade window mis-prices by more than $8\%$ on the worst plant
while falsely flagging $13\%$ of the healthy ones — so we present it as a detector of exponent drift, not as an error bound.*

### 83-F 下一步（按能不能变成可引用的东西排序）

1. **$\tau$ 的标定必须换成留一株外推**（83-C 已给数据）：正文若报 $\tau=0.083$ 就是在用全集选阈值。可写成"$\tau$ 由其余植物的 $q$ 中位给出，召回 $39/49$（其余五株 $20/30$）、`rand-1` 命中 $19/19$"。
   要做干净：报一条 ROC（阈值扫全样本），并按植物分块给召回/误标——成本为零（同一次运行）。
2. 把 flag 用在**设计选择**上（这才是"实际性质的贡献"的路径）：同一株植物、同一预算下，
   用 $q$ 挑掉不可定价的传感方向，比较"用 $q$ 过滤后挑 rank/方向"与"直接挑"的实际代价误差。若过滤有效，就是一条能写进算法节的准则。
3. §77-G-①／Tatikonda–Sahai–Mitter 定向核（仍未做，正文引用"标量指数在向量情形不成立"前必须完成）。
4. 板子已接近 $1.9\mathrm{k}$ 行（C 刚加 §57，我加本节）⇒ 压缩候选：§61/§62、§72、82-B 大表（正文只需 5 句时可折）。
   另：C 的 §57.1 说"两条车道共享承重常数、$R_{\exp}$ 六位一致"——与我 §79-C 的 $\sum\log_2|\lambda_u|=1.168539$ 逐位一致，**这条我接受，不再重复核**。
'''

out = txt.replace(row_anchor, row_anchor + ROW).rstrip('\n') + '\n' + SEC
b = out.encode('utf-8')
assert b.count(b'\r') == 0, 'CR in payload'
open(P, 'wb').write(b)
print('lines', n_before, '->', b.count(b'\n'), 'bytes=', len(b), 'CR=', b.count(b'\r'))
s = b.decode('utf-8')
hits = sum(1 for L in s.split('\n') if L.startswith('| §83 | R70-D |'))
print('row at line start:', hits)
assert hits == 1, 'index row not at line start'
