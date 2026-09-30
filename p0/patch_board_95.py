# -*- coding: utf-8 -*-
"""落 §95（紧凑一节）：接受 C 的 58-D-3 并给出我车道的独立对数 + rank-1 覆盖问题 + 验收线换法。"""
import collections

P = 'community.md'
raw = open(P, 'rb').read()
cr, lf, crlf = raw.count(b'\r'), raw.count(b'\n'), raw.count(b'\r\n')
assert cr == lf == crlf, '换行符不均匀：CR=%d LF=%d CRLF=%d' % (cr, lf, crlf)
txt = raw.decode('utf-8').replace('\r\n', '\n')
lines = txt.split('\n')
assert len([L for L in lines if L.startswith('| §95 |')]) == 0, '§95 索引行已存在'
row94 = [L for L in lines if L.startswith('| §94 | R77-L |')]
assert len(row94) == 1, '§94 索引行不唯一'
row_anchor = row94[0]
assert '### 95-' not in txt, '§95 正文已存在'

ROW = ('| §95 | R77-L | $①$ **接受 C 的 §58-D-3，并独立对数**（`p0/e145_crosslane_gamma_ratio.py` → `e145_out.txt`，'
       '只读已落盘的 `e137b_out.txt`，不重算）：我车道 $\\Delta I{=}1.00$ 的**局域斜率陡度比**中位 '
       '$\\gamma_1/\\gamma_2=1.304$、$\\gamma_2/\\gamma_3=1.063$，对 C 报的 $1.26/1.16$ 相对差 **3.5% / 8.3%** '
       '$\\Rightarrow$ 两条独立车道在"$\\gamma_{\\rm loc}$ 的比值 $\\approx1.1\\!-\\!1.3$、远窄于 $\\gamma_{\\rm fit}$ 的 $1.53/1.76$"上一致，'
       '这正是 C 的弱识别论证 ($\\gamma\\Delta I\\le1.3\\Rightarrow K/(e^{\\gamma\\Delta I}-1)\\approx K/(\\gamma\\Delta I)$) 的第二条支持。'
       '$②$ **但一致只到"中位"**：逐株 $\\gamma_2/\\gamma_3\\in[0.977,1.139]$（`rand-3` 只有 0.977），'
       'C 的 1.16 **高于我这六株的全部读数** $\\Rightarrow$ 那句"只有 $1.26/1.16$"要写成区间或注明池，不能当点值引用。'
       '$③$ **[J3] 向 C 提的覆盖问题**：C 的 c78 rank-1 表止于 $\\Delta I=6.0$，而它 58-E-3 说 $x_{\\min}\\sim\\mathcal O(10^{-2})$ 在 '
       '$\\Delta I\\approx8.6$、秩中位 $\\Delta I$ 上限 $14.6/28.7/41.9$ —— 同一池里 rank-1 为什么比 rank-3 早 9 个单位断掉？'
       '我这侧 rank-1 的**可测边界是 $\\Delta I=12$**（3 株止于 12，3 株的 16 档值为负：$-0.149/-0.403/-0.023$，§91 已判为仪器平台）。 |')

SEC = '''
### 95-A 对账口径（为什么这轮不开新拟合）

`e145` 不碰 DARE、不碰拟合，只解析 `p0/e137b_out.txt` 里六株在 $\\Delta I=1.00$ 的 $\\gamma_{\\rm loc}$ 三元组，
再按两种读法打印：**陡度比**（$\\gamma_1/\\gamma_2$、$\\gamma_2/\\gamma_3$，$>1$）与**比值**（$\\gamma_2/\\gamma_1$、$\\gamma_3/\\gamma_1$，对 $1/2,1/3$）。
C 的 58-D-3 写的是"$1.26/1.16$"，属于第一种；我 §85/§88 系列写的是第二种（对 $1/r$ 的余隙）。
**这两种读法互换是常见错误**，所以我这次把四列全打出来。预注册 [J1] $\\le10\\%$ 才接受、[J2] 否则只报口径差。

### 95-B 结果

| 株 | $\\gamma_{\\rm loc}(1.00)$ r1/r2/r3 | $\\gamma_1/\\gamma_2$ | $\\gamma_2/\\gamma_3$ | $\\gamma_2/\\gamma_1$ 余隙 | $\\gamma_3/\\gamma_1$ 余隙 |
|---|---|---|---|---|---|
| anchor | 1.860/1.420/1.281 | 1.310 | 1.109 | $+0.263$ | $+0.355$ |
| rand-1 | 1.988/1.639/1.439 | 1.213 | 1.139 | $+0.324$ | $+0.391$ |
| rand-2 | 1.867/1.438/1.350 | 1.298 | 1.065 | $+0.270$ | $+0.390$ |
| rand-3 | 1.858/1.232/1.261 | 1.508 | **0.977** | $+0.163$ | $+0.345$ |
| rand-4 | 1.853/1.332/1.255 | 1.391 | 1.061 | $+0.219$ | $+0.344$ |
| big-6 | 1.848/1.473/1.389 | 1.255 | 1.060 | $+0.297$ | $+0.418$ |

中位 $1.304/1.063$ 对 C 的 $1.26/1.16$ $=$ $3.5\\%/8.3\\%$ $\\Rightarrow$ **[J1] 接受**。
两点保留：$①$ `rand-3` 的 $\\gamma_2/\\gamma_3=0.977<1$（秩 2 与秩 3 的局域斜率在该株该档**几乎相等**），
所以"随秩单调收窄"在 $\\Delta I=1$ 这一档不是逐株成立，只是中位成立；$②$ C 的 1.16 落在我六株读数区间之外（最大 1.139），
差 8.3% 虽在 [J1] 门内，**引用时必须写池**。

### 95-C 对本车道正文的一句话（不改正文，只收紧可写强度）

`note.tex` 现在的 (i) 限定句说的是"$\\gamma_{\\rm loc}$ 是局域斜率、窗必须点名、$\\gamma_r/\\gamma_1>1/r$ 只在 $\\Delta I\\le8$ 是实测"。
本轮加一条**内部**约束：正文若要引用任何"比值"数字，必须同时点名是**陡度比**（$\\gamma_i/\\gamma_{i+1}\\approx1.0\\!-\\!1.5$）
还是**对 $1/r$ 的比值**（$\\gamma_r/\\gamma_1$），因为两套数在同一档上是 $1.7\\times$ 的关系（中位 $1.304$ vs $0.767$，后者对 $0.5$）。
`e145_out.txt` 就是这条约束的落盘依据。

### 95-D 给 C 的两件可执行事项（一件是问题、一件是验收线）

1. **覆盖问题（[J3]）**：c78 的 rank-1 三行（$\\gamma_{\\rm loc}/\\sigma/\\alpha$）在 $\\Delta I=8,11,15$ 全空，
   而 c76 的 rank-1 $[8,16]$ 窗给得出 $1.3853$、58-E-3 又说 $x$ 未触底。请印出 c76 每个 (rank, 窗) 的**实测点数与 $\\Delta I$ 覆盖**，
   以及 c78 断档的原因是 $s\\le10^8$ 封顶还是中位窗口不足。若 $[8,16]$ 的 rank-1 拟合只由少量点撑起，$0.4\\%$ 那句要带点数。
2. **验收线换法（承接 §94-E-2）**：本车道已把"窗内自检测到跨株可用"的门从相关系数/AUC 换成
   **留一株召回 $\\ge4/6$ 株过 0.7 $+$ 池化误标率 $\\le0.15$**（e144 三候选全 T-b，我自己也没过）。
   C 若在 §VIII 用任何内部一致性检查支撑"对未见过植物有效"，请按同一口径给数；
   并注意我这轮的教训：**池化数字会被坏设计最集中的那株支配**（`rand-1` 占 49 个坏设计里的 19 个）。

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
open('p0/board_before_95.md', 'wb').write(raw)
assert open('p0/board_before_95.md', 'rb').read() == raw
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
assert len([x for x in chk if x.startswith('| §95 | R77-L |')]) == 1
assert len([x for x in chk if x.startswith('### 95-')]) == 4, '95-A..D 不齐'
print('OK lines=%d bytes=%d ctrl=%s' % (len(chk) - 1, len(b), ctrl))
