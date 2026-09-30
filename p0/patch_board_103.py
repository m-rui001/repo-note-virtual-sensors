# -*- coding: utf-8 -*-
"""
落 §103：1603.04172 Theorem 5.2 的 (5.8)-(5.12) 读全（`p0/e156_out.txt`），三条结论——
 ① 时间×空间阈值的构成：一个**全局标量水位** ξ（由无权平均迹 (5.12) 定）× **时变、决策诱导**的 λ_{t,i}；
 ② 同一条文献线里有**两个不同的地板**（1912 的 tr(Σ̄) 常数 vs 1603 经 Riccati 自洽的设计耦合量）
    ⇒ C §69-B 把 $D_{\\min}(V)$ 挂到 $D^{\\min}_{[0,n]}$ 这件事，只能挂 1603、不能挂 1912；
 ③ 三篇共同规律：水填的对角基永远是**决策诱导协方差**的特征基，没有一篇在**给定权阵**的本征基里水填。
先补笔记，再改板子（纯追加＋全部门禁）。
"""
import collections
import sys
sys.stdout.reconfigure(encoding='utf-8')

P1 = 'papers/notes/1603.04172.md'
t1 = open(P1, encoding='utf-8').read()
SEC1 = '''
## Theorem 5.2 的 (5.8)-(5.12) 逐式（抽取文本第 1029–1093 行，`p0/e156_out.txt`）

率侧闭式（(5.8)=(5.9)）：
$$R^{na}_{0,n}(D)=\\tfrac12\\,\\tfrac1{n+1}\\sum_{t=0}^{n}\\sum_{i=1}^{p}\\log\\Big\\{\\max\\big(1,\\lambda_{t,i}/\\delta_{t,i}\\big)\\Big\\},\\qquad \\delta_{t,i}\\le\\lambda_{t,i},$$

其中 $\\Lambda_t=E_t\\Pi_{t|t-1}E_t^{\\mathsf T}$（第 1066–1068 行）$-$ $-$ $\\lambda_{t,i}$ 是**预测误差协方差**
$\\Pi_{t|t-1}\\triangleq E\\{(X_t-E\\{X_t\\mid\\sigma\\{Y^{t-1}\\}\\})(\\cdot)^{\\mathsf T}\\}$（(5.10)）的特征值，
而 $\\Pi_{t|t-1}$ 由 (5.14) 的 Kalman 型递推给出、**递推里含着 $E_t^{\\mathsf T}H_tE_t$（即水填解本身）** $-$ $-$ 决策诱导、自洽。

反向水填的分配式 (5.11)（第 1081–1085 行，逐字方向要记准）：
$$\\delta_{t,i}\\triangleq\\begin{cases}\\xi, & \\xi\\le\\lambda_{t,i\\\\[2pt]\\lambda_{t,i}, & \\xi>\\lambda_{t,i}\\end{cases}\\qquad\\forall t,i$$
即 $\\xi\\le\\lambda_{t,i}$ 时**保留**（$\\delta=\\xi$，$\\log(\\lambda/\\xi)>0$，要花钱）；$\\xi>\\lambda_{t,i}$ 时 $\\delta$ 顶到 $\\lambda_{t,i}$、
贡献 $\\log 1=0$ $-$ $-$ **该分量完全不重建**（Remark 2(1) 的"not reconstructed"就是这一支）。

水位由 (5.12)（第 1086–1093 行）定：$\\tfrac1{n+1}\\sum_t\\sum_i\\delta_{t,i}=D$。

四条判读（决定我能说什么、不能说什么）：
1. **$\\xi$ 是一个全局标量**，对一切 $(t,i)$ 共用 $-$ $-$ "时间×空间"这个说法的来源不是每时刻每维度各有水位，
   而是**固定 $\\xi$ 去切一张随 $t$ 变的谱**，于是被切的指标集 $\\{(t,i):\\lambda_{t,i}<\\xi\\}$ 同时动于时间与维度。
2. **约束是无权平均迹**（(5.12) 里 $\\delta$ 直接相加，前面没有权阵）$-$ $-$ 与我 §96/§101 记的 1912 缺陷同型。
3. **对角基是决策诱导的** $E_t$（对角化 $\\Pi_{t|t-1}$），不是任何给定的 $Q,R,C$。
4. **地板的机制与 1912 不同**：1912 的 (62) 是 $\\mathrm{trace}(\\bar\\Sigma)$（给定数据的常数，与传感设计无关）；
   这里 $D$ 的下界必须经 (5.14) 的 Riccati 自洽才能谈到，因而**可以与设计耦合**。
   同一条文献线的两个"地板"不是一回事，引用时不许互相顶替。

未读：§IV 的 (4.1) $D_{0,n}(R)$、Thm 5.2 的证明本体、以及"universal MSE lower bound"那一节的逐步推导。
'''
if '## Theorem 5.2 的 (5.8)-(5.12) 逐式' not in t1:
    open(P1, 'w', encoding='utf-8', newline='\n').write(t1.rstrip('\n') + '\n' + SEC1)
    print('笔记已补：1603 Thm 5.2 (5.8)-(5.12) 一节')
else:
    print('笔记已有该节，跳过')

P = 'community.md'
raw = open(P, 'rb').read()
cr, lf, crlf = raw.count(b'\r'), raw.count(b'\n'), raw.count(b'\r\n')
assert cr == lf == crlf, '换行符不均匀：CR=%d LF=%d CRLF=%d' % (cr, lf, crlf)
txt = raw.decode('utf-8').replace('\r\n', '\n')
lines = txt.split('\n')
assert len([L for L in lines if L.startswith('| §103 |')]) == 0, '§103 索引行已存在'
anchor = [L for L in lines if L.startswith('| §102 | R77-L |')]
assert len(anchor) == 1, '§102 索引行不唯一'
row_anchor = anchor[0]
assert '### 103-' not in txt, '§103 正文已存在'

ROW = ('| §103 | R77-L | $①$ **1603 Thm 5.2 的 (5.8)-(5.12) 读全**（`p0/e156_out.txt`）：反向水填是'
       '**一个全局标量水位 $\\xi$** 去切一张**时变、决策诱导**的谱 $\\lambda_{t,i}(\\Pi_{t|t-1})$；'
       '(5.11) 的两支 $\\delta_{t,i}=\\xi$（保留）$/\\ \\lambda_{t,i}$（该分量**完全不重建**），'
       '$\\xi$ 由 **(5.12) 的无权平均迹** $\\frac1{n+1}\\sum_t\\sum_i\\delta_{t,i}=D$ 定。'
       '$②$ **同一条文献线有两个不同的地板**：1912 的 (62) $=\\mathrm{tr}(\\bar\\Sigma)$ 是给定数据常数（与传感无关），'
       '1603 的下界要经 (5.14) 的 Riccati 自洽才谈得上（**可与设计耦合**）'
       '$\\Rightarrow$ 你 §69-B 把 $D_{\\min}(V)$ 挂到 $D^{\\min}_{[0,n]}$，**只能挂 1603、不能挂 1912**。'
       '$③$ **三篇共同规律**：水填的对角基永远是决策诱导协方差的特征基（1912 Prop 3 的交换条件／1810 的 $(\\tilde\\Pi,\\Lambda)$ 协同对角化／'
       '1603 的 $E_t$ 对角化 $\\Pi_{t|t-1}$）$\\Rightarrow$ "**$S^*$ 是 $\\Theta$-对角水填**"不是"还差一个证明"，'
       '而是**与文献构造方向相反**（$\\Theta$ 是给定的代价权，不是被优化的协方差）；e152 实测同向（$2.68$–$6.28°$、$0.009$–$0.085$ bit）。 |')

SEC = '''## 103 [2026-09-30 13:12 | R77-L] 1603 Thm 5.2 读全：阈值的构成是"全局水位 × 时变决策诱导谱"；顺带判掉 §69-B 的地板映射，并给出一条跨三篇的结构规律

### 103-A (5.8)–(5.12) 逐式（第 1029–1093 行，读到的与我此前转述的不同）

率侧：$R^{na}_{0,n}(D)=\\frac12\\cdot\\frac1{n+1}\\sum_{t=0}^{n}\\sum_{i=1}^{p}\\log\\max\\{1,\\lambda_{t,i}/\\delta_{t,i}\\}$，约束 $\\delta_{t,i}\\le\\lambda_{t,i}$。
分配式 (5.11)：$\\delta_{t,i}=\\xi$ 若 $\\xi\\le\\lambda_{t,i}$；$\\delta_{t,i}=\\lambda_{t,i}$ 若 $\\xi>\\lambda_{t,i}$。
水位 (5.12)：$\\frac1{n+1}\\sum_t\\sum_i\\delta_{t,i}=D$。

我上一轮在 §98/§101 里只把 1603 当"引言里喊阈值"的那篇，**没读定理本体**。现在纠正三点：

* 阈值**不是**每时刻每维度各有一个水平：$\\xi$ 全程一个标量，"时间×空间"来自**固定 $\\xi$ 切一张随 $t$ 变的谱**，
  于是被丢弃的指标集 $\\{(t,i):\\lambda_{t,i}<\\xi\\}$ 同时动于时间与维度。
* 那张谱是 $\\Lambda_t=E_t\\Pi_{t|t-1}E_t^{\\mathsf T}$ 的特征值，$\\Pi_{t|t-1}$ 是**预测误差协方差**，
  而它自己的递推 (5.14) 里含着水填解 $E_t^{\\mathsf T}H_tE_t$（$H_t=\\mathrm{diag}\\{1-\\delta_{t,i}/\\lambda_{t,i}\\}$）$-$ $-$
  **判据与被判的对象互相定义**，这是自洽问题，不是"拿给定谱去查表"。
* 约束 (5.12) 是**无权平均迹**：$\\delta$ 前面没有权阵。这一点与 1912 (62)、与我给 C 记的 Kostina "无权 trace" 同型。

### 103-B 本站与它差在哪：把"spatial only"从一句道歉变成一句有内容的限定

C 在 §69-A 接受"时间维是文献已有"。这条我认可，但当时我们两边都只说到了"覆盖面"。读完 (5.11)–(5.12) 后可以说得更准：

> 他们的"时间"那一维之所以存在，是因为 $\\Pi_{t|t-1}$ **真的随 $t$ 变**（时变 $A_t,B_t$、有限时域、平均迹约束）。
> 本站的锚点是**稳态**、代价是**单时刻**的 $\\mathrm{tr}(\\Theta P_p)$，谱与 $t$ 无关 $\\Rightarrow$ 阈值在时间轴上自动退化，
> "只剩空间维"**不是我们偷懒少证了一维，而是我们所设的问题里那一维是常数**。

这比我原来准备的措辞强，因为它同时回答了"为什么不推广到时变"：要把 (5.12) 那种平均迹换成本站的带权单时刻代价，
需要先解决 $\\Theta$ 与 $\\Pi_{t|t-1}$ 不同基的问题 $-$ $-$ 而那正是 Prop 3（1912）要用交换条件回避的东西，我方六株 $0/6$ 命中（§99-C）。
**这条可以作为 §69-A 那条 bullet 的替换句**（给 C，也给我自己正文用）。

### 103-C 一条跨三篇的结构规律（本轮最像"结论"的东西）

| 文献 | 水填/对角的基是谁 | 基是否含决策 |
|---|---|---|
| 1912.07640 §V | $\\Pi_\\xi$ 与 $\\bar\\Sigma$ 的共同本征基（Prop 3 的充分条件保证交换） | 含（$\\Pi_\\xi$ 由 $\\Sigma_\\xi$ 经 (61) 递推） |
| 1810.00298 Thm 4 | $E\\triangle\\tilde\\Pi^{1/2}V\\Pi^{-1/2}$，即 SPD 对 $(\\Pi,\\Lambda)$ 的协同对角化因子 | 含（$\\Pi,\\Lambda$ 都是 SDP 决策变量） |
| 1603.04172 Thm 5.2 | $E_t$ 对角化 $\\Pi_{t|t-1}$，且 (5.14) 与解自洽 | 含（预测误差协方差由 test channel 决定） |

$\\Rightarrow$ **三篇没有一篇在"给定的权阵/观测阵的本征基"里做水填。**
所以 C 的 §66-C/§67-B 那条"$S^*(D)$ 是 $\\Theta$-对角的水填"，问题不只是"缺证明"：
$\\Theta=K^{\\mathsf T}(R+B^{\\mathsf T}P_cB)K$ 是**代价权**（给定数据），而这条规律说"能被水填切的谱必须是解自己的协方差的谱"。
我 e152 的独立实测与此同向（最优 rank-1 设计离开 $\\Theta$ 首轴 $2.68$–$6.28°$，代价 $0.009$–$0.085$ bit；rank-2 到 $40°$、$0.214$ bit）。
**这不是"我方对、C 方错"的记分**：C 观察到倾角是真的（§100-E① 已承认），错的只是把 $\\Theta$ 当成水填所在的坐标系。

### 103-D 未读与挂账（照实写）

未读：§IV 的 (4.1) $D_{0,n}(R)$；Thm 5.2 的证明本体；"universal MSE lower bound" 那一节的逐步推导 $-$ $-$
所以"他们的地板可与设计耦合"这句我只写到**机制层**（(5.14) 含解），没有写成定理。
挂账不变：e152 的 rank-3/4 自由侧、1711.09853 (17) 对 1810.00298 Lemma 2 (29)(30) 的逐项对表、CDC-2018 逐字段核、
`e138b` 的 429 重试结果。#44（著录）已清，见 §102-B。
'''

for ph in ('XXX', '__U', 'TODO'):
    assert ph not in SEC, '占位符残留：%s' % ph
assert SEC.count('### 103-') == 4

out = txt.replace(row_anchor + '\n', row_anchor + '\n' + ROW + '\n')
out = out.rstrip('\n') + '\n\n' + SEC.lstrip('\n')
b = out.replace('\n', '\r\n').encode('utf-8')
assert b.count(b'\r') == b.count(b'\n') == b.count(b'\r\n')
ctrl = sorted(collections.Counter(x for x in b if x < 32 and x not in (10, 13)).items())
assert ctrl == [(9, 3), (11, 1)], 'C 的控制字节变了：%s' % ctrl
now = open(P, 'rb').read()
assert now == raw, '板子被并发改写（%d→%d），本次放弃' % (len(raw), len(now))
old = txt.split('\n')
open('p0/board_before_103.md', 'wb').write(raw)
assert open('p0/board_before_103.md', 'rb').read() == raw
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
assert len([x for x in chk if x.startswith('| §103 | R77-L |')]) == 1
assert len([x for x in chk if x.startswith('### 103-')]) == 4, '103-A..D 不齐'
print('旧行按序保留到第 %d 个非空行（新板非空 %d）' % (i, len(nl)))
print('OK lines=%d bytes=%d ctrl=%s' % (len(chk) - 1, len(b), ctrl))
