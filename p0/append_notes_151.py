# -*- coding: utf-8 -*-
"""把 e151 的文献侧结论写进两篇笔记：1912 的 Prop 3 三条件（原先列为"没读完"）与 1810 的 Thm 4 对角化基。"""
import sys
sys.stdout.reconfigure(encoding='utf-8')

# ---------- 1912.07640.md ----------
P1 = 'papers/notes/1912.07640.md'
t1 = open(P1, encoding='utf-8').read()
old_line = 'Proposition 3(i) 的"强结构性质"具体是什么条件没读完（这决定 Thm 5 的闭式在我锚株上是否适用）；'
assert t1.count(old_line) == 1, '1912 笔记的存疑行不唯一/不存在'
new_line = 'Proposition 3 的"强结构性质"**已读全**（见下节），其在本站六株上的适用性已判定为 $0/6$；'
t1 = t1.replace(old_line, new_line)

SEC1 = '''
## Proposition 3 / Theorem 5 的门（本轮读全：抽取文本第 2073–2092 行）

这条线的"结构简化到逆水填"**不是无条件的**。原文：

> Proposition 3 (Strong structural properties on (60))：在 (60) 的刻画里，若 $(A,\\bar\\Sigma)$ 之间满足三者之一 $\\Rightarrow$
> (i) $A=\\alpha I_p$（标量阵）且 $\\bar\\Sigma\\succeq0$；(ii) $A$ 实对称且 $\\bar\\Sigma=\\sigma^2_{\\bar\\Sigma}I_p$（标量阵）；
> (iii) $A=\\bar\\Sigma\\succ0$。则 $(A,\\Sigma_\\xi,\\bar\\Sigma)$ **两两对易（commute by pairs）**，进而 $(\\Sigma_\\xi,\\Pi_\\xi)$ 对易。

> Theorem 5 (Optimal numerical solution of (60))：**Suppose that one of the conditions of Proposition 3 hold.**

三点必须记准，否则会误引：

1. $\\bar\\Sigma$ **不是初值协方差**。式 (61) 是 $\\Pi_\\xi=A\\Sigma_\\xi A^{\\mathsf T}+\\bar\\Sigma$（预报协方差的 Lyapunov 形），
   且第 2889 行明写 "the solution is $\\bar\\Sigma=\\lim_{t\\to\\infty}\\Sigma_n$" $\\Rightarrow$ 它是**噪声协方差的稳态值**。
   本站植物是时不变噪声，所以 $\\bar\\Sigma=W$ 逐字成立，不需要取极限。
2. 被水填的对象是 **$(\\Sigma_\\xi,\\Pi_\\xi)$ 这一对**，而 $\\Sigma_\\xi$ 正是 (60) 的决策变量
   （$\\inf_{0\\prec\\Sigma_\\xi\\preceq\\Pi_\\xi,\\ \\mathrm{trace}(\\Sigma_\\xi)\\le D-D^{min}_{[0,\\infty]}}\\tfrac12\\log|\\Pi_\\xi|/|\\Sigma_\\xi|$，
   式 (60) 原文第 1989–1996 行）$\\Rightarrow$ **"在哪个基里对角"由解自己决定**。
   这与本站 `E83` 的机制是同一句话：被水填的那组数本身随设计移动。
3. 第 2140 行补了逃生门："However, if for instance in Proposition 3, (i), $\\bar\\Sigma\\succeq0$ we use $\\dots$"
   $\\Rightarrow$ 条件 (i) 只钉 $A$ 为标量阵、$\\bar\\Sigma$ 可以一般；也说明作者把这三条当**充分而非必要**，
   所以"条件不成立"**不等于**"闭式不成立"，只能说"此处无文献依据可借"。

**本站判定**（`p0/e151_waterfill_basis_audit.py` $\\to$ `e151_out.txt` [W2]，六株逐一打印
$|A-\\alpha I|/|A|$、$|A-A^{\\mathsf T}|/|A|$、$|W-\\sigma I|/|W|$、$|\\Theta-\\tau I|/|\\Theta|$）：
锚点 $+$ `rand-1..4` $+$ `big-6` 里 **(i)(ii)(iii) 全 False，命中 $0/6$**。
锚点 $|A-A^{\\mathsf T}|/|A|=1.136$（$A$ 远离对称）、$|W-\\sigma I|/|W|=0.535$、$|\\Theta-\\tau I|/|\\Theta|=0.808$。
$\\Rightarrow$ 引 Thm 5 的闭式必须带条件号；在我这六株上它不适用，只能写成"文献在 $A=\\alpha I$ 等特殊结构下有闭式"。
'''
open(P1, 'w', encoding='utf-8', newline='\n').write(t1.rstrip('\n') + '\n' + SEC1)
chk = open(P1, encoding='utf-8').read()
assert 'Proposition 3 / Theorem 5 的门' in chk and old_line not in chk
print('1912 笔记 OK，字符 %d → %d' % (len(t1), len(chk)))

# ---------- 1810.00298.md ----------
P2 = 'papers/notes/1810.00298.md'
t2 = open(P2, encoding='utf-8').read()
SEC2 = '''
## Thm 4 的"对角"是在**谁的**基里（第 774–790、1530–1566 行）——这条决定能不能用它支撑"固定基对角"

Theorem 4 给的是 (27) 的**等价实现方案**：$y_t=E^{-1}\\tilde H E x_t+(I-E^{-1}\\tilde HE)Ay_{t-1}+E^{-1}\\Theta v_t$ (31)，
其中 "$E\\in\\mathbb R^{p\\times p}$ is a **non-singular matrix that simultaneously diagonalizes $\\Pi\\succ0,\\Lambda\\succ0$**" (32)，
水填设计阵 $\\tilde H=I-\\tilde\\Pi\\tilde\\Lambda^{-1}\\equiv\\Theta\\Phi$ (33a)。附录 B（第 1532–1561 行）说清了这个基怎么来的：
取 $\\Pi=U^{\\mathsf T}\\tilde\\Pi U$ 与 $\\Pi^{-1/2}\\Lambda\\Pi^{-1/2}=V^{\\mathsf T}SV$，则 **$E\\triangleq\\tilde\\Pi^{1/2}V\\Pi^{-1/2}$**，
并且作者注明这是 "a version of the cogredient diagonalization approach derived in [47, Theorem 8.3.1]"
（其 [47] $=$ D. S. Bernstein, *Matrix Mathematics: Theory, Facts, and Formulas*, 2nd ed., Princeton, 2011 $-$ $-$ **不是 Horn–Johnson**，
抄 bibitem 时别按惯例猜）。

三条判读：

1. $\\Pi,\\Lambda$ 都是 Lemma 2 的半定表示 (29)(30) 里的**决策变量** $\\Rightarrow$ 对角化的基是**由解导出的**，
   不是任何给定数据阵（源协方差、任务加权、$A$ 的对称化）的本征基。
   两阵 SPD 的合同同时对角化**总是存在**（Bernstein 8.3.1 那一类事实），所以 Thm 4 是无条件的——
   但正因为它无条件，它**不可能**给出"在某个事先指定的基里对角"这种带内容的结构陈述。
1b. **水填只在秩亏档发生**（第 795–798 行原文）：
    "a) Full-rank $\\tilde H$: If $\\tilde H\\succ0$ $\\dots$ then, **no reverse-waterfilling occurs (in dimension)** and $\\Theta\\succ0,\\Phi\\succ0$；
    b) Rank-deficient $\\tilde H$: $\\dots$ then, **the reverse-waterfilling kicks in** and $\\Theta\\succeq0,\\Phi\\succeq0$"，
    其中 (33b) $\\Theta=\\tilde\\Sigma_v^{1/2}$、$\\tilde\\Sigma_v=\\tilde\\Pi\\tilde H$，(33c) $\\Phi=(\\tilde H\\tilde\\Pi^{-1})^{1/2}$
    $\\Rightarrow$ 他们的"水填"是**按维数关掉若干维**这一件事，且 $r=p$ 时压根不叫水填。
    这条与我 Prop 2 的开关律同型（开关由水位定），可作为"秩亏才谈水填"的文献先例引用。
2. 该文第 135 行把这类基称为 "non-singular joint diagonalizers (**KLT matrices**)"，第 660 行说方案 "makes use of **joint diagonalization matrices**"
   $\\Rightarrow$ 作者自己的词是"同时对角化阵"，读者不能替换成"源的特征基"。
3. 注意符号撞车：这里 (33a)–(33b) 的 $\\Theta\\triangleq\\tilde\\Sigma_v^{1/2}$（$\\tilde\\Sigma_v=\\tilde\\Pi\\tilde H$，对角 PSD 缩放阵），
   与本站代价侧的任务加权 $\\Theta=K^{\\mathsf T}(R+B^{\\mathsf T}P_cB)K$ **同名不同物**；
   该篇另一处 (31) 的 $E^{-1}\\Theta v_t$ 用的还是这个对角阵。引用时若不点名，必然写成错误归因。

$\\Rightarrow$ 用法限定：想引 1810.00298 说"逆水填在**给定**矩阵的本征基里做"，本篇**不支持**；
它支持的是"存在一个由 $(\\Pi,\\Lambda)$ 决定的基，使解在该基里对角、并可读成逐维逆水填"。
'''
open(P2, 'w', encoding='utf-8', newline='\n').write(t2.rstrip('\n') + '\n' + SEC2)
chk2 = open(P2, encoding='utf-8').read()
assert 'Thm 4 的"对角"是在**谁的**基里' in chk2
print('1810 笔记 OK，字符 %d → %d' % (len(t2), len(chk2)))
