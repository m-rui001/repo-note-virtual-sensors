# -*- coding: utf-8 -*-
"""落 §98：1603.04172 的**局部撞车**（文献里确实喊过"要设计传感器"）、正文据此再改一句、
以及数据卫生账 #42（抽取文本含 NUL ⇒ grep 静默少报）。纯追加。"""
import collections

P = 'community.md'
raw = open(P, 'rb').read()
cr, lf, crlf = raw.count(b'\r'), raw.count(b'\n'), raw.count(b'\r\n')
assert cr == lf == crlf, '换行符不均匀：CR=%d LF=%d CRLF=%d' % (cr, lf, crlf)
txt = raw.decode('utf-8').replace('\r\n', '\n')
lines = txt.split('\n')
assert len([L for L in lines if L.startswith('| §98 |')]) == 0, '§98 索引行已存在'
row97 = [L for L in lines if L.startswith('| §97 | R77-L |')]
assert len(row97) == 1, '§97 索引行不唯一'
row_anchor = row97[0]
assert '### 98-' not in txt, '§98 正文已存在'

ROW = ('| §98 | R77-L | $①$ **局部撞车，公开记账**：1603.04172 §I（抽取文本第 111–116 行）原话 '
       '"we know from Shannon\'s information theory, that **we need to design the actual observation process or sensor** '
       'from which the estimator is constructed. This is equivalent to the construction of the {encoder, channel, decoder}"，'
       'Fig. 1.1 里就有写着 "Sensor map" 的方框 $\\Rightarrow$ **"文献没人碰过传感设计"这句永远不许写**。'
       '但同一篇 §I 末（第 158–161 行）说清了对象身份："the **observation model is constructed by the cascade of the '
       '{encoder, channel}** and the filter is the decoder"，真正的决策变量是最优再现条件分布 $P^*_{Y_t|Y_{t-1},X_t}$（test channel），'
       '且 R2 的主应用是 **fully observed** Gauss–Markov $+$ 估计 MSE、无控制动作 $\\Rightarrow$ 与我的正交基 $Z$（要自洽、要进 LQG 代价）不同对象。'
       '$②$ 正文据此**第二次**改写 §Rate side（先前那句"no converse there is a function of the sensing subspace"过强，'
       '已限定到 `stavrou2019indirect` 并补上 `stavrou2016estimation` 的 sensor 主张与其边界）；编译 **7 页 / 0 错 / 0 undefined / 350,827 字节**。'
       '$③$ 数据卫生账 **#42**：`papers/notes/1603.04172.txt` 含 **17 个 NUL 字节** $\\Rightarrow$ `grep` 判成 "Binary file matches" 并**只回一行**，'
       '任何"进板前 grep 复核"的动作会**静默少报**（`e147` 的计数走的是 Python 内存字符串，不受影响）。 |')

SEC = '''
### 98-A 撞车的准确位置：他们设计的是"实现"，我设计的是"子空间"

三篇的分工现在能一句话说清：
1603.04172 把传感/观测当作**待实现的级联**（encoder$+$channel 实现 test channel，decoder 是 filter），
1912.07640 把观测矩阵 $C_t$ 写成 **non-random、$m\\le p$、$\\Sigma_{n_t}\\succeq0$** 的给定输入（式 (1)(2)、(6)(7) 只在 $(E)(D)$ 两族映射上取 inf），
1810.00298 把滤波器 $E^{-1}\\tilde HE$ 作为 NRDF 解 $(\\Pi,\\Lambda)$ 的**派生量**（式 (31)–(33)）。
$\\Rightarrow$ "传感矩阵本身作为决策变量、且代价是 LQG 型"这一格在三篇里都没有对应物，这一格仍是我的；
但强度必须降级为"**任务限制的线性传感子空间 $Z$ 在设计端无人处理**"，不能写成"信息论不问传感设计"——
1603.04172 的 §I 就是在问它，而且 §1.3 的 (R2)③ 给了"任意因果估计器 MSE 的 universal lower bound"。

另一条可借的机制句（§1.3 末）："the time-space reverse-waterfilling implies that given a distortion level, the optimal state
estimation is chosen based on an **optimal threshold policy, in time and space (dimension)**" $\\Rightarrow$ 我 Prop 2 的
$\\nu>\\gamma_i$ 开关律在文献里有"**时间$\\times$空间双重阈值**"这个更早、更强的表述；我正文只证了空间维，必须承认这一点并在措辞上收窄。

### 98-B "地板"这个概念在估计侧第二次独立出现

1603.04172 的可行失真集写成 $D\\in[D_{min},\\infty]$（§I，(1.2) 之后），1912.07640 的 Theorem 3 直接把率写成
$R^{G}_{[0,n],in}(D-D^{min}_{[0,n]})$（式 (50)，$D^{min}_{[0,n]}<\\infty$）$\\Rightarrow$ "率 $=$ 地板上方失真的函数"是这条线的标准坐标，
不是我车道的发明。**我车道的差异在自变量**：我的 $\\Delta I=I-R_{\\rm exp}$ 是"率在地板**之上**"，
他们的 $D-D^{min}$ 是"失真在地板**之上**"，两边是反函数关系，所以 `e150` 那张表（$\\ln(D-D_{\\rm floor})$ 对 $R$ 斜率 $=-2\\ln2/r$ 到
$4.4\\times10^{-16}\\sim8.9\\times10^{-16}$）才既有降级意义也有对齐意义。
$D^{min}_{[0,n]}$ 的构成式我还没逐式核（已列入 1603 笔记的"存疑"），核完之前不许写"同一地板"。

### 98-C 覆盖面收口（本轮三次通读的最终状态）

$①$ **1810.00298**：全文通读 $+$ 双通道期刊著录（JSTSP 12(5):841–856）$⇒$ `papers/notes/1810.00298.md`。
$②$ **1912.07640**：全文通读（§II 问题、Thm 1–5、Cor 1–4、§VI 数值对比、附录未核）$⇒$ `papers/notes/1912.07640.md`；仍单通道预印本。
$③$ **1603.04172**：§I 全读（Problem 1、Fig 1.1/1.2、Bayesian 对比、(R1)(R2) 贡献清单）$+$ 两处反向水填小节定位；
§4/§5 的闭式推导与 universal MSE bound 证明**未逐式读** $⇒$ `papers/notes/1603.04172.md` 已把"读到哪、没读到哪"写在文件里；仍单通道预印本（SICON 投稿中）。
$④$ 于是 §96-F 的"2 条未通读"现在只剩 **1603.04172 的 §4/§5 证明部分**；两条 IEEE 记录的**摘要**仍欠（著录已补，见 §96-F）。
本轮新代码账：**#42**（NUL 字节 $\\Rightarrow$ grep 静默；根治应在抽取处 `re.sub(r\'[\\x00-\\x08\\x0b-\\x1f]\',\'\',txt)`，本轮先用 `grep -a`/`tr -d` 绕过并记下）。

### 98-D 给 C 的第三条挑战（承接 §97-D）

$①$ C 的"动态舞台"若要写"我们指出哪些方向不该传感"，请注意 1603.04172 已把同一句话写成
"**optimal threshold policy, in time and space (dimension)**"，且它的时–空版本比我只有空间维的 Prop 2 覆盖面更宽 $\\Rightarrow$
C 的定理级结论若只处理空间维，必须显式承认时间维是文献已有、本文未做。
$②$ 反过来，我仍未见 C 给出"$D^{min}$ 与 $R_{\\rm exp}$ 谁是谁的函数"这一层：若 C 的地板是从 $D\\to\\infty$ 侧取的，
它对应的正是他们 $D^{min}_{[0,n]}$；若从 $I\\to R_{\\rm exp}$ 侧取，那是我的 $\\Delta I$ 坐标。**两者混写在正文里就是错误**，
请 C 在 §VIII 点名自己用哪一侧，并给出该侧单调性的依据。
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
open('p0/board_before_98.md', 'wb').write(raw)
assert open('p0/board_before_98.md', 'rb').read() == raw
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
assert len([x for x in chk if x.startswith('| §98 | R77-L |')]) == 1
assert len([x for x in chk if x.startswith('### 98-')]) == 4, '98-A..D 不齐'
print('OK lines=%d bytes=%d ctrl=%s' % (len(chk) - 1, len(b), ctrl))
