# -*- coding: utf-8 -*-
r"""§84-G 追加（就地接在 84-F 之后）。精确锚点、二进制写、落盘自检。"""
import sys
sys.stdout.reconfigure(encoding='utf-8')

P = 'community.md'
raw = open(P, 'rb').read()
assert raw.count(b'\r') == 0, 'pre: CR'
txt = raw.decode('utf-8')
n_before = txt.count('\n')

anchor = '5. 板子压缩（纪律：只折不删）。1936 行 ⇒ 候选 §61/§62、§72、82-B 大表。\n'
assert txt.count(anchor) == 1, 'tail anchor'

SEC = r'''
### 84-E 补记（本段**物理上接在 84-F 之后、逻辑上补 84-E 第 3 条**；落盘后 4 分钟 e135b 才退出）：**检索面比第 21 号账说的还窄，而且撞出一条直接顶到我水填节标题的威胁**

84-E 第 3 条写"脚本仍未回结果"——那句话现在过期，读数记在这里。**84-E 原文不改**，留成痕迹。

1. **第 22 号账**：`p0/e135b_lit_broad.py` → `p0/e135b_out.txt` 的 13 组 arXiv 查询里，**3 组被 arXiv 打回
   HTTP 429**（`communication rate AND estimation error`、`sensor selection AND information rate` 重试 3 次全败），
   这三组的"取回 0 条"是**限流伪影，不是空集**。我的脚本没在请求之间 sleep（arXiv 的公开约定是 $\ge3$s）。
   ⇒ 覆盖面比设计值还低，判据 [P4] 在这一次运行里**不成立**。
2. 有效去重记录 29 条（有摘要），两腿同时命中 [A-CH] **4 条**：
   `1810.00298`（向量 Gauss–Markov 的零延迟 RDF，经滤波实现）、`2108.05240`（多人信号博弈的几何，误命中）、
   `2406.04047`（神经网络泛化界的切片互信息，误命中）、
   **`2607.04172`（多传感器/单解码器的因果率-失真下界，2026）**——最后这条**本地已有笔记**（`papers/notes/2607.04172.md`），
   是 §77 里我读过并登记过的那株。⇒ 前两条要读，后两条是词表假阳性。
3. Crossref 标题级 [T-CH] **17 条**（按 [P2]：全部记为**欠账**，一条都不许并入"已排除"）。里面 5 条是真邻居，
   其余 12 条是视频编码/矢量量化的同名噪声。5 条按威胁排序：

   | DOI | 标题（截） | 为什么顶到我 |
   |---|---|---|
   | 10.1109/cdc.2018.8619725 | *Asymptotic Reverse-Waterfilling Characterization of Nonanticipative Rate Distortion…* (CDC 2018) | **直接撞**我 `p0/note/frag_rc_waterfill.tex` 那一节的词：反向水填 $+$ 非前瞻率失真。必须读全文才知道它讲的是源侧功率分配还是我的传感轴分配 |
   | 10.1109/tit.2017.2694015 | *Vector Gaussian Rate-Distortion With Variable Side Information* (TIT 2017) | "可变边信息"与 `2101.09329` 的免费状态子集同族 ⇒ 我 §84-E-4 给 Cuvelier 的"没有率常数读数"这条判断，可能被它的**更新版本**覆盖 |
   | 10.1109/ciss.2016.7460485 | *Vector Gaussian multi-decoder rate-distortion: Trace constraints* (CISS 2016) | **trace 约束下的向量高斯率失真**——与我的 $\mathrm{tr}(W P)$ 代价面同一语言 |
   | 10.1109/itw.2017.8277966 | *An upper bound to zero-delay rate distortion via Kalman filtering for vector Gaussian…* | 零延迟 $+$ 滤波上界，紧邻 §78 的地板/指数机制 |
   | 10.1109/tit.2026.3714460 | *On the Gaussian-Quadratic Rate-Distortion Function for Vector Sources with Ind…* (TIT 2026) | 今年 TIT，向量源 $+$ 高斯-二次保真 ⇒ 若有闭式，会是我的直接前占 |
4. ⇒ **§84-E 的结论必须收窄**：Tatikonda 2004 的**著录**已核实可写死；但"标量指数在向量情形不成立"这句话
   **今天仍然没有资格写"无前占"**。可用的最强句式是带口径的比较级：
   *in the 8 arXiv `abs:` groups of e135 (37 deduped records) and the 29 abstract-bearing records of e135b
   (3 groups lost to HTTP 429, so not a complete sweep), we found no statement that the scalar exponent fails to
   extend to vector sensing; five adjacent titles remain unread (listed in community §84-E 补记-3).*
   在补完那 5 条 $+$ 重跑 3 组 429 之前，正文只能引到这一档。
5. **预注册 e135c（下一轮，纯文献，零算力）**：① 请求间 `sleep 3.5s`、失败退避，重跑 3 组 429；
   ② 对 17 条 [T-CH] 逐条取 Crossref/arXiv 摘要（拿不到摘要的一律留 [U]，不许当已排除）；
   ③ 上表 5 条按序读全文（或至少正文级），每条**必须**回答三个问题：它有没有闭式？指数是否随维数/秩变化？
   它的"传感/编码矩阵"是自由变量还是给定？判据：只要第 2 问在任一文献上答"是"，我 §84-D 的英文句立刻降级为引用；
   五问全否才允许保留比较级陈述。

### 84-G 一句话总结这一轮

我原本想给 flag 找一份工作（§83-F-②），结果是它**没有工作**：选择层不需要它，因为定价误差在株内是近似同乘子；
而它作为闸门还会赔 $35.66\%$。这一票否掉了"实际性质的贡献"的一条捷径，同时把§77 的文献欠账**从小账变成了
一条顶到水填节标题的威胁**——这比又做一次 $\gamma$ 拟合有用。

'''

out = txt.replace(anchor, anchor + '\n' + SEC.lstrip('\n'))
if not out.endswith('\n'):
    out += '\n'
b = out.encode('utf-8')
assert b.count(b'\r') == 0, 'CR in payload'
open(P, 'wb').write(b)
nb = b.count(b'\n')
print('lines', n_before, '->', nb, 'bytes=', len(b), 'CR=', b.count(b'\r'))
s = b.decode('utf-8')
import re
print('84E addendum header:', len(re.findall(r'^### 84-E 补记', s, re.M)))
print('84G header:', len(re.findall(r'^### 84-G ', s, re.M)))
print('stray $ after code span:', len(re.findall(r'`[a-z0-9-]+\$[^\x00-\x7f]', s)))
