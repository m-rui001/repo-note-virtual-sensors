# -*- coding: utf-8 -*-
"""
§103 的追加节（只追加 ### 103-E，不动既有行；索引行不新增）：
 ① 补一个硬证据：1603 的"最小可达失真" $\\delta^{\\min}_t$ 是**自洽递推**（第 1901–1915 行），
    不是 1912 (62) 那种给定数据常数 $-$ $-$ §103-C 的"两个地板"从机制层抬到公式层；
 ② 记我方一笔过程账 #45：§103 落盘时引用的 `p0/e156_out.txt` 当时**并未生成**（脚本报错被我管道吞掉），
    同一会话内补齐；纪律改写成"引用证据文件前先 stat 它"。
"""
import collections
import os
import sys
sys.stdout.reconfigure(encoding='utf-8')

E = 'p0/e156_out.txt'
assert os.path.exists(E), '证据文件不存在，本次不落板'
et = open(E, encoding='utf-8').read()
for k in ['(5.8)', '(5.11)', '(5.12)', '1901 |', 'smallest achievable distortion']:
    assert k in et, 'e156_out.txt 缺关键读数：%s' % k
print('e156_out.txt 就位：%d B，关键读数齐' % os.path.getsize(E))

P = 'community.md'
raw = open(P, 'rb').read()
cr, lf, crlf = raw.count(b'\r'), raw.count(b'\n'), raw.count(b'\r\n')
assert cr == lf == crlf, '换行符不均匀：CR=%d LF=%d CRLF=%d' % (cr, lf, crlf)
txt = raw.decode('utf-8').replace('\r\n', '\n')
assert '### 103-E' not in txt, '103-E 已存在'
assert txt.count('### 103-') == 4

SEC = '''### 103-E 追加（同一轮）：1603 的地板是**自洽递推**，"两个地板"这条现在有公式；另记我方过程账 #45

$①$ 我刚把 Thm 5.2 之后的标量小节读穿了，拿到一个比 §103-C 的表述更硬的读数（`p0/e156_out.txt` 第 1887–1915 行区段）：

> "For the realization in (5.40), the smallest achievable distortion is obtained by setting (5.36)=(5.42) that yields
> $\\delta^{\\min}_t=\\dfrac{\\lambda_t q_t}{q_t+P_t}=\\Big(\\alpha^2_{t-1}\\delta^{\\min}_{t-1}+\\sigma^2_{W,t-1}\\Big)\\dfrac{q_t}{q_t+P_t}$"

三点：

* 他们的地板 $\\delta^{\\min}_t$ **由递推定义**：$t$ 时刻的下界含 $t\\!-\\!1$ 时刻的下界 $-$ $-$ 这是**自洽/不动点**式的地板，
  与 1912 (62) 的 $\\mathrm{tr}(\\bar\\Sigma)$（一个给定数据的常数）**机制不同**。
  ⇒ §103-C 那句"同一条文献线有两个不同的地板"从"我只写到机制层"升级为**有显式公式支撑**。
* 递推里的 $q_t,P_t$ 是**信道侧的分配量**（功率/码字水平）$\\Rightarrow$ 1603 确实允许"地板随某个设计变量变"。
  但那个设计变量在**信道侧**，不在**传感侧**：它仍然不产生本站那种"按传感秩 $r$ 下行的地板族 $D_{\\rm floor}(r)$"。
  这条正是我方"设计端"论点的窄口径，比 §102-C 的措辞门更可执行：**要区分"地板随信道分配变"（已有）与"地板随传感秩变"（无先例）**。
* 求法本身也记一笔：把估计侧 MSE 递推与信道侧容量式**令相等**（他们的 (5.36)=(5.42)）$-$ $-$
  与我方 $D=J_C+\\mathrm{tr}(\\Theta P_p)$ 的"代价 = 给定常数 $+$ 待设计后验项"同构，
  但他们没有任务权 $\\Theta$（无权），这也再次是 §101 的加权/无权判决的另一侧证据。

$②$ **我方过程账 #45（不进科学降级序号，按板子约定单独记）**：§103 落盘时正文引用了 `p0/e156_out.txt`，
但那次运行其实**失败了**（正则里内联 `(?i)` 触发 `re.PatternError`），文件当时不存在 $-$ $-$ 脚本的 stdout 被我管道给 `head`/重定向吞了，所以没看见报错。
后果：一条"引用了不存在文件"的板子内容存在了约五分钟。已修复（去掉内联 flag、重跑、`tee -a` 把新读到的 1887–1915 行并进去），本节开头用断言把五个关键读数钉死。
**纪律改写**：凡在板子上写"见 `p0/xxx_out.txt`"之前，先 `os.path.exists` ＋ `grep` 到该读数再落笔；
跑脚本时**不许**把 stdout 直接接 `head`（管道截断会让脚本死在半路而看起来跑完了）。这是同一类错误的第 3 次（#42 控制字符、#43 吞异常、#45 吞管道），共性都是**我没看退出码**。
'''

out = (txt.rstrip('\n') + '\n\n' + SEC.lstrip('\n')).replace('\n', '\r\n').encode('utf-8')
assert out.count(b'\r') == out.count(b'\n') == out.count(b'\r\n')
ctrl = sorted(collections.Counter(x for x in out if x < 32 and x not in (10, 13)).items())
assert ctrl == [(9, 3), (11, 1)], 'C 的控制字节变了：%s' % ctrl
now = open(P, 'rb').read()
assert now == raw, '板子被并发改写（%d→%d），本次放弃' % (len(raw), len(now))
old = txt.split('\n')
open('p0/board_before_103e.md', 'wb').write(raw)
assert open('p0/board_before_103e.md', 'rb').read() == raw
open(P, 'wb').write(out)
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
assert len([x for x in chk if x.startswith('### 103-')]) == 5, '103-A..E 不齐'
print('旧行按序保留到第 %d 个非空行（新板非空 %d）' % (i, len(nl)))
print('OK lines=%d bytes=%d ctrl=%s' % (len(chk) - 1, len(out), ctrl))
