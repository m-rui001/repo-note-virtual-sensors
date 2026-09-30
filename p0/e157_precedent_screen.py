#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
e157：把"无先例"这句话的证据补硬。走 arXiv 通道（export API），绕开 Semantic Scholar 的 429。
本轮要回答的是本站唯一那句强陈述的前提：
  有没有一篇公开工作，同时做到 (a) 把**传感/观测矩阵本身**当决策变量，
  (b) 失真用**控制代价权**（加权 trace，不是无权 trace），
  (c) 给出**随传感秩下行的地板族**？
判据在本文件开头写死，不看结果再挑：
  [Q0] 我方正文里现在有哪些"无先例/无人处理"级句子（逐行打印，含文件名行号）。
  [Q1] ITW-2017（arXiv 1701.06368）到底是谁、讲什么：取回摘要，按三问逐项标。
  [Q2] 五组关键词检索，命中条目逐条标 (a)(b)(c) 三个关键词族是否出现；
       命中里有 **三项同时出现** 的条目 => [Q-NEG] 触发，"无先例"当场作废、改写为"与 X 的差别是 ..."。
       没有任何条目三项同现 => "无先例"降级为"**关键词层面未见**"，仍不许写成定理级"无人处理"。
"""
import re, sys, time, urllib.parse, urllib.request
import xml.etree.ElementTree as ET
sys.stdout.reconfigure(encoding='utf-8')

NS = {'a': 'http://www.w3.org/2005/Atom'}
OUT = []


def p(*a):
    s = ' '.join(str(x) for x in a)
    OUT.append(s)
    print(s)


def atom(url):
    req = urllib.request.Request(url, headers={'User-Agent': 'rc-lane-note/1.0'})
    return ET.fromstring(urllib.request.urlopen(req, timeout=60).read())


def entries(root):
    out = []
    for e in root.findall('a:entry', NS):
        t = re.sub(r'\s+', ' ', (e.findtext('a:title', '', NS) or '')).strip()
        ab = re.sub(r'\s+', ' ', (e.findtext('a:summary', '', NS) or '')).strip()
        aid = (e.findtext('a:id', '', NS) or '').strip()
        auth = [re.sub(r'\s+', ' ', (n.findtext('a:name', '', NS) or '')).strip()
                for n in e.findall('a:author', NS)]
        pub = (e.findtext('a:published', '', NS) or '')[:10]
        out.append((aid, t, auth, pub, ab))
    return out


A = re.compile(r'sensor|sensing|observation matrix|measurement matrix|encoding matrix|\bC_k\b|probe', re.I)
B = re.compile(r'linear[- ]quadratic|LQG|control cost|cost of control|certainty equivalence|weighted trace|task[- ]driven|cost function', re.I)
C = re.compile(r'rank[- ]constrained|rank of the (sensor|measurement|observation)|dimension of the (sensor|measurement)|reduced[- ]rank|low[- ]rank (sensor|measurement|design)', re.I)


def flag(txt):
    return ('a' if A.search(txt) else '-') + ('b' if B.search(txt) else '-') + ('c' if C.search(txt) else '-')


p('== [Q0] 我方正文里的"无先例/无人处理"级句子 ==')
for f in ('p0/note/note.tex', 'p0/note/frag_rc_waterfill.tex', 'p0/note/frag_D_repro_infeasible_nofactor.tex'):
    src = open(f, encoding='utf-8', errors='replace').read().split('\n')
    PAT = re.compile(r'no prior|no precedent|unaddressed|nobody|no published|not (?:been )?(?:treated|studied|addressed)|'
                     r'no counterpart|absent from|is the first|only .{0,20}who', re.I)
    hits = [(i + 1, l.strip()) for i, l in enumerate(src) if PAT.search(l)]
    p('--- %s：%d 处' % (f, len(hits)))
    for n, l in hits:
        ctx = ' '.join(x.strip() for x in src[max(0, n - 3):n + 2])
        p('  %d| %s' % (n, re.sub(r'\s+', ' ', ctx)[:330]))

p('')
p('== [Q1] ITW-2017 = arXiv 1701.06368（arXiv 通道，绕开 429）==')
try:
    r = atom('http://export.arxiv.org/api/query?id_list=1701.06368')
    for aid, t, au, pub, ab in entries(r):
        p('id      : %s' % aid)
        p('标题    : %s' % t)
        p('作者    : %s' % ', '.join(au))
        p('发表    : %s' % pub)
        p('摘要    : %s' % ab)
        p('[三问标] %s  （a=传感矩阵是决策变量 b=控制代价/加权失真 c=随秩地板）—— 关键词级筛查，不是通读' % flag(t + ' ' + ab))
except Exception as e:
    p('[Q1-FAIL] %s：%s' % (type(e).__name__, str(e)[:160]))

p('')
p('== [Q2] 五组检索（arXiv 通道；判据见文件头）==')
Q = [
    'all:"reverse waterfilling" AND all:sensor',
    'all:"rate-distortion" AND all:"sensor design"',
    'all:"linear quadratic" AND all:"rate distortion" AND all:"limited sensing"',
    'abs:"water-filling" AND abs:"measurement matrix"',
    'all:"task-oriented sensing" AND all:"water filling"',
]
tot = 0
triple = []
for q in Q:
    url = 'http://export.arxiv.org/api/query?' + urllib.parse.urlencode(
        {'search_query': q, 'start': 0, 'max_results': 8})
    try:
        r = atom(url)
        es = entries(r)
    except Exception as e:
        p('-- %s\n   [FAIL] %s：%s' % (q, type(e).__name__, str(e)[:140]))
        continue
    p('-- %s  命中 %d 条' % (q, len(es)))
    for aid, t, au, pub, ab in es:
        tot += 1
        fl = flag(t + ' ' + ab)
        p('   [%s] %s (%s) %s' % (fl, t[:96], pub, aid.split('/abs/')[-1]))
        if fl == 'abc':
            triple.append((t, aid, ab[:300]))
    time.sleep(3.0)

p('')
p('== [Q-判决行]（判据：三项同现 => "无先例"作废；否则只许写"关键词层面未见"）==')
p('[Q2] 累计列出 %d 条（含跨查询重复）；三项同现条目 %d 条' % (tot, len(triple)))
for t, aid, ab in triple:
    p('   !!! 撞车候选：%s  %s' % (t, aid))
    p('       %s' % ab)
if not triple:
    p('[结论] 未见三项同现 $-$ $-$ "无先例"只能写成 **"关键词检索未见同现"**，不许写"无人处理"。')
    p('[遗留] 关键词筛查不等于通读；撞车风险最高的仍是 (a) 单问那一类（sensor/observation 设计），'
      '因为 (b)(c) 两族词在本领域摘要里出现率本就低。')
else:
    p('[结论] 有撞车候选 => 正文那句"无先例"必须改写为与该条目的具体差别陈述。')

open('p0/e157_out.txt', 'w', encoding='utf-8').write('\n'.join(OUT) + '\n')
