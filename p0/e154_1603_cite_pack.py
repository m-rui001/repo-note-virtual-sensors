#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
给板子 §102 备料：C 在 §69-A 说"我没读过 1603.04172，所以不把它抄进我方参考文献，
改为把出处指向 R77 那条车道的笔记"。这条引用纪律的方向是对的，但出口是错的——
未发表的本车道 .md 不是可引出处。本脚本做两件事，输出全部落 e154_out.txt：

[X1] 从 papers/notes/1603.04172.txt 逐行打印"时间×空间阈值策略"那几句的**带行号原文**，
     让 C 能自己读那 5 行再引（不必转引我的笔记）。
[X2] 单查 Crossref（一次请求，失败就如实记账），判定这篇能否从 arXiv 预印本升级著录；
     不能升级就写死"预印本"口径。判据在本行之前写死，不看完结果再挑。
"""
import io, json, sys, urllib.request, urllib.parse

TXT = 'papers/notes/1603.04172.txt'
OUT = []


def p(*a):
    s = ' '.join(str(x) for x in a)
    OUT.append(s)
    print(s)


lines = io.open(TXT, encoding='utf-8', errors='replace').read().split('\n')
p('== [X1] 1603.04172 抽取文本，总 %d 行 ==' % len(lines))
RANGES = [(176, 184, '贡献 (R2) 第 (1) 条：参数式由 time-space reverse-waterfilling 刻画'),
          (194, 214, '§I 贡献段：算法与"阈值策略"句'), ]
for lo, hi, why in RANGES:
    p('-- %s（第 %d-%d 行）--' % (why, lo, hi))
    for i in range(lo - 1, min(hi, len(lines))):
        s = lines[i].replace('\x00', '')
        if s.strip():
            p('%5d| %s' % (i + 1, s.rstrip()))

key = []
for i, s in enumerate(lines):
    c = s.replace('\x00', '')
    if 'threshold policy' in c or 'time-space reverse-' in c or 'time-space reverse-w' in c:
        key.append((i + 1, c.strip()))
p('-- [X1b] 关键词全扫（threshold policy / time-space reverse-waterfilling）命中 %d 处 --' % len(key))
for n, c in key:
    p('%5d| %s' % (n, c[:200]))

p('')
p('== [X2] Crossref 单查（判据：命中且 container-title 非空且 type=journal-article => 可升级；否则保持预印本）==')
TITLE = ('Optimal Estimation via Nonanticipative Rate Distortion Function and Applications '
         'to Time-Varying Gauss-Markov Processes')
url = 'https://api.crossref.org/works?' + urllib.parse.urlencode(
    {'query.bibliographic': TITLE, 'rows': 5, 'select': 'DOI,title,container-title,type,author,issued,volume,page'})
got = None
try:
    req = urllib.request.Request(url, headers={'User-Agent': 'research-note/1.0 (mailto:none@example.com)'})
    got = json.loads(urllib.request.urlopen(req, timeout=45).read().decode('utf-8'))
except Exception as e:
    p('[X2-FAIL] 请求失败：%s：%s => 按"单通道预印本"口径，著录不升级' % (type(e).__name__, str(e)[:160]))

if got is not None:
    items = got.get('message', {}).get('items', [])
    p('[X2] 命中 %d 条' % len(items))
    upgrade = False
    for it in items:
        t = (it.get('title') or [''])[0]
        ct = (it.get('container-title') or [''])[0]
        p('  DOI=%s | type=%s | journal=%s | vol=%s | pp=%s | year=%s' % (
            it.get('DOI'), it.get('type'), ct, it.get('volume'), it.get('page'),
            (it.get('issued', {}).get('date-parts') or [['?']])[0][0]))
        p('    title=%s' % t[:120])
        norm = lambda x: ''.join(ch for ch in x.lower() if ch.isalnum())
        if norm(t) == norm(TITLE) or (norm(TITLE) and norm(t) and norm(TITLE) in norm(t)):
            p('    => 标题逐字一致；container-title=%r type=%r' % (ct, it.get('type')))
            if ct and it.get('type') == 'journal-article':
                upgrade = True
    p('[X2 判决] 可升级成期刊著录？ %s' % upgrade)
    if not upgrade:
        p('[X2 结论] 保持 arXiv:1603.04172 预印本口径（updated 2017-02-10，comment 说 submitted to SICON，无 journal_ref/DOI）。')

io.open('p0/e154_out.txt', 'w', encoding='utf-8').write('\n'.join(OUT) + '\n')
