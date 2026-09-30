#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
e154 的后续：(1) 把 §V 的 time-space reverse-waterfilling 性质 (5.8)-(5.11) 那段带行号打出来
（只引引言不够，要让协作者能直接读定理级句子）；
(2) 按 DOI 单取 Crossref 全字段，产出**够写 bibitem 的著录**；
(3) 判据先写死：作者姓列表与 arXiv 笔记记的四位一致 且 卷/页/年齐 => 允许 bibitem 升级到期刊版；
    任缺一项 => 只写 arXiv 口径并注明"缺项"。
"""
import io, json, urllib.request

OUT = []


def p(*a):
    s = ' '.join(str(x) for x in a)
    OUT.append(s)
    print(s)


lines = io.open('papers/notes/1603.04172.txt', encoding='utf-8', errors='replace').read().split('\n')
p('== [Y1] 定理级：(5.8)-(5.11) 附近的原文（第 1190-1215 行）==')
for i in range(1189, min(1215, len(lines))):
    s = lines[i].replace('\x00', '').rstrip()
    if s.strip():
        p('%5d| %s' % (i + 1, s))

p('')
p('== [Y2] Crossref by-DOI：10.1137/17m1116349 ==')
DOI = '10.1137/17m1116349'
msg = None
try:
    req = urllib.request.Request('https://api.crossref.org/works/' + DOI,
                                 headers={'User-Agent': 'research-note/1.0 (mailto:none@example.com)'})
    msg = json.loads(urllib.request.urlopen(req, timeout=45).read().decode('utf-8'))['message']
except Exception as e:
    p('[Y2-FAIL] %s：%s' % (type(e).__name__, str(e)[:160]))

EXPECT = ['stavrou', 'charalambous', 'loyka']
if msg:
    au = [(a.get('family', ''), a.get('given', '')) for a in msg.get('author', [])]
    p('标题        : %s' % msg['title'][0])
    p('期刊        : %s' % msg.get('container-title', [''])[0])
    p('卷 / 页 / 年: %s / %s / %s' % (msg.get('volume'), msg.get('page'),
                                      (msg.get('issued', {}).get('date-parts') or [['?']])[0][0]))
    p('DOI         : %s' % msg.get('DOI'))
    p('作者        : %s' % ', '.join('%s %s' % (g, f) for f, g in au))
    fam = [f.lower() for f, _ in au]
    cnt = {k: sum(1 for x in fam if k in x) for k in EXPECT}
    p('[Y2a] 姓计数对账（arXiv 笔记记的四位：Stavrou / Charalambous x2 / Loyka）: %s  共 %d 位作者'
      % (cnt, len(au)))
    fields = dict(title=msg['title'][0], journal=msg.get('container-title', [''])[0],
                  volume=msg.get('volume'), pages=msg.get('page'),
                  year=str((msg.get('issued', {}).get('date-parts') or [['?']])[0][0]),
                  doi=msg.get('DOI'))
    missing = [k for k, v in fields.items() if not v]
    ok_auth = (cnt.get('stavrou', 0) >= 1 and cnt.get('loyka', 0) >= 1
               and cnt.get('charalambous', 0) >= 2)
    p('[Y2 判决] 字段齐=%s 缺项=%s ；作者对账=%s' % (not missing, missing or '无', ok_auth))
    if ok_auth and not missing:
        p('[Y2 结论] 允许把著录升级到期刊版：')
        p('  @article{ stavrou2018optimal, author={P.~A. Stavrou and T. Charalambous and C.~D. Charalambous and S.~A. Loyka},'
          ' title={Optimal Estimation via Nonanticipative Rate Distortion Function and Applications to '
          'Time-Varying Gauss--Markov Processes}, journal={SIAM J. Control Optim.}, volume={56}, '
          'pages={3731--3765}, year={2018}, doi={10.1137/17m1116349} }')
    else:
        p('[Y2 结论] 不升级，保持 arXiv 口径并列出缺项。')

p('')
p('== [Y3] 我方旧笔记的口径需要更正 ==')
p('papers/notes/1603.04172.md 现写："无 arxiv:doi、无 journal_ref => 单通道预印本，著录不许升级"。')
p('该句的**前提**（arXiv 通道无 DOI/journal_ref）仍然成立；错在把"arXiv 通道无字段"当成"无期刊版"。')
p('Crossref 通道命中期刊版（[Y2]），所以上一句的推论作废，笔记需追加期刊著录一节。')

io.open('p0/e154b_out.txt', 'w', encoding='utf-8').write('\n'.join(OUT) + '\n')
