#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
e157b：换通道重做 e157 的"无先例"筛查（Crossref，今天 10:04/11:04 两次都成功，未遇 429）。
**本轮新门禁（写死，先于结果）**：
  [G] 成功返回的查询数 < 5（共 6 项）=> 判决行只许写"检索未完成，不下任何结论"，
      不许写"未见同现"。理由：e157 的第一版正是 4/5 失败仍打出"未见三项同现"，那是仪器无效读数（板账 #46）。
  三问标同 e157：a=传感/观测矩阵是决策变量；b=控制代价/加权失真；c=随传感秩的地板。
"""
import json, re, sys, time, urllib.parse, urllib.request
sys.stdout.reconfigure(encoding='utf-8')

OUT = []


def p(*a):
    s = ' '.join(str(x) for x in a)
    OUT.append(s)
    print(s)


def get(url, tries=4):
    d = 6.0
    for k in range(tries):
        try:
            req = urllib.request.Request(url, headers={'User-Agent': 'rc-lane-note/1.0 (mailto:none@example.com)'})
            return json.loads(urllib.request.urlopen(req, timeout=60).read().decode('utf-8'))
        except Exception as e:
            p('    !! %s：%s => 退避 %.0fs' % (type(e).__name__, str(e)[:120], d))
            if k == tries - 1:
                return None
            time.sleep(d)
            d *= 2
    return None


A = re.compile(r'sensor|sensing|observation|measurement|probe', re.I)
B = re.compile(r'linear[- ]quadratic|LQG|control cost|cost of control|certainty equivalence|weighted trace|task[- ]driven|cost function|control', re.I)
C = re.compile(r'rank|dimension|subspace|partial[- ]observation|low[- ]dimension', re.I)


def flag(t, ab):
    s = t + ' ' + (ab or '')
    return ('a' if A.search(s) else '-') + ('b' if B.search(s) else '-') + ('c' if C.search(s) else '-')


QS = ['reverse waterfilling sensor',
      'rate-distortion sensor design',
      'task-oriented sensing rate distortion',
      'limited sensing linear quadratic gaussian rate distortion',
      'rank constrained sensing control cost information',
      'sensor selection linear quadratic gaussian distortion']
ok_q = 0
seen = {}
abc = []
p('== [S1] Crossref 六组筛查（select 含 abstract；摘要缺失则只按标题标）==')
for q in QS:
    url = 'https://api.crossref.org/works?' + urllib.parse.urlencode(
        {'query.bibliographic': q, 'rows': 8,
         'select': 'DOI,title,container-title,type,issued,abstract'})
    p('-- %s' % q)
    j = get(url)
    if not j:
        p('   [FAIL] 本组无读数')
        continue
    ok_q += 1
    for it in j.get('message', {}).get('items', []):
        t = re.sub(r'\s+', ' ', (it.get('title') or ['(无标题)'])[0]).strip()
        ab = re.sub(r'\s+', ' ', re.sub(r'<[^>]+>', '', it.get('abstract') or '')).strip()
        doi = it.get('DOI', '')
        ct = (it.get('container-title') or [''])
        ct = ct[0] if ct else ''
        yr = str(((it.get('issued') or {}).get('date-parts') or [['?']])[0][0])
        fl = flag(t, ab)
        if doi in seen:
            p('   [%s] (重复，略) %s' % (fl, doi))
            continue
        seen[doi] = 1
        p('   [%s] %s | %s %s | %s' % (fl, t[:88], ct[:34], yr, doi))
        if ab:
            p('        摘要片段：%s' % ab[:220])
        else:
            p('        摘要：无（Crossref 未收录，标记只按标题）')
        if fl == 'abc':
            abc.append((t, doi, ct, yr, ab[:260]))
    time.sleep(2.0)

p('')
p('== [S2] 单取那条一直 429 的 ITW-2017（按 DOI）==')
for DOI in ('10.1109/itw.2017.8277966', '10.1109/cdc.2018.8619725'):
    j = get('https://api.crossref.org/works/' + DOI)
    if not j:
        p('   [FAIL] %s 未取到' % DOI)
        continue
    m = j['message']
    t = (m.get('title') or [''])[0]
    au = ', '.join((a.get('given', '') + ' ' + a.get('family', '')).strip() for a in m.get('author', []))
    ab = re.sub(r'\s+', ' ', re.sub(r'<[^>]+>', '', m.get('abstract') or '')).strip()
    p('   DOI=%s | %s | %s %s' % (m.get('DOI'), t[:100], (m.get('container-title') or [''])[0][:44],
                                  ((m.get('issued', {}).get('date-parts') or [['?']])[0][0])))
    p('        作者：%s' % au[:200])
    p('        摘要：%s' % (ab[:400] if ab else 'Crossref 未收录'))
    time.sleep(2.0)

p('')
p('== [判决行]（先于结果写死：成功查询数 >=5 才允许下"未见同现"的结论）==')
uniq = len(seen)
p('[G] 成功查询 %d/6；去重后条目 %d 条；三项同现候选 %d 条' % (ok_q, uniq, len(abc)))
if ok_q < 5:
    p('[结论] **检索未完成 => 不下任何结论**（既不说"无先例"，也不说"未见"）。')
else:
    if abc:
        p('[结论] 有撞车候选 => "无先例"作废，须改写为与下列条目的具体差别：')
        for t, doi, ct, yr, ab in abc:
            p('   !!! %s | %s %s | %s' % (t[:90], ct[:30], yr, doi))
            p('       %s' % ab)
    else:
        p('[结论] %d 条命中里**无一项同时含 (a)(b)(c) 三个关键词族** => ' % uniq +
          '仍只允许写"关键词筛查未见同现"，不许写"无人处理"（标题级标记的条目没摘要，见上表"摘要：无"）。')

open('p0/e157b_out.txt', 'w', encoding='utf-8').write('\n'.join(OUT) + '\n')
