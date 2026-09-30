# -*- coding: utf-8 -*-
# e182b：补齐"条目原文里没有 DOI／卷／期／页"的那几条（Crossref 标题等值 + 作者重叠双门槛）。
# #48：升版常常**改标题**，所以标题等值检索有假阴性；本节的"没查到"一律写成"本次未命中"。
# 双门槛：标题等值 **且** 首位作者姓氏出现在我条目里 ⇒ 才允许写"可补 DOI"。
import io
import json
import re
import sys
import time
import urllib.parse
import urllib.request
import urllib.error

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')
OUT = []


def p(s):
    OUT.append(s)
    print(s, flush=True)


HDR = {'User-Agent': 'rc-lane/1.0 (research note; mailto:m-rui001@example.com)'}


def get(url, tries=3, sleep=5.0):
    last = ''
    for k in range(tries):
        try:
            req = urllib.request.Request(url, headers=HDR)
            with urllib.request.urlopen(req, timeout=60) as r:
                return r.read(), 'ok'
        except urllib.error.HTTPError as e:
            last = 'HTTP %d' % e.code
            if e.code in (429, 503):
                time.sleep(sleep * (k + 1))
                continue
            return None, last
        except Exception as e:  # noqa: BLE001
            last = type(e).__name__ + ':' + str(e)[:70]
            time.sleep(sleep * (k + 1))
    return None, last


def norm_t(s):
    return re.sub(r'[^a-z0-9]+', ' ', (s or '').lower()).strip()


BIB = open('preprint/sec_bib.tex', encoding='utf-8').read()
ENTRIES = {}
for m in re.finditer(r'(?m)^\\bibitem\{([^}]+)\}\s*(.*?)(?=^\\bibitem\{|^\\end\{thebibliography\})',
                     BIB, re.S):
    body = re.sub(r'\s+', ' ', m.group(2)).strip()
    q = re.search(r"``(.+?)''", body, re.S)
    ENTRIES[m.group(1)] = {'body': body,
                           'title': (q.group(1).rstrip(',') if q else ''),
                           'doi': re.findall(r'10\.\d{4,9}/[^\s,;]+', body)}

# 待补：原文没有 DOI，或明显缺卷／期／页的期刊条目
TODO = ['tanaka2018di', 'tanaka2017srd', 'tanaka2016prefix', 'tzoumas2020co',
        'nair2004stabilizability', 'fox2016a', 'fox2016b', 'cuvelier2021side',
        'sabag2021reducing', 'stavrou2017revisited', 'stavrou2019indirect',
        'kostina2019rate']
# 手工升级项：e181 双通道已核的 TAC 期刊版（改过标题，标题等值检索命中不了）
MANUAL = ('stavrou2018asymptotic',
          'Asymptotic Reverse Waterfilling Algorithm of NRDF for Certain Classes of '
          'Vector Gauss-Markov Processes',
          '10.1109/tac.2021.3099444')

p('== [U0] 待补条数 ==')
p('  待核 %d 条；其中原文无 DOI 的 %d 条' % (len(TODO),
   sum(1 for k in TODO if not ENTRIES[k]['doi'])))
assert all(k in ENTRIES for k in TODO), '有条目名对不上 sec_bib.tex'

p('')
p('== [U1] Crossref 标题等值 + 作者重叠双门槛 ==')
FILL = []
miss = []
for k in TODO:
    e = ENTRIES[k]
    q = urllib.parse.quote(e['title'][:120])
    url = ('https://api.crossref.org/works?query.bibliographic=%s&rows=6'
           '&select=title,author,container-title,volume,issue,page,issued,DOI,type' % q)
    obj, st = get(url)
    time.sleep(1.2)
    p('  --- %-24s %s' % (k, e['title'][:66]))
    seg0 = e['body'].split('``')[0]
    p('      作者段=%s' % seg0[:80])
    if obj is None:
        p('      查询失败 %s ⇒ 不判' % st)
        miss.append(k)
        continue
    items = json.loads(obj.decode('utf-8', 'replace'))['message'].get('items', [])
    # 姓氏抽取：作者块是 "R.~Fox and N.~Tishby," 这种形式，只按逗号切会漏掉第一个作者，
    # 所以先取 `` 之前的作者段，再按 ',' 与 ' and ' 分段、各取末尾词。
    seg = e['body'].split('``')[0]
    mine_fam = set()
    for piece in re.split(r',|\band\b|et al\.', seg):
        toks = re.findall(r"[A-Za-z\'\u00c0-\u024f]+", piece.replace('~', ' '))
        if toks:
            mine_fam.add(toks[-1].lower())
    shown = 0
    for it in items:
        t = (it.get('title') or [''])[0]
        if norm_t(t) != norm_t(e['title']):
            continue
        auths = it.get('author') or []
        fams = {(a.get('family') or '').lower() for a in auths if a.get('family')}
        first = (auths[0].get('family') or '').lower() if auths else ''
        overlap = len(fams & mine_fam)
        c = (it.get('container-title') or [''])[0]
        y = ((it.get('issued') or {}).get('date-parts') or [[None]])[0][0]
        d = it.get('DOI', '')
        gate = (first in mine_fam) and overlap >= 1 and d.lower() not in \
            [x.lower() for x in e['doi']]
        shown += 1
        p('      标题等值：%s [type=%s] %s %s vol=%s iss=%s pp=%s' %
          (d, it.get('type', ''), c, y, it.get('volume', '') or '-',
           it.get('issue', '') or '-', it.get('page', '') or '-'))
        p('        首作者=%s 姓氏交集=%d/%d ⇒ %s' %
          (first or '?', overlap, len(fams), '双门槛过（可补）' if gate else '门槛不过'))
        if gate:
            FILL.append((k, c, y, it.get('volume', ''), it.get('issue', ''),
                         it.get('page', ''), d))
    if shown == 0:
        p('      本次未命中标题等值记录（**不等于**不存在；升版会改标题，#48）')
        miss.append(k)
p('  可补条目 %d；未命中／失败 %d（未命中清单：%s）' % (len(FILL), len(miss), miss))

p('')
p('== [U2] 按 key 汇总（同一 DOI 只留第一次命中）==')
seen = set()
for k, c, y, v, i, pg, d in FILL:
    if d.lower() in seen:
        continue
    seen.add(d.lower())
    p('  %-24s 可补 doi: %s  （%s, %s, vol=%s iss=%s pp=%s）' % (k, d, c, y, v or '-', i or '-', pg or '-'))
p('  去重后可补 DOI %d 个' % len(seen))

p('')
p('== [U3] 手工升级项（e181 双通道核过的改标题期刊版）==')
obj, st = get('https://api.crossref.org/works/' + MANUAL[2])
time.sleep(1.0)
assert obj is not None, '手工项通道失败，不许写升级建议'
j = json.loads(obj.decode('utf-8', 'replace'))['message']
p('  %s 的 CDC-2018 版已被期刊版覆盖：' % MANUAL[0])
p('      期刊标题：%s' % (j.get('title') or [''])[0])
p('      出处：%s vol=%s iss=%s pp=%s year=%s doi=%s' %
  ((j.get('container-title') or [''])[0], j.get('volume', ''), j.get('issue', ''),
   j.get('page', ''), ((j.get('published-print') or j.get('published-online') or {})
                       .get('date-parts') or [[None]])[0][0], MANUAL[2]))
p('      我的标题=%s' % ENTRIES[MANUAL[0]]['title'][:60])
p('      两标题等值=%s ⇒ 本条是 #48 说的**改标题升版**，标题检索必然假阴性' %
  (norm_t((j.get('title') or [''])[0]) == norm_t(ENTRIES[MANUAL[0]]['title'])))
p('      作者交集：CDC 版 Stavrou/Charalambous/Charalambous/Loyka/Skoglund；'
  '期刊版 Stavrou/Skoglund ⇒ 同一组的收紧结果（按两通道著录＋作者列表核过）')

p('')
p('== [U4] 条数与退出（#46）==')
p('  待核 %d；可补 DOI %d（去重后 %d）；未命中 %d；手工升级项 1' %
  (len(TODO), len(FILL), len(seen), len(miss)))
p('  本节只出**著录补全**建议，不改任何数字或科学判决。')

open('p0/e182b_out.txt', 'w', encoding='utf-8').write('\r\n'.join(OUT) + '\r\n')
