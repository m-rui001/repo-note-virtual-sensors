# -*- coding: utf-8 -*-
# e182：preprint/sec_bib.tex 的**逐条著录核对**（Crossref 解 DOI、arXiv API 解编号、会议→期刊升版探测）。
# 纪律：#46 每类都要打印可解析条数；#48 不许用"标题像"判同一文献，只报**标题等值/不等值**；
#       读不到全文就不写"读过证明"。本节只出**书目层**的错与缺，不碰科学判决。
import io
import json
import re
import sys
import time
import urllib.request
import urllib.error
import urllib.parse
import xml.etree.ElementTree as ET

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')
OUT = []


def p(s):
    OUT.append(s)
    print(s, flush=True)


HDR = {'User-Agent': 'rc-lane/1.0 (research note; mailto:m-rui001@example.com)'}


def get(url, tries=3, sleep=4.0):
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


# ---------------- [T0] 解析书目条目 ----------------
BIB = open('preprint/sec_bib.tex', encoding='utf-8').read()
ENTRIES = []
for m in re.finditer(r'(?m)^\\bibitem\{([^}]+)\}\s*(.*?)(?=^\\bibitem\{|^\\end\{thebibliography\})',
                     BIB, re.S):
    key, body = m.group(1), m.group(2)
    body = re.sub(r'\s+', ' ', body).strip()
    q = re.search(r"``(.+?)''", body, re.S)
    dois = re.findall(r'10\.\d{4,9}/[^\s,;]+', body)
    axids = re.findall(r'(?i)arxiv:?\s*(\d{4}\.\d{4,5})', body)
    VEN = r'(Proc\.|Conf\.|Workshop|CDC|ACC|ITW|CISS|ISIT)'
    ENTRIES.append({'key': key, 'body': body, 'title': (q.group(1) if q else ''),
                    'doi': [d.rstrip('.') for d in dois], 'ax': axids,
                    'conf': bool(re.search(VEN, body))})
N_DOI = sum(len(e['doi']) for e in ENTRIES)
N_AX = sum(len(e['ax']) for e in ENTRIES)
NOLOC = [e['key'] for e in ENTRIES if not e['doi'] and not e['ax']]
p('== [T0] 解析 sec_bib.tex ==')
p('  条目 %d；DOI %d 个；arXiv 编号 %d 个；无任何定位符 %d 条 %s'
  % (len(ENTRIES), N_DOI, N_AX, len(NOLOC), NOLOC))
p('  标题抽取成功 %d／%d（无 ``...'' 引对的条目列在下面）'
  % (sum(1 for e in ENTRIES if e['title']), len(ENTRIES)))
for e in ENTRIES:
    if not e['title']:
        p('    无引号标题：%s' % e['key'])
assert len(ENTRIES) == 24, '条目数变了，先查解析器再谈判决'
assert N_DOI == 9 and N_AX == 15, '定位符条数与上一轮不一致'

# ---------------- [T1] DOI -> Crossref ----------------
p('')
p('== [T1] DOI 逐条 Crossref 解析（返回标题／出处／卷期页年，并回填条目原文核对）==')
DOIROWS = []
ok_doi = 0
for e in ENTRIES:
    for d in e['doi']:
        obj, st = get('https://api.crossref.org/works/' + d.lower())
        time.sleep(1.0)
        if obj is None:
            p('  %-24s %s  ->  %s（**解析失败**）' % (e['key'], d, st))
            DOIROWS.append((e['key'], d, None, st))
            continue
        j = json.loads(obj.decode('utf-8', 'replace'))['message']
        rt = (j.get('title') or [''])[0]
        rc = (j.get('container-title') or [''])[0]
        rv, rip, rp = j.get('volume', ''), j.get('issue', ''), j.get('page', '')
        ry = ((j.get('published-print') or j.get('published-online') or {})
              .get('date-parts') or [[None]])[0][0]
        ty = j.get('type', '')
        same = norm_t(rt) == norm_t(e['title'])
        year_in = bool(ry) and (str(ry) in e['body'])
        page_in = bool(rp) and (rp in e['body'].replace('--', '-'))
        p('  %-24s %s' % (e['key'], d))
        p('      CR 标题：%s' % rt)
        p('      标题等值：%s   我的标题：%s' % (same, e['title']))
        p('      CR 出处：%s [type=%s] vol=%s iss=%s page=%s year=%s'
          % (rc, ty, rv or '-', rip or '-', rp or '-', ry))
        p('      条目原文含该年：%s；含该页码区间：%s' % (year_in, page_in))
        DOIROWS.append((e['key'], d, j, 'ok'))
        ok_doi += 1
p('  DOI 可解析 %d／%d' % (ok_doi, N_DOI))
assert ok_doi > 0

# ---------------- [T2] arXiv 批量 ----------------
p('')
p('== [T2] arXiv 编号批量解析（一次 id_list；返回条数必须对上，未返回的显式列）==')
ALLAX = []
for e in ENTRIES:
    for a in e['ax']:
        if a not in ALLAX:
            ALLAX.append(a)
obj, st = get('http://export.arxiv.org/api/query?id_list=%s&max_results=%d'
              % (','.join(ALLAX), len(ALLAX) + 5), tries=4, sleep=8.0)
GOT = {}
if obj is None:
    p('  arXiv 批量请求失败：%s ⇒ 本节所有编号判"未核"' % st)
else:
    root = ET.fromstring(obj.decode('utf-8', 'replace'))
    ns = {'a': 'http://www.w3.org/2005/Atom'}
    for en in root.findall('a:entry', ns):
        eid = (en.findtext('a:id', default='', namespaces=ns) or '')
        base = re.sub(r'v\d+$', '', eid.rsplit('/', 1)[-1])
        ttl = re.sub(r'\s+', ' ', (en.findtext('a:title', default='', namespaces=ns) or '')).strip()
        pub = (en.findtext('a:published', default='', namespaces=ns) or '')[:4]
        GOT[base] = (ttl, pub)
p('  请求 %d 个编号，返回 %d 条记录' % (len(ALLAX), len(GOT)))
bad_ax, ok_ax = [], 0
for e in ENTRIES:
    for a in e['ax']:
        if a not in GOT:
            p('  %-24s arXiv:%s  ->  **无记录**（编号不存在或格式不符）' % (e['key'], a))
            bad_ax.append((e['key'], a))
            continue
        ttl, pub = GOT[a]
        same = norm_t(ttl) == norm_t(e['title'])
        p('  %-24s arXiv:%s 年=%s 标题等值=%s' % (e['key'], a, pub, same))
        p('      arXiv 标题：%s' % ttl)
        ok_ax += 1
p('  arXiv 可解析 %d／%d；**不可解析 %d**' % (ok_ax, N_AX, len(bad_ax)))
for k, a in bad_ax:
    p('    待处理：%s 的 arXiv:%s' % (k, a))

# ---------------- [T3] 会议条目 -> 是否存在同名期刊版 ----------------
p('')
p('== [T3] 会议／workshop 条目的"期刊升版"探测（Crossref query.bibliographic；只报标题等值项）==')
UPG = []
for e in ENTRIES:
    if not e['conf'] or not e['title']:
        continue
    q = urllib.parse.quote(e['title'][:120])
    obj, st = get('https://api.crossref.org/works?query.bibliographic=%s&rows=8'
                  '&select=title,container-title,issued,DOI,volume,page,type' % q)
    time.sleep(1.2)
    p('  --- %s（会议条目，标题：%s）' % (e['key'], e['title'][:70]))
    if obj is None:
        p('      查询失败 %s ⇒ 不判' % st)
        continue
    items = json.loads(obj.decode('utf-8', 'replace'))['message'].get('items', [])
    hits = []
    for it in items:
        t = (it.get('title') or [''])[0]
        if norm_t(t) != norm_t(e['title']):
            continue
        c = (it.get('container-title') or [''])[0]
        y = ((it.get('issued') or {}).get('date-parts') or [[None]])[0][0]
        if it.get('DOI', '').lower() == (e['doi'][0].lower() if e['doi'] else 'x'):
            continue  # 就是它自己
        hits.append((c, y, it.get('DOI'), it.get('volume', ''), it.get('page', '')))
    if hits:
        for c, y, d, v, pg in hits:
            p('      同名异出处：%s %s vol=%s pp=%s doi=%s' % (c, y, v or '-', pg or '-', d))
        UPG.append((e['key'], hits))
    else:
        p('      无标题等值的异出处记录 ⇒ 本条不提级（**不等于**不存在；只说本次没查到）')
p('  探测到候选升版 %d 条' % len(UPG))

# ---------------- [T4] 判决表 ----------------
p('')
p('== [T4] 判决（只书目层；每条都带条数，#46）==')
tmis = [(k, d) for k, d, j, s in DOIROWS if j and norm_t((j.get('title') or [''])[0]) != norm_t(
    next(e['title'] for e in ENTRIES if e['key'] == k))]
unres = [(k, d) for k, d, j, s in DOIROWS if j is None]
p('  条目 24；DOI 解析 %d／%d；arXiv 解析 %d／%d；编号不可解析 %d 个；'
  '无定位符 %d 条；标题不等值 %d 条；候选升版 %d 条'
  % (ok_doi, N_DOI, ok_ax, N_AX, len(bad_ax), len(NOLOC), len(tmis), len(UPG)))
for k, d in unres:
    p('    DOI 不可解析：%s %s' % (k, d))
for k, d in tmis:
    p('    标题不等值：%s %s' % (k, d))
p('  本节**只**核对著录元数据；三篇最近邻（e181）无开放全文 ⇒ 不得写"读过其证明"。')

open('p0/e182_out.txt', 'w', encoding='utf-8').write('\r\n'.join(OUT) + '\r\n')
