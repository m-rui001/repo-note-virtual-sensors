# -*- coding: utf-8 -*-
# e181：三条最近邻的**著录 + 开放全文**取证（Crossref 与 OpenAlex 双通道；OpenAlex 还带摘要与 OA 链接）。
# 纪律：#46 每个通道打印成功条数；#48 "A 的期刊版是 B"只许用标题等值或全文定理号互涉判同一。
import io
import json
import re
import sys
import time
import urllib.request
import urllib.error

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')
OUT = []
def p(s):
    OUT.append(s)
    print(s, flush=True)

HDR = {'User-Agent': 'rc-lane/1.0 (research note; mailto:m-rui001@example.com)'}
TARGETS = [
    ('TAC-2022', '10.1109/tac.2021.3099444'),
    ('ACC-2012', '10.1109/acc.2012.6314650'),
    ('ACC-2023', '10.23919/acc55779.2023.10156247'),
]


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
            last = type(e).__name__ + ':' + str(e)[:60]
            time.sleep(sleep * (k + 1))
    return None, last


def norm_title(s):
    return re.sub(r'[^a-z0-9]+', ' ', (s or '').lower()).strip()


def cr_query(doi):
    obj, st = get('https://api.crossref.org/works/' + doi)
    if obj is None:
        return None, st
    j = json.loads(obj.decode('utf-8', 'replace'))['message']
    return {
        'title': (j.get('title') or [''])[0],
        'authors': [('%s %s' % (a.get('given') or '', a.get('family') or '')).strip() for a in j.get('author', [])],
        'container': (j.get('container-title') or [''])[0],
        'volume': j.get('volume', ''),
        'issue': j.get('issue', ''),
        'page': j.get('page', ''),
        'year': ((j.get('published-print') or j.get('published-online') or {}).get('date-parts') or [[None]])[0][0],
        'abstract': re.sub(r'<[^>]+>', ' ', j.get('abstract', '') or ''),
        'cited': j.get('is-referenced-by-count', ''),
        'doi': doi,
    }, 'ok'


def oa_query(doi):
    obj, st = get('https://api.openalex.org/works/doi:' + doi)
    if obj is None:
        return None, st
    j = json.loads(obj.decode('utf-8', 'replace'))
    inv = j.get('abstract_inverted_index') or {}
    abstract = ''
    if inv:
        pos = {}
        for w, idxs in inv.items():
            for i in idxs:
                pos[i] = w
        abstract = ' '.join(pos[i] for i in sorted(pos))
    oa = j.get('open_access') or {}
    pdf = ''
    # OpenAlex 新 schema 里 locations 是 list，旧的是 dict；两种都要读，
    # 否则 pdf_url 静默为空会被误判成"无开放全文"。
    locs = j.get('locations')
    src = []
    if isinstance(locs, dict):
        src = locs.get('sources', []) or []
    elif isinstance(locs, list):
        src = locs
    bl = [b for b in src if isinstance(b, dict) and b.get('pdf_url')]
    if bl:
        pdf = bl[0].get('pdf_url', '')
    if not pdf:
        rl = (j.get('best_oa_location') or {}).get('pdf_url')
        pdf = rl or ''
    return {
        'title': j.get('display_name') or j.get('title') or '',
        'authors': [((a.get('author') or {}).get('display_name') or '') for a in (j.get('authorships') or [])],
        'container': ((j.get('primary_location') or {}).get('source') or {}).get('display_name', ''),
        'year': j.get('publication_year'),
        'oa': oa.get('is_oa'),
        'oa_url': oa.get('oa_url') or '',
        'pdf_url': pdf,
        'cited_by': j.get('cited_by_count', ''),
        'work_id': (j.get('id') or '').rsplit('/', 1)[-1],
        'abstract': abstract,
    }, 'ok'


CR, OA = {}, {}
p('== [S0] 双通道取数（Crossref + OpenAlex；每条都记成功/失败）==')
for tag, doi in TARGETS:
    c, sc = cr_query(doi)
    time.sleep(1.5)
    o, so = oa_query(doi)
    time.sleep(1.5)
    CR[tag], OA[tag] = c, o
    p('  %-9s crossref=%-9s openalex=%-9s' % (tag, sc, so))
    p('    CR 标题  : %s' % (c['title'] if c else '（无）'))
    p('    CR 作者  : %s（%d 位）' % ('; '.join(c['authors']) if c else '（无）', len(c['authors']) if c else 0))
    p('    CR 出处  : %s vol=%s iss=%s page=%s year=%s cited=%s' %
      (c['container'] if c else '', c['volume'] if c else '', c['issue'] if c else '',
       c['page'] if c else '', c['year'] if c else '', c['cited'] if c else ''))
    p('    CR 摘要  : %s' % ('有（%d 字）' % len(c['abstract']) if c and c['abstract'].strip() else '无'))
    p('    OA 标题  : %s' % (o['title'] if o else '（无）'))
    p('    OA 作者  : %s（%d 位）' % ('; '.join(o['authors']) if o else '（无）', len(o['authors']) if o else 0))
    p('    OA 出处  : %s year=%s is_oa=%s cited_by=%s' % (o['container'] if o else '', o['year'] if o else '', o['oa'] if o else '', o['cited_by'] if o else ''))
    p('    OA 链接  : oa_url=%s' % (o['oa_url'] if o else ''))
    p('    OA pdf   : %s' % ((o['pdf_url'] or '（无）') if o else '（无）'))
    p('    OA 摘要  : %s' % ('有（%d 字）' % len(o['abstract']) if o and o['abstract'].strip() else '无'))
    p('    标题两通道一致：%s' % (bool(c and o and norm_title(c['title']) == norm_title(o['title']))))

ok_cr = sum(1 for t, _ in TARGETS if CR[t])
ok_oa = sum(1 for t, _ in TARGETS if OA[t])
p('')
p('  通道成功条数：crossref %d/3，openalex %d/3' % (ok_cr, ok_oa))
assert ok_cr == 3 and ok_oa == 3, '双通道没取满，不许出判决'

p('')
p('== [S1] 摘要全文（OpenAlex 重建的 inverted index，逐条打印，供通读比对）==')
for tag, doi in TARGETS:
    ab = (OA[tag]['abstract'] or '').strip()
    p('  --- %s (%s)：' % (tag, doi))
    p('      %s' % (ab if ab else '（OpenAlex 无摘要 ⇒ 只能声称著录，不能声称读过）'))

p('')
p('== [S2] 下载开放全文（只有 OA 链接存在才下；下完立刻转文本并打印首行）==')
import os
os.makedirs('papers/parts', exist_ok=True)
DLED = []
for tag, doi in TARGETS:
    url = (OA[tag]['pdf_url'] or OA[tag]['oa_url'] or '') if OA[tag] else ''
    if not url:
        p('  %-9s 无 OA 链接 ⇒ 不下（不编造全文）' % tag)
        continue
    raw, st = get(url)
    if raw is None:
        p('  %-9s 下载失败 %s（%s）' % (tag, st, url[:70]))
        continue
    fn = 'papers/parts/e181_%s' % tag
    open(fn + ('.pdf' if url.lower().endswith('.pdf') or raw[:4] == b'%PDF' else '.html'), 'wb').write(raw)
    DLED.append((tag, fn, len(raw)))
    p('  %-9s 落盘 %s（%d B，%s）' % (tag, fn + ('.pdf' if raw[:4] == b'%PDF' else '.html'), len(raw), url[:60]))

p('')
p('== [S3] 条数与退出（#46）==')
p('  目标 3 条；crossref 成功 %d；openalex 成功 %d；摘要可得 %d；OA 全文落盘 %d' %
  (ok_cr, ok_oa, sum(1 for t, _ in TARGETS if OA[t] and OA[t]['abstract'].strip()), len(DLED)))
p('  本节**不**判"A 的期刊版是 B"（#48：需标题等值或全文定理号互涉）；只出著录＋可读全文。')

open('p0/e181_out.txt', 'w', encoding='utf-8').write('\r\n'.join(OUT) + '\r\n')
