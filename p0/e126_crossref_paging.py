r"""E126 = 结清 §65-G/§66-G 的 ①：Crossref 翻页 + 精筛 35 条定档。

跑前写死的判据（跑后不许改）：
 [M1] 先实测 Crossref `query.bibliographic` 的 total-results 量级。它是相关性排序的模糊检索，
      total 是**语料库级**而非命中集级 ⇒ "第二通道全集"这个说法作废，改写成可执行的：
      按相关性翻页到**饱和深度**（连续 2 页不新增相关条目即停），将 D 作为有记录的覆盖口径。
 [M2] ≥3 组互为同义扩展的查询（实际 8 组 bibliographic + 6 组 title），每组落盘：
      查询串 / total-results / 抓取数 / 相关数 / 停止原因（纪律 26）。
 [M3] 真正**有界**的语料用第二条通道：**引文邻域** —— 承重论文的 Crossref `reference` 字段
      并集。它是有限集，所以"在这个集合里没人做过 X"这种句子只有在这里才允许写。
 [M4] 精筛 35 条逐条定档，代码只有四种：
      [R]     本地有全文笔记 ⇒ 已读
      [A-IR]  摘要含控制轴但缺"传感约束词 或 反证/对偶词"之一 ⇒ 摘要级判无关（只用于覆盖统计）
      [A-CH]  摘要**同时**含传感约束词与反证/对偶词 ⇒ 威胁项，必须读全文才能进任何新意陈述
      [U]     其余 ⇒ 未读、不作前占判断（纪律 24③）
      不许由 [A-IR]/[U] 外推成"没人做过"。
"""
import sys, os, re, json, time, urllib.parse, urllib.request
sys.stdout.reconfigure(encoding='utf-8')

UA = {'User-Agent': 'prior-art-census/1.0 (mailto:p0-lane@example.org)'}


def gj(url, tries=4):
    last = None
    for a in range(tries):
        try:
            req = urllib.request.Request(url, headers=UA)
            with urllib.request.urlopen(req, timeout=90) as f:
                d = json.loads(f.read().decode('utf-8', 'replace'))
            return d
        except Exception as e:
            last = e
            time.sleep(5.0)
    raise RuntimeError(f'failed {tries}x: {type(e).__name__} {str(last)[:90]}')


def norm(s):
    return re.sub(r'[^a-z0-9]+', ' ', (s or '').lower()).strip()


CTRL = ['lqg', 'lqr', 'linear quadratic', 'quadratic cost', 'control cost', 'controller',
        'control', 'stabiliz', 'estimation', 'kalman', 'lateral control']
COMM = ['directed information', 'directed mutual information', 'causal rate', 'rate distortion',
        'rate-distortion', 'data rate', 'communication rate', 'mutual information', 'entropy',
        'quantiz', 'coding', 'bandwidth', 'communication']
SENSC = ['sensor', 'measurement', 'observation', 'actuator', 'selection', 'placement', 'scheduling',
         'partial', 'limited', 'task driven', 'task based', 'subspace', 'rank']
DUALC = ['converse', 'lower bound', 'bound', 'dual', 'lagrang', 'semidefinite', 'sdp', 'exact',
         'optimal', 'achievable', 'approximation ratio', 'np hard', 'inapproximab']

print('== (0) [M1] 这条通道能不能承载"全集"这个词：实测 total-results 量级 ==')
probe = [('directed information LQG', 'directed+information+LQG'),
         ('data rate theorem linear systems', 'data+rate+theorem+linear+systems'),
         ('sensor selection control', 'sensor+selection+control')]
for lab, q in probe:
    d = gj(f'https://api.crossref.org/works?query.bibliographic={q}&rows=1')
    print(f'  {lab:<34} total-results = {d["message"]["total-results"]:,}  （语料库 1.6e8 量级）')

print('\n== (1) [M2] 8 组同义扩展：按相关性翻页到饱和深度 ==')
BX = [
    ('X1', 'directed information LQG control', 'directed+information+lqg+control'),
    ('X2', 'directed mutual information linear quadratic', 'directed+mutual+information+linear+quadratic'),
    ('X3', 'causal rate distortion control', 'causal+rate+distortion+control'),
    ('X4', 'sensor selection LQR control cost', 'sensor+selection+lqr+control+cost'),
    ('X5', 'task based sensing rate distortion control', 'task+based+sensing+rate+distortion+control'),
    ('X6', 'data rate theorem communication rate stabilization', 'data+rate+theorem+communication+rate+stabilization'),
    ('X7', 'sensor scheduling LQG converse bound', 'sensor+scheduling+lqg+converse+bound'),
    ('X8', 'actuator sensor placement control cost semidefinite', 'actuator+sensor+placement+control+cost+semidefinite'),
]
PAGE, MAXOFF = 1000, 4000


def relevant(title):
    t = ' ' + norm(title) + ' '
    return any(w in t for w in CTRL) and any(w in t for w in COMM)


union, log = {}, []
for tag, lab, q in BX:
    off, newrel, stop, tot = 0, 0, 'max-offset', None
    deadpages = 0
    while off <= MAXOFF:
        u = (f'https://api.crossref.org/works?query.bibliographic={q}'
             f'&rows={PAGE}&offset={off}&select=DOI,title,container-title,issued,is-referenced-by-count')
        d = gj(u)
        m = d['message']
        if tot is None:
            tot = m['total-results']
        items = m['items']
        hit = 0
        for it in items:
            ti = (it.get('title') or [''])[0]
            doi = (it.get('DOI') or '').lower()
            cur = union.setdefault(doi, dict(title=ti, cont=(it.get('container-title') or [''])[0],
                                             yr=(it.get('issued', {}).get('date-parts') or [[None]])[0][0],
                                             cited=it.get('is-referenced-by-count'), frm=[], rel=False))
            if tag not in cur['frm']:
                cur['frm'].append(tag)
            if relevant(ti):
                cur['rel'] = True
                hit += 1
        newrel += hit
        deadpages = deadpages + 1 if hit == 0 else 0
        print(f'   [{tag}] offset={off:<5} 抓到={len(items):<5} 相关={hit:<4} 累计相关={newrel:<5} 死页={deadpages}', flush=True)
        if deadpages >= 2:
            stop = 'saturated(连续2页0相关)'
            break
        if len(items) < PAGE:
            stop = 'end-of-list'
            break
        off += PAGE
        time.sleep(1.2)
    log.append((tag, q, tot, min(off + PAGE, MAXOFF + PAGE), newrel, stop))

print('\n  逐查询串台账（纪律 26：查询串 + 原始命中数一起落盘）')
print(f'  {"tag":<4}{"total-results":>16}{"翻到":>8}{"相关":>7}  停止原因 / query')
for tag, q, tot, depth, nr, stop in log:
    print(f'  {tag:<4}{tot:>16,}{depth:>8,}{nr:>7,}  {stop:<26} {q}')
rel = {k: v for k, v in union.items() if v['rel']}
print(f'\n [M2] 抓取并集 {len(union):,} 条 DOI → 相关性过滤 {len(rel):,} 条；'
      f'覆盖口径写死为"每组查询前 {MAXOFF+PAGE} 条相关性排序位"，不是全集')

print('\n== (2) [M3] 有界通道：承重论文的引文邻域（这才是能写"在这个集合里"的地方）==')
ANCH = ['1510.04214', '1612.02126', '2109.12246', '2003.11951', '1811.11792', '2606.31396']


def arxiv_title(aid):
    import xml.etree.ElementTree as ET
    NS = {'a': 'http://www.w3.org/2005/Atom'}
    u = 'https://export.arxiv.org/api/query?id_list=' + aid
    req = urllib.request.Request(u, headers=UA)
    with urllib.request.urlopen(req, timeout=90) as f:
        root = ET.fromstring(f.read().decode('utf-8', 'replace'))
    return ' '.join(root.find('.//a:entry/a:title', NS).text.split())


cited, anchors = set(), []
for aid in ANCH:
    try:
        ti = arxiv_title(aid)
    except Exception as e:
        print(f'  {aid} arXiv 取标题失败 {type(e).__name__}')
        continue
    time.sleep(2.5)
    d = gj('https://api.crossref.org/works?query.title=' + urllib.parse.quote_plus(ti) + '&rows=3&select=DOI,title')
    cand = [(it['DOI'], norm((it.get('title') or [''])[0])) for it in d['message']['items']]
    doi = next((x for x, n in cand if norm(ti)[:24] in n or n[:24] in norm(ti)[:24]), None)
    print(f'\n  [{aid}] "{ti[:64]}" → DOI={doi}')
    anchors.append((aid, ti, doi))
    if not doi:
        continue
    time.sleep(1.5)
    try:
        w = gj('https://api.crossref.org/works/' + urllib.parse.quote(doi))
    except Exception as e:
        print(f'    works/ 失败 {e}')
        continue
    refs = w['message'].get('reference') or []
    n_unres = sum(1 for r in refs if not r.get('DOI'))
    print(f'    reference 字段 {len(refs)} 条（DOI 未解析 {n_unres}）；is-referenced-by-count='
          f'{w["message"].get("is-referenced-by-count")}')
    for r in refs:
        rt = (r.get('unstructured') or '') + ' ' + (r.get('article-title') or '')
        if relevant(rt) or any(k in norm(rt) for k in ['directed information', 'causal rate',
                                                       'rate distortion', 'sensor selection']):
            cited.add((aid, norm(rt)[:120]))
    time.sleep(1.5)

print(f'\n  引文邻域里"信息轴∩控制轴"字样的被引条目（去重后）：{len(cited)} 条')
for aid, rt in sorted(cited):
    print(f'    <-{aid}  {rt[:104]}')

print('\n== (3) [M4] 精筛 35 条逐条定档 ==')
raw = json.load(open('p0/e97_census_raw.json', encoding='utf-8'))
union97, pids = raw['union'], raw['precise']
NOTES = {m.group(1) for p in os.listdir('papers/notes')
         if (m := re.match(r'(\d{4}\.\d{4,5})', p))}
PDFS = {m.group(1) for p in os.listdir('papers') if (m := re.match(r'(\d{4}\.\d{4,5})', p))}
rows = []
for i in pids:
    h = union97[i]
    t = norm(h['title'] + ' ' + h['abstract'])
    sens = [w for w in SENSC if w in t]
    dual = [w for w in DUALC if w in t]
    if i in NOTES:
        st = '[R]'
    elif sens and dual:
        st = '[A-CH]'
    elif (sens or dual):
        st = '[A-IR]'
    else:
        st = '[U]'
    rows.append((st, i, h['pub'][:10], bool(sens), bool(dual), h['title'][:66], ','.join((sens + dual)[:5])))
order = {'[A-CH]': 0, '[R]': 1, '[A-IR]': 2, '[U]': 3}
for st, i, pub, s, du, ti, kw in sorted(rows, key=lambda r: (order[r[0]], r[1])):
    print(f' {st:<7}{i:<12}{pub}{"传" if s else " "}"{"双" if du else " "}  {ti:<66} {kw[:34]}')
from collections import Counter
cnt = Counter(r[0] for r in rows)
print(f'\n  定档统计：{dict(cnt)}  合计 {len(rows)}')
print(f'  [A-CH] 威胁清单：{[r[1] for r in sorted(rows) if r[0]=="[A-CH]"]}')
print(f'  本地全文覆盖：{cnt.get("[R]",0)}/{len(rows)}；旧说法"3/35"复核：'
      f'{[i for i in pids if i in NOTES]}')
ch = [r[1] for r in sorted(rows) if r[0] == '[A-CH]']
print(f'\n  [M4] 结论口径：{cnt.get("[U]",0)}+[A-IR] 条只作覆盖统计，'
      f'任何"无已发表对偶界"式句子必须先读完 {len(ch)} 条 [A-CH]')

json.dump(dict(crossref_log=log, crossref_rel_total=len(rel),
               crossref_rel={k: union[k] for k in sorted(rel)},
               anchors=anchors, cited_n=len(cited), status=rows),
          open('p0/e126_out.json', 'w', encoding='utf-8'), ensure_ascii=False, indent=1)
print('\n落盘 p0/e126_out.json')
