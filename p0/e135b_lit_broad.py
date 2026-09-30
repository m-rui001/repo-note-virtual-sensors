# -*- coding: utf-8 -*-
r"""实验 135b（接 e135）：**把 e135 那条"检索"从 8 组扩到 20 组，并明写它够不够格叫"无前占"**。

e135 第 1 版的账：8 组 `abs:` 严格短语查询，取回 37 条，其中 3 组返回 **0 条** ⇒ 这个覆盖面**不足以**
支撑正文里任何"没人写过"的句子。本脚本按同一判据扩面，并且**先声明哪条通道有判决权**：

 [P1] 只有**能拿到摘要**的通道才允许做 [A-CH]/[A-IR] 分类（arXiv 有 summary；Crossref 的 abstract 字段
      在 IEEE 记录上普遍缺失）。Crossref 在此只承担两件事：DOI 权威著录、以及**标题级**候选登记。
 [P2] 标题级命中两腿的 Crossref 条目一律标 [T-CH]，并明写：**[T-CH] 不许当"已读"、也不许当"无关"**（纪律 24③），
      它只是"待取全文"清单。
 [P3] 三档词表与 e135 逐字相同（$S$ 传感/向量腿、$L$ 律/闭式腿），判据也相同：两腿同时命中才算 [A-CH]。
 [P4] 判决口径：$\rho_{\rm cov}=$（有摘要且命中两腿的**去重**条数）。$0$ ⇒ 正文可写"在 X 组查询、Y 条去重记录里没有前占"，
      并**必须**带上 X/Y 两个数。任何 [T-CH] 条目单独列为欠账，不许并入"已排除"。
"""
import sys, os, re, json, time, glob, urllib.parse, urllib.request
import xml.etree.ElementTree as ET
sys.stdout.reconfigure(encoding='utf-8')

UA = {'User-Agent': 'prior-art-check/1.0 (mailto:p0-lane@example.org)'}
NS = {'a': 'http://www.w3.org/2005/Atom'}


def get(url, tries=3):
    last = None
    for _ in range(tries):
        try:
            req = urllib.request.Request(url, headers=UA)
            with urllib.request.urlopen(req, timeout=90) as f:
                return f.read().decode('utf-8', 'replace')
        except Exception as e:
            last = e
            time.sleep(3.0)
    print('    !! 取不到 %s: %s %s' % (url[:80], type(last).__name__, str(last)[:70]))
    return None


def norm(s):
    return re.sub(r'[^a-z0-9]+', ' ', (s or '').lower())


S = ['multidimensional', 'multi-dimensional', 'vector', 'sensing rank', 'rank', 'sensor selection',
     'measurement matrix', 'sensing matrix', 'partial observation', 'subspace', 'multivariate',
     'sensor capacity', 'multiple sensor', 'measurement dimension']
L = ['exponent', 'closed form', 'closed-form', 'rate distortion', 'rate-distortion', 'data rate',
     'capacity', 'logarithmic', 'exponential', 'rate-cost', 'trade-off', 'tradeoff', 'slope',
     'characterization', 'exact']

local = set()
for p in glob.glob('papers/notes/*.md'):
    m = re.search(r'(\d{4}\.\d{4,5})', os.path.basename(p))
    if m:
        local.add(m.group(1))

AQ = ['all:"rate distortion" AND all:"control cost"',
      'all:"data rate" AND all:"LQG"',
      'all:"communication rate" AND all:"estimation error"',
      'all:"sensor selection" AND all:"information rate"',
      'all:"measurement rank" AND all:"control cost"',
      'all:"Kalman filtering" AND all:"bandwidth constraint"',
      'all:"remote state estimation" AND all:"rate distortion"',
      'all:"minimal data rate" AND all:"stabilization"',
      'all:"directed information" AND all:"linear quadratic"',
      'all:"cost of communication" AND all:"linear systems"',
      'all:"dimension" AND all:"rate-distortion" AND all:"control"',
      'all:"sensor scheduling" AND all:"estimation" AND all:"rate"',
      'all:"Gaussian source" AND all:"remote estimation" AND all:"multiple"']
CQ = ['sensor+selection+control+cost+rate', 'measurement+rank+kalman+rate+distortion',
      'data+rate+LQG+cost+tradeoff', 'communication+rate+control+closed+form',
      'remote+estimation+rate+distortion+vector']

seen, CH, TCH = {}, [], []
print('== [P1][P3] arXiv %d 组（全字段）==' % len(AQ))
for q in AQ:
    hits = []
    for st in (0, 100, 200):
        u = ('https://export.arxiv.org/api/query?search_query=' + urllib.parse.quote_plus(q)
             + '&start=%d&max_results=100&sortBy=relevance' % st)
        s = get(u)
        if s is None:
            continue
        try:
            root = ET.fromstring(s)
        except Exception:
            continue
        for e in root.findall('a:entry', NS):
            aid = (e.findtext('a:id', '', NS) or '').split('/abs/')[-1]
            hits.append((aid, (e.findtext('a:published', '', NS) or '')[:4],
                         re.sub(r'\s+', ' ', e.findtext('a:title', '', NS) or '').strip(),
                         re.sub(r'\s+', ' ', e.findtext('a:summary', '', NS) or '').strip()))
        if len(hits) < 100 * (st // 100 + 1):
            break
    n = 0
    for aid, yr, ttl, bs in hits:
        key = re.sub(r'v\d+$', '', aid.split('/abs/')[-1]) if '/abs/' in aid else re.sub(r'v\d+$', '', aid)
        key = key.split(':')[0] if ':' in key and not key[:4].isdigit() else key
        rec = seen.setdefault(key, dict(yr=yr, title=ttl, abs=bs, qs=[]))
        if q not in rec['qs']:
            rec['qs'].append(q)
    print('  %-58s 取回 %3d' % (q[:58], len(hits)))

for key, v in seen.items():
    txt = norm(v['title'] + ' ' + v['abs'])
    ms = [w for w in S if w in txt]
    ml = [w for w in L if w in txt]
    if ms and ml:
        CH.append((key, v['yr'], v['title'], ms, ml, v['abs'], len(v['qs'])))
CH.sort(key=lambda c: (-c[6], c[1]))
print('\n== [P4] 去重后 %d 条有摘要记录，其中 [A-CH]（两腿同时命中）= %d 条，按命中查询组数排序 ==' % (len(seen), len(CH)))
for i, (aid, yr, ttl, ms, ml, bs, nq) in enumerate(CH, 1):
    print('%2d. %s (%s) 命中 %d 组  S=%s L=%s\n    %s\n    %s' % (i, aid, yr, nq, ms[:4], ml[:4], ttl, bs[:260]))

print('\n== [P2] Crossref：标题级候选（无摘要 ⇒ 只登记，不判决）==')
for q in CQ:
    s = get('https://api.crossref.org/works?query.bibliographic=%s&rows=60&select=DOI,title,container-title,issued' % q)
    if s is None:
        continue
    try:
        msg = json.loads(s)['message']
    except Exception:
        continue
    k = 0
    for it in msg.get('items', []):
        ttl = (it.get('title') or [''])[0]
        txt = norm(ttl)
        ms = [w for w in S if w in txt]
        ml = [w for w in L if w in txt]
        if ms and ml:
            k += 1
            doi = it.get('DOI', '?')
            if doi.lower() not in [c[0].lower() for c in TCH]:
                TCH.append((doi, (it.get('issued', {}).get('date-parts') or [[None]])[0][0],
                            ttl, (it.get('container-title') or [''])[0]))
    print('  %-42s total=%12s 取 60  标题级两腿命中=%d' % (q[:42], format(msg.get('total-results', 0), ','), k))
print('\n  [T-CH] 清单（%d 条，全部是**欠账**，不许并入"已排除"）：' % len(TCH))
for doi, yr, ttl, cn in TCH[:25]:
    print('   · %s (%s) %s — %s' % (doi, yr, ttl[:78], cn[:34]))
Rloc = sorted(k for k in seen if any(k.startswith(x) or x.startswith(k) for x in local))
print('\n== 汇总 ==')
print('  arXiv 去重记录 %d；[A-CH] %d；本地已通读且命中两腿 %d 篇：%s' % (len(seen), len(CH), len(Rloc), Rloc))
print('  Crossref 标题级候选 %d 条（[P2]：未读全文 ⇒ 不构成排除）' % len(TCH))
print('  判据 [P4]：ρ_cov=%d ⇒ %s' % (len(CH), '仍须逐条读完才能定档' if CH else '带口径可写"无前占"'))
