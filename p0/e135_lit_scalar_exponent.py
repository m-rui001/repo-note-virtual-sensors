# -*- coding: utf-8 -*-
r"""实验 135（§83-F-③ 预注册，结清 §77-G-① 的文献欠账）：
**"标量闭式 $\gamma=2\ln2$ 在向量传感情形不成立"这句话，有没有人已经写过？**

正文要引这条（§81-F/§82-G 的限定句都靠它），所以我必须先假设**它已经被写过**，把它当成威胁项去找，
而不是按我自己的措辞去找。通道两条，分工写死：

 [N1] **Crossref**：按 DOI 直取 Tatikonda–Sahai–Mitter TAC 49(9):1549–1561 (2004), doi 10.1109/TAC.2004.834430
      的记录（标题/期刊/年份/被引数），再用 `query.bibliographic` 补 4 组近义查询。
      按 §66-M1 的结论，**不许**把这里的 total-results 当"全集"。
 [N2] **arXiv API**：8 组 `abs:` 查询，每组取回 $\le200$ 条，按相关性排序翻页到该组上限。
 [N3] 三档分类（判据跑前写死，词表在下面），只有 [A-CH] 才允许威胁我的陈述：
      传感/向量腿词 $S$ = {rank, multidimensional, vector, sensor selection, measurement matrix, sensing matrix, partial observation, subspace}
      律/闭式腿词 $L$ = {exponent, closed form, rate distortion, data rate, capacity, logarithmic, exponential, rate-cost, tradeoff}
      [A-CH] = 标题或摘要同时命中 $S$ 与 $L$（两腿都要）；[A-IR] = 只命中一腿；[U] = 都不命中。
      另有 [R] = 本地 `papers/notes/` 已有该 arXiv 号 ⇒ 已通读，须回笔记核。
 [N4] 判决口径：$k=$ [A-CH] 条数。$k=0$ ⇒ 这句话在我检索到的语料里无前占，正文可写，但要带上检索口径。
      $k\ge1$ ⇒ **逐条读完摘要并列出**，凡真的写了"$\gamma$（或等价指数）随传感维数/秩变化"或"标量闭式在向量情形失效"者，
      我的陈述立刻降级为引用；剩下的记为"无关但同名"。不许用"它没给证明"糊过去。
"""
import sys, os, re, json, time, glob, urllib.parse, urllib.request
import xml.etree.ElementTree as ET
sys.stdout.reconfigure(encoding='utf-8')

UA = {'User-Agent': 'prior-art-check/1.0 (mailto:p0-lane@example.org)'}
NS = {'a': 'http://www.w3.org/2005/Atom'}


def get(url, tries=4):
    last = None
    for _ in range(tries):
        try:
            req = urllib.request.Request(url, headers=UA)
            with urllib.request.urlopen(req, timeout=90) as f:
                return f.read().decode('utf-8', 'replace')
        except Exception as e:
            last = e
            time.sleep(4.0)
    print('    !! 取不到 %s: %s %s' % (url[:90], type(last).__name__, str(last)[:80]))
    return None


def gj(url):
    s = get(url)
    if s is None:
        return None
    try:
        return json.loads(s)
    except Exception:
        return None


def norm(s):
    return re.sub(r'[^a-z0-9]+', ' ', (s or '').lower())


S = ['multidimensional', 'multi-dimensional', 'vector', 'sensing rank', 'rank', 'sensor selection',
     'measurement matrix', 'sensing matrix', 'partial observation', 'subspace', 'multivariate']
L = ['exponent', 'closed form', 'closed-form', 'rate distortion', 'rate-distortion', 'data rate',
     'capacity', 'logarithmic', 'exponential', 'rate-cost', 'trade-off', 'tradeoff']

local = set()
for p in glob.glob('papers/notes/*.md'):
    m = re.search(r'(\d{4}\.\d{4,5})', os.path.basename(p))
    if m:
        local.add(m.group(1))

print('== [N1] Crossref：Tatikonda–Sahai–Mitter 2004 的权威记录 ==')
d = gj('https://api.crossref.org/works/10.1109/TAC.2004.834430')
if d:
    it = d['message']
    print('  title   = %s' % (it.get('title') or ['?'])[0])
    print('  container=%s  year=%s  cited-by=%s  pages=%s'
          % ((it.get('container-title') or ['?'])[0],
             it.get('issued', {}).get('date-parts', [[None]])[0][0], it.get('is-referenced-by-count'),
             it.get('page')))
    print('  authors =', ', '.join((a.get('family') or '?') for a in it.get('author', [])))
    ab = it.get('abstract')
    print('  abstract= %s' % (re.sub(r'<[^>]+>', ' ', ab)[:600] if ab else '（Crossref 无摘要）'))
    print('  命中律腿词 %s / 传感腿词 %s' % ([w for w in L if w in norm((it.get('title') or [''])[0] + ' ' + (ab or ''))],
                                             [w for w in S if w in norm((it.get('title') or [''])[0] + ' ' + (ab or ''))]))
    refs = gj('https://api.crossref.org/works/10.1109/TAC.2004.834430/references')
    if refs:
        rl = refs['message'].get('reference', []) or []
        print('  它的参考文献 %d 条（下面只列同时命中两腿的）' % len(rl))
        n = 0
        for r in rl:
            t = norm(' '.join([str(r.get('unstructured') or ''), str((r.get('articleTitle') or '')),
                               str((r.get('series-title') or ''))]))
            if any(w in t for w in S) and any(w in t for w in L):
                n += 1
                print('    · %s' % (r.get('unstructured') or r.get('articleTitle'))[:150])
        if n == 0:
            print('    （0 条）')
else:
    print('  !! DOI 记录未取到 ⇒ 本项欠账仍未结，不许写"已核")')

print('\n== [N2][N3] arXiv 8 组查询 + 三档分类（本地笔记 %d 篇标 [R]）==' % len(local))
QS = ['abs:"rate distortion" AND abs:"linear quadratic"',
      'abs:"directed information" AND abs:"LQG"',
      'abs:"sensor selection" AND abs:"Kalman" AND abs:"rate"',
      'abs:"data rate" AND abs:"estimation error" AND abs:"Riccati"',
      'abs:"cost of communication" AND abs:"control"',
      'abs:"remote estimation" AND abs:"sensing" AND abs:"rank"',
      'abs:"information-based cost" AND abs:"communication"',
      'abs:"exponent" AND abs:"data rate" AND abs:"Kalman"']
seen, CH = {}, []
for q in QS:
    hits = []
    for st in (0, 100):
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
            ttl = re.sub(r'\s+', ' ', e.findtext('a:title', '', NS) or '').strip()
            bs = re.sub(r'\s+', ' ', e.findtext('a:summary', '', NS) or '').strip()
            yr = (e.findtext('a:published', '', NS) or '')[:4]
            hits.append((aid, yr, ttl, bs))
        if st == 0 and len(hits) < 100:
            break
    tot_txt = len(hits)
    nch = nir = nu = nR = 0
    for aid, yr, ttl, bs in hits:
        key = aid.split('v')[0]
        rec = seen.setdefault(key, dict(yr=yr, title=ttl, abs=bs, qs=[]))
        rec['qs'].append(q)
        txt = norm(ttl + ' ' + bs)
        ms = [w for w in S if w in txt]
        ml = [w for w in L if w in txt]
        isR = any(key.startswith(x) or x.startswith(key) for x in local)
        if ms and ml:
            nch += 1
            if not isR:
                CH.append((key, yr, ttl, ms, ml, bs))
        elif ms or ml:
            nir += 1
        else:
            nu += 1
        if isR:
            nR += 1
    print('  %-58s 取回 %3d  [A-CH]=%3d [R]=%2d [A-IR]=%3d [U]=%2d'
          % (q[:58], tot_txt, nch, nR, nir, nu))

CH2 = sorted({c[0]: c for c in CH}.values(), key=lambda c: c[1], reverse=True)
print('\n== [N4] 合并去重后 [A-CH] 共 %d 条（按年份倒序，全部列出摘要前 320 字）==' % len(CH2))
for i, (aid, yr, ttl, ms, ml, bs) in enumerate(CH2, 1):
    print('%2d. %s (%s)  S=%s L=%s\n    %s\n    %s' % (i, aid, yr, ms, ml, ttl, bs[:320]))
print('\n  本地已通读且命中两腿的：', sorted(
    k for k, v in seen.items() if any(k.startswith(x) or x.startswith(k) for x in local)))
print('  判据 [N4]：k=%d ⇒ %s' % (len(CH2), '无前占（检索口径见上），正文可写这条 + 检索口径'
                                  if len(CH2) == 0 else '有候选威胁项，必须逐条读完再定档'))
