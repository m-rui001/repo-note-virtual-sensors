# -*- coding: utf-8 -*-
r"""实验 135c（§84-E 补记-5 预注册）：**把 e135b 的三笔缺陷补掉，并把 5 条邻居文献读到正文级**。

跑前写死（照 §84-E 补记-5 的三条，不扩权）：
 [Q1] arXiv 请求间 `sleep 3.5s`、429 退避 $20s$ 重试 3 次；重跑 e135b 里被限流打回的 3 组，
      并把这 3 组与 e135b 已成功的 10 组**合并**成一个真正的去重记录池再分类（词表与 e135/$S$-$L$ 逐字相同）。
 [Q2] 17 条 [T-CH] 的 DOI 从 `p0/e135b_out.txt` 现场解析（不抄数），逐条向 Crossref 要**摘要**：
      拿到摘要 ⇒ 用 $S$/$L$ 两腿重判 $=$ [A-CH]/[A-IR]/[U]；拿不到 ⇒ **留 [U]**，不许并入"已排除"（纪律 24③）。
 [Q3] §84-E 补记-3 表里的 5 条邻居，逐条额外做**开放获取查找**（arXiv `ti:` 近似标题检索），
      命中则取回完整摘要 $+$ 年份 $+$ 作者，并按预注册的三问逐条回答：
        (a) 它有没有闭式（还是只有界/算法/渐近式）？
        (b) 它的指数/斜率**是否随维数或秩变化**（这是顶我 §84-D 那句话的那一问）？
        (c) 它的"传感/编码矩阵"是自由变量还是给定？
      判据：任一文献在 (b) 答"是" ⇒ 84-D 英文句降级为引用；全部答"否/未答" ⇒ 保持比较级陈述。
 [Q4] 本脚本**不做**任何数值实验，产出只有：合并覆盖面计数、逐条定档、5 张三问卡。
      任何"读不到全文"的条目一律写明"摘要级"，不许用摘要级读数为定理级陈述背书。
"""
import sys, re, json, time, urllib.parse, urllib.request, urllib.error
import xml.etree.ElementTree as ET
sys.stdout.reconfigure(encoding='utf-8')

UA = {'User-Agent': 'prior-art-check/1.0 (mailto:p0-lane@example.org)'}
NS = {'a': 'http://www.w3.org/2005/Atom'}
SLEEP = 3.5

S = ['multidimensional', 'multi-dimensional', 'vector', 'sensing rank', 'rank', 'sensor selection',
     'measurement matrix', 'sensing matrix', 'partial observation', 'subspace', 'multivariate',
     'sensor capacity', 'multiple sensor', 'measurement dimension']
L = ['exponent', 'closed form', 'closed-form', 'rate distortion', 'rate-distortion', 'data rate',
     'capacity', 'logarithmic', 'exponential', 'rate-cost', 'trade-off', 'tradeoff', 'slope',
     'characterization', 'exact']


def norm(s):
    return re.sub(r'[^a-z0-9]+', ' ', (s or '').lower())


def legs(t):
    x = norm(t)
    return [w for w in S if w in x], [w for w in L if w in x]


def fetch(url, tries=3, backoff=20.0):
    for a in range(tries):
        try:
            req = urllib.request.Request(url, headers=UA)
            with urllib.request.urlopen(req, timeout=90) as f:
                return f.read().decode('utf-8', 'replace')
        except urllib.error.HTTPError as e:
            if e.code == 429:
                print('      (429，退避 %.0fs 后重试)' % backoff)
                time.sleep(backoff)
                backoff *= 1.8
                continue
            print('      !! HTTP %s %s' % (e.code, url[:80]))
            return None
        except Exception as e:
            print('      !! %s %s' % (type(e).__name__, url[:80]))
            return None
    print('      !! 重试耗尽 %s' % url[:80])
    return None


def arxiv(q, cap=300):
    out = []
    st = 0
    while st < cap:
        s = fetch('https://export.arxiv.org/api/query?search_query=' + urllib.parse.quote_plus(q)
                  + '&start=%d&max_results=100&sortBy=relevance' % st)
        time.sleep(SLEEP)
        if s is None:
            break
        try:
            root = ET.fromstring(s)
        except Exception:
            break
        got = 0
        for e in root.findall('a:entry', NS):
            aid = re.sub(r'v\d+$', '', (e.findtext('a:id', '', NS) or '').split('/abs/')[-1])
            out.append((aid, (e.findtext('a:published', '', NS) or '')[:4],
                        re.sub(r'\s+', ' ', e.findtext('a:title', '', NS) or '').strip(),
                        re.sub(r'\s+', ' ', e.findtext('a:summary', '', NS) or '').strip()))
            got += 1
        if got < 100:
            break
        st += 100
    return out


RETRY = ['all:"communication rate" AND all:"estimation error"',
         'all:"sensor selection" AND all:"information rate"',
         'all:"Gaussian source" AND all:"remote estimation" AND all:"multiple"']
DONE10 = ['all:"rate distortion" AND all:"control cost"', 'all:"data rate" AND all:"LQG"',
          'all:"measurement rank" AND all:"control cost"', 'all:"Kalman filtering" AND all:"bandwidth constraint"',
          'all:"remote state estimation" AND all:"rate distortion"', 'all:"minimal data rate" AND all:"stabilization"',
          'all:"directed information" AND all:"linear quadratic"', 'all:"cost of communication" AND all:"linear systems"',
          'all:"dimension" AND all:"rate-distortion" AND all:"control"',
          'all:"sensor scheduling" AND all:"estimation" AND all:"rate"']

print('== [Q1] 重跑 3 组被限流的查询（sleep %.1fs）==' % SLEEP)
pool = {}
for q in RETRY + DONE10:
    tag = 'RETRY' if q in RETRY else 'rescan'
    hits = arxiv(q)
    n0 = len(pool)
    for aid, yr, ttl, bs in hits:
        pool.setdefault(aid, dict(yr=yr, title=ttl, abs=bs, qs=[]))
        if q not in pool[aid]['qs']:
            pool[aid]['qs'].append(q)
    print('  [%-6s] %-58s 取回 %3d  新增 %d' % (tag, q[:58], len(hits), len(pool) - n0))

CH = []
for aid, v in pool.items():
    ms, ml = legs(v['title'] + ' ' + v['abs'])
    if ms and ml:
        CH.append((aid, v['yr'], v['title'], ms, ml, v['abs']))
print('\n  合并去重池 %d 条；[A-CH]（两腿同时命中）= %d 条' % (len(pool), len(CH)))
for aid, yr, ttl, ms, ml, bs in sorted(CH, key=lambda c: c[1]):
    inpool135b = ''
    print('   · %s (%s) S=%s L=%s  %s' % (aid, yr, ms[:3], ml[:3], ttl[:88]))

print('\n== [Q2] 从 e135b_out.txt 现场解析 [T-CH] DOI，逐条向 Crossref 要摘要 ==')
log = open('p0/e135b_out.txt', encoding='utf-8').read()
DOIS = re.findall(r'·\s+(10\.[^\s]+)\s', log)
DOIS = list(dict.fromkeys(DOIS))
print('  解析到 %d 条 DOI' % len(DOIS))
ABS, NOABS = [], []
for doi in DOIS:
    s = fetch('https://api.crossref.org/works/' + urllib.parse.quote(doi, safe=''))
    time.sleep(0.6)
    if s is None:
        NOABS.append((doi, '?', '取记录失败'))
        continue
    try:
        it = json.loads(s)['message']
    except Exception:
        NOABS.append((doi, '?', 'JSON 失败'))
        continue
    ttl = (it.get('title') or [''])[0]
    ab = re.sub(r'<[^>]+>', ' ', it.get('abstract') or '')
    ab = re.sub(r'\s+', ' ', ab).strip()
    yr = (it.get('issued', {}).get('date-parts') or [[None]])[0][0]
    if ab:
        ms, ml = legs(ttl + ' ' + ab)
        grade = 'A-CH' if (ms and ml) else ('A-IR' if (ms or ml) else 'U')
        ABS.append((doi, yr, ttl, grade, ms, ml, ab))
        print('  [%s] %s (%s) %s' % (grade, doi, yr, ttl[:74]))
    else:
        NOABS.append((doi, yr, ttl))
        print('  [无摘要→留U] %s (%s) %s' % (doi, yr, ttl[:74]))
print('  小结：拿到摘要 %d 条（其中 [A-CH] %d）；仍无摘要 ⇒ 留 [U] %d 条'
      % (len(ABS), sum(1 for a in ABS if a[3] == 'A-CH'), len(NOABS)))
for doi, yr, ttl in NOABS:
    print('    [U] %s (%s) %s' % (doi, yr, ttl[:76]))

print('\n== [Q3] 5 条邻居：开放获取查找 + 三问卡（a 闭式？b 指数随维数/秩？c 传感矩阵自由？）==')
NEIGH = [
    ('10.1109/cdc.2018.8619725', 'Asymptotic Reverse Waterfilling Characterization Nonanticipative Rate Distortion'),
    ('10.1109/tit.2017.2694015', 'Vector Gaussian Rate-Distortion With Variable Side Information'),
    ('10.1109/ciss.2016.7460485', 'Vector Gaussian multi-decoder rate-distortion Trace constraints'),
    ('10.1109/itw.2017.8277966', 'upper bound to zero-delay rate distortion via Kalman filtering vector Gaussian'),
    ('10.1109/tit.2026.3714460', 'Gaussian-Quadratic Rate-Distortion Function for Vector Sources'),
]
for doi, t in NEIGH:
    print('\n  ── %s' % doi)
    hits = arxiv('ti:"%s"' % t) or arxiv('all:"%s"' % t)
    if not hits:
        hits = arxiv('all:"%s"' % ' '.join(t.split()[:6]))
    if not hits:
        print('     arXiv 无命中 ⇒ 只有 IEEE 付费版本，本轮读不到全文（记 [U]，摘要级判据都不给）')
        continue
    for aid, yr, ttl, bs in hits[:3]:
        print('     %s (%s) %s' % (aid, yr, ttl[:86]))
        print('       摘要：%s' % bs[:520])
        print('       两腿命中：S=%s L=%s' % (legs(ttl + ' ' + bs)))
