# -*- coding: utf-8 -*-
r"""实验 138（§86-G-1 预注册）：两条**公开通道拿不到摘要**的欠账
（$10.1109$/cdc.2018.8619725 反向水填、$10.1109$/ciss.2016.7460485 trace 约束多译码器）
——换通道再试一次，并把 16 条 [U] 一并过一遍。**不许猜 URL**，只允许 API 返回的字段。

 [S1] 端点：Semantic Scholar `graph/v1/paper/DOI:{doi}`，
      `fields=title,abstract,year,venue,externalIds,openAccessPdf,citationCount,authors`。
      请求间 `sleep 1.2s`，429/403 退避 $10s\times1.6^k$ 重试 3 次；**只用返回值**，不拼接任何 URL。
 [S2] 判据（跑前写死）：
      若某条 `openAccessPdf.url` 非空 $\Rightarrow$ 记 [OA]，并把该 url 原样打印（是否可读由下一条决定）；
      若 `externalIds.ArXiv` 存在 $\Rightarrow$ 记 [AX]，随后向 arXiv `abs:` 取**完整摘要**；
      两者都无 $\Rightarrow$ 维持 [U]，并在正文那句 "unread through public channels" 里点名，**不许并入"已排除"**。
 [S3] 拿到摘要的每条回答三问（同 §86-C）：
      (a) 是否给出**控制代价**对信息率的闭式？(b) 是否出现随**秩/维数**变化的指数？(c) 传感矩阵是否为决策变量？
      三问全 "否" $\Rightarrow$ 定 [A-IR]；(a)(b) 有 "是" $\Rightarrow$ 立刻标 **[COLLIDE]** 并停下手头的其它活。
 [S4] 引用计数只作旁证（cdc.2018 若 citationCount 很小，说明它不是主流引用点，前占风险低）。
 [S5] 本脚本不做数值实验。
"""
import sys
import time
import json
import urllib.parse
import urllib.request
sys.stdout.reconfigure(encoding='utf-8')

UA = {'User-Agent': 'p0-e138 literature check (academic use)'}
SLEEP = 1.2
FIELDS = 'title,abstract,year,venue,externalIds,openAccessPdf,citationCount,authors'
S2 = 'https://api.semanticscholar.org/graph/v1/paper/DOI:'


def get(url, tries=3):
    for k in range(tries):
        try:
            req = urllib.request.Request(url, headers=UA)
            with urllib.request.urlopen(req, timeout=40) as r:
                return r.read().decode('utf-8', 'replace'), None
        except Exception as e:
            code = getattr(e, 'code', None)
            if code in (429, 403, 500, 502, 503):
                w = 10.0 * (1.6 ** k)
                print('    !! %s %s ⇒ 退避 %.0fs' % (type(e).__name__, code, w))
                time.sleep(w)
                continue
            return None, '%s %s' % (type(e).__name__, e)
    return None, 'RETRY-EXHAUSTED'


DOIS = ['10.1109/cdc.2018.8619725', '10.1109/ciss.2016.7460485',
        '10.1109/itw.2017.8277966', '10.1109/tit.2017.2694015', '10.1109/tit.2026.3714460',
        '10.1109/isit.2014.6874973', '10.1109/ssp.2005.1628751', '10.1109/isit63088.2025.11195430',
        '10.1109/icip.1995.529053', '10.1109/icip.1997.638663', '10.1109/icip.2001.958522',
        '10.1109/icip.1998.727409', '10.1109/icassp.2000.859200', '10.1109/ijcnn.2001.938831',
        '10.1049/el.2011.1734', '10.1007/978-3-642-55753-8_42']

LIMB_L = ['rate distortion', 'closed form', 'closed-form', 'exact', 'characterization', 'capacity',
          'law', 'watermark', 'water-filling', 'waterfilling', 'reverse water']
LIMB_S = ['sensor', 'vector', 'subspace', 'multiple sensor', 'dimension', 'rank', 'matrix']

pool = {}
print('== [S1][S2] Semantic Scholar 逐条（只用返回字段，不拼 URL）==')
for doi in DOIS:
    u = S2 + urllib.parse.quote(doi, safe='') + '?' + urllib.parse.urlencode({'fields': FIELDS})
    time.sleep(SLEEP)
    body, err = get(u)
    if err:
        print('  !! %-40s %s' % (doi, err))
        continue
    try:
        j = json.loads(body)
    except Exception:
        print('  !! %-40s JSON 解析失败' % doi)
        continue
    if not j.get('title'):
        print('  [404 ] %-40s S2 无此记录' % doi)
        continue
    oa = (j.get('openAccessPdf') or {}).get('url') or ''
    ax = (j.get('externalIds') or {}).get('ArXiv') or ''
    ab = j.get('abstract') or ''
    tags = []
    if oa:
        tags.append('[OA]')
    if ax:
        tags.append('[AX %s]' % ax)
    if not (oa or ax or ab):
        tags.append('[U]')
    print('  %-6s %-40s %s  cit=%s  %s' % ((''.join(tags))[:14], doi, j.get('year'),
                                            j.get('citationCount'), (j.get('title') or '')[:58]))
    if oa:
        print('        OA url（原样）: %s' % oa)
    if ab:
        print('        摘要：%s' % ab[:400].replace('\n', ' '))
    rec = pool.setdefault(doi, {})
    rec.update(title=j.get('title'), year=j.get('year'), cit=j.get('citationCount'),
               oa=oa, ax=ax, ab=ab, venue=j.get('venue'))

print('\n== [S2 续] 有 ArXiv id 但 S2 没给摘要的，去 arXiv abs 补摘要 ==')
for doi, r in pool.items():
    if r.get('ax') and not r.get('ab'):
        time.sleep(3.5)
        u = 'https://export.arxiv.org/api/query?id_list=' + r['ax']
        body, err = get(u)
        if err:
            print('  !! %s arXiv %s' % (doi, err))
            continue
        import re as _re
        m = _re.search(r'<summary>(.*?)</summary>', body, _re.S)
        if m:
            r['ab'] = ' '.join(m.group(1).split())
            print('  [补摘要 OK] %s (%s)：%s' % (doi, r['ax'], r['ab'][:200]))
        else:
            print('  [补摘要失败] %s' % doi)

print('\n== [S3] 三问卡（a 控制代价闭式？ b 指数随秩/维数？ c 传感矩阵为决策变量？）==')
KEY_A = ['control cost', 'lq cost', 'linear quadratic', 'estimation error cost', 'mse cost', 'cost function']
KEY_B = ['exponent', 'converges to', 'asymptotic rate', 'decreases with dimension', 'per dimension', 'rank-dependent']
KEY_C = ['sensor selection', 'sensor placement', 'design the sensor', 'measurement matrix', 'sensing matrix',
         'jointly optimize', 'encoder design']
n_collide = 0
for doi, r in sorted(pool.items()):
    ab = (r.get('ab') or '').lower()
    if not ab:
        print('  %-40s 无摘要 ⇒ 维持 [U]，三问不判（**不许当已排除**）' % doi)
        continue
    a = [k for k in KEY_A if k in ab]
    b = [k for k in KEY_B if k in ab]
    c = [k for k in KEY_C if k in ab]
    s = [k for k in LIMB_S if k in ab]
    l = [k for k in LIMB_L if k in ab]
    collide = bool(a) and bool(b)
    n_collide += collide
    print('  %-40s a=%-2s b=%-2s c=%-2s   %s' % (doi, '是' if a else '否', '是' if b else '否', '是' if c else '否',
                                                 '**[COLLIDE]**' if collide else ''))
    if a or b or c:
        print('      命中词：a=%s b=%s c=%s' % (a, b, c))
    print('      %s' % (r['title'] or '')[:80])

n_u = sum(1 for r in pool.values() if not r.get('ab'))
n_oa = sum(1 for r in pool.values() if r.get('oa'))
n_ax = sum(1 for r in pool.values() if r.get('ax'))
print('\n== [S2][S4] 汇总 ==')
print('  查询 %d 条：S2 命中 %d；有 OA url %d；有 ArXiv id %d；拿到摘要 %d；仍 [U] %d'
      % (len(DOIS), len(pool), n_oa, n_ax, len(pool) - n_u, n_u))
print('  [COLLIDE]（a$\wedge$b 同时为是）= %d ⇒ %s'
      % (n_collide, '0 ⇒ §86-E 的覆盖面陈述可维持（未读者仍点名）' if n_collide == 0
         else '**>0 ⇒ 立刻停手，先把它读完再谈任何新颖性**'))
print('  未读的必须点名：' + ', '.join(sorted(d for d, r in pool.items() if not r.get('ab'))))
print('EXIT=0')
