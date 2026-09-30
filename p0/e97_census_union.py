r"""E97 = §64-F-2 欠的普查：纪律 26 要求的"≥3 组同义扩展取并"，产出真正的并集全集表。

背景（为什么这条是欠账而不是重复劳动）：`p0/results/novelty_waterfill_20260930.md:35` 自称
"普查整个对撞线：12 条全集"，它来自单条布尔查询 `all:"directed information" AND all:LQG`。
§64-A 判定这条查询结构性失效两重（短语打不中 *directed mutual information*；`all:LQG` 打不中 LQR 口径）。
本轮先实测一条**更准的漏检机制**：arXiv 的 `all:` 只检索**元数据**（标题/摘要/评论/作者），不检索正文，
所以 §64-A 在 `papers/txt/1612.02126.txt` 上做的"全文各 6 次"核对**不能**直接用来解释命中与否。
实测（`p0/e97_out.txt` §0）：
  `all:"directed information" AND all:LQG` → 12 条，**不含** 1612.02126；
  只把**信息轴**扩成 `("directed information" OR "directed mutual information" OR "causal rate-distortion")`
  再接控制轴 → 14 条，**含** 1612.02126。
⇒ 漏检的真凶是**信息词形**（directed mutual information），控制轴（LQG→LQR）单独并不够。
这条修正本身要上报，因为它改写 §64-A 的归因。

普查要回答的实质问题只有一个（不是收集编号）：
  **受限传感（给定 range(S)/平面/秩 r）下，率-代价前沿的对偶界/反证界是否已有人发表？**

== 判据（跑前写死，跑后不许改）==
 [N1] 并集**必须**包含 `1612.02126`（Kostina–Hassibi）。不在 ⇒ 不许写"全集"，只许写"并集样本"。
 [N2] ≥3 组互为同义扩展的查询，每组单独落盘查询串 + totalResults + 抓取数；
      抓取数 < totalResults 时必须分页抓全，否则明写截断。
 [N3] 并集每一条要么本地已有笔记（视为已读），要么在表里明写"未读、不作前占判断"（纪律 24③）。
 [N4] 前占判断只看摘要层面的关键词证据（converse / lower bound / dual / Lagrang / SDP / exact / rank），
      命中者逐条人工读摘要后才进正文；未命中者记为"该查询未给出对偶界证据"，
      **不许**由"未命中"外推成"没人做过"。
 [N5] 在并集、与对偶界相关、但本地无 PDF = 新增义务清单，按相关度排序，不许只列编号。
 [N6] 精筛规则（跑前写死）：一条命中算"在对撞线上" ⟺ 标题+摘要**同时**含
      (i) 信息词形之一 {directed information, directed mutual information, causal rate-distortion}
      (ii) 控制/代价词形之一 {LQG, LQR, linear quadratic, quadratic cost, control cost, controller,
            linear system, stabilization, estimation}。
      并报"宽口径并集 → 精筛后"各层数量；baseline 的 12 条与 1612.02126 必须都存活，否则精筛规则作废。
 [N7] 第二数据库通道（Crossref，只认 DOI/元数据，不依赖 arXiv）：≥2 组同义扩展查询，
      记录查询串与 `items-per-page` 上限下的命中数；凡精筛存活且本地无 PDF 者进 N5 义务清单。
"""
import sys, os, re, time, json, urllib.parse, urllib.request
import xml.etree.ElementTree as ET
sys.stdout.reconfigure(encoding='utf-8')

NS = {'a': 'http://www.w3.org/2005/Atom', 'os': 'http://a9.com/-/spec/opensearch/1.1/'}
INFO = '(all:"directed information" OR all:"directed mutual information" OR all:"causal rate-distortion")'
CTRL = ('(all:LQG OR all:LQR OR all:"linear quadratic" OR all:"quadratic cost" OR all:"control cost")')
SENS = ('(all:"sensor" OR all:"actuator" OR all:"sensor selection" OR all:"partial observation" '
        'OR all:"task-based sensing" OR all:"measurement")')
DUAL = '(all:"converse" OR all:"lower bound" OR all:"dual" OR all:"Lagrangian" OR all:"semidefinite")'
RATE = '(all:"rate-distortion" OR all:"rate distortion" OR all:"data rate" OR all:"communication rate")'

QUERIES = [
    ('baseline', '旧那条（已知会漏，留作差集与 [N6] 存活检查）', f'all:"directed information" AND all:LQG'),
    ('G1', '信息轴同义扩展 × 控制轴同义扩展（两轴同时放开）', f'{INFO} AND {CTRL}'),
    ('G2', '信息轴 × 传感/观测对象轴（受限传感这一支）', f'{INFO} AND {SENS}'),
    ('G3', '信息轴 × 对偶/下界轴（直接问"反证界是否已知"）', f'{INFO} AND {DUAL}'),
    ('G4', '信息轴 × 率失真轴（rate-distortion 口径的 DI 控制论文）', f'{INFO} AND {RATE}'),
]
INFO_WORDS = ['directed information', 'directed mutual information', 'causal rate-distortion',
              'causal rate distortion']
CTRL_WORDS = ['lQG', 'lQR', 'linear quadratic', 'quadratic cost', 'control cost', 'controller',
              'linear system', 'stabilization', 'estimation']
CTRL_WORDS = [w.lower() for w in CTRL_WORDS]
KW = dict(converse=r'\bconverse\b', lowerbound=r'lower bound', dual=r'\bdual|lagrang|kkt',
          sdp=r'\bsdp|semidefinite|convex', exact=r'\bexact\b|achievable|optimal',
          rank=r'\brank|dimension|subset|selection')


def get(url):
    for a in range(4):
        try:
            with urllib.request.urlopen(url, timeout=90) as f:
                d = f.read().decode('utf-8', 'replace')
            if len(d) > 400:
                return d
        except Exception as e:
            print(f'   retry {a+1}: {type(e).__name__} {str(e)[:60]}')
        time.sleep(6.0)
    raise RuntimeError('failed 4x: ' + url)


def arxiv(q, cap=1200):
    hits, tot = [], None
    start = 0
    while True:
        u = ('https://export.arxiv.org/api/query?search_query=' + urllib.parse.quote_plus(q)
             + f'&start={start}&max_results=200&sortBy=relevance&sortOrder=descending')
        root = ET.fromstring(get(u))
        if tot is None:
            tot = int(root.find('os:totalResults', NS).text)
        es = root.findall('a:entry', NS)
        for e in es:
            aid = re.sub(r'v\d+$', '', e.find('a:id', NS).text.split('/abs/')[-1])
            hits.append(dict(id=aid,
                             title=' '.join(e.find('a:title', NS).text.split()),
                             pub=e.find('a:published', NS).text[:10],
                             abstract=' '.join(e.find('a:summary', NS).text.split())))
        start += len(es)
        print(f'   页 start={start-len(es)} 抓到 {len(es)}，累计 {len(hits)}/{tot}', flush=True)
        time.sleep(3.5)
        if not es or start >= tot or start >= cap:
            break
    return tot, hits


print('== (0) 漏检机制的实测分解（先归因，再普查）==')
for lab, q in [('只放开信息轴（控制轴仍写死 LQG）', 'all:"directed mutual information" AND all:LQG'),
               ('只放开控制轴（信息轴仍写死短语）', 'all:"directed information" AND (all:LQR OR all:"linear quadratic")'),
               ('两轴都放开', f'{INFO} AND {CTRL}')]:
    tot, hits = arxiv(q)
    ids = {h['id'] for h in hits}
    kh = '1612.02126' in ids
    print(f'  {lab}\n    query={q}\n    total={tot}  KH(1612.02126) 命中={kh}')
    time.sleep(3.0)

all_hits, per_query = {}, []
for tag, why, q in QUERIES:
    print(f'\n== [{tag}] {why} ==\n  query: {q}', flush=True)
    try:
        tot, hits = arxiv(q)
    except RuntimeError as e:
        print(f'  !! 失败，本轮并集不含它：{e}')
        per_query.append((tag, q, -1, -1, []))
        continue
    print(f'  totalResults={tot}  抓取={len(hits)}  [N2] {"完整" if len(hits) >= tot else "截断"}')
    per_query.append((tag, q, tot, len(hits), [h['id'] for h in hits]))
    for h in hits:
        cur = all_hits.setdefault(h['id'], h)
        cur.setdefault('from', [])
        if tag not in cur['from']:
            cur['from'].append(tag)

ids = sorted(all_hits, key=lambda k: (all_hits[k]['pub'], k))
base = set(per_query[0][4])
print(f'\n== (1) [N1/N2] 宽口径并集：{len(ids)} 条（逐条列举，方为"全集"；纪律 26）==')
for i in ids:
    h = all_hits[i]
    print(f'{i:<12} {h["pub"]} [{">".join(h["from"]):<28}] {h["title"][:86]}')
print(f'\n 1612.02126 在并集？{"在" if "1612.02126" in all_hits else "不在——禁写全集"}'
      f'  来自 {all_hits.get("1612.02126", {}).get("from")}   旧 baseline 含它？{"含" if "1612.02126" in base else "不含"}')
print(f' 相对旧单条查询净新增 {len(set(ids) - base)} 条；旧 12 条全在并集内？'
      f'{"是" if base <= set(ids) else "否（说明分页/上限有问题）"}')

print('\n== (2) [N6] 精筛：宽口径并集 → 对撞线 ==')
def precise(h):
    t = (h['title'] + ' ' + h['abstract']).lower()
    return any(w in t for w in INFO_WORDS) and any(w in t for w in CTRL_WORDS)
pids = [i for i in ids if precise(all_hits[i])]
print(f' 宽口径 {len(ids)} → 精筛 {len(pids)}')
surv = base <= set(pids)
print(f' [N6] baseline 12 条全部存活？{"是" if surv else "否——精筛规则作废，改用宽口径"}')
print(f' [N6] 1612.02126 精筛存活？{"是" if "1612.02126" in pids else "否"}')

PDFS, NOTES = set(), set()
for p in os.listdir('papers'):
    m = re.match(r'(\d{4}\.\d{4,5})', p)
    if m: PDFS.add(m.group(1))
for p in os.listdir('papers/notes'):
    m = re.match(r'(\d{4}\.\d{4,5})', p)
    if m: NOTES.add(m.group(1))

print('\n== (3) [N3/N4] 精筛清单：本地覆盖 + 摘要关键词证据 ==')
rows = []
for i in pids:
    h = all_hits[i]
    t = (h['title'] + ' ' + h['abstract']).lower()
    flags = [k for k, rx in KW.items() if re.search(rx, t)]
    local = 'PDF+note' if i in NOTES else ('PDF only' if i in PDFS else '未入库')
    rows.append((i, h, flags, local))
    star = '★' if ({'converse', 'lowerbound', 'dual'} & set(flags)) and ({'rank', 'sdp'} & set(flags)) else ' '
    print(f' {star}{i:<12}{local:<9}kw={",".join(flags) if flags else "-":<36}{h["title"][:72]}')

print(f'\n== (4) [N5] 精筛且"★"（对偶/下界 ∩ 秩/凸）但本地无 PDF —— 逐条读摘要 ==')
cand = [r for r in rows if r[3] == '未入库' and ({'converse', 'lowerbound', 'dual'} & set(r[2]))
        and ({'rank', 'sdp', 'exact'} & set(r[2]))]
print(f' 共 {len(cand)} 条')
for i, h, flags, _ in sorted(cand, key=lambda r: -len(r[2])):
    print(f'\n--- {i} kw={flags}  {h["title"]}  ({h["pub"]})')
    print('    ' + h['abstract'][:1400])

print('\n== (5) 旧 baseline 的入库状态（§53 失误复现检查）==')
for i in sorted(base):
    print(f'   {i:<12}{"已入库" if i in (PDFS | NOTES) else "**未入库**":<8}{all_hits[i]["title"][:74]}')

print('\n== (6) [N7] 第二通道 Crossref（不依赖 arXiv 元数据）==')
CR = [('CR1', 'directed information LQG', 'directed+information+LQG'),
      ('CR2', 'directed mutual information linear quadratic', 'directed+mutual+information+linear+quadratic'),
      ('CR3', 'task-based sensing rate distortion control cost', 'task-based+sensing+rate+distortion+control+cost')]
for tag, lab, q in CR:
    u = (f'https://api.crossref.org/works?query.bibliographic={q}&rows=100&select=DOI,title,container-title,issued')
    try:
        d = json.loads(get(u))
        items = d['message']['items']
        print(f'\n [{tag}] {lab}\n   query={u[:110]}...\n   命中 items={len(items)} (items-per-page=100 上限)')
        keep = []
        for it in items:
            tt = (it.get('title') or [''])[0].lower()
            if any(w in tt for w in ['directed', 'rate-distortion', 'sensor', 'estimation', 'lqr', 'lqg']):
                keep.append((it.get('DOI'), tt[:74], (it.get('container-title') or [''])[0][:28]))
        for doi, tt, ct in keep[:14]:
            print(f'     {doi:<34}{ct:<28}{tt}')
        print(f'   标题相关 {len(keep)} 条')
    except Exception as e:
        print(f' [{tag}] 失败 {type(e).__name__} {str(e)[:80]}')
    time.sleep(2.0)

json.dump(dict(per_query=per_query, union=all_hits, precise=pids),
          open('p0/e97_census_raw.json', 'w', encoding='utf-8'), ensure_ascii=False, indent=1)
print('\n原始数据落盘 p0/e97_census_raw.json')
