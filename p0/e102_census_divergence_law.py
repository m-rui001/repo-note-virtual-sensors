r"""E102 = e101b 的 [U3] 文献义务：我量到"率贴地板时代价按 1/(I−R_exp) 发散、指数对一切可检测
观测子空间等于 1"，必须先查这条**是否已发表**（纪律 26 + e96f 自己写死的 [U3]）。
e96f 查的是另一头（D→D_min 时 DI 对 log(1/x) 的斜率 = r/2），本轮查的是**地板那一头**。

普查协议（跑前写死，跑后不许改）：
 [V1] ≥3 组互为同义扩展的查询；每条落盘查询串 + totalResults + 实抓数；未抓全必须分页抓全，
      否则明写截断。`all:` 只检索元数据（标题/摘要/评论/作者），不检索正文 ⇒ 结论只到摘要级。
 [V2] 一条命中算"已发表近地板发散律" ⟺ 标题+摘要**同时**含
      (a) 地板词形 {critical data rate, minimum data rate, infimum rate, rate asymptote, horizontal asymptote}
      (b) 发散词形 {cost diverges, unbounded cost, blows up, infinite cost, goes to infinity,
                    tends to infinity, no finite cost}。
      任一命中 ⇒ 我的指数 1 律**降级为复现**，并当场在板上自首。
 [V3] 只含 (a) 的记"邻居"（地板本身确实早被发表：1510.04214 Cor.1 的 1.169、2606.31396 Th.1 的 R_exp），
      不许把它当成"发散速率也发表了"。
 [V4] 不许由"未命中"外推"没人做过"；只报本次并集的覆盖边界。
"""
import json
import re
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET

OUT = open('p0/e102_out.txt', 'w', encoding='utf-8')


def p(*a):
    s = ' '.join(str(x) for x in a)
    OUT.write(s + '\n')
    OUT.flush()


NS = {'a': 'http://www.w3.org/2005/Atom', 'os': 'http://a9.com/-/spec/opensearch/1.1/'}
DI = ('(all:"directed information" OR all:"directed mutual information" '
      'OR all:"causal rate distortion" OR all:"causal rate-distortion")')
RATE = ('(all:"rate-distortion" OR all:"rate distortion" OR all:"data rate" '
        'OR all:"communication rate" OR all:"asymptote")')
CTRL = ('(all:LQG OR all:LQR OR all:"linear quadratic" OR all:"quadratic cost" '
        'OR all:"mean-square stability")')

QUERIES = [
    ('A1', '地板词形直查（最便宜的反占探针）',
     'all:"critical data rate" OR all:"minimum data rate" OR all:"rate asymptote"'),
    ('A2', '地板 × 代价/控制（近地板发散若发表应在此）',
     '(all:"critical data rate" OR all:"minimum data rate" OR all:"horizontal asymptote") '
     'AND (all:cost OR all:LQG OR all:"linear quadratic" OR all:"control")'),
    ('A3', '有向信息 × 控制轴（我车道的主对撞线，换发散口径）',
     f'{DI} AND {CTRL}'),
    ('A4', '有向信息 × 渐近/发散词形',
     f'{DI} AND {RATE}'),
    ('A5', '受限传感/观测 × 率失真 × 不稳植物（组合支的另一种写法）',
     '(all:"sensor selection" OR all:"partial observation" OR all:"task-based sensing" '
     'OR all:"measurement vector") AND (all:"rate-distortion" OR all:"rate distortion" '
     'OR all:"directed information") AND (all:unstable OR all:"mean-square")'),
]
FLOOR_WORDS = ['critical data rate', 'minimum data rate', 'infimum rate', 'rate asymptote',
               'horizontal asymptote', 'lowest data rate', 'fundamental lower bound on the rate']
DIV_WORDS = ['cost diverges', 'diverges', 'unbounded cost', 'blow up', 'blows up', 'infinite cost',
             'goes to infinity', 'tends to infinity', 'no finite cost', 'becomes arbitrarily large',
             'arbitrarily large cost', 'inverse of the']


def get(url):
    for a in range(4):
        try:
            with urllib.request.urlopen(url, timeout=90) as f:
                d = f.read().decode('utf-8', 'replace')
            if len(d) > 400:
                return d
        except Exception as e:
            p(f'   retry {a+1}: {type(e).__name__} {str(e)[:60]}')
        time.sleep(5.0)
    raise RuntimeError('failed 4x')


def arxiv(q, cap=400):
    hits, tot, start = [], None, 0
    while True:
        u = ('https://export.arxiv.org/api/query?search_query=' + urllib.parse.quote_plus(q)
             + f'&start={start}&max_results=200&sortBy=relevance&sortOrder=descending')
        root = ET.fromstring(get(u))
        if tot is None:
            tot = int(root.find('os:totalResults', NS).text)
        es = root.findall('a:entry', NS)
        for e in es:
            aid = re.sub(r'v\d+$', '', e.find('a:id', NS).text.split('/abs/')[-1])
            hits.append(dict(id=aid, pub=e.find('a:published', NS).text[:10],
                             title=' '.join(e.find('a:title', NS).text.split()),
                             abstract=' '.join(e.find('a:summary', NS).text.split())))
        start += len(es)
        p(f'   页 start={start-len(es)} 抓 {len(es)}，累计 {len(hits)}/{tot}')
        time.sleep(3.0)
        if not es or start >= tot or start >= cap:
            break
    return tot, hits


p('== (0) 逐查询命中（[V1] 落盘查询串与抓取完整性） ==')
union = {}
for lab, why, q in QUERIES:
    tot, hits = arxiv(q)
    ids = {h['id'] for h in hits}
    p(f'\n {lab}  {why}')
    p(f'   query = {q}')
    p(f'   totalResults = {tot}   实抓 = {len(hits)}   抓全={"是" if len(hits) >= tot else "否（截断，须声明）"}')
    p(f'   唯一 id 数 = {len(ids)}')
    for h in hits:
        if h['id'] not in union:
            union[h['id']] = dict(h, _from=[lab])
        elif lab not in union[h['id']]['_from']:
            union[h['id']]['_from'].append(lab)
    time.sleep(2.0)

p('\n== (1) 并集与本地语料的重叠 ==')
import glob
import os
local = {os.path.splitext(os.path.basename(f))[0] for f in glob.glob('papers/txt/*.txt')}
p(f' 并集唯一条数 = {len(union)}')
loc = sorted(set(union) & local)
p(f' 本地已有全文的 = {len(loc)}：{loc}')
p(f' 只到摘要级（本地无全文，[V4] 不许外推） = {len(union)-len(loc)}')

p('\n== (2) [V2] 判定：谁同时含地板词形与发散词形 ==')
v2, v3 = [], []
for hid, h in sorted(union.items()):
    t = (h['title'] + ' || ' + h['abstract']).lower()
    a = [w for w in FLOOR_WORDS if w in t]
    b = [w for w in DIV_WORDS if w in t]
    if a and b:
        v2.append((hid, h, a, b))
    elif a:
        v3.append((hid, h, a, b))
p(f' [V2] 同时命中（=已发表近地板发散律）：{len(v2)} 条')
for hid, h, a, b in v2:
    p(f'   {hid} {h["pub"]}  floor={a} diverg={b}  《{h["title"][:90]}》 来源{h["_from"]}')
p(f' [V3] 只命中地板词形（邻居，不等于发散率已发表）：{len(v3)} 条')
for hid, h, a, b in v3:
    p(f'   {hid} {h["pub"]}  floor={a}  《{h["title"][:90]}》 来源{h["_from"]}'
      f'{"  [本地有全文]" if hid in local else ""}')
if v2:
    p('\n ⇒ [V2] 非空：我的指数-1 律必须降级为复现，正文与板都要自首。')
else:
    p('\n ⇒ [V2] 空：在本次并集的**摘要级**覆盖内，没有任何一条同时陈述"地板"与"代价发散"。'
      '\n   按 [V4]：这只说明"该查询集未给出前占证据"，不说明没人做过。'
      '\n   已发表的只有地板本身：1510.04214 Cor.1（1.169 bit/sample）与 2606.31396 Th.1（R_exp，任意观测律）。')

p('\n== (3) 本地语料的定点复核（摘要级检索打不中正文，必须回全文找） ==')
p(' 在 papers/txt/1510.04214.txt / 1612.02126.txt / 2109.12246.txt / 2606.31396.txt 里找 '
  '"as rate approaches / diverges as / 1/(R−Rmin)" 型陈述：')
for f in ['1510.04214', '1612.02126', '2109.12246', '2606.31396', '2601.12782', '1606.01946']:
    path = f'papers/txt/{f}.txt'
    if not os.path.exists(path):
        p(f'   {f}: 本地无')
        continue
    t = re.sub(r'\s+', ' ', re.sub(r'[\x00-\x08\x0b-\x1f]', ' ', open(path, encoding='utf-8', errors='replace').read()))
    pat = re.compile(r'(diverg\w+|arbitrar\w+ large|blow\w*\s*up|tends to infinity|'
                     r'unbounded\w*\s+cost|as the rate|as the data rate|approaches the (?:critical|minimum))')
    ct = re.compile(r'rate|cost|distortion|information')
    n = 0
    for m in pat.finditer(t):
        w = t[max(0, m.start()-160):m.end()+160]
        if ct.search(w) and ('cost' in w or 'distortion' in w):
            n += 1
            if n <= 4:
                p(f'   {f}#{n}: …{w}…')
    p(f'   {f}: 相关句 {n} 处')
OUT.close()
print('written p0/e102_out.txt')
