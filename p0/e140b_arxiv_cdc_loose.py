# -*- coding: utf-8 -*-
"""[Z5] CDC 2018 那条的 arXiv 侧兜底检索：e140 的 ti: 精确短语命中 0 条，
**0 条不等于不存在**（大小写、连字符、短语切分都能让 ti: 落空）。这里换三条更松的查询，
并把每条的原始命中数印出来；只有三条都空，才允许写"arXiv 侧无对应预印本（三条查询）"。"""
import re
import time
import urllib.parse
import urllib.request

UA = {'User-Agent': 'rc-lane-e140b (research citation check)'}


def get(url, tries=(0, 10, 16, 26, 42)):
    for k, wait in enumerate(tries):
        if wait:
            print('    !! 退避 %ds (第 %d 次)' % (wait, k))
            time.sleep(wait)
        try:
            return urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=45).read().decode('utf-8', 'replace')
        except Exception as ex:
            print('    !! %s %s' % (type(ex).__name__, str(ex)[:70]))
    return None


def parse(xml):
    tot = re.search(r'<opensearch:totalResults[^>]*>(\d+)</opensearch:totalResults>', xml)
    rows = []
    for e in re.findall(r'<entry>(.*?)</entry>', xml, re.S):
        t = re.search(r'<title>(.*?)</title>', e, re.S)
        a = re.findall(r'<name>(.*?)</name>', e)
        idd = re.search(r'<id>https?://arxiv.org/abs/(.*?)</id>', e, re.S)
        pub = re.search(r'<published>(.*?)</published>', e, re.S)
        rows.append((idd.group(1) if idd else '?', ' '.join(t.group(1).split()) if t else '?',
                     a, pub.group(1)[:10] if pub else '?'))
    return (int(tot.group(1)) if tot else -1), rows


QS = [
    ('au:Stavrou AND au:Loyka', '五人组合里两个独特姓氏'),
    ('au:Stavrou AND all:"nonanticipative rate distortion"', '同作者 + 主题词'),
    ('all:"reverse waterfilling" AND all:"Gauss-Markov" AND all:"nonanticipative"', '纯主题词，不含作者'),
]
print('== [Z5] 三条松查询（每条都印 totalResults，不只印前几条）==')
hit_any = 0
for q, why in QS:
    x = get('http://export.arxiv.org/api/query?search_query=%s&max_results=8' % urllib.parse.quote(q))
    print('\n  query: %s   [%s]' % (q, why))
    if not x:
        print('    RETRY-EXHAUSTED ⇒ 该条**未结**，不许计入"空"')
        hit_any += 1
        continue
    tot, rows = parse(x)
    print('    totalResults=%d  （返回 %d 条）' % (tot, len(rows)))
    for a, t, au, p in rows:
        print('      arxivid=%s  %s\n        作者=%s  提交日=%s' % (a, t[:100], '; '.join(au), p))
    time.sleep(6)

print('\n== [Z5] 判决 ==')
print('  若上面出现 Stavrou/Loyka/Charalambous/Skoglund 的组合 ⇒ arXiv 侧**有**预印本，'
      '把它的作者顺序与 Crossref 五人逐位对表。')
print('  若三条都 totalResults=0 ⇒ 只能写"arXiv 侧未检索到（三条查询，见 e140b stdout）"，'
      '署名口径维持 Crossref 单通道＋明写未结。')
print('EXIT=0')
