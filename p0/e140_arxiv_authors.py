# -*- coding: utf-8 -*-
"""[Z4] arXiv 侧作者顺序交叉核对：把 §90-B 的 Crossref 单通道署名变成双通道。
只用 arXiv Atom API 返回的字段，不拼 URL、不猜 id；id 来自 e135c 已落盘的 externalIds。
限流纪律：请求间 sleep 6s，429/超时按 10/16/26/42s 退避，耗尽就明写"未结"，不许把空结果读成否证。"""
import re
import time
import urllib.parse
import urllib.request

UA = {'User-Agent': 'rc-lane-e140 (research citation check)'}


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
    entries = re.findall(r'<entry>(.*?)</entry>', xml, re.S)
    out = []
    for e in entries:
        t = re.search(r'<title>(.*?)</title>', e, re.S)
        a = re.findall(r'<name>(.*?)</name>', e)
        idd = re.search(r'<id>https?://arxiv.org/abs/(.*?)</id>', e, re.S)
        pub = re.search(r'<published>(.*?)</published>', e, re.S)
        out.append((idd.group(1) if idd else '?',
                    ' '.join(t.group(1).split()) if t else '?',
                    a,
                    pub.group(1)[:10] if pub else '?'))
    return out


# 1) 已知预印本 id 的：Stavrou–Ostergaard–Charalambous–Derpich, ITW 2017
IDS = ['1701.06368']
print('== [Z4a] id_list 直查 ==')
for i in IDS:
    x = get('http://export.arxiv.org/api/query?id_list=%s&max_results=1' % i)
    print('  id %s :' % i, end='')
    if not x:
        print(' RETRY-EXHAUSTED ⇒ 未结')
        continue
    for a, t, au, p in parse(x):
        print('\n    arxivid=%s  %s\n    作者(arXiv 顺序)=%s\n    提交日=%s' % (a, t, '; '.join(au), p))
    time.sleep(6)

# 2) CDC 2018 那条没有 arXiv id（OA 是 Aalto 仓库 PDF）⇒ 用标题检索看有没有预印本对应
Q = 'ti:"Asymptotic Reverse-Waterfilling Characterization of Nonanticipative Rate Distortion"'
print('\n== [Z4b] 标题检索（CDC 2018 Stavrou 等五人）==')
x = get('http://export.arxiv.org/api/query?search_query=%s&max_results=5' % urllib.parse.quote(Q))
if not x:
    print('  RETRY-EXHAUSTED ⇒ 未结')
else:
    rows = parse(x)
    print('  命中 %d 条' % len(rows))
    for a, t, au, p in rows:
        print('    arxivid=%s  %s\n      作者=%s  提交日=%s' % (a, t[:96], '; '.join(au), p))
time.sleep(6)

print('\n== [Z4c] 对账口径 ==')
print('  Crossref 侧（e139/e139b 已落盘）：')
print('    ITW 2017 = P. A. Stavrou; J. Ostergaard; C. D. Charalambous; M. Derpich（4 人）')
print('    CDC 2018 = P. A. Stavrou; T. Charalambous; C. D. Charalambous; S. Loyka; M. Skoglund（5 人）')
print('  ⇒ 只有 arXiv 侧署名与上面**逐位一致**才算双通道确认；顺序不同也要明写，不许静默取其一。')
print('EXIT=0')
