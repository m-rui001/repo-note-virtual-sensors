# -*- coding: utf-8 -*-
r"""实验 139（§89-A 的落地前置）：给两篇"机制邻居"取**权威著录字段**，供 `note.tex` 加 bibitem。

 [Z1] $10.1109$/cdc.2018.8619725（渐近反向水填 $+$ RAE）与 $10.1109$/itw.2017.8277966
      （向量高斯零延迟 RD 的上界，arXiv $1701.06368$）。
 [Z2] 只允许两个通道：Crossref `works/{DOI}`（期刊/会议记录的权威作者与卷期页）与
      arXiv `id_list=`（预印本作者与日期）。**不许猜、不许拼 URL**。
 [Z3] 判据：作者/标题/venue/年/页 任一缺失 ⇒ 该条 bibitem 不写进正文，改在板上点名"著录不全"。
 [Z4] 顺带把 e135c 已读到的 $1612.03455$（同工作 ISIT2014/TIT2017）一并取全，供后续引用。
"""
import sys
import time
import json
import urllib.parse
import urllib.request
sys.stdout.reconfigure(encoding='utf-8')

UA = {'User-Agent': 'p0-e139 bibliographic check (academic use)'}
SLEEP = 4.0


def get(url, tries=4):
    for k in range(tries):
        try:
            req = urllib.request.Request(url, headers=UA)
            with urllib.request.urlopen(req, timeout=40) as r:
                return r.read().decode('utf-8', 'replace'), None
        except Exception as e:
            code = getattr(e, 'code', None)
            if code in (429, 403, 500, 502, 503) and k < tries - 1:
                w = 10 * (1.6 ** k)
                print('    !! HTTPError %s ⇒ 退避 %ds' % (code, w))
                time.sleep(w)
                continue
            return None, '%s %s' % (type(e).__name__, code or e)
    return None, 'exhausted'


print('== [Z2] Crossref works/{DOI} ==')
CR = {}
for doi in ('10.1109/cdc.2018.8619725', '10.1109/itw.2017.8277966',
            '10.1109/isit.2014.6874973', '10.1109/tit.2017.2694015'):
    body, err = get('https://api.crossref.org/works/' + urllib.parse.quote(doi))
    time.sleep(SLEEP)
    if err:
        print('  !! %-36s %s' % (doi, err))
        continue
    j = json.loads(body)['message']
    au = j.get('author', [])
    fam = lambda a: a.get('family', '?')
    giv = lambda a: (a.get('given', '') or '').replace(' ', '')
    CR[doi] = j
    print('  %-36s %s' % (doi, (j.get('title') or ['?'])[0][:78]))
    print('     作者 %d 人: %s' % (len(au), '; '.join('%s %s' % (giv(a), fam(a)) for a in au)))
    print('     venue=%s  年=%s  页=%s  出版者=%s' % (
        (j.get('container-title') or ['—'])[0][:52],
        (j.get('published-print') or j.get('published-online') or {}).get('date-parts', [['?']])[0][0],
        j.get('page', '—'), j.get('publisher', '?')))
    print('     event=%s  DOI=%s' % ((j.get('event') or {}).get('name', '—'), j.get('DOI')))

print('\n== [Z2] arXid id_list 取作者（预印本侧，用于交叉核对署名顺序）==')
ARX = {}
for aid, lab in (('1701.06368', 'itw2017/upper-bound'), ('1612.03455', 'sidinfo'),
                 ('1802.08376', 'sanity: tzoumas 已知')):
    body, err = get('http://export.arxiv.org/api/query?id_list=%s&max_results=1' % aid)
    time.sleep(SLEEP)
    if err:
        print('  !! %-12s %s' % (aid, err))
        continue
    import re
    t = re.search(r'<title>(.*?)</title>', body.split('<entry>')[1], re.S).group(1).strip()
    aus = re.findall(r'<name>(.*?)</name>', body.split('<entry>')[1], re.S)
    pub = re.search(r'<published>(\d{4})', body.split('<entry>')[1]).group(1)
    ARX[aid] = (t, aus, pub)
    print('  %s (%s, %s)  %s' % (aid, pub, lab, t[:70]))
    print('     作者: %s' % '; '.join(aus))

print('\n== [Z3] 著录完整性判决 ==')
need = ('10.1109/cdc.2018.8619725', '10.1109/itw.2017.8277966')
ok = []
for doi in need:
    j = CR.get(doi)
    if not j:
        print('  %-36s 通道未返回 ⇒ **bibitem 不写**，板上点名' % doi)
        continue
    have = bool((j.get('title') or [''])[0]) and bool(j.get('author')) and bool(j.get('DOI'))
    print('  %-36s 标题√作者√DOI√ = %s ⇒ %s' % (doi, have, '可写' if have else '著录不全'))
    if have:
        ok.append(doi)
print('  可写条数 = %d / %d' % (len(ok), len(need)))
print('  arXiv 侧作者与 Crossref 侧是否一致：')
for doi, aid in (('10.1109/itw.2017.8277966', '1701.06368'),):
    j = CR.get(doi)
    if j and aid in ARX:
        crf = [a.get('family', '') for a in j.get('author', [])]
        axf = [s.split()[-1] for s in ARX[aid][1]]
        print('    %s vs %s: %s（CR=%s / AX=%s）' % (doi, aid, crf == axf, crf, axf))
print('EXIT=0')
