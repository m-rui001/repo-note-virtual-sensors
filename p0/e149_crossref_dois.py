# -*- coding: utf-8 -*-
r"""实验 149：Crossref 逐字段现取三条 DOI，把 §96 要用的期刊元数据**落盘**（不许只在终端里过一遍）。

三条：
 [X1] 10.1109/JSTSP.2018.2855046 —— e146 从 arXiv `arxiv:doi` 发现的期刊版；要印 vol/issue/page/作者序/ISSN，
      与 `note.tex` 里 `stavrou2018zerodelay` 新 bibitem 逐字段对账。
 [X2] 10.1109/ciss.2016.7460485、[X3] 10.1109/ssp.2005.1628751 —— §92 记为"无公开摘要的真欠"两条；
      Crossref 的标题/作者/容器即使没有摘要也算**部分清偿**，必须分开写"有元数据 / 无摘要"。

退避是通道正确性条件（纪律 12）：SLEEP=2.0s，重试 4 次，每次印 HTTP 状态；429/超时**绝不写成"查无此项"**。
判定分支预注册：
 [B1] 200 且有 title ⇒ 打印全部字段（"元数据可得"）。
 [B2] 404 ⇒ 打印"Crossref 无此 DOI"（这才允许写"无记录"）。
 [B3] 其它（429/5xx/超时/异常）⇒ 打印状态并标 **通道未完**，本轮不得据此写任何"无先例/无记录"句子。
"""
import json
import sys
import time
import urllib.error
import urllib.request

sys.stdout.reconfigure(encoding='utf-8')
SLEEP = 2.0
HDR = {'User-Agent': 'p0-e149/1.0 (mailto:research.audit@example.org)',
       'Accept': 'application/json'}
DOIS = [('X1', '10.1109/JSTSP.2018.2855046'),
        ('X2', '10.1109/ciss.2016.7460485'),
        ('X3', '10.1109/ssp.2005.1628751')]

for tag, doi in DOIS:
    url = 'https://api.crossref.org/works/' + doi
    got, status, tries = None, None, 0
    for tries in range(1, 5):
        try:
            req = urllib.request.Request(url, headers=HDR)
            with urllib.request.urlopen(req, timeout=40) as r:
                status = r.getcode()
                got = json.loads(r.read().decode('utf-8'))
            break
        except urllib.error.HTTPError as e:
            status = e.code
            if e.code in (404,):
                got = None
                break
            time.sleep(SLEEP * (2 ** (tries - 1)))
        except Exception as e:
            status = 'ERR:%s' % type(e).__name__
            time.sleep(SLEEP * (2 ** (tries - 1)))
    time.sleep(SLEEP)
    print('== [%s] %s  tries=%d  HTTP=%s' % (tag, doi, tries, status))
    if status == 200 and got and got.get('message', {}).get('title'):
        m = got['message']
        print('  标题：%s' % m['title'][0])
        print('  容器：%s' % (m.get('container-title') or ['(none)'])[0])
        print('  事件：%s' % json.dumps(m.get('event', {}), ensure_ascii=False))
        print('  volume=%s  issue=%s  page=%s' % (m.get('volume'), m.get('issue'), m.get('page')))
        for k in ('published-print', 'published-online', 'issued', 'created'):
            if k in m:
                print('  %s：%s' % (k, m[k].get('date-parts')))
        print('  ISSN=%s' % m.get('ISSN'))
        print('  作者序：%s' % '; '.join('%s %s' % (a.get('given', '?'), a.get('family', '?'))
                                          for a in m.get('author', [])))
        ab = m.get('abstract')
        print('  摘要：%s' % ('(Crossref 未提供摘要)' if not ab else ab[:300].replace('\n', ' ')))
        print('  引用数=%s  被引列表条数=%s' % (m.get('is-referenced-by-count'), len(m.get('reference', []))))
        print('  ⇒ 分支 [B1] 元数据可得')
    elif status == 404:
        print('  ⇒ 分支 [B2] Crossref 无此 DOI（仅此情形允许写"无记录"）')
    else:
        print('  ⇒ 分支 [B3] **通道未完**（HTTP=%s）：本轮不得据此写"无先例/无记录"' % status)
print('EXIT=0')
