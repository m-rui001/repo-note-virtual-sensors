# -*- coding: utf-8 -*-
"""e139b（§89-A 落地前置）：只取 CDC-2018 / ITW-2017 两条的**完整标题**。
e139 的打印把标题截到 78 字符，不许据此写 `bibitem`；本脚本落盘完整值。
通道：Crossref `works/{DOI}`（版本记录），单请求间隔 4s，429 退避 8/16/24/32s。"""
import sys
import time
import json
import urllib.request
sys.stdout.reconfigure(encoding='utf-8')
UA = {'User-Agent': 'p0-e139b full-title fetch (academic use)'}
for doi in ('10.1109/cdc.2018.8619725', '10.1109/itw.2017.8277966'):
    j = None
    for k in range(4):
        try:
            with urllib.request.urlopen(urllib.request.Request(
                    'https://api.crossref.org/works/' + doi, headers=UA), timeout=40) as r:
                j = json.loads(r.read().decode('utf-8', 'replace'))['message']
                break
        except Exception as e:
            w = 8 * (k + 1)
            print('  !! %s 退避 %ds' % (getattr(e, 'code', e), w))
            time.sleep(w)
    time.sleep(4)
    if j:
        print('%s\n  完整标题: %s\n  页 %s / 年 %s / venue %s' % (
            doi, j['title'][0], j.get('page', '—'),
            (j.get('published-print') or j.get('published-online') or {}).get('date-parts', [['?']])[0][0],
            (j.get('container-title') or ['—'])[0]))
print('EXIT=0')
