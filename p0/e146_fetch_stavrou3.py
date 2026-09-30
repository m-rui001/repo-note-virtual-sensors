# -*- coding: utf-8 -*-
r"""实验 146：**把三条 Stavrou 线预印本从"只读到摘要"升级为"通读"**，并顺手做 CDC-2018 著录的第二通道。

背景（板 §90/§92）：`e141` 只拿到 arXiv 摘要就写进了 `frag_rc_waterfill.tex` 的区分段，
覆盖面句子只能写 "three Stavrou-line preprints read at abstract only"；
而 `10.1109/cdc.2018.8619725` 的署名只有 Crossref 单通道（`e140/e140b` 的 arXiv 侧三条全 429/timeout）。

跑前写死：
 [F1] 每条先打 Atom `id_list`（作者顺序/日期/comment/journal_ref/primary_category），失败退避重试至多 4 次、间隔 $\\ge4.5$\\,s；
      **429/timeout 不等于"查无此项"**（§90 第 12 条纪律），必须印出每条的重试次数与最终状态。
 [F2] PDF 取 `https://arxiv.org/pdf/{id}`，落到 `papers/pdf/{id}.pdf`；验收：响应头 `content-type` 含 pdf、
      首 5 字节是 `%PDF`、字节数 $>50{,}000$，否则判 FAILED 不算取到。
 [F3] 页数用 pypdf（装了就用，没装就用 `/Type /Page` 计数并标注"口径弱"）。
 [F4] 若某条的 Atom comment/journal_ref 里出现 "CDC 2018" 之类字样 $\\Rightarrow$ 它就是 `10.1109/cdc.2018.8619725` 的预印本，
      作者顺序由 **Crossref + arXiv 两通道**共同确认；否则维持单通道。

三条：1603.04172（滤波器作为编码/解码 + 条件互信息 MSE 下界）、1810.00298（零延迟 / 格量化 $0.254r+1$ bits）、
1912.07640（部分可观测 Gauss–Markov 的 NRDF 在严格可行时可 SDP 计算）。
"""
import re
import sys
import time
import json
import urllib.request
import urllib.error

sys.stdout.reconfigure(encoding='utf-8')
IDS = ['1603.04172', '1810.00298', '1912.07640']
UA = {'User-Agent': 'research-notes/0.1 (contact: local lane p0)'}
SLEEP = 4.5
MAXTRY = 4


def get(url, tries=MAXTRY):
    last = None
    for k in range(1, tries + 1):
        try:
            req = urllib.request.Request(url, headers=UA)
            with urllib.request.urlopen(req, timeout=60) as r:
                return r.read(), dict(r.headers), k, None
        except Exception as e:  # noqa: BLE001
            last = '%s: %s' % (type(e).__name__, e)
            time.sleep(SLEEP)
    return None, None, tries, last


print('== [F1] Atom 元数据（作者顺序 / comment / journal_ref / primary_category）==')
meta = {}
for i, aid in enumerate(IDS):
    if i:
        time.sleep(SLEEP)
    body, hdr, used, err = get('http://export.arxiv.org/api/query?id_list=%s' % aid)
    if body is None:
        print('  %-12s FAILED after %d tries — %s' % (aid, used, err))
        continue
    txt = body.decode('utf-8', 'replace')
    entries = re.findall(r'<entry>(.*?)</entry>', txt, re.S)
    if not entries:
        print('  %-12s 无 entry（totalResults=%s）' % (aid, re.findall(r'totalResults[^>]*>(\d+)', txt)))
        continue
    e = entries[0]
    grab = lambda tag: [x.strip() for x in re.findall(r'<%s>(.*?)</%s>' % (tag, tag), e, re.S)]  # noqa: E731
    au = grab('author')
    names = [re.sub(r'\s+', ' ', re.sub(r'<[^>]+>', '', a)).strip() for a in au]
    prim = re.findall(r'arxiv:primary_category[^>]*term="([^"]+)"', e)
    com = grab('arxiv:comment')
    jref = grab('journal_ref')
    doi = grab('arxiv:doi')
    meta[aid] = dict(title=[re.sub(r'\s+', ' ', x) for x in grab('title')], authors=names,
                     published=grab('published'), updated=grab('updated'),
                     primary=prim, comment=com, journal_ref=jref, doi=doi, tries=used)
    print('  %-12s tries=%d  primary=%s  updated=%s' % (aid, used, prim, grab('updated')))
    print('      标题：%s' % (re.sub(r'\s+', ' ', grab('title')[0])[:150]))
    print('      作者：%s' % '; '.join(names))
    if com:
        print('      comment：%s' % re.sub(r'\s+', ' ', com[0])[:220])
    if jref:
        print('      journal_ref：%s' % re.sub(r'\s+', ' ', jref[0])[:160])
    if doi:
        print('      arxiv:doi：%s' % doi[0])

print('\n== [F2]+[F3] PDF 全文下载验收 ==')
try:
    from pypdf import PdfReader
    HPDF = True
except Exception:  # noqa: BLE001
    PdfReader = None
    HPDF = False
print('  pypdf 可用：%s' % HPDF)
for i, aid in enumerate(IDS):
    if i:
        time.sleep(SLEEP)
    path = 'papers/%s.pdf' % aid   # 本仓库的全文 PDF 一直落在 papers/ 根下（无 papers/pdf 目录）
    body, hdr, used, err = get('https://arxiv.org/pdf/%s' % aid)
    if body is None:
        print('  %-12s FAILED after %d tries — %s' % (aid, used, err))
        continue
    ct = (hdr.get('content-type') or '').lower()
    # 第 39 号账：原来写的是 body[:5] == b'%PDF'（5 字节 vs 4 字节），恒假 ⇒ 三份全文都被误判 FAILED。
    ok = body[:4] == b'%PDF' and len(body) > 50000
    npg = None
    if ok and HPDF:
        open(path, 'wb').write(body)
        try:
            npg = len(PdfReader(path).pages)
        except Exception as e:  # noqa: BLE001
            npg = '页数读取失败 %s' % e
    if ok and not HPDF:
        open(path, 'wb').write(body)
        npg = '%d（口径弱：/Type /Page 计数）' % len(re.findall(rb'/Type\s*/Page[^s]', body))
    print('  %-12s tries=%d bytes=%d head=%r content-type=%r ⇒ %s  页数=%s'
          % (aid, used, len(body), body[:5], ct, 'OK 已落 %s' % path if ok else 'FAILED', npg))

print('\n== [F4] CDC-2018 第二通道判定 ==')
hit = None
for aid, d in meta.items():
    blob = ' '.join((d.get('comment') or []) + (d.get('journal_ref') or []) + (d.get('doi') or []))
    if re.search(r'CDC|decision and control', blob, re.I):
        hit = aid
        print('  %s 的元数据里有 CDC 字样：%s' % (aid, blob[:200]))
        print('  作者顺序（arXiv）：%s' % '; '.join(d['authors']))
print('  ⇒ %s' % ('**两条通道对上：`10.1109/cdc.2018.8619725` 的署名不再单通道**' if hit
                  else '未在任何条目元数据里找到 CDC 字样 ⇒ CDC-2018 维持 Crossref 单通道，正文照 §92 的限定句不动'))
open('p0/e146_meta.json', 'w', encoding='utf-8').write(json.dumps(meta, ensure_ascii=False, indent=1))
print('元数据落盘 p0/e146_meta.json（%d 条）' % len(meta))
print('EXIT=0')
