# -*- coding: utf-8 -*-
"""[Z6] e140b 的松查询捞出三条**类内新邻居**（都在 Stavrou 这条线上），逐条取摘要＋三问卡。
只报 arXiv Atom API 原样返回的字段；三问是关键词级**候选**判读，最终措辞由我按摘要原文定。
三问卡（与 e138/e138b 同一张，便于横向比）：
  a 控制代价（LQG/估计误差型泛函）是否给出闭式？
  b 指数/斜率是否随秩或维数变化？
  c 传感/观测矩阵是否是决策变量？"""
import re
import time
import urllib.request

UA = {'User-Agent': 'rc-lane-e141 (research citation check)'}
IDS = ['1810.00298', '1912.07640', '1603.04172']

KW = {
    'a': ['closed-form', 'closed form', 'explicit', 'computable', 'characterization', 'single-letter',
          'reverse-waterfilling', 'reverse waterfilling', 'riccati'],
    'b': ['dimension', 'vector-valued', 'rank', 'number of unstable', 'state dimension', 'per-mode',
          'multiplicity', 'eigen'],
    'c': ['sensor', 'measurement matrix', 'sensing', 'output matrix', 'design of the', 'joint source',
          'encoder', 'estimator', 'controller'],
}


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


xml = get('http://export.arxiv.org/api/query?id_list=%s&max_results=5' % ','.join(IDS))
if not xml:
    print('  RETRY-EXHAUSTED ⇒ 三条全未结，不许当已排除')
    raise SystemExit('EXIT=2')

entries = re.findall(r'<entry>(.*?)</entry>', xml, re.S)
print('== [Z6] 取回 %d/%d 条 ==' % (len(entries), len(IDS)))
if len(entries) != len(IDS):
    print('  !! 返回 %d 条 != 请求 %d 条 ⇒ 缺的条目算未结，逐条点名' % (len(entries), len(IDS)))
for e in entries:
    idd = re.search(r'<id>https?://arxiv.org/abs/(.*?)</id>', e, re.S).group(1)
    t = ' '.join(re.search(r'<title>(.*?)</title>', e, re.S).group(1).split())
    au = [x for x in re.findall(r'<name>(.*?)</name>', e)]
    ab = ' '.join(re.search(r'<summary>(.*?)</summary>', e, re.S).group(1).split())
    print('\n  [%s] %s' % (idd, t))
    print('    作者：%s' % '; '.join(au))
    print('    摘要全文（原样）：%s' % ab)
    low = ab.lower()
    hits = {k: [w for w in v if w in low] for k, v in KW.items()}
    print('    三问候选：a=%s  b=%s  c=%s' % (
        '有' if hits['a'] else '无', '有' if hits['b'] else '无', '有' if hits['c'] else '无'))
    for k in ('a', 'b', 'c'):
        print('      %s 命中词：%s' % (k, hits[k] if hits[k] else '（无）'))
    time.sleep(6)

print('\n== [Z6] 汇总 ==')
print('  这三条是 e138 那 16 条之外的**新增类内候选**（松查询捞到的），'
      '进未读账必须点名；三问里 a$\\wedge$b 同时为"有"才构成撞车，逐条人工复核后再定性。')
print('EXIT=0')
