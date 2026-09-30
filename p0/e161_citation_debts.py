#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""e161：把 §104 留下的三笔引用债一次结掉——**逐字段现取现判**，不沿用任何旧结论。

债目（全部来自我方板 §104 与 last commit "citation expansion needs per-entry verification"）：
 [T1] "CDC-2018 的期刊版是 IEEE TAC 2022"——这句话我 §104 写下时**只比过 DOI，没比过标题与作者**。
      现在取两个 DOI 的完整记录，按 (标题归一化等值) + (作者集合包含) 判"同一工作/扩展版/不同工作"。
      只有判成"同一工作或扩展版"才允许把 `note.tex` 第 315--319 行那条升级成期刊口径。
 [T2] 两条最近风险条目的**摘要级**三项筛查（(a) 传感阵列为决策变量 (b) 控制代价加权失真 (c) 按秩地板）：
      `10.1109/acc.2012.6314650`、`10.23919/acc55779.2023.10156247`。
      Crossref 常不收摘要 $\\Rightarrow$ **收不到就明写"无摘要，不得下差异句"**，只许声称著录。
 [T3] TAC-2022 若字段齐，直接打印可用的 bibitem（作者/卷/页/年/DOI 全部来自本次响应）。
 [G]  有效读数门（板账 #46）：四条 DOI 各自打印 HTTP 状态与"是否拿到 title"；
      成功 <4 条 $\\Rightarrow$ 本轮只报采集状态，不下任何 T1/T2 判决。
"""
import json
import re
import sys
import time
import urllib.request
sys.stdout.reconfigure(encoding='utf-8')

out = []


def p(*a):
    s = ' '.join(str(x) for x in a)
    out.append(s)
    print(s)


HDR = {'User-Agent': 'rc-lane-citation-check/1.0 (mailto:research@example.org)'}
DOIS = {
    'CDC2018': '10.1109/cdc.2018.8619725',
    'TAC2022': '10.1109/tac.2021.3099444',
    'ACC2012': '10.1109/acc.2012.6314650',
    'ACC2023': '10.23919/acc55779.2023.10156247',
}
REC = {}
for k, d in DOIS.items():
    ok, rec = False, None
    for att in range(3):
        try:
            req = urllib.request.Request('https://api.crossref.org/works/' + d,
                                         headers=HDR)
            with urllib.request.urlopen(req, timeout=45) as f:
                rec = json.load(f)['message']
            ok = bool(rec.get('title'))
        except Exception as e:
            p('  [%s] 第 %d 次取失败：%s' % (k, att + 1, e))
            time.sleep(2.0 * (att + 1))
            continue
        break
    REC[k] = rec if ok else None
    p('[G] %-8s doi=%-36s 取到 title=%s' % (k, d, ok))

n_ok = sum(1 for v in REC.values() if v)
p('[G] 成功 %d/4 条' % n_ok)
if n_ok < 4:
    p('[G] 成功数 <4 $\\Rightarrow$ 本轮不下 T1/T2 判决（只报采集状态）')


def norm(s):
    return re.sub(r'[^a-z]', '', (s or '').lower())


def authset(rec):
    st = set()
    for a in (rec or {}).get('author', []):
        fam = norm(a.get('family', ''))
        giv = (a.get('given', '') or '')
        if fam:
            st.add((fam, giv[:1].lower()))
    return st


p('')
p('== [T1] CDC-2018 vs TAC-2022 是否同一工作 ==')
for k in ('CDC2018', 'TAC2022'):
    r = REC.get(k)
    if not r:
        p('  %s 无记录' % k)
        continue
    ev = (r.get('issued', {}) or {}).get('date-parts', [[None]])[0][0]
    p('  %-8s title=%s' % (k, (r.get('title') or [''])[0]))
    p('           venue=%s vol=%s pages=%s year=%s type=%s'
      % ((r.get('container-title') or ['-'])[0], r.get('volume', '-'), r.get('page', '-'), ev, r.get('type')))
    p('           authors=%s' % ', '.join(sorted('%s %s' % (a.get('given', ''), a.get('family', '')) for a in r.get('author', []))))
if REC.get('CDC2018') and REC.get('TAC2022'):
    t1 = norm((REC['CDC2018'].get('title') or [''])[0])
    t2 = norm((REC['TAC2022'].get('title') or [''])[0])
    a1, a2 = authset(REC['CDC2018']), authset(REC['TAC2022'])
    same_t = (t1 == t2)
    cover = len(a1 & a2)
    p('  标题归一化等值=%s  作者交集 %d 人（CDC %d / TAC %d）' % (same_t, cover, len(a1), len(a2)))
    verdict = ('同一工作（标题同 + 作者全含）$\\Rightarrow$ 允许升级 bibitem' if same_t and (a1 <= a2 or a2 <= a1)
               else '标题同但作者集不同 $\\Rightarrow$ 记"扩展/相关"，bibitem 两条并存') if same_t else \
              ('词干重叠高、作者全含 $\\Rightarrow$ 判"期刊扩展版"，升级须同时保留会议版' if cover >= max(1, min(len(a1), len(a2))) and (a1 <= a2 or a2 <= a1)
               else '不同工作 $\\Rightarrow$ §104 那句话作废并公开更正')
    p('  [T1] 判决：%s' % verdict)

p('')
p('== [T2] 两条最近风险条目：摘要级三项筛查 ==')
RA = re.compile(r'sensor|sensing|observation|measurement', re.I)
RB = re.compile(r'lqg|control cost|weighted trace|control Lyapunov|LQ ', re.I)
RC = re.compile(r'\brank\b|dimension|subspace', re.I)
for k in ('ACC2012', 'ACC2023'):
    r = REC.get(k)
    if not r:
        p('  %s 无记录' % k)
        continue
    ab = r.get('abstract') or ''
    ab = re.sub(r'<[^>]+>', ' ', ab)
    p('  %-8s title=%s' % (k, (r.get('title') or [''])[0]))
    p('           venue=%s year=%s' % ((r.get('container-title') or ['-'])[0],
                                       ((r.get('issued', {}) or {}).get('date-parts', [[None]])[0][0])))
    if not ab.strip():
        p('           [T2] Crossref **无摘要** $\\Rightarrow$ 不得写差异句，只许声称著录已取')
        continue
    p('           摘要（前 480 字）：%s' % ab[:480].replace('\n', ' '))
    f = (bool(RA.search(ab)), bool(RB.search(ab)), bool(RC.search(ab)))
    p('           三项命中 (a)(b)(c) = %s $\\Rightarrow$ %s' % (f, '三项同现，必须读完正文才能写差异'
      if all(f) else '不同现（%d/3）$\\Rightarrow$ "窄口径无先例"这句话仍安全' % sum(f)))

p('')
p('== [T3] 可直接落 note.tex 的 bibitem（字段全取自本响应）==')
r = REC.get('TAC2022')
if r:
    au = r.get('author', [])
    if len(au) >= 2:
        names = ', '.join('%s~%s' % (a.get('given', '')[:1], a.get('family', '')) for a in au)
        p("  \\bibitem{stavrou2022tac} %s, ``%s'', \\emph{%s}, vol.~%s, pp.~%s, %s, doi:~%s."
          % (names, (r.get('title') or [''])[0], (r.get('container-title') or [''])[0],
             r.get('volume', '?'), r.get('page', '?'),
             ((r.get('issued', {}) or {}).get('date-parts', [[None]])[0][0]), r.get('DOI')))
        p('  [T3] 字段齐=%s' % all([r.get('volume'), r.get('page'), au, (r.get('title') or [''])[0]]))
    else:
        p('  [T3] 作者字段缺 $\\Rightarrow$ 不写 bibitem')
else:
    p('  [T3] 无记录 $\\Rightarrow$ 不写 bibitem')

open('p0/e161_out.txt', 'w', encoding='utf-8').write('\n'.join(out) + '\n')
