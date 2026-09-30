#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
e155：核我自己正文的两处改动（bibitem 升级 + 那句"他们的阈值断在什么对象上"），
并把 note.tex 重新编译的读数固定到 p0/e155_out.txt，供板 §102 引用。
判据先写死：
  [V1] note.tex 里不再出现 arXiv:1603.04172 作为出版口径；期刊卷/页/年/doi 四项齐。
  [V2] 新增限定句存在，且包含 "eigendirections"（阈值断在误差协方差的特征方向上，不是给定权阵）。
  [V3] 编译门槛：错误行数 == 0 且页数 <= 7。
"""
import io, re, subprocess, sys
sys.stdout.reconfigure(encoding='utf-8')

T = 'p0/note/note.tex'
src = io.open(T, encoding='utf-8', errors='replace').read()
out = []


def p(*a):
    s = ' '.join(str(x) for x in a)
    out.append(s)
    print(s)


bib = src[src.find('\\bibitem{stavrou2016estimation}'):src.find('\\bibitem{stavrou2018zerodelay}')]
p('== [V1] stavrou2016estimation 的 bibitem 现状 ==')
p(bib.strip())
need = ['SIAM J. Control Optim.', 'vol.~56', 'pp.~3731--3765', '2018', 'doi:10.1137/17m1116349']
miss = [k for k in need if k not in bib]
p('[V1] 缺项=%s ；仍写 arXiv 口径？ %s' % (miss or '无', 'arXiv:1603.04172 [cs.IT]' in bib))

p('')
p('== [V2] 新增限定句 ==')
m = re.search(r'pricing no control action;.{0,400}', src, re.S)
p((m.group(0).replace('\n', ' ')[:400] if m else '未找到'))
p('[V2] 含 "eigendirections"？ %s ；含 "fully observed"？ %s'
  % (bool(m and 'eigendirections' in m.group(0)), bool(m and 'fully observed' in src)))

p('')
p('== [V3] 编译（pdflatex 两遍）==')
for i in (1, 2):
    subprocess.run(['pdflatex', '-interaction=nonstopmode', 'note.tex'], cwd='p0/note', capture_output=True)
log = io.open('p0/note/note.log', encoding='utf-8', errors='replace').read()
errs = [l for l in log.split('\n') if l.startswith('!')]
mm = re.search(r'Output written on note.pdf \((\d+) pages, (\d+) bytes\)', log)
p('错误行数 = %d' % len(errs))
p('产出      = %s' % (mm.group(0) if mm else '未产出'))
pages = int(mm.group(1)) if mm else -1
size = int(mm.group(2)) if mm else -1
p('[V3] 0 错误且页数 <=7 ？ %s（页 %d，字节 %d）' % (len(errs) == 0 and 0 < pages <= 7, pages, size))

io.open('p0/e155_out.txt', 'w', encoding='utf-8').write('\n'.join(out) + '\n')
