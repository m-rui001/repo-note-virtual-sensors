#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
e156：读 1603.04172 的 Theorem 5.2 与 time-space reverse-waterfilling (5.8)-(5.11) 全块。
本轮只问三件有明确判据的事（判据在打印前写死）：
 [T1] 水位 $\\xi$ 是否带时间下标（即逐时刻各自的水平，还是一个全局水平）；
 [T2] 该式的"地板"是否由 $\\mathrm{tr}(\\Pi_{t|t-1})$ 这类**给定/递推数据**给出，
      还是含传感设计变量的极值（本站 $D_{\\rm floor}(r)$ 是后者）；
 [T3] 目标泛函是估计 MSE（无权迹）还是带权代价 $\\mathrm{tr}(\\Theta\\cdot)$。
"""
import io, re, sys
sys.stdout.reconfigure(encoding='utf-8')

L = io.open('papers/notes/1603.04172.txt', encoding='utf-8', errors='replace').read().replace('\x00', '').split('\n')
out = []


def p(*a):
    s = ' '.join(str(x) for x in a)
    out.append(s)
    print(s)


p('== 定位 ==')
for pat in (r'Theorem 5\.2', r'\(5\.8\)', r'\(5\.9\)', r'\(5\.10\)', r'\(5\.11\)', r'reverse-water'):
    ls = [i + 1 for i, s in enumerate(L) if re.search(pat, s, re.I)]
    p('%-18s 命中行：%s' % (pat, ls[:14]))

i0 = [i for i, s in enumerate(L) if 'Theorem 5.2' in s]
if i0:
    a = i0[0]
    p('')
    p('== [T0] Theorem 5.2 起（第 %d 行）往后 75 行 ==' % (a + 1))
    for i in range(a, min(a + 75, len(L))):
        s = L[i].rstrip()
        if s.strip():
            p('%5d| %s' % (i + 1, s[:170]))

p('')
p('== [T1] 出现 "xi"（水位）的行，看是否带时间下标 ==')
for i, s in enumerate(L):
    if re.search(r'\bxi\b|ξ', s) and 500 < i < len(L):
        p('%5d| %s' % (i + 1, s.strip()[:170]))

p('')
p('== [T2]/[T3] 目标泛函与地板：搜 trace/MSE/水位相关句 ==')
PAT = re.compile(r'(min distortion|achievable distortion|trace of|is chosen such that|reverse-water|unweighted)', re.I)
for i, s in enumerate(L):
    if PAT.search(s):
        p('%5d| %s' % (i + 1, s.strip()[:170]))

io.open('p0/e156_out.txt', 'w', encoding='utf-8').write('\n'.join(out) + '\n')
