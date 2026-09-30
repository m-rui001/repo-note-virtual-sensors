# -*- coding: utf-8 -*-
import sys, re
sys.stdout.reconfigure(encoding='utf-8')
p = sys.argv[1]
s = open(p, encoding='utf-8').read()
for i, l in enumerate(s.split('\n'), 1):
    for m in re.finditer(r'(?<!\\)\\[^\\\'"`abfnrtv0-9xuNN]', l):
        print(i, repr(l[max(0, m.start()-40):m.start()+40]))
