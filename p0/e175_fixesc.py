# -*- coding: utf-8 -*-
"""把 e175 代码区（docstring 之外）里 LaTeX 的单反斜杠补成双反斜杠。
   保护：已双写的 `\\\\xx` 不动；字符串字面量 `'\\n'`（真换行）不动；含 r'...' 的整行跳过。
   先 compile 通过再写盘。"""
import re
import sys

path = sys.argv[1]
src = open(path, encoding='utf-8').read()
parts = src.split('"""', 2)
assert len(parts) == 3, 'docstring 分隔符数目不对，不写'
head, doc, rest = parts
TOK = '\x00NL\x00'
pat = re.compile(r'(?<!\\)\\(?=[A-Za-z])')
fixed = []
tot = 0
for ln in rest.split('\n'):
    if "r'" in ln or 'r"' in ln:
        fixed.append(ln)
        continue
    pr = ln.replace("'\\n'", TOK)
    new, k = pat.subn(lambda m: '\\\\', pr)
    tot += k
    fixed.append(new.replace(TOK, "'\\n'"))
new_src = head + '"""' + doc + '"""' + '\n'.join(fixed)
import warnings
with warnings.catch_warnings():
    warnings.simplefilter('error')
    compile(new_src, path, 'exec')
open(path, 'w', encoding='utf-8').write(new_src)
print('补了 %d 处单反斜杠；clean compile（-W error）通过' % tot)
