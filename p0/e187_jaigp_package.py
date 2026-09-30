"""Flatten the manuscript into one self-contained .tex for the JAIGP LaTeX
upload, then reconcile every cross-reference against the compiled aux/bbl state.
Nothing here re-derives a scientific number; it audits the artefact."""
import io
import os
import re
import subprocess
import sys
import zipfile

HERE = os.path.dirname(os.path.abspath(__file__))
PRE = os.path.join(os.path.dirname(HERE), 'preprint')
OUT = []


def p(line=''):
    OUT.append(line)
    io.open(os.path.join(HERE, 'e187_out.txt'), 'w', encoding='utf-8',
            newline='\r\n').write('\n'.join(OUT) + '\n')


def read(path):
    return io.open(path, encoding='utf-8').read()


def count(text, pat):
    return len(re.findall(pat, text))


p('== [P1] 扁平化（JAIGP 可选 LaTeX 源：单文件最稳） ==')
main = read(os.path.join(PRE, 'preprint.tex'))
parts = []
for line in main.split('\n'):
    m = re.match(r'\s*\\input\{([^}]+)\}\s*$', line)
    if m:
        frag = read(os.path.join(PRE, m.group(1) + '.tex'))
        parts.append('%% ---- begin %s.tex ----' % m.group(1))
        parts.append(frag.rstrip('\n'))
        parts.append('%% ---- end %s.tex ----' % m.group(1))
    else:
        parts.append(line)
flat = '\n'.join(parts)
assert '\\input{' not in flat and '\\include{' not in flat, '扁平化后仍有 \\input'
assert flat.count('\\begin{document}') == 1 and flat.count('\\end{document}') == 1
assert flat.count('\\begin{abstract}') == 1 and flat.count('\\maketitle') == 1
flat_tex = os.path.join(PRE, 'jaigp_flat.tex')
io.open(flat_tex, 'w', encoding='utf-8', newline='\n').write(flat)
p('   源件 %d 个（主文件 + \\input 目标），扁平文件 %s = %d B / %d 行'
  % (1 + count(main, r'\\input\{'), 'jaigp_flat.tex',
     os.path.getsize(flat_tex), flat.count('\n') + 1))

# ---------------------------------------------------------------- 编译扁平版
env = dict(os.environ)
r = subprocess.run(['pdflatex', '-interaction=nonstopmode', 'jaigp_flat'],
                   cwd=PRE, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT,
                   env=env)
r2 = subprocess.run(['pdflatex', '-interaction=nonstopmode', 'jaigp_flat'],
                    cwd=PRE, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT,
                    env=env)
lg = read(os.path.join(PRE, 'jaigp_flat.log'))
p('   pdflatex 两次退出码 %s / %s；日志中 "^!" 错误行数 = %d'
  % (r.returncode, r2.returncode, count(lg, r'(?m)^!')))
assert count(lg, r'(?m)^!') == 0, '扁平版编译有错误'
assert 'undefined' not in lg.lower(), '扁平版有未定义引用'
pgs = re.search(r'Output written on jaigp_flat\.pdf \((\d+) pages, (\d+) bytes\)', lg)
assert pgs, '读不到扁平版页数'
pg_main = re.search(r'Output written on preprint\.pdf \((\d+) pages, (\d+) bytes\)',
                    read(os.path.join(PRE, 'preprint.log')))
assert pg_main, '读不到主版页数'
p('   扁平版 %s 页 / %s B；多文件版 %s 页 / %s B ⇒ 页数 %s'
  % (pgs.group(1), pgs.group(2), pg_main.group(1), pg_main.group(2),
     '一致' if pgs.group(1) == pg_main.group(1) else '不一致'))
assert pgs.group(1) == pg_main.group(1), '扁平版与多文件版页数不一致'

# ---------------------------------------------------------------- [P2] 引用对账
p('')
p('== [P2] 交叉引用逐条对账（扁平版 = 将要上传的那一份） ==')
aux = read(os.path.join(PRE, 'jaigp_flat.aux'))
defined_lbl = set(re.findall(r'\\newlabel\{([^}]+)\}', aux))
used_lbl = set(re.findall(r'\\(?:ref|eqref|S|Cref|autoref)\{([^}]+)\}', flat))
bib_used = set()
for grp in re.findall(r'\\cite[tp]?\*?(?:\[[^\]]*\])*\{([^}]*)\}', flat):
    for k in grp.split(','):
        if k.strip():
            bib_used.add(k.strip())
bib_def = set(re.findall(r'\\bibitem\{([^}]+)\}', flat))
p('   label：定义 %d 个，正文引用 %d 个，引用但未定义 %s'
  % (len(defined_lbl), len(used_lbl), sorted(used_lbl - defined_lbl) or '0 个'))
p('   bibitem：定义 %d 条，正文引用 %d 条，引用但无条目 %s'
  % (len(bib_def), len(bib_used), sorted(bib_used - bib_def) or '0 个'))
p('   未被引用的条目 %d 条：%s' % (len(bib_def - bib_used), sorted(bib_def - bib_used) or '无'))
assert not (used_lbl - defined_lbl), '存在引用不到的 label'
assert not (bib_used - bib_def), '存在没有著录的 \\cite'

# ---------------------------------------------------------------- [P3] 草稿痕迹
p('')
p('== [P3] 投稿不该留下的痕迹（扁平版正文计数，含注释行，全部须为 0） ==')
body = flat
for m in re.finditer(r'(?m)^%.*$', body):
    body = body.replace(m.group(0), '')
checks = [('preprint（不分大小写）', r'(?i)preprint'), ('draft（不分大小写）', r'(?i)draft'),
          ('to be settled', r'to be settled'), ('\\TODO 使用', r'\\TODO\{'),
          ('\\TODO 宏定义', r'\\newcommand\{\\TODO'), ('\\Pv 使用', r'\\Pv\{'),
          ('\\today', r'\\today'), ('XXX/TBD/FIXME', r'(?i)\b(XXX|TBD|FIXME)\b'),
          ('匿名/占位作者', r'(?i)(anonymous|placeholder|TBC)')]
for name, pat in checks:
    n_all, n_body = count(flat, pat), count(body, pat)
    if name.startswith('draft'):
        flag = '逐条见下'
    else:
        flag = 'OK' if n_body == 0 else '须处理'
    p('   %-22s 全文 %2d 处 / 去注释后 %2d 处  %s' % (name, n_all, n_body, flag))
assert count(body, r'(?i)preprint') == 0, '正文仍含 preprint 字样'
assert count(body, r'to be settled') == 0, '作者块仍未定'
assert count(body, r'\\TODO\{|\\Pv\{|\\today') == 0, '仍有草稿宏或 \\today'
# "draft" 允许留下的唯一用法：正文自述本 note 自己的早期版本（撤回账的一部分），
# 逐条列出并核对，不允许再出现在标题/作者块/文件名里。
keep = [ln for ln in body.split('\n') if re.search(r'(?i)\bdraft\b', ln)]
p('   正文 "draft" %d 处，逐条如下（均为对自家早期版本的自述，属撤回账，不删）：' % len(keep))
for ln in keep:
    p('      | %s' % ' '.join(ln.split())[:112])
assert all('earlier draft' in ln.lower() or 'off the draft' in ln.lower()
           for ln in keep), '出现了未被授权的 draft 用法'

# ---------------------------------------------------------------- [P4] 摘要单段
p('')
p('== [P4] JAIGP 表单字段对账 ==')
ab = re.search(r'\\begin\{abstract\}(.*?)\\end\{abstract\}', flat, re.S).group(1)
paras = [x for x in re.split(r'\n\s*\n', ab.strip()) if x.strip()]
plain = re.sub(r'\$[^$]*\$', ' ', ab)          # 先摘掉行内数学式
words_prose = len(re.sub(r'\\[a-zA-Z]+|[{}&\\]', ' ', plain).split())
words = len(re.sub(r'\\[a-zA-Z]+|[{}$&\\]', ' ', ab).split())
title = re.search(r'\\title\{(.*?)\}', flat, re.S).group(1)
title_one = ' '.join(re.sub(r'\\\\|\n|~', ' ', title).split())
assert 'preprint' not in title_one.lower(), '标题含 preprint'
kw = re.search(r'\\begin\{IEEEkeywords\}(.*?)\\end\{IEEEkeywords\}', flat, re.S).group(1)
p('   Title（单行，可直接粘表单）：%s' % title_one)
p('   Abstract：%d 段（表单要求单段）⇒ %s；散文 %d 词（不计行内数学式；含式记 %d 词）'
  % (len(paras), '合' if len(paras) == 1 else '不合', words_prose, words))
p('   Keywords：%s' % ' '.join(kw.split()))
byline = re.search(r'\\author\{(.*?)\\thanks', flat, re.S).group(1)
byline_one = ' '.join(re.sub(r'\\\\|~|%', ' ', byline).split())
p('   Byline：%s' % byline_one)
assert len(paras) == 1, '摘要不是单段'
assert 120 <= words_prose <= 260, '摘要散文 %d 词不在期刊常规区间' % words_prose


def prose(text):
    return len(re.sub(r'\\[a-zA-Z]+|[{}&\\]', ' ', re.sub(r'\$[^$]*\$', ' ', text)).split())


# 同一口径量一次改动前的摘要（从 git 取，不硬编码）
old_words = None
try:
    g = subprocess.run(['git', '-C', os.path.dirname(os.path.dirname(PRE)),
                        'show', 'HEAD:out/rc/preprint/preprint.tex'],
                       stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    if g.returncode == 0:
        old_abs = re.search(r'\\begin\{abstract\}(.*?)\\end\{abstract\}',
                            g.stdout.decode('utf-8', 'replace'), re.S)
        if old_abs:
            old_words = prose(old_abs.group(1))
except OSError:
    pass
p('   摘要散文词数：改前 %s → 改后 %d（同一口径：摘掉行内数学式后计）'
  % (str(old_words) if old_words else 'git 无旧版', words_prose))
# 表单可直接粘贴的纯文本（摘要去掉 \emph 包装与 --- 破折号）
abs_txt = re.sub(r'\\emph\{([^{}]*)\}', r'\1', ab)
abs_txt = re.sub(r'\\,', ' ', abs_txt)
abs_txt = ' '.join(abs_txt.replace('---', '—').replace('--', '–').split())
io.open(os.path.join(HERE, 'jaigp_title.txt'), 'w', encoding='utf-8',
        newline='\r\n').write(title_one + '\n')
io.open(os.path.join(HERE, 'jaigp_abstract.txt'), 'w', encoding='utf-8',
        newline='\r\n').write(abs_txt + '\n')
p('   已写出可直接粘贴的纯文本：p0/jaigp_title.txt、p0/jaigp_abstract.txt'
  '（摘要纯文本 %d 词）' % len(abs_txt.split()))

# ---------------------------------------------------------------- [P5] 打包
p('')
p('== [P5] 上传包（JAIGP：PDF ≤20 MB，LaTeX 源 ≤30 MB） ==')
pdf = os.path.join(PRE, 'preprint.pdf')
zip_path = os.path.join(PRE, 'jaigp_submission.zip')
if os.path.exists(zip_path):
    os.remove(zip_path)
with zipfile.ZipFile(zip_path, 'w', zipfile.ZIP_DEFLATED) as z:
    z.write(os.path.join(PRE, 'fig_wall.pdf'), 'fig_wall.pdf')
    # 单一入口：包内只有一个 .tex，且它不 \input 任何东西，转换管线不会选错主文件
    z.write(flat_tex, 'preprint.tex')
with zipfile.ZipFile(zip_path) as z:
    names = z.namelist()
    bad = z.testzip()
mb = lambda b: '%.2f MB' % (b / 1048576.0)
p('   PDF  %s = %s（上限 20 MB）' % ('preprint.pdf', mb(os.path.getsize(pdf))))
p('   源包 %s = %s（上限 30 MB），条目 %s，完整性 testzip %s'
  % ('jaigp_submission.zip', mb(os.path.getsize(zip_path)), names, '通过' if bad is None else '失败 ' + str(bad)))
assert bad is None and os.path.getsize(pdf) < 20 * 1048576 and os.path.getsize(zip_path) < 30 * 1048576
p('')
p('   [P6] 本清单只审 artefact，不产新科学结论；页数/字节/计数均可由 '
  'preprint/jaigp_flat.log 复核。')
p('')
p('退出码 0 = 扁平化成功且两版页数一致；引用 0 缺失；摘要单段且在词数区间内；'
  '正文无 preprint/to be settled/\\TODO/\\Pv/\\today（"draft" 仅允许出现在逐条列出的'
  '自家早期版本自述句里）；PDF 与源包均在表单限额内。')
print('OK')
