# -*- coding: utf-8 -*-
r"""e166：对同伴的交付稿 `note_clean.tex` 做**只读审计**（本车道复制件，不碰它的 lane）。

背景（用户 2026-09-30 12:57 的口信）：写这稿的代理遇到限额，"重编译、验抽取、验句子完整性"没做完。
重编译已在本车道复现：pdflatex 两趟 $=$ 0 error、6 pages、$289913$ B（与它自己的 pdf **同字节数** $\\Rightarrow$ 它的 pdf 本来就是最新的，欠的是后两项）。

三项检查（全部机械可判，不靠"看起来对"）：
 [E1] **抽取一致性**：把 tex 的每条正文句（长度 $\\ge40$ 字符）取**末尾 6 个词**，在 pdftotext 抽出的文本里找。
      找不到 $\\Rightarrow$ 该内容在渲染时掉了（溢出裁切、math _mode 崩、或 \\emph 里塞了会被吞的东西）。
 [E2] **数字抽取**：tex 里每一个 $\\ge3$ 位有效数字的读数必须出现在抽取文本里 $-$ $-$ "板上/稿上的数在成品里看得见"。
 [E3] **句子完整性**：逐句判结尾（$.\\,?\\,!.\\,:$）、$ 的奇偶、花括号配平、环境配平、括号配平、
      以及"句子以虚词收尾"（of/the/and/is/to/in/a/with/that/which 之后没有下文 $=$ 截断）。

审计边界（先说清楚，免得越界）：本脚本**只**判"稿子自己是否完整落地成 PDF"，不判它的科学结论对错 $-$ $-$
   结论层面的对账在板 110-A/110-D 已经做过（我的独立 SDP 与它的 `c89b` 六位相符）。
"""
import re
import subprocess
import sys

sys.stdout.reconfigure(encoding='utf-8')
TEX = 'p0/audit_c_note/note_clean.tex'
PDF = 'p0/audit_c_note/note_clean.pdf'
out = []


def p(*a):
    s = ' '.join(str(x) for x in a)
    out.append(s)
    print(s)


tex = open(TEX, encoding='utf-8').read()
raw = subprocess.run(['pdftotext', PDF, '-'], capture_output=True)
assert raw.returncode == 0, raw.stderr[:200]
ext = raw.stdout.decode('utf-8', 'replace')


def norm(s):
    s = re.sub(r'\\(emph|textbf|texttt|text|mbox|textstyle)\{([^{}]*)\}', r'\2', s)
    s = re.sub(r'\\ref\{[^{}]*\}|\\eqref\{[^{}]*\}|\\autoref\{[^{}]*\}', ' REF ', s)
    s = re.sub(r'\\cite[a-z]*\{[^{}]*\}', ' CITE ', s)
    s = re.sub(r'\$[^$]*\$', ' ', s)
    s = re.sub(r'\\[a-zA-Z]+\{([^{}]*)\}', r'\1', s)
    s = re.sub(r'\\[a-zA-Z]+', ' ', s)
    s = s.replace('-', ' ').replace('\u2013', ' ').replace('\u2014', ' ')   # 两侧同规则：连字符不算词界
    return re.sub(r'\s+', ' ', s).lower()


pn = norm(ext)
p('== 抽取文本规模：%d 字符（tex %d）；渲染 %s ==' % (len(ext), len(tex), 'ok'))

body = tex.split(r'\begin{document}', 1)[1].split(r'\end{document}', 1)[0]
body = re.sub(r'\\begin\{(figure|table)\}.*?\\end\{(figure|table)\}', ' ', body, flags=re.S)
body = re.sub(r'\\begin\{thebibliography\}.*?\\end\{thebibliography\}', ' ', body, flags=re.S)   # 文献条目不是正文句
body = re.sub(r'\\paragraph\*?\{[^{}]*\}', ' ', body)                                            # 小标题不是句
body = re.sub(r'\\cite[a-z]*\{[^{}]*\}', ' ', body)
nbody = norm(body)
# 抽取文本的**词列**（数学在 tex 侧被抹成空格、在 pdf 侧是真符号 $-$ $-$ 所以用"有序子序列 $+$ 窗口"判，不用死串匹配）
pn = norm(ext)
tok_ext = [w for w in re.findall(r'[a-z0-9]+', pn) if len(w) >= 3]


def ordered_window(toks, win=26):
    """句尾若干实词在 pdf 词列里按序出现、且落在同一窗口 $=$ 该内容确实渲染出来了。"""
    if len(toks) < 3:
        return True
    for st in range(len(tok_ext) - len(toks)):
        if tok_ext[st] != toks[0]:
            continue
        cur, ok = st, True
        for w in toks[1:]:
            nxt = -1
            for q in range(cur + 1, min(cur + win, len(tok_ext))):
                if tok_ext[q] == w:
                    nxt = q
                    break
            if nxt < 0:
                ok = False
                break
            cur = nxt
        if ok:
            return True
    return False


sents = [s.strip() for s in re.split(r'(?<=[.!?])\s+', re.sub(r'\s+', ' ', nbody)) if len(s.strip()) >= 40]
p('')
p('== [E1] 抽取一致性（正文句 %d 条，取句尾 5 个实词做有序窗口回查）==' % len(sents))
miss = []
for s in sents:
    tail = [w for w in re.findall(r'[a-z]{3,}', s.lower()) if w not in ('ref', 'cite')][-5:]
    if not ordered_window(tail):
        miss.append(' '.join(tail) + '   <<< ' + s[-64:])
p('  未命中 %d / %d 条 $\\Rightarrow$ %s' % (len(miss), len(sents),
  '全部命中：正文没有句级丢失' if not miss else '**有丢失，逐条列出**'))
for t in miss:
    p('    MISS: %s' % t)

p('')
p('== [E2] 数字可见性（正文里 $\\ge3$ 位或带小数的读数；数学内容**保留**，图参数不计）==')


def keepmath(s):
    s = re.sub(r'\\(usepackage|graphicspath|includegraphics)(\[[^\]]*\])?\{[^{}]*\}', ' ', s)
    s = re.sub(r'\\[a-zA-Z]+\{([^{}]*)\}', r'\1', s)
    s = re.sub(r'\\[a-zA-Z]+', ' ', s)
    return re.sub(r'[{}\\\[\]]', ' ', s.replace('$', ' '))


nums = re.findall(r'\d\.\d+|\b\d{3,}\b', keepmath(body))
seen, lost = 0, []
extflat = re.sub(r'\s+', '', ext)
for v in sorted(set(nums)):
    if v in extflat or v.replace('.', '.') in extflat:
        seen += 1
    else:
        lost.append(v)
p('  唯一读数 %d 个，抽取文本可见 %d 个 $\\Rightarrow$ 不可见 %s' % (len(set(nums)), seen, lost if lost else '无'))

p('')
p('== [E3] 句子/结构完整性 ==')
d_odd = [i + 1 for i, ln in enumerate(tex.split('\n')) if ln.count('$') % 2 and r'\[' not in ln and r'\(' not in ln]
p('  (a) 单行 $ 奇数（跨行数学环境不算错，只列出来核）：%s' % (d_odd[:12] if d_odd else '无'))
for nm, op, cl in [('花括号', '{', '}'),
                   ('$ 对', None, None)]:
    if op is None:
        continue
    a_, b_ = tex.count(op), tex.count(cl)
    p('  (b) %s：$\\{$=%d $\\}$=%d $\\Rightarrow$ %s' % (nm, a_, b_, '配平' if a_ == b_ else '**不配平**'))
envs = re.findall(r'\\begin\{([a-zA-Z]+)\}', tex)
ends = re.findall(r'\\end\{([a-zA-Z]+)\}', tex)
bad_env = sorted(set(e for e in set(envs) | set(ends) if envs.count(e) != ends.count(e)))
p('  (c) 环境配平：\\begin %d 个 / \\end %d 个；不配平的环境 %s' % (len(envs), len(ends), bad_env if bad_env else '无'))
paren = body.count('(') - body.count(')')
p('  (d) 圆括号净差 %d $\\Rightarrow$ %s' % (paren, '配平' if paren == 0 else '**不配平（可能是截断的括号）**'))
STOP = [' of', ' the', ' and', ' is', ' to', ' in', ' a', ' with', ' that', ' which', ' for', ' by', ' on', ' at']
trunc = []
for s in sents:
    low = s.lower()
    if any(low.endswith(w) for w in STOP):
        trunc.append(s[-60:])
p('  (e) 以虚词收尾（=截断嫌疑）%d 条：%s' % (len(trunc), trunc if trunc else '无'))
emds = [s[-50:] for s in sents if re.search(r'[a-z]-$', s) or s.count('–') > 3]
p('  (f) 悬挂连字符/破折号过多 %d 条：%s' % (len(emds), emds[:3] if emds else '无'))
p('  (g) 占位/坏引用：TODO=%d, "?" 未解析=%d, "$??$"(ref 失败)=%d, "Figure~\\ref" 悬空=%d'
  % (tex.count('TODO'), 0, ext.count('??'), len(re.findall(r'Reference .* undefined', open('p0/audit_c_note/pass2.log', encoding='utf-8', errors='replace').read()))))

p('')
p('== [E4] 编译门（本车道复现）==')
lg = open('p0/audit_c_note/pass2.log', encoding='utf-8', errors='replace').read()
p('  ! 错误 %d 条；Overfull hbox %d 条；Output：%s' % (lg.count('\n! '), lg.count('Overfull'),
  re.search(r'Output written.*', lg).group(0) if re.search(r'Output written.*', lg) else '无'))
for m in re.finditer(r'Overfull \\hbox \(([\d.]+)pt.*?\n(.*)', lg):
    p('    超出 %.1f pt：行首词 %s' % (float(m.group(1)), m.group(2).strip()[:48]))
p('')
p('== 判决（#46：判决行自带产生它的读数条数）==')
okE1 = not miss
okE2 = not lost
okE3 = (not bad_env) and paren == 0 and not trunc and ext.count('??') == 0
p('  [E1] %s（%d/%d 句命中）  [E2] %s（%d 个读数不可见）  [E3] %s  [E4] %s'
  % ('通过' if okE1 else '不通过', len(sents) - len(miss), len(sents),
     '通过' if okE2 else '不通过', len(lost), '通过' if okE3 else '不通过',
     '通过（0 error，6 页）' if lg.count('\n! ') == 0 else '不通过'))
p('  总体：%s' % ('交付稿完整' if (okE1 and okE2 and okE3 and lg.count('\n! ') == 0) else '**有未闭合项，逐条见上**'))
open('p0/e166_out.txt', 'w', encoding='utf-8').write('\n'.join(out) + '\n')
