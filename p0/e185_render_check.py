# Rendered-text integrity check for the two files this lane edited.
#
# Context: the manuscript is IEEEtran [journal] = TWO-COLUMN. `pdftotext -layout`
# reconstructs by absolute position, so it merges the left- and right-column
# fragments that share a y-coordinate into one output line (observed: "withdraw
# it -- together with              tificate only where explicitly said"). That
# artifact -- not a truncated sentence -- is what made six phrases look missing in
# the earlier spot check. So this script audits both extraction modes and stops
# spot-checking: it verifies EVERY prose sentence of the two edited files against
# the rendered document.
#
# Why the gate is letter-level, not token-level ([W1] is kept as a diagnostic only):
# the extractor's spacing is not reproducible -- "active-mode" becomes "active- mode"
# while "phase-transition" becomes "phasetransition", `Remark` splits as "re mark" at
# two of its fourteen occurrences -- and the render INSERTS text the stripped source
# does not have (a `Section~\ref` prints the roman numeral "iv", which lowercases into
# prose). Letters are insensitive to all of that. So a sentence passes when each part
# of it is covered by a run of >= 6 letters that occurs in the render, in order, with
# block jumps allowed for the float/page-break reflow that puts one sentence of a
# paragraph in another block. Insertions cost nothing; dropped letters cannot hide --
# and [W1e] proves that by deleting 12 letters from the render and requiring a flag.
#
# Known extraction limitation, recorded rather than hidden: pdftotext does not map
# newtxMathItalic glyphs, so `$C(V)$` reaches the text layer as " ()" ("the
# constant  () was left unspecified"). Math is invisible to this check by
# construction; it verifies PROSE survival, and takes equation NUMBERING from the
# .aux, which is what LaTeX actually shipped to the PDF.
#
# Compile evidence is read from the build, not re-run here (the build is shared with
# the other lane): preprint/preprint.log says "Output written on preprint.pdf
# (12 pages, 353945 bytes)" with zero "^!" lines, and no .tex file is newer than the
# .aux, so preprint/preprint.pdf is the current sources' output.
#
# Verdicts:
#  [W0] extraction-mode audit: the same phrases against both modes.
#  [W0b] those phrases must be verbatim source text (self-check of the checker --
#        a phrase list typed from memory proves nothing).
#  [W1] diagnostic only: token-order alignment per sentence/paragraph.
#  [W1c] THE GATE: every prose sentence's letters present, in order, in the render.
#  [W1e] negative + positive control on that gate, so a pass means something.
#  [W1b] float captions under the same letter rule.
#  [W2] equation numbering from preprint.aux: eq:wf -> (8), eq:cond -> (9), and no
#        label referenced by the text is missing from the aux.
#  [W3] citation graph closed: every \cite key resolves, no orphan bibitem.
#  [W3b] every bibliography entry's letters present in the render.
import os
import re
import subprocess
import sys
import difflib
import unicodedata

os.chdir(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
out = []


def p(*a):
    s = ' '.join(str(x) for x in a)
    out.append(s)
    print(s)
    sys.stdout.flush()


def norm(t):
    for a, b in [('\ufb00', 'ff'), ('\ufb01', 'fi'), ('\ufb02', 'fl'),
                 ('\ufb03', 'ffi'), ('\ufb04', 'ffl')]:
        t = t.replace(a, b)
    for a in ['\u2013', '\u2014', '\u2212', '\u2010', '\u2011']:
        t = t.replace(a, '-')
    for a in ['\u2019', '\u2018']:
        t = t.replace(a, '')          # U+2019 is an apostrophe, not a dash:
    for a in ['\u201c', '\u201d']:    # "project's" must stay one token
        t = t.replace(a, '"')
    return t.lower().replace('\u00a0', ' ')


def toks(s):
    """Pure-letter tokens, one recipe for both sides.

    Hyphens are dropped on both sides because a line break inside
    `rank-deficient` reaches the text layer as "rank-" + newline + "deficient",
    which no dash rule can tell apart from a real dash; quotes because the source
    has ``'' and the rendering has U+201D; accents because the source has L\\"owner
    and the rendering has an o-umlaut; digits because math-adjacent numerals do
    not extract. What survives is prose, and prose is what must be present.
    """
    s = unicodedata.normalize('NFD', norm(s))
    s = ''.join(c for c in s if not unicodedata.combining(c))
    s = s.replace('-', '')          # pdftotext is inconsistent at line breaks:
    s = s.replace("'", '').replace('"', '')   # "project's"->projects on both sides
    return re.sub(r'[^a-z]+', ' ', s).split()   # but "active-mode"->"active- mode"


def chars(s):
    """The letter string of a run of text: the spacing-proof form of it.

    This is what makes the sentence check sound. pdftotext renders the same
    source word as "active- mode", "phasetransition" or "phase transition"
    depending only on where the line break falls, and it splits `Remark` into
    "re mark" at two of its fourteen occurrences; every one of those has the
    same letter string as the source. Digits and math are excluded on both
    sides (newtxMathItalic glyphs do not map, and math-adjacent numerals do not
    extract), so the surviving letters are pure prose and can be compared
    verbatim instead of approximately.
    """
    s = unicodedata.normalize('NFD', norm(s))
    s = ''.join(c for c in s if not unicodedata.combining(c))
    return ''.join(re.findall(r'[a-z]+', s.replace('-', '').replace("'", '')))


def strip_tex(s):
    s = s.replace('\\%', '@@PCT@@')
    s = re.sub(r'\\{1,2}\[[^\]]*\]', ' ', s)              # \\[1pt] is a length, not prose
    s = re.sub(r'%[^\n]*', '', s)
    s = re.sub(r'\\(["^~=])', '', s)                     # L\"owner -> Lowner
    s = re.sub(r'\$[^$]*\$', ' ', s)                     # inline math
    s = re.sub(r'(?s)\\begin\{(?:equation|align\*?|gather\*?|eqnarray|array|'
               r'tabular\w*|thebibliography)\}.*?\\end\{(?:equation|align\*?|'
               r'gather\*?|eqnarray|array|tabular\w*|thebibliography)\}', ' ', s)
    s = re.sub(r'\\(?:begin|end)\{[^{}]*\}', ' ', s)     # env tags are not prose;
    # \begin{proof}[Proof sketch] keeps its bracket argument, which the render prints
    s = re.sub(r'\\cite\w*\{[^{}]*\}', ' ', s)          # keys are not prose
    s = re.sub(r'\\(?:ref|eqref|label|pageref|index)\{[^{}]*\}', ' ', s)
    s = re.sub(r'\\(?:emph|texttt|textmd|textsc|textbf|textit|mbox|textrm|'
               r'caption|item|prop|rm|bf)\b\*?\{?', ' ', s)
    s = re.sub(r'\\[a-zA-Z@]+\*?', ' ', s)              # remaining commands
    s = re.sub(r'[{}&|~]', ' ', s)
    return s.replace('@@PCT@@', '%')


def inslice(seq, sub):
    m = len(sub)
    return m == 0 or any(seq[k:k + m] == sub for k in range(len(seq) - m + 1))


def align_from(S, R, start, slack=300):
    """Count the source tokens absent from the rendered region of one paragraph.

    Two hand-rolled greedy aligners gave false failures here: first-occurrence
    scanning lets a stale common word (`section` recurs every few lines) leapfrog
    the pointer past the paragraph, and an expectation-anchored variant is worse.
    difflib's matching blocks take insertions (cross-reference numerals, which the
    PDF prints and the stripper removes from the source) and deletions without
    heuristics, so a paragraph printed verbatim aligns completely.
    """
    n = len(S)
    stop = min(len(R), start + n + slack)
    region = R[max(0, start - 6):stop]
    sm = difflib.SequenceMatcher(None, S, region, autojunk=False)
    matched = set()
    for a, b, size in sm.get_matching_blocks():
        matched.update(range(a, a + size))
    unm = [w for t, w in enumerate(S) if t not in matched]
    return unm, [size for _, _, size in sm.get_matching_blocks()]


def anchors(sub, R, n=8):
    """Positions where the first n source tokens sit contiguously.

    Anchoring first is what makes the forward window meaningful: this file's prose
    sits 2/3 of the way through the document, so a window measured from position 0
    can never find it (that was the 315-token failure of the previous draft).
    """
    return [k for k in range(len(R) - n) if R[k:k + n] == sub[:n]]


FILES = ['preprint/frag_rc_waterfill.tex', 'preprint/sec_related.tex']
SRC = {f: open(f, encoding='utf-8').read() for f in FILES}
BIBTXT = open('preprint/sec_bib.tex', encoding='utf-8').read()
PRE = open('preprint/preprint.tex', encoding='utf-8').read()
AUX = open('preprint/preprint.aux', encoding='utf-8', errors='replace').read()

for args, dst in [(['-enc', 'UTF-8'], 'p0/e185_pdf_default.txt'),
                  (['-enc', 'UTF-8', '-layout'], 'p0/e185_pdf_layout.txt')]:
    r = subprocess.run(['pdftotext'] + args + ['preprint/preprint.pdf', dst],
                       capture_output=True, text=True)
    assert r.returncode == 0, args
DEF = open('p0/e185_pdf_default.txt', encoding='utf-8', errors='replace').read()
LAY = open('p0/e185_pdf_layout.txt', encoding='utf-8', errors='replace').read()
RT, LT = toks(DEF), toks(LAY)
CD, CL = chars(DEF), chars(LAY)
p('[W0] 渲染 prose token 数  default(栏序)=%d  layout(按行)=%d' % (len(RT), len(LT)))

PROSE = [
    'their theorem does not price the object of this section',
    'needs no commutativity at all',
    'no determinant appears',
    'a different statement from the design-dependent monotonicity',
    'satisfies none of the three',
    'neither corroborates nor contradicts it',
    'closing a question left open by',
    'which holds in the three structural cases',
    "what makes hadamard's inequality tight",
    'move with the design',
    'is entirely in the requirement that',
    'the sensor map enters one level up through',
]
BIPH = [
    'asymptotic reverse waterfilling algorithm of',
    'optimal sensor selection of linear',
    'sensor data scheduling for linear',
]
for q in PROSE + BIPH:
    qt = toks(q)
    p('[W0] %-56s default=%-4s layout=%-4s'
      % (q[:56], 'OK' if inslice(RT, qt) else 'MISS', 'OK' if inslice(LT, qt) else 'MISS'))

STOK = {f: toks(strip_tex(SRC[f])) for f in FILES}
BTOK = toks(BIBTXT)
badph = [q for q in PROSE if not any(inslice(s, toks(q)) for s in STOK.values())]
badph += [q for q in BIPH if not inslice(BTOK, toks(q))]
p('[W0b] 抽查短语逐字可溯  prose=%d  bib=%d  不合格=%s'
  % (len(PROSE), len(BIPH), badph or 'none'))
assert not badph, 'spot-check phrase not verbatim in source: %s' % badph

# A whole-section monotone alignment is the wrong model for an IEEE two-column
# paper: floats (tables, figures) are pulled to the top of a page, so their
# captions reach the text layer out of source order, and one early divergence
# cascades (that is how the previous draft reported 320 "missing" tokens while
# the sentences were plainly intact). Paragraphs, in contrast, do not move. So
# each source paragraph is anchored independently and aligned forward.
def locate(S):
    """Best (start, unmatched-tokens) for one token run, over the whole text layer.

    Prefix probes come first because a heading run-in or a marginal note can sit
    in front of the paragraph; the suffix probes matter when a cross-reference
    numeral is inserted near the start and breaks every prefix.
    """
    best = None
    probes = []
    for off in (0, 1, 2, 3, 6, 10, 16):          # debris can head a paragraph
        for n in (8, 6, 5, 4):
            if off + n <= len(S):
                probes.append((off, S[off:off + n]))
    for n in (8, 6, 5):                          # a cross-ref inserted near the
        if len(S) >= n:                          # start can defeat every prefix
            probes.append((len(S) - n, S[-n:]))
    for off, sub in probes:
        for a in anchors(sub, RT, n=len(sub))[:6]:
            start = a - off
            if start < 0:
                continue
            unm, jumps = align_from(S, RT, start)
            score = (len(unm), -sum(jumps))
            if best is None or score < best[0]:
                best = (score, start, unm)
        if best and best[0][0] <= 0.03 * len(S):
            break
    return None if best is None else (best[1], best[2])


n_bad = 0
GATE = 0.95
for f in FILES:
    body = re.sub(r'(?s)\\begin\{(table|tabular\w*|figure)\}.*?\\end\{(table|tabular\w*|figure)\}',
                  ' ', SRC[f])
    paras = re.split(r'\n\s*\n', body)
    para_rows = []
    sent_rows = []
    kept = 0
    for pi, para in enumerate(paras):
        plain = strip_tex(para)
        S = toks(plain)
        if len(S) < 10:
            continue
        kept += 1
        pb = locate(S)
        para_rows.append((1.0 - (len(pb[1]) / len(S) if pb else 1.0), pi, len(S),
                          ' '.join(S[:8])))
        # The sentence is the gated unit, not the paragraph. A paragraph can be
        # split in the text layer without a word being lost: a page break inside a
        # float-heavy section sends one sentence to a different block (Remark 5
        # does exactly that -- its last sentence is extracted 113 lines before its
        # first), so no single window covers a whole paragraph. Sentences do not
        # straddle blocks, so each is anchored independently and must align alone.
        for sn in re.split(r'(?<=[.;:!?])\s+', plain):
            T = toks(sn)
            if len(T) < 6:
                continue
            hb = locate(T)
            fr = 1.0 - (len(hb[1]) / len(T) if hb else 1.0)
            sent_rows.append((fr, pi, len(T), ' '.join(T[:7]),
                              'NO-ANCHOR' if hb is None else ' '.join(hb[1][:6])))
    bad = [r for r in sent_rows if r[0] < GATE]
    p('[W1] %-26s 段=%d 句=%d  句级最低=%.3f  低于门限(%.2f)的句=%d  段级最低=%.3f'
      % (os.path.basename(f), kept, len(sent_rows),
         min([r[0] for r in sent_rows] or [1.0]), GATE, len(bad),
         min([r[0] for r in para_rows] or [1.0])))
    for fr, pi, ln, head, unm in sorted(sent_rows)[:6]:
        p('[W1]   句@段%-3d 对齐=%.3f (%2d tok) 未对齐=%-38s | 句首=%s'
          % (pi, fr, ln, unm[:38], head[:44]))
    for fr, pi, ln, head in sorted(para_rows)[:4]:
        p('[W1]   段%-3d 对齐=%.3f (%d token) | 段首=%s' % (pi, fr, ln, head[:52]))
    for fr, pi, ln, head, unm in bad:
        p('[W1]   >>> 段%d 句子未对齐(%d tok, %.3f): %s || 句首=%s'
          % (pi, ln, fr, unm, head))
    n_bad += sum(1 for r in bad if r[0] == 0.0) * 1000 + len(bad)

# [W1c] THE GATE: sentence integrity, letter by letter.
# Why not token alignment ([W1], kept above only as a diagnostic) and not plain
# substring search either:
#   * letters are insensitive to every spacing/hyphen artifact of the extraction,
#     and to `Remark` splitting as "re mark";
#   * substring search is too strict, because the render INSERTS text the stripped
#     source does not have -- `Section~\ref{...}` prints the roman numeral "iv"
#     (lowercase, so it looks like prose) and math-adjacent numerals -- while it
#     never deletes a source letter;
#   * a global window is too narrow, because floats and column overflow send one
#     sentence of a paragraph to another block (Remark 5's last sentence extracts
#     113 lines before its first), so no contiguous window covers it.
# The rule used here: walk the sentence's letter string and demand that every part
# of it is explained by a contiguous run of >= MINRUN letters that really occurs in
# the render. Insertions cost nothing; a deletion leaves a run too short to pass.
MINRUN = 6


def gaps_of(cs, hay, minrun=MINRUN):
    """Letters of `cs` that no run of >= minrun letters of `hay` explains, in order.

    Returns (gap runs, block jumps). The walk consumes `cs` left to right and always
    looks for the longest run starting at the current letter at or after the position
    reached so far, so order is enforced; when a run exists only BEHIND that position
    the walk retries globally and counts a block jump, which is what a float or a
    column overflow does to a paragraph in the text layer (Remark 5's last sentence
    is extracted 113 lines before its first). A run shorter than minrun that occurs
    nowhere is a real letter loss and is reported: insertions in the render cost
    nothing, deletions from it cannot be hidden.
    """
    miss, jumps, i, k, pos = [], 0, 0, len(cs), 0
    while i < k:
        lo, hi, best = minrun, k - i, 0
        while lo <= hi:
            mid = (lo + hi) // 2
            if hay.find(cs[i:i + mid], pos) >= 0:
                best, lo = mid, mid + 1
            else:
                hi = mid - 1
        if best < minrun:
            lo, hi, best2 = minrun, k - i, 0
            while lo <= hi:
                mid = (lo + hi) // 2
                if hay.find(cs[i:i + mid]) >= 0:
                    best2, lo = mid, mid + 1
                else:
                    hi = mid - 1
            if best2 >= minrun:
                jumps += 1
                pos = hay.find(cs[i:i + best2]) + best2
                i += best2
                continue
            miss.append(i)
            i += 1
        else:
            pos = hay.find(cs[i:i + best], pos) + best
            i += best
    runs = []
    for x in miss:
        if runs and x == runs[-1][1]:
            runs[-1][1] = x + 1
        else:
            runs.append([x, x + 1])
    return [(cs[a:b], cs[max(0, a - 20):a], cs[b:b + 20]) for a, b in runs], jumps


n_sent = n_sent_bad = n_jump = n_strict = n_strict_lay = 0
SENT = []
for f in FILES:
    body = re.sub(r'(?s)\\begin\{(table|tabular\w*|figure)\}.*?\\end\{(table|tabular\w*|figure)\}',
                  ' ', SRC[f])
    for para in re.split(r'\n\s*\n', body):
        for sn in re.split(r'(?<=[.;:!?])\s+', strip_tex(para)):
            cs = chars(sn)
            if len(cs) < 14:
                continue
            n_sent += 1
            n_strict += cs in CD
            n_strict_lay += cs in CL
            g, jumps = gaps_of(cs, CD)
            n_jump += jumps
            if g:
                n_sent_bad += 1
                SENT.append((cs, g))
p('[W1c] 句子=%d  字母级缺口句=%d  跨块跳转=%d  (字母串原样连续: 栏序=%d / 按行=%d, '
  '差值即双栏同行拼接造成的抽取假象)' % (n_sent, n_sent_bad, n_jump, n_strict, n_strict_lay))
for cs, g in SENT[:10]:
    p('[W1c]   >>> 句=%s' % cs[:60])
    for miss, pre, post in g[:3]:
        p('[W1c]       缺口="%s"  前=%s  后=%s' % (miss, pre, post))

# [W1e] negative control: the gate must be able to fail, in the right direction. The
# risk being tested for is a sentence the printer dropped part of, so the control
# deletes letters from the RENDER (a page-budget cut is exactly that) and requires
# gaps_of to flag them. Deleting from the source side would prove nothing.
probe_sent = next((c for c in (chars(sn) for f in FILES
                               for para in re.split(r'\n\s*\n', SRC[f])
                               for sn in re.split(r'(?<=[.;:!?])\s+', strip_tex(para)))
                   if 60 <= len(c) <= 200 and c in CD), None)
assert probe_sent, 'no clean sentence found for the negative control'
q = CD.find(probe_sent)
w = q + len(probe_sent) // 2
hay2 = CD[:w] + CD[w + 12:]                      # 12 printed letters go missing
g2, _ = gaps_of(probe_sent, hay2)
p('[W1e] 阴性对照: 从渲染文本中删去该句中段 12 个字母 -> gaps_of 报缺失=%s '
  '(若为 False 则本门限形同虚设); 报出的缺口=%s'
  % (bool(g2), [m for m, _a, _b in g2][:3]))
assert g2, 'the sentence gate cannot detect letters dropped from the render -- it is vacuous'
g3, _ = gaps_of(probe_sent, CD)
assert not g3
p('[W1e] 阳性对照: 同一句对未改动的渲染 -> 无缺口 %s ; 对照有效' % (not g3,))

# [W1b] float captions: excluded from [W1] because floats move, checked here with
# the same letter rule because their text must still survive verbatim.
n_cap = n_cap_bad = 0
for f in FILES:
    for m in re.finditer(r'(?s)\\caption\{(.*?)\}\s*\n?\s*\\label', SRC[f]):
        cs = chars(strip_tex(m.group(1)))
        if len(cs) < 14:
            continue
        n_cap += 1
        g, j = gaps_of(cs, CD)
        n_cap_bad += len(g)
        p('[W1b] caption %2d (%3d 字母) 字母缺口=%d 跨块跳转=%d%s'
          % (n_cap, len(cs), len(g), j,
             '' if not g else '  >> ' + ' '.join(x[0] for x in g[:4])))

# --- [W2] numbering straight from the .aux --------------------------------
num = {}
for m in re.finditer(r'\\newlabel\{([^}]+)\}\{\{([^}]*?)\}', AUX):
    num.setdefault(m.group(1), m.group(2))
ALL = PRE + ''.join(SRC.values())
for base in sorted(set(re.findall(r'\\input\{([^}]+)\}', PRE))):
    c = os.path.join('preprint', base + '.tex')
    if os.path.exists(c):
        ALL += open(c, encoding='utf-8').read()
refs = set(re.findall(r'\\(?:eq)?ref\{([^}]+)\}', ALL))
missing = sorted(refs - set(num))
p('[W2] aux 里有编号的 label=%d  正文引用 label=%d  引用但 aux 无编号=%s'
  % (len(num), len(refs), missing or 'none'))
for lb in ['eq:cost', 'eq:decay', 'eq:TG', 'eq:red', 'eq:wf', 'eq:cond']:
    p('[W2] %-9s -> aux 记录 (%s)  被引用 %2d 次'
      % (lb, num.get(lb, '?'), len(re.findall(r'\\(?:eq)?ref\{%s\}' % re.escape(lb), ALL))))
assert num.get('eq:wf') == '8', num.get('eq:wf')
assert num.get('eq:cond') == '9', num.get('eq:cond')
assert not missing, missing

# --- [W3] citation graph --------------------------------------------------
keys = set(re.findall(r'\\bibitem(?:\[[^\]]*\])?\{([^}]+)\}', BIBTXT))
cited = set()
for base in set(re.findall(r'\\input\{([^}]+)\}', PRE)):
    c = os.path.join('preprint', base + '.tex')
    if os.path.exists(c):
        for m in re.finditer(r'\\cite\w*\{([^}]*)\}', open(c, encoding='utf-8').read()):
            cited.update(x.strip() for x in m.group(1).split(','))
used_here = set()
for f in FILES:
    for m in re.finditer(r'\\cite\w*\{([^}]*)\}', SRC[f]):
        used_here.update(x.strip() for x in m.group(1).split(','))
p('[W3] bibitem=%d  被引用=%d  从未被引用=%s'
  % (len(keys), len(cited & keys), sorted(keys - cited) or 'none'))
p('[W3] 两个编辑文件的 \\cite key 全部可解析: %s  缺失=%s'
  % (used_here <= keys, sorted(used_here - keys)))
assert used_here <= keys
assert not (keys - cited)

# [W3b] the bibliography is prose too: every entry's letters must reach the render.
n_t = n_t_bad = 0
for m in re.finditer(r'(?s)\\bibitem(?:\[[^\]]*\])?\{[^}]+\}(.*?)(?=\\bibitem|\\end\{thebibliography\}|\Z)',
                     BIBTXT):
    cs = chars(strip_tex(m.group(1)))
    if len(cs) < 30:
        continue
    n_t += 1
    g, j = gaps_of(cs, CD)
    if g:
        n_t_bad += 1
        p('[W3b]   >>> 条目缺字母: %s || 上下文=%s…%s'
          % (' '.join(x[0] for x in g[:3]), g[0][1][-24:], g[0][2][:24]))
p('[W3b] bibitem 条目=%d  字母级缺口条目=%d' % (n_t, n_t_bad))

p('[W1] 参考: token 级对齐低于门限的句=%d（NO-ANCHOR 按 1000 计）。这些缺口是 Remark/'
  'Proposition 被拆成 "re mark"、连字符断行把 "arithmetic-geometric" 切成两段、'
  '\\ref 数字插入所致，逐句由 [W1c] 的字母级检查解释掉，不是丢句。'
  '   [W1b] caption=%d 字母缺口=%d   [W1c] 句子字母缺口=%d   [W3b] 书目缺口=%d'
  % (n_bad, n_cap, n_cap_bad, n_sent_bad, n_t_bad))
assert n_sent_bad == 0 and n_cap_bad == 0 and n_t_bad == 0, \
    'a sentence is missing from the rendered PDF'
p('[W4] 通过：%d 个 prose 句子的每个字母都在 12 页渲染文本里按序出现（run 级覆盖：允许渲染端'
  '插入 \\ref 编号、允许浮动体把句子拆到别的块，共 %d 次跨块跳转；不允许丢字母）；'
  'caption %d/%d 完整；书目 %d 条完整；eq:wf/eq:cond 编号由 aux 确认为 (8)/(9)；'
  '\\cite key 全部可解析、无孤儿条目；[W1e] 的阴性对照证明该门限不是空转。'
  % (n_sent, n_jump, n_cap, n_cap, n_t))

with open('p0/e185_out.txt', 'w', encoding='utf-8', newline='\r\n') as fh:
    fh.write('\n'.join(out) + '\n')
p('[done] p0/e185_out.txt 行数=%d' % len(out))
