import re
import sys

sys.stdout.reconfigure(encoding='utf-8', errors='replace')
raw = open('preprint/frag_rc_waterfill.tex', encoding='utf-8').read()
body = re.sub(r'(?s)\\begin\{(table|tabular\w*|figure)\}.*?\\end\{(table|tabular\w*|figure)\}',
              ' ', raw)
ps = re.split(r'\n\s*\n', body)
for pi in (16, 13):
    print('===== 源文段 %d（%d 行）=====' % (pi, ps[pi].count(chr(10)) + 1))
    print(ps[pi][:1500])
    print()

T = open('p0/e185_pdf_default.txt', encoding='utf-8', errors='replace').read()
TL = re.sub(r'\s+', ' ', T.replace('\ufb01', 'fi').lower())
for pat in ['boxing eigenvalues', 'what is not claimed', 'order phase',
            'on the sensing side', 'spectral thresholding', 'branch and bound',
            'branch-bound', 'a closed-form rank-one counterexample']:
    m = re.search(re.escape(pat), TL)
    print('%-38s %s' % (pat, 'FOUND' if m else '*** ABSENT ***'))
    if m:
        print('    …%s…' % TL[max(0, m.start() - 130):m.start() + 190])
