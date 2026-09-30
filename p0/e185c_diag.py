import re
import sys
sys.path.insert(0, 'p0')

sys.stdout.reconfigure(encoding='utf-8', errors='replace')
import importlib.util
spec = importlib.util.spec_from_file_location('e185', 'p0/e185_render_check.py')
# do not execute the whole checker; re-implement the two helpers by importing text
src = open('p0/e185_render_check.py', encoding='utf-8').read()
ns = {'__file__': 'p0/e185_render_check.py'}
exec(compile(src.split("FILES = [")[0], 'e185helpers', 'exec'), ns)
toks, strip_tex, anchors = ns['toks'], ns['strip_tex'], ns['anchors']

raw = open('preprint/frag_rc_waterfill.tex', encoding='utf-8').read()
body = re.sub(r'(?s)\\begin\{(table|tabular\w*|figure)\}.*?\\end\{(table|tabular\w*|figure)\}',
              ' ', raw)
ps = re.split(r'\n\s*\n', body)
S = toks(strip_tex(ps[2]))
print('源文段2 tokens(%d):' % len(S))
print(' '.join(S))
print()
print('原文:'); print(ps[2][:600])
print()
RT = toks(open('p0/e185_pdf_default.txt', encoding='utf-8', errors='replace').read())
i = RT.index('semilog') if 'semilog' in RT else -1
print('渲染侧 semilog 附近:')
print(' '.join(RT[i - 20:i + 40]) if i >= 0 else 'ABSENT')
print()
for n in (8, 6, 5, 4):
    for off in (0, 1, 2, 3):
        h = anchors(S[off:off + n], RT, n=n)
        if h:
            print('off=%d n=%d -> %d 个锚点 首个@%d' % (off, n, len(h), h[0]))
print('（无输出 = 没有任何前缀命中）')

align_from = ns['align_from']
start = 3719 - 1
print('start=%d  RT[start]=%r  RT[start-2:start+4]=%s'
      % (start, RT[start], ' '.join(RT[start - 2:start + 4])))
unm, jumps = align_from(S, RT, start, hi=start + len(S) + 90)
print('align_from -> 未对齐=%d  样本=%s' % (len(unm), ' '.join(unm[:10])))
for off in (0, 1, 2, 3):
    for n in (8, 6, 5, 4):
        h = anchors(S[off:off + n], RT, n=n)
        if not h:
            continue
        st = h[0] - off
        u, j = align_from(S, RT, st, hi=st + len(S) + 90)
        print('probe off=%d n=%d start=%d -> 未对齐=%d' % (off, n, st, len(u)))
        break
    else:
        continue
    break

