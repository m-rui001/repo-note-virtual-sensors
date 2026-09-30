r"""E68b = 注记 Fig.2 的图：从 wall_approach.csv 画"墙是水平渐近线 + 原文窗口读不到墙"。"""
import sys
sys.stdout.reconfigure(encoding='utf-8')
import csv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

rows = list(csv.DictReader(open('p0/note/data/wall_approach.csv')))
fam = {}
for x in rows:
    fam.setdefault(x['family'], []).append((float(x['I']), float(x['rel_gap_pct'])))

fig, ax = plt.subplots(figsize=(3.4, 2.3))
sty = {'schur_tail_2': ('C0', '-'), 'schur_tail_3': ('C1', '-'),
       'best_support_r1': ('C2', '--'), 'best_support_r2': ('C3', '--'),
       'best_support_r3': ('C4', '-.')}
for k, pts in sorted(fam.items()):
    pts.sort()
    c, ls = sty.get(k, ('k', ':'))
    lab = k.replace('best_support_', 'optimal support, r=').replace('schur_tail_', 'Schur tail ')
    ax.semilogy([p[0] for p in pts], [max(p[1], 1e-4) for p in pts], ls,
                color=c, lw=1.1, label=lab)
ax.axvspan(1.0, 5.0, color='0.85', zorder=0)
ax.axhline(0.28, color='k', lw=0.7, ls=(0, (2, 2)))
ax.text(19.4, 0.32, r'digitisation band ($\pm0.09$ in $D$)', ha='right', fontsize=6)
ax.set_xlabel(r'rate $I$ [bit/sample]', fontsize=8)
ax.set_ylabel(r'$(D-D_{\min})/D_{\min}$ [\%]', fontsize=8)
ax.tick_params(labelsize=7)
ax.set_xlim(1, 20)
ax.set_ylim(1e-3, 1.4e2)
ax.legend(fontsize=5.6, loc='upper right', framealpha=0.95)
fig.tight_layout(pad=0.3)
fig.savefig('p0/note/fig_wall.pdf')
print('written p0/note/fig_wall.pdf')
