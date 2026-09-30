# -*- coding: utf-8 -*-
r"""e178：把 e177 的 M3（$2/3$ 支，$\Delta_2\propto(\lambda_3/\lambda_1)^2$）**拿到 C 的 $1/2$ 支上对账**，
并交 C 在 81-A 点名要的那张表（"best known upper bound 逐格取两车道最小值"，板上写的是"**待你的表**"）。

数据来源（全部只读落盘日志，#45/#51：不转录板、不引用对方手写读数）：
  我方 $\lambda$ 与 $\Delta_2$：`p0/e177_out.txt`（9 格，seed 20260930）、`p0/e175_out.txt`（5 格，seed 777）、
        `p0/e176_out.txt`（[Y1] 11 格：C 的读数 + 我方原变量 $S$ 与 KKT 阵 $\Lambda$ 的 $\lambda_2/\lambda_1$）
  C 方 $\Delta$：`.work3/c90b_out.txt`、`.work3/c91_out.txt`、`.work3/c92_out.txt`、`.work3/c94_out.txt`（只读）

预注册判据（跑前写死）：
 [V0] **我方支内 M3 的稳健性**：照 C 补记二十一 82-B 新采纳的第 6 条（$\delta_{\rm sw}=0$ 的格不许充当判据证据），
      把 $34.41$ 剔出后重算 M3 的极差；判法：剔后仍 $\mathrm{spr}\le 2N$ $\\Rightarrow$ M3 **不靠阈值格**；否则 M3 降级。
 [V1] **有效指数**：$p=\Delta\ln\Delta_2/\Delta\ln(\lambda_3/\lambda_1)$ 的（i）全程端点值（剔阈值格）与（ii）末端 3 段中位。
      判法：$|p_{\rm near}-2|\le0.3$ 且 $|p_{\rm all}-2|>0.3$ $\\Rightarrow$ "平方律是**阈值端渐近**形式，全域有效指数偏大"。
 [V2] **跨支对账（本轮的实质问题）**：在 $1/2$ 支上用**我方** e176 的 $\lambda_2/\lambda_1(\Lambda)$ 列（不用 C 转录）
      与 C 每个搜索器自己的 $\Delta$ 算 $c=\Delta/(\lambda_2/\lambda_1)^2$。取 $c_{\min}^{\rm mine}$（我方剔阈值格的最小 prefactor）
      为**保守**下沿，对每格算 $\hat\Delta=c_{\min}^{\rm mine}(\lambda_2/\lambda_1)^2$：
      该格 C 的**全部**读数 $<\hat\Delta$ $\\Rightarrow$ "**低于保守预测**"（要么 prefactor 支依赖，要么该读数受 79-A 非稳定根影响）；
      读数把 $\hat\Delta$ **夹住** $\\Rightarrow$ "**该格不可判**"，并印 C 自身的分歧倍数。
 [V3] 取最小表：$r{=}2$ 每格列两车道全部读数（带搜索器名与 $k/n$ 覆盖），逐格取最小并标出处。
 [V4] 条数门（#46）：每条判决行印有效格数与退出格数。
"""
import re
import numpy as np
import sys
sys.stdout.reconfigure(encoding='utf-8')

out = []


def p(*a):
    s = ' '.join(str(x) for x in a)
    out.append(s)
    print(s)
    sys.stdout.flush()


def rd(path):
    with open(path, 'rb') as f:
        return f.read().decode('utf-8')


E177 = rd('p0/e177_out.txt')
E176 = rd('p0/e176_out.txt')
E175 = rd('p0/e175_out.txt')

# ---- 我方 9 格（e177）----
MINE = {}
for m in re.finditer(r'^\s*(\d+\.\d+)\s+(\d+\.\d+)\s+(\d)\s+([\d.]+)\s+([\d.]+)\s+([\d.eE+-]+)\s+([\d.eE+-]+)'
                     r'\s+([+-]?[\d.eE+-]+)\s+(\S+)\s+(\d+)/25\s+([\d.eE+-]+)\s*$', E177, re.M):
    MINE[float(m.group(1))] = dict(ds=float(m.group(2)), rstar=int(m.group(3)),
                                   lam31=float(m.group(7)), d2=float(m.group(8)),
                                   q=m.group(9), k=int(m.group(10)))
assert len(MINE) == 9, ('e177 应为 9 格，实际 %d' % len(MINE))
# ---- 我方 seed 777（e175，秩-2 行）----
MINE175 = {}
for m in re.finditer(r'^\s*(\d+\.\d+)\s+(\d)\s+(\d)\s+([\d.]+)\s+([+-][\d.eE+-]+)\s+', E175, re.M):
    if int(m.group(3)) == 2:
        MINE175[float(m.group(1))] = float(m.group(5))
assert len(MINE175) == 5, ('e175 秩-2 行应为 5，实际 %d' % len(MINE175))
# ---- e176 [Y1]：11 格的 $\lambda_2/\lambda_1$（C 读数 / 我方 S / 我方 Lam）----
Y1 = {}
for m in re.finditer(r'^\s*(\d+\.\d+)\s+([\d.eE+-]+)\s+([\d.eE+-]+)\s+([\d.eE+-]+)\s+([\d.]+)\s+([\d.]+)\s+(\S.*)$',
                     E176, re.M):
    Y1[float(m.group(1))] = dict(c=float(m.group(2)), s=float(m.group(3)), lam=float(m.group(4)), obj=m.group(7).strip())
assert len(Y1) == 11, ('e176 [Y1] 应为 11 格，实际 %d' % len(Y1))
assert all(v['obj'].startswith('- Lam') for v in Y1.values()), 'e176 相符列不是 Lam，对象识别要重跑'

# ---- C 的读数：逐文件按各自列序解析；用 I_unc 与我方 SDP 交叉核；剔非有限 ----
IUNC = {33.00: 6.172695, 33.50: 5.593895, 34.00: 5.169971, 34.41: 4.893329,
        42.00: 3.004696, 44.00: 2.795086, 45.00: 2.705659, 46.00: 2.624326, 48.00: 2.481442,
        50.00: 2.359474, 52.00: 2.253741, 54.00: 2.160951, 55.00: 2.118635, 56.00: 2.078709}
CD = {}
CHECKED = []


def grab(tag, path, pat, rpos, ipos, dpos, npos, r_fixed=None):
    txt = rd(path)
    for m in re.finditer(pat, txt, re.M):
        g = m.groups()
        D = float(g[0])
        r = r_fixed if r_fixed is not None else int(g[rpos])
        iu = float(g[ipos])
        if D in IUNC:
            CHECKED.append((tag, D, iu))
            assert abs(iu - IUNC[D]) < 5e-4, ('I_unc 与我方 SDP 不符 $\\Rightarrow$ 列序看错了', tag, D, iu, IUNC[D])
        dl = float(g[dpos].replace('+', ''))
        if not np.isfinite(dl) or dl <= 0:
            continue
        CD.setdefault((D, r), []).append((tag, dl, g[npos].replace(' ', '')))


grab('c90b', '.work3/c90b_out.txt',
     r'^\s*(\d+\.\d+)\s+(\d)\s+(\d)\s+([\d.]+)\s+(\S+)\s+([+-][\d.eE+-]+)\s+(\d+/\s*\d+)', 1, 3, 5, 6)
grab('c91', '.work3/c91_out.txt',
     r'^\s*(\d+\.\d+)\s+(\d)\s+(\d)\s+([\d.]+)\s+([\d.]+)\s+([+-][\d.eE+-]+)\s+(\d+/\s*\d+)', 2, 3, 5, 6)
grab('c92', '.work3/c92_out.txt',
     r'^\s*(\d+\.\d+)\s+(\d)\s+(\d)\s+([\d.]+)\s+(\S+)\s+([+-][\d.eE+-]+)\s+(\d+/\s*\d+)', 2, 3, 5, 6)
# c94 的行没有 $r$ 列（C 自述其覆盖 $r{=}1$、$D\in[42,56]$ 一支）
grab('c94', '.work3/c94_out.txt',
     r'^\s*(\d+\.\d+)\s+([\d.]+)\s+([+-][\d.eE+-]+)\s+(\d+/\s*\d+)', None, 1, 2, 3, r_fixed=1)
assert len(CHECKED) >= 20, ('I_unc 交叉核太少：%d' % len(CHECKED))

p('== 解析条数（#46）：我方 e177 %d 格 / e175 %d 格 / e176 [Y1] %d 格；C 方 %d 个 $(D,r)$ 格、有限读数 %d 条；'
  'I_unc 与我方 SDP 交叉核 %d 次（差 $<5\\times10^{-4}$）=='
  % (len(MINE), len(MINE175), len(Y1), len(CD), sum(len(v) for v in CD.values()), len(CHECKED)))
p('   C 的逐格读数（只读日志、非板转录）：')
for (D, r), lst in sorted(CD.items()):
    p('   %6.2f r=%d  %s' % (D, r, '  '.join('%s %+.3e (%s)' % (t, d, c) for t, d, c in lst)))

# ================= [V0] =================
p('')
p('== [V0] 我方 $2/3$ 支 M3 的 prefactor（$c=\\Delta_2/(\\lambda_3/\\lambda_1)^2$）与剔阈值格的稳健性 ==')
allc = {t: v['d2'] / v['lam31'] ** 2 for t, v in sorted(MINE.items())}
for t in sorted(allc):
    p('   $D=%5.2f$  $\\lambda_3/\\lambda_1=%9.3e$  $\\Delta_2=%9.3e$  $c=%8.1f$  %s  $\\delta_{sw}=%.2f$'
      % (t, MINE[t]['lam31'], MINE[t]['d2'], allc[t],
         '(阈值格，按 C 第 6 条剔出)' if MINE[t]['ds'] == 0.0 else ' ' * 22, MINE[t]['ds']))
spr_all = max(allc.values()) / min(allc.values())
drop0 = {k: v for k, v in allc.items() if MINE[k]['ds'] > 0.0}
spr_drop = max(drop0.values()) / min(drop0.values())
N177 = float(re.search(r'\$N=([\d.]+)\$（我方', E177).group(1))
assert N177 == 1.53, ('[Z0] 的 $N$ 不是 1.53', N177)
p('   含阈值格 $\\mathrm{spr}=%.2f$；**剔除 $\\delta_{sw}=0$ 的 $34.41$** 后 $\\mathrm{spr}=%.2f$（$2N=%.2f$）$\\Rightarrow$ %s'
  % (spr_all, spr_drop, 2 * N177,
     '**M3 不靠阈值格**（照 C 补记 82-B 第 6 条仍过线）' if spr_drop <= 2 * N177 else 'M3 **降级**（过线只靠阈值格）'))
near = {k: v for k, v in drop0.items() if MINE[k]['ds'] <= 0.41}
p('   近阈值子集（$\\delta_{sw}\\le0.41$，%d 格）$c\\in[%.0f,%.0f]$ 带宽 $%.2f\\times$；远端（$\\delta_{sw}\\ge0.66$）最大 $c=%.0f$'
  % (len(near), min(near.values()), max(near.values()), max(near.values()) / min(near.values()),
     max(v for k, v in drop0.items() if MINE[k]['ds'] > 0.41)))

# ================= [V1] =================
p('')
p('== [V1] 有效指数 $p=\\Delta\\ln\\Delta_2/\\Delta\\ln(\\lambda_3/\\lambda_1)$ ==')
ks = sorted(drop0)


def slope(a, b):
    return (np.log(MINE[b]['d2']) - np.log(MINE[a]['d2'])) / (np.log(MINE[b]['lam31']) - np.log(MINE[a]['lam31']))


p_all = slope(ks[0], ks[-1])
segs = [slope(a, b) for a, b in zip(ks, ks[1:])]
p_near = float(np.median(segs[-3:]))
p('   全程（$%.2f\\to%.2f$，两端都剔阈值格）$p=%.2f$；逐段（%d 段，沿 $D\\uparrow$）%s；末端 3 段中位 $p=%.2f$'
  % (ks[0], ks[-1], p_all, len(segs), ' '.join('%.2f' % v for v in segs), p_near))
p('   [V1] 判法 $|p_{\\rm near}-2|\\le0.3$ 且 $|p_{\\rm all}-2|>0.3$：分别 %.2f / %.2f $\\Rightarrow$ %s'
  % (abs(p_near - 2), abs(p_all - 2),
     '**平方律是阈值端的渐近形式，全域有效指数偏大（$p_{\\rm all}\\approx%.2f$）**' % p_all
     if abs(p_near - 2) <= .3 and abs(p_all - 2) > .3
     else ('**全域平方律**' if abs(p_all - 2) <= .3 else '**两都不像 2，M3 只算噪声内过线、不升格为指数结论）**')))

# ================= [V2] =================
p('')
p('== [V2] 把 M3 拿到 C 的 $1/2$ 支：$c=\\Delta/(\\lambda_2/\\lambda_1(\\Lambda))^2$，$\\lambda$ 列用我方 e176 读数 ==')
cm_min = min(drop0.values())
p('   保守下沿 = 我方（剔阈值格）最小 prefactor $c_{\\min}^{\\rm mine}=%.1f$（$\\Rightarrow$ $\\hat\\Delta$ 是**下沿预测**）' % cm_min)
rows = []
CC = {}
for (D, r), lst in sorted(CD.items()):
    if r != 1 or D not in Y1:
        continue
    lam = Y1[D]['lam']
    dh = cm_min * lam ** 2
    CC[D] = [(t, d / lam ** 2) for t, d, c in lst]
    below = [x for x in lst if x[1] < dh]
    above = [x for x in lst if x[1] >= dh]
    spread = max(x[1] for x in lst) / min(x[1] for x in lst)
    verdict = ('**全部读数低于下沿预测**' if not above else
               ('读数把预测夹住 $\\Rightarrow$ **该格不可判**' if below else '全部高于预测 $\\Rightarrow$ 与下沿相容'))
    rows.append((D, lam, dh, spread, verdict, below, above))
    p('   $D=%5.2f$  $\\lambda_2/\\lambda_1=%.3e$  $\\hat\\Delta=%.3e$  C：%s  自身分歧 $%.1f\\times$  $\\Rightarrow$ %s'
      % (D, lam, dh, '  '.join('%s %+.3e ($c=%.0f$, %s)' % (t, d, d / lam ** 2, c) for t, d, c in lst), spread, verdict))
nhard = len([x for x in rows if x[5] and not x[6]])
Dall = sorted(CC)
cflat = [v for D in Dall for _, v in CC[D]]
c_lo, c_hi = min(cflat), max(cflat)
p('   [V2] 进检格 %d、"全部读数低于我方下沿"的格 %d、被同格两族分歧夹住的格 %d（$D=55.00$：$52.9\\times$）'
  % (len(rows), nhard, len([x for x in rows if x[5] and x[6]])))
p('   **[V2-更正（该我认，跑后自查发现判据方向写反）]**：$c=\\Delta_{\\rm ub}/\\lambda^2$ 里的 $\\Delta$ 是**上界** $\\Rightarrow$')
p('   "$c$ 比我方下沿小"**只表示那一格的上界更紧**，**不构成对平方律的反证**（上界不能否证一个上界型律）。')
lowD = [x[0] for x in rows if x[5] and not x[6]]
lowc = [v for D in lowD for _, v in CC[D]]
p('   所以那 %d 格（$D\\in\\{%s\\}$）不作反证用，只登记为"$1/2$ 支这些格的 true prefactor 的**上界** $c\\in[%.0f,%.0f]$，'
  '比我方 $2/3$ 支的 $[%.0f,%.0f]$ 低一个档"'
  % (nhard, ', '.join('%.2f' % v for v in lowD), min(lowc), max(lowc), min(drop0.values()), max(drop0.values())))
p('   真正可判的跨支比较必须**同族工具**（同一搜索器、同一种子流），否则比的是搜索质量而不是律 $\\Rightarrow$ 交 e179 预注册。')
p('   对上界合法的那一条仍然成立：C 方 $1/2$ 支的 $c$ 全体极差 $=%.0f\\times$（同格两族 $52.9\\times$）$\\Rightarrow$ 常数性在 C 的数据层面不成，'
  % (c_hi / c_lo))
p('   但这只说明"C 的上界还没紧到能测指数"，**不说明平方律被否证**（与上面同一条逻辑）。')
c94 = [v for D in Dall for t, v in CC[D] if t == 'c94']
p('   另记一条可直接给 C 的核对：c94 的 $c$ 沿 $D\\uparrow$ 为 %s $\\Rightarrow$ **在阈值端反而涨**，'
  % ' '.join('%.0f' % v for v in c94))
p('   与我方支的单调降（$333\\to117$）反向 $-$ $-$ 按 C 79-A 的守卫（负/异常 $\\Delta$ 作废）方向上这**不是**同一类错误（$\\Delta>0$），')
p('   但它是"覆盖 $1$–$5/84$ 的上界随阈值升高"这种**采样太少**的签名，与我 e177 [Z0] 的 $N=1.53$ 对照即知差多少倍。')

# ================= [V3] =================
p('')
p('== [V3] 81-A 要的"逐格取最小"表（$r{=}2$；$\\Delta$ 一律是**上界**；两车道 + 我方两种子并列）==')
p('   %6s %11s %11s %16s %20s %11s   %s' % ('D', 'e175(777)', 'e177(0930)', 'c90b', 'c91 / c92', '本格最小', '出处（最松者的倍数）'))
win = 0
for t in sorted(MINE):
    cvals = {}
    for (D, r), lst in CD.items():
        if D == t and r == 2:
            for tag, dl, cov in lst:
                cvals[tag] = (dl, cov)
    allv = []
    if t in MINE175:
        allv.append(('e175', MINE175[t], '25 起点'))
    allv.append(('e177', MINE[t]['d2'], '%d/25' % MINE[t]['k']))
    allv += [(k, v[0], v[1]) for k, v in sorted(cvals.items())]
    bt, bv, bc = min(allv, key=lambda z: z[1])
    wt, wv, wc = max(allv, key=lambda z: z[1])
    if bt.startswith('e1'):
        win += 1
    p('   %6.2f %11s %11s %16s %20s %11.3e   %s %s（%s 的 %.1f$\\times$）'
      % (t,
         '%.3e' % MINE175[t] if t in MINE175 else '——',
         '%.3e' % MINE[t]['d2'],
         '%.3e (%s)' % cvals['c90b'] if 'c90b' in cvals else '——',
         ' / '.join('%.3e' % cvals[k][0] for k in ('c91', 'c92') if k in cvals) or '——',
         bv, bt, bc, wt, wv / bv))
p('   有效格 %d，其中最小读数出自我方 %d 格 $-$ $-$ #46 条数；$r{=}1$ 支不入本表（我方本轮没测 $\\Delta_1$）' % (len(MINE), win))

p('')
p('== [V4] 条数与退出 ==')
p('   我方 %d 格全进 [V0]/[V1]（[V0] 再剔 1 个阈值格）；C 方 $(D,r)$ 格 %d、有限读数 %d 条；进 [V2] 的 $r{=}1$ 格 %d'
  % (len(MINE), len(CD), sum(len(v) for v in CD.values()), len(rows)))
open('p0/e178_out.txt', 'w', encoding='utf-8').write('\n'.join(out) + '\n')
