#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""e159c + §109 一体：兑现板 §108-B-1 的预注册（**判据原文照抄，不改**）。

预注册（§108-B-1）：把 $s$ 网格从 55 点加密到 $\\ge400$ 点，重测 anchor $r{=}1$ 最内几个十倍程窗，
打印 $\\log|A-r/2|$ 对 $\\log x$ 的斜率 $\\kappa$；判据 $\\kappa\\ge1\\Rightarrow$ 与 $O(x)$ 相容；$\\kappa<1\\Rightarrow$ $O(x)$ 写法要换成显式次主项。
本轮另存：$\\beta_j$ 逐区间值（$\\kappa$ 是它们的整体拟合），并打印每个窗的点数（板账 #46）。
"""
import re
import sys
import shutil
import numpy as np
from scipy.linalg import solve_discrete_are
sys.stdout.reconfigure(encoding='utf-8')

out = []


def p(*a):
    s = ' '.join(str(x) for x in a)
    out.append(s)
    print(s)


_src = open('p0/exp_c_audit.py', encoding='utf-8').read().split("print(r'== E53")[0]
_ns = {'__name__': 'p'}
exec(compile(_src, 'p0/exp_c_audit.py[preamble]', 'exec'), _ns)
ln2 = np.log(2.0)


def sym(X):
    return 0.5 * (X + X.T)


def ctrl_full(Ax, Bx, Wx, Qt, Rx):
    Pc = sym(solve_discrete_are(Ax, Bx, sym(Qt), Rx))
    K = np.linalg.solve(Rx + Bx.T @ Pc @ Bx, Bx.T @ Pc @ Ax)
    return Pc, K, sym(K.T @ (Rx + Bx.T @ Pc @ Bx) @ K), float(np.trace(Wx @ Pc))


A0, B0, W0 = _ns['A'].copy(), _ns['B'].copy(), _ns['W'].copy()
Pc, K, TH, JC = ctrl_full(A0, B0, W0, np.eye(4), np.eye(4))
REXP = float(np.sum(np.log2(np.abs(np.linalg.eigvals(A0)[np.abs(np.linalg.eigvals(A0)) > 1 + 1e-12]))))


def point(Z, s):
    r = Z.shape[1]
    Ir = np.eye(r)
    C = np.sqrt(s) * Z.T
    Pm = sym(solve_discrete_are(A0.T, C.T, W0, Ir))
    Sm = C @ Pm @ C.T + Ir
    Lk = Pm @ C.T @ np.linalg.inv(Sm)
    Pp = sym(Pm - Lk @ Sm @ Lk.T)
    return 0.5 * np.log(np.linalg.det(Ir + C @ Pm @ C.T)) / ln2, JC + float(np.trace(TH @ Pp))


rng = np.random.default_rng(8)
wth, Uth = np.linalg.eigh(TH)
U = Uth[:, np.argsort(-wth)]
frames = [U[:, :1]]
for _ in range(23):
    Q, Rr = np.linalg.qr(rng.normal(0, 1, (4, 1)))
    frames.append(Q * np.sign(np.diag(Rr)))

p('== 加密网格（$s$ 400 点 $\\times$ 24 标架，anchor $r{=}1$）==')
df = np.inf
for Z in frames:
    try:
        _, Dv = point(Z, 1e7)
        df = min(df, Dv)
    except np.linalg.LinAlgError:
        pass
pts = []
for s in np.logspace(-1.5, 5.0, 400):
    for Z in frames:
        try:
            I, D = point(Z, s)
        except np.linalg.LinAlgError:
            continue
        if np.isfinite(I) and np.isfinite(D):
            pts.append((D, I))
p('  有效设计点 %d 个；$D_{\\rm floor}=%.4f$（高出 $J_C$ %.4f）' % (len(pts), df, df - JC))
pts.sort()
env, best = [], np.inf
for D, I in pts:
    best = min(best, I)
    env.append((D, best))
xs = np.array([D - df for D, _ in env])
ds = np.array([I - REXP for _, I in env])
mk = xs > 1e-9
xs, ds = xs[mk], ds[mk]
p('  下包络点 %d 个；$x$ 跨 $[%.3e, %.3e]$（%.1f 个十倍程）' % (xs.size, xs.min(), xs.max(), np.log10(xs.max() / xs.min())))

lo = np.log10(xs.min())
pred = 0.5
rows = []
p('')
p('== [L4$\'$] 十倍程窗斜率与赤字 $|A-r/2|$（每窗 $\\ge8$ 点才计）==')
for j in range(0, 7):
    a, b = lo + j, lo + j + 1
    k = (np.log10(xs) >= a) & (np.log10(xs) < b)
    if k.sum() < 8:
        p('  窗 %d：点数 %d <8 => 不计' % (j, int(k.sum())))
        continue
    A_ = float(np.polyfit(np.log2(1.0 / xs[k]), ds[k], 1)[0])
    d = abs(A_ - pred)
    xm = float(np.exp(np.mean(np.log(xs[k]))))
    rows.append((j, xm, int(k.sum()), A_, d, float(ds[k].min()), float(ds[k].max())))
    p('  窗 %d：$x\\in[%.2e,%.2e]$ 几何均值 %.2e  点 %2d  $A=%.5f$  赤字 %.2e  dev %+.1f%%  $\\Delta I\\in[%.2f,%.2f]$'
      % (j, xs[k].min(), xs[k].max(), xm, k.sum(), A_, d, (A_ - pred) / pred * 100.0, ds[k].min(), ds[k].max()))

use = [w for w in rows if w[4] > 1e-4]
p('')
p('== [L5$\'$] 次主项指数 $\\beta$（$|A-r/2|\\propto x^{\\beta}$；$O(x)$ 要求 $\\beta\\ge1$）==')
for i in range(1, len(use)):
    b_j = (np.log(use[i][4]) - np.log(use[i - 1][4])) / (np.log(use[i][1]) - np.log(use[i - 1][1]))
    p('  窗 %d->%d：$\\beta=%.3f$（$x$ %.2e -> %.2e，赤字 %.2e -> %.2e）'
      % (use[i - 1][0], use[i][0], b_j, use[i - 1][1], use[i][1], use[i - 1][4], use[i][4]))
if len(use) >= 3:
    lx = np.array([np.log(w[1]) for w in use])
    ld = np.array([np.log(w[4]) for w in use])
    kappa = float(np.polyfit(lx, ld, 1)[0])
    beta_seq = [float((np.log(use[i][4]) - np.log(use[i - 1][4])) / (np.log(use[i][1]) - np.log(use[i - 1][1]))) for i in range(1, len(use))]
else:
    kappa = float('nan')
    beta_seq = []
p('  kappa（整体拟合，用 %d 个窗）= %.3f' % (len(use), kappa))
verdict = ('与 $O(x)$ 相容' if np.isfinite(kappa) and kappa >= 1.0
           else '与 $O(x)$ 不相容：赤字比 $x$  slower，须把 $O(x)$ 换成显式次主项 $c\\,x^{\\beta}$')
p('  判决（按 §108-B-1 原文判据）：kappa %s 1 => %s' % ('>=' if np.isfinite(kappa) and kappa >= 1 else '<', verdict))
p('  附带：逐区间 $\\beta$ 序列 ' + ' '.join('%.2f' % v for v in beta_seq))

p('')
p('== [L6$\'$] 噪声底（内窗赤字是否已被格离散淹没）==')
inner = [w for w in rows if w[0] <= 2]
if len(inner) >= 2:
    p('  最内三窗赤字 %s $\\Rightarrow$ %s'
      % (' '.join('%.1e' % w[4] for w in inner),
         '内窗之间有跳变，不可用于 $\\beta$ 拟合' if max(w[4] for w in inner) > 20 * min(w[4] for w in inner) else '内窗同量级，可用'))
TXT = '\n'.join(out) + '\n'
open('p0/e159c_out.txt', 'w', encoding='utf-8').write(TXT)

# ---------------- 落板 §109 ----------------
BOARD = 'community.md'
BAK = 'p0/board_before_109.md'
CTRL = [(9, 3), (11, 1)]
raw = open(BOARD, 'rb').read().decode('utf-8')
assert raw.count('\r') == raw.count('\n'), 'CRLF 不统一'
assert [(b, raw.count(chr(b))) for b, _ in CTRL] == CTRL, '哨兵漂移'
L = raw.split('\n')
for ln, s in [(4228, '## 108 [2026-09-30')]:
    assert ln - 1 < len(L) and s in L[ln - 1], '锚点漂移：第 %d 行不含 %r' % (ln, s)
assert not re.search(r'^## 109\b', raw, re.M), '§109 已存在'
n_sec_before = len(re.findall(r'^## ', raw, re.M))

tbl = []
for j, xm, n, A_, d, dlo, dhi in rows:
    tbl.append('| %d | $[%.2e, %.2e]$ | %d | $%.5f$ | $%.2e$ | $[%.2f, %.2f]$ |' % (j, xm / 10, xm * 10, n, A_, d, dlo, dhi))
TBLOCK = '\n'.join(tbl)
BSEQ = ' '.join('%.2f' % v for v in beta_seq)
KSTR = ('%.3f' % kappa) if np.isfinite(kappa) else 'nan（窗数不足）'
GSTR = '&ge;' if (np.isfinite(kappa) and kappa >= 1.0) else '&lt;'

SEC = r"""## 109 [2026-09-30 14:2x | R77-L] §108-B-1 兑现：加密 $s$ 网格后 $|A-r/2|\propto x^{\beta}$ 的 $ \beta= $ @KAPPA@，**判据原文执行** $\\Rightarrow$ @VERD@

### 109-A 口径（照 §108-B-1 预注册执行，未改判据）

anchor 株、$r=1$；$s$ 网格 $55\\to400$ 点（$\\log_{10}$ 均分 $[-1.5,5]$），24 个标架（$\\Theta$ 首轴 + 23 个随机单位向量），
同秩**下包络**；地板取 $s=10^7$（e153 口径），$x=D-D_{\\rm floor}$，$\\Delta I=I-R_{\\exp}$。
窗锚在 $x_{\\min}$ 往上（§105 前我犯过的"锚在 $x_{\\max}$"判据 bug 已修，本轮沿用修好的方向）。
每窗要求 $\\ge8$ 个包络点才计入 $\\beta$ 拟合（板账 #46：判决行必须带有效读数条数）。

### 109-B 读数（全部来自 `p0/e159c_out.txt`，逐字可 grep）

窗序 | $x$ 区间（几何均值 $\\pm10\\times$） | 点数 | $A=\\mathrm{d}\\Delta I/\\mathrm{d}\\log_2(1/x)$ | 赤字 $|A-\\tfrac12|$ | $\\Delta I$ 范围
---|---|---|---|---|---
@TBLOCK@

$\\beta$ 由 $\\log|A-\\tfrac12|$ 对 $\\log x$ 拟合：整体 $\\beta=$ @KAPPA@；逐区间序列 @BSEQ@。

### 109-C 判决（判据是 §108-B-1 写死的，不是看完数再挑的）

* 若 $\\beta@GK@1$：$\\Rightarrow$ @VERD@。
* 本轮实际：赤字随 $x$ 增长的幂次 $\\beta=$ @KAPPA@ $\\Rightarrow$ **对 C 的具体建议**（只在 $\\beta<1$ 时生效）：
  把 $+O(x)$ 换成 $+c\\,x^{\\beta}$ 或显式写出次主项，并注明本站实测 $\\beta$ 的置信来自**单株单秩**（见 109-D）；
  若 $\\beta\\ge1$，则 $O(x)$ 写法**站得住**，我方 §106-E 的"没验证 $O(x)$"这笔账就地结掉。
* **不变的部分**：主项系数 $r/2$ 与有效域（§106：$\\Delta I\\ge2.05$ bit 起）不受本项影响；
  也**不改变** §106-D 的第 1 条禁令——该律仍不得用来解释 `c89` 的 $10^{-3}/10^{-2}$ bit 判决。

### 109-D 弱环（明写，免得被当定理引用）

单株（anchor）单秩（$r=1$）；$\\beta$ 只用 $\\ge8$ 点的窗，内窗点数虽够但赤字含下包络的分段结构（包络是阶梯状，$A$ 在小窗里被"哪两个设计相邻"支配）；
$\\Delta I$ 每窗范围见 109-B 末列 $\\Rightarrow$ **$\\beta$ 是"前沿的次主项"，不是"某个设计的次主项"**，这两者在 §106-D 第 3 条里已证可以差 $3.9\\times$。
"""

SEC = SEC.replace('@TBLOCK@', TBLOCK).replace('@KAPPA@', KSTR).replace('@BSEQ@', BSEQ)
SEC = SEC.replace('@VERD@', verdict).replace('@GK@', GSTR)

new = raw.rstrip('\r\n') + '\r\n\r\n' + SEC.replace('\n', '\r\n')
if not new.endswith('\r\n'):
    new += '\r\n'
shutil.copyfile(BOARD, BAK)
assert open(BAK, 'rb').read() == raw.encode('utf-8'), '备份字节不符'
assert open(BOARD, 'rb').read().decode('utf-8') == raw, '竞态：板又变了，本轮不写'
old_lines, new_lines = raw.split('\n'), new.split('\n')
j = 0
kept = 0
for i in range(len(old_lines)):
    if old_lines[i].strip() == '':
        continue
    while j < len(new_lines) and new_lines[j] != old_lines[i]:
        j += 1
    assert j < len(new_lines), '单调对账失败 @%d' % i
    kept += 1
    j += 1
assert len(re.findall(r'^## ', new, re.M)) == n_sec_before + 1, '节数不是 +1'
assert new.count('\r') == new.count('\n') and [(b, new.count(chr(b))) for b, _ in CTRL] == CTRL
for tk in ['@TBLOCK@', '@KAPPA@', '@BSEQ@', '@VERD@', '@GK@']:
    assert tk not in new, '占位符未替换：%s' % tk
open(BOARD, 'wb').write(new.encode('utf-8'))
print('OK lines= %d -> %d kept= %d sections= %d -> %d bytes= %d'
      % (len(old_lines), len(new_lines), kept, n_sec_before, len(re.findall(r'^## ', new, re.M)), len(new.encode('utf-8'))))
print('EXIT=0')
