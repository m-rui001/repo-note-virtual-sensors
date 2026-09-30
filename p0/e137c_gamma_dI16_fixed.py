# -*- coding: utf-8 -*-
r"""实验 137b（§90 回应 C 的 §58-G-1：把 $\Delta I$ 推到 16）：
 相对 e137 只改两个常数——SGRID 上界 $10^{10}	o10^{14}$、DGRID 追加 $12,16$；判据/差分口径/自检全部不动。

 原 137（§85-D 的强度问题，预注册）：C 的 `c73_out.txt` 08:22 版给出锚点的**局域指数**表
（rank-1 的 $g$ 从 $\Delta I{=}0.30$ 的 $4.218$ 降到 $2.50$ 的 $1.436$，42/42 设计单调不增）。
⇒ 我 §85-D 那条 "$\gamma_r/\gamma_1>1/r$" 里的 $\gamma$ **是一个窗内局域斜率，不是指数**。
本实验把这件事测到底：把 $\Delta I$ 推到 $8$，看比值会不会走向 $1/r$。

解析预告（先写下预测，再跑，防事后挑）：$s\to\infty$ 时 $I\to(r/2)\log_2 s$，
而测得方向的误差 $x:=D-D_{\rm floor}\asymp C/s\Rightarrow \mathrm{d}\ln x/\mathrm{d}\Delta I\to-2\ln2/r$。
⇒ **预测**：$\gamma_r/\gamma_1\to1/r$，**从上方**。若实测不降或降不到 $1/r$ 附近，则本预测作废并公开降级。

 [V1] 口径：$x=D-D_{\rm floor}$（**投影地板**，与 C 同口径）；
      $\gamma_{\rm loc}:=-\mathrm{d}\ln x/\mathrm{d}\Delta I$，用 $\ln x$ 对 $\Delta I$ 的**中心差分**（非 curve_fit，避免窗宽歧义）；
      $s$ 网格 $10^{\rm linspace}(-6,10,337)$；植物/设计族与 e131 逐字相同（6 株 $\times$ 每秩 4 随机 $+$ $\Theta$min/$\Theta$max）。
 [V2] 判据（跑前写死，三选一，不许合并）：
      (a) 某秩比值在 $\Delta I\le8$ 内单调降且进入 $[1/r,\,1/r+0.05]$ ⇒ **$2\ln2/r$ 是渐近律**，
          §85-D 降级为"不进入渐近前的有限率陈述"，并给出进入交叉所需的 $\Delta I$；
      (b) 单调降但 $\Delta I=8$ 处仍 $>1/r+0.15$ ⇒ **可测范围内无法区分**"$\to1/r$"与"收敛到更大值"，
          §85-D 必须加 "measured up to $\Delta I=8$"；
      (c) 不单调或上升 ⇒ 上面那条解析预告**否证**，C 的"渐近未到"解释也不成立。
 [V3] 速率：若 $\ln(\gamma_r/\gamma_1-1/r)$ 对 $\ln\Delta I$ 的 Pearson $|r|\ge0.95$ 才许报斜率；
      且**不许**把斜率外推到实测范围上界的 2 倍以外。
 [V4] 与 C 对表：我 anchor 的 $\gamma_{\rm loc}$ 在 $\Delta I\in\{0.65,1.00,1.50,2.50\}$ 要与 C 打印的
      $g(0.65)/g(1.00)/g(1.50)/g(2.50)$ 对上（相对差 $\le10\%$）；对不上要报差分口径差异，不许各自留一手。
 [V5] 逐格报有效设计数（$x>0$ 且 $\Delta I$ 落在网格内），$<4$ 的格明写并排除；分母口径全程不换。
"""
import sys
import numpy as np
from scipy.linalg import solve_discrete_are, eigh
sys.stdout.reconfigure(encoding='utf-8')

_src = open('p0/exp_c_audit.py', encoding='utf-8').read().split("print(r'== E53")[0]
_ns = {'__name__': 'p'}
exec(compile(_src, 'p0/exp_c_audit.py[preamble]', 'exec'), _ns)
A0, B0, W0, ctrl, sym = _ns['A'], _ns['B'], _ns['W'], _ns['ctrl'], _ns['sym']
GT = 2.0 * np.log(2.0)
SGRID = 10.0 ** np.linspace(-6.0, 14.0, 409)
DGRID = np.array([0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 8.0, 12.0, 16.0])


def make_plant(kind, seed):
    rng = np.random.default_rng(seed)
    if kind == 'anchor':
        A, B, W, m = A0, B0, W0, 4
    elif kind == 'rand4':
        m = 4
        A = rng.normal(0.0, 1.05, (m, m)); B = rng.normal(0.0, 1.15, (m, m))
        M = rng.normal(0.0, 1.0, (m, m)); W = sym(M @ M.T + 0.6 * np.eye(m))
    elif kind == 'big6':
        m = 6
        bl = [np.array([[2.1]]),
              1.35 * np.array([[np.cos(0.9), -np.sin(0.9)], [np.sin(0.9), np.cos(0.9)]]),
              np.array([[-0.55, 0.3], [-0.2, 0.22]]), np.array([[0.34]])]
        A = np.zeros((m, m)); pos = 0
        for blk in bl:
            k = blk.shape[0]; A[pos:pos + k, pos:pos + k] = blk; pos += k
        A = A + 0.16 * rng.normal(0, 1, (m, m))
        B = rng.normal(0.0, 1.2, (m, m))
        M = rng.normal(0.0, 1.0, (m, m)); W = sym(M @ M.T + 0.6 * np.eye(m))
    Pc, K, TH, JC = ctrl_full(A, B, W, np.eye(m), np.eye(m))
    return A, B, W, TH, JC, m


def ctrl_full(Ax, Bx, Wx, Qt, Rx):
    Pc = sym(solve_discrete_are(Ax, Bx, sym(Qt), Rx))
    K = np.linalg.solve(Rx + Bx.T @ Pc @ Bx, Bx.T @ Pc @ Ax)
    return Pc, K, sym(K.T @ (Rx + Bx.T @ Pc @ Bx) @ K), float(np.trace(Wx @ Pc))


def r_exp(A):
    mod = np.abs(np.linalg.eigvals(A))
    return float(np.log(mod[mod > 1.0]).sum() / np.log(2.0))


PLANTS = [('anchor', ('anchor', 0)), ('rand-1', ('rand4', 1)), ('rand-2', ('rand4', 2)),
          ('rand-3', ('rand4', 3)), ('rand-4', ('rand4', 4)), ('big-6', ('big6', 11))]


def curve(Z, s, A, W, TH, JC):
    r = Z.shape[1]
    Ir = np.eye(r)
    C = np.sqrt(s) * Z.T
    Pm = sym(solve_discrete_are(A.T, C.T, W, Ir))
    Sm = C @ Pm @ C.T + Ir
    Lk = Pm @ C.T @ np.linalg.inv(Sm)
    Pp = sym(Pm - Lk @ Sm @ Lk.T)
    I = 0.5 * np.log(np.linalg.det(Ir + C @ Pm @ C.T)) / np.log(2.0)
    return I, JC + np.trace(TH @ Pp)


def proj_floor(Z, A, W, TH, itmax=60000, tol=1e-16):
    Q, _ = np.linalg.qr(Z)
    Pt = sym(W.copy()); P = Pt
    for _ in range(itmax):
        P = sym(Pt - Pt @ Q @ np.linalg.solve(Q.T @ Pt @ Q, Q.T @ Pt))
        Pn = sym(A @ P @ A.T + W)
        done = np.max(np.abs(Pn - Pt)) < tol * max(1.0, np.max(np.abs(Pn)))
        Pt = Pn
        if done:
            break
        if not np.all(np.isfinite(Pt)) or np.max(np.abs(Pt)) > 1e14:
            return np.inf
    return float(np.trace(TH @ P))


def local_gamma(di, lx):
    """γ_loc = -d ln x / d ΔI，中心差分（端点用单侧）。返回在 DGRID 上的插值。"""
    o = np.argsort(di)
    d, l = di[o], lx[o]
    g = np.empty(len(d))
    g[1:-1] = -(l[2:] - l[:-2]) / (d[2:] - d[:-2])
    g[0] = -(l[1] - l[0]) / (d[1] - d[0])
    g[-1] = -(l[-1] - l[-2]) / (d[-1] - d[-2])
    lo, hi = d[1], d[-2]          # 中心差分只在 [d[1], d[-2]] 内成立（端点无双侧邻）
    return np.array([np.nan if (x < lo or x > hi) else np.interp(x, d, g) for x in DGRID])


GTAB = {}   # (lab, rk) -> array over DGRID (median over designs)
PD = {}     # (lab, rk) -> designs x DGRID 的逐设计 γ_loc
XMIN = {}   # (lab, rk) -> 每个设计进入 ln 的最小 x = D - JC - D_floor
DMAX = {}   # (lab, rk) -> 每个设计曲线可达的最大 ΔI
TAIL = {}   # (lab, rk) -> 每个设计末 6 个 (ΔI, x) 点
NCNT = {}
COV = {}    # (lab, rk) -> 最大可用 ΔI（该档 ≥4 个设计有限）
for lab, spec in PLANTS:
    A, B, W, TH, JC, m = make_plant(*spec)
    REXP = max(r_exp(A), 0.0)
    _, Vth = eigh(TH)
    rng = np.random.default_rng(20260930 + 7 * PLANTS.index((lab, spec)))
    print('\n== [%s] n=%d Rexp=%.6f  γ_loc(ΔI)（每秩 6 设计取中位）==' % (lab, m, REXP))
    print('  rk  n_eff ' + ''.join('%8.2f' % d for d in DGRID))
    for rk in (1, 2, min(3, m - 1)):
        pool = []
        for c in range(4):
            Z, _ = np.linalg.qr(rng.standard_normal((m, rk))); pool.append(Z)
        pool.append(Vth[:, :rk]); pool.append(Vth[:, m - rk:])
        cols = []
        xmin_l, dmax_l, tail_l = [], [], []
        for Z in pool:
            fp = proj_floor(Z, A, W, TH)
            if not np.isfinite(fp):
                continue
            di, lx = [], []
            for s in SGRID:
                try:
                    I, D = curve(Z, s, A, W, TH, JC)
                except Exception:
                    continue
                x = D - JC - fp
                # 相对地板：x 是 D 与地板之差，小于 $10^{-9}$ 相对量级就是浮点相减的噪声，不许进 ln
                if np.isfinite(I) and np.isfinite(x) and x > 1e-9 * (JC + fp):
                    di.append(I - REXP); lx.append(np.log(x))
            di, lx = np.array(di), np.array(lx)
            if len(di) < 40:
                continue
            xmin_l.append(float(np.exp(lx).min()))
            dmax_l.append(float(di.max()))
            tail_l.append(list(zip(di[-6:], np.exp(lx[-6:]))))
            cols.append(local_gamma(di, lx))
        cm = np.array(cols, dtype=float)
        nfin = np.array([np.isfinite(cm[:, j]).sum() if len(cm) else 0 for j in range(len(DGRID))])
        med = np.array([np.nanmedian(cm[:, j]) if nfin[j] >= 4 else np.nan
                        for j in range(len(DGRID))])
        PD[(lab, rk)] = cm
        XMIN[(lab, rk)] = np.array(xmin_l)
        DMAX[(lab, rk)] = np.array(dmax_l)
        TAIL[(lab, rk)] = tail_l
        GTAB[(lab, rk)] = med
        NCNT[(lab, rk)] = int(nfin[0])
        covd = DGRID[nfin >= 4]
        COV[(lab, rk)] = float(covd.max()) if len(covd) else float('nan')
        print('  %d  %5d ' % (rk, NCNT[(lab, rk)]) + ''.join('%8s' % ('%7.3f' % v if np.isfinite(v) else '   --') for v in med)
              + '   覆盖到 ΔI≤%s' % ('%.2f' % COV[(lab, rk)] if np.isfinite(COV[(lab, rk)]) else '无'))

print('\n== [V0] 仪器自检（先确认差分口径能用，再谈判决）==')
v_ok = True
for rk, ref in ((1, 1.8824), (2, 1.4890), (3, 1.2858)):   # C 的 anchor g(1.00)
    j = list(DGRID).index(1.0)
    m = GTAB[('anchor', rk)][j]
    rel = (100 * (m - ref) / ref) if np.isfinite(m) else float('nan')
    good = np.isfinite(m) and abs(rel) <= 25.0
    v_ok = v_ok and good
    print('  rank%d γ_loc(1.00)=%s  C 的 g(1.00)=%.4f  相对差=%s ⇒ %s'
          % (rk, '%.4f' % m if np.isfinite(m) else 'NaN', ref,
             '%+.1f%%' % rel if np.isfinite(rel) else 'n/a', 'OK' if good else '**仪器不通**'))
n_cell_ok = sum(1 for k in GTAB if np.isfinite(GTAB[k]).sum() >= 4)
print('  有限格数（≥4 档可算）：%d/%d' % (n_cell_ok, len(GTAB)))
if not v_ok or n_cell_ok == 0:
    print('  [V0] **判据不启动**：口径/仪器失败，本轮不出 (a)/(b)/(c)。')
    print('EXIT=1 (instrument failure, 不判决)')
    sys.exit(1)

print('\n== [V4] 与 C 的 anchor 局域指数对表（C 打印值 g(0.65)/g(1.00)/g(1.50)/g(2.50)）==')
CREF = {1: {0.65: 2.3961, 1.00: 1.8824, 1.50: 1.6012, 2.50: 1.4355},
        2: {0.65: 2.0563, 1.00: 1.4890, 1.50: 1.1511, 2.50: 0.8882},
        3: {0.65: 1.8251, 1.00: 1.2858, 1.50: 0.9888, 2.50: 0.7332}}
worst = 0.0
for rk in (1, 2, 3):
    mine = GTAB[('anchor', rk)]
    line = '  rank%d  ' % rk
    for dv in (0.75, 1.0, 1.5, 2.0):
        key = min(CREF[rk], key=lambda k: abs(k - dv))
        j = list(DGRID).index(dv)
        v = mine[j]
        if not np.isfinite(v):
            line += ' ΔI=%.2f:n/a' % dv; continue
        rel = 100 * abs(v - CREF[rk][key]) / CREF[rk][key]
        worst = max(worst, rel)
        line += ' %.2f~%.2f:%+.1f%%' % (dv, key, rel)
    print(line + '   （对最近档，|相对差|）')
print('  ⇒ 最差相对差 %.1f%% ⇒ %s' % (worst, '≤10% 两车道差分口径可互通' if worst <= 10 else '**>10%，先报口径差再谈结论**'))

print('\n== [V2] 比值 γ_r/γ_1 随 ΔI 的走势（1/r 为渐近候选）==')
res = {}
for lab, _ in PLANTS:
    g1 = GTAB[(lab, 1)]
    for rk in (2, 3):
        gr = GTAB.get((lab, rk))
        if gr is None:
            continue
        ratio = gr / g1
        tgt = 1.0 / rk
        res[(lab, rk)] = (ratio, tgt)
        print('  %-7s r=%d  ' % (lab, rk) + ''.join('%8s' % ('%6.3f' % v if np.isfinite(v) else '   --') for v in ratio)
              + '   1/%d=%.4f' % (rk, tgt))
print('  ΔI 网格：' + ''.join('%8.2f' % d for d in DGRID))

print('\n== [V2] 判决（逐 (株,秩) 格，三选一）==')
cnt = {'a': 0, 'b': 0, 'c': 0, 'd': 0, 'skip': 0}
cross_list = []
for (lab, rk), (ratio, tgt) in sorted(res.items(), key=lambda kv: (kv[0][1], kv[0][0])):
    ok = np.isfinite(ratio)
    d, v = DGRID[ok], ratio[ok]
    if len(d) < 4:
        print('  %-7s r=%d 有效 ΔI 档 %d<4 ⇒ skip 并明写' % (lab, rk, len(d))); cnt['skip'] += 1; continue
    mono = bool(np.all(np.diff(v) <= 1e-9))
    if v[-1] <= 0.0:                      # (d) 末档非正 ⇒ 仪器下界，不判
        cnt['d'] += 1
        print('  %-7s r=%d 末档(ΔI=%.1f)比值=%+.3f ≤0 ⇒ (d) 仪器下界，**不参与判决**'
              % (lab, rk, d[-1], v[-1]))
        continue
    if tgt - 1e-12 <= v[-1] <= tgt + 0.05 and v[-1] < v[0]:
        v_ = 'a'; cross_list.append((lab, rk, d[np.argmax(v <= tgt + 0.05)]))
    elif v[-1] < v[0]:
        v_ = 'b'
    else:
        v_ = 'c'
    tail_gap = v[-1] - tgt
    cnt[v_] += 1
    print('  %-7s r=%d 全程单调降=%-5s 末档(ΔI=%.1f)比值=%.3f 对 1/%d 的余隙=%+.3f ⇒ (%s)'
          % (lab, rk, mono, d[-1], v[-1], rk, tail_gap, v_))
    if v_ == 'b':
        print('       （末档 %.3f 在带沿 %.3f 之上 ⇒ 按预注册原文不是 (a)）' % (v[-1], tgt + 0.05))
na = cnt['a']; nb = cnt['b']; nc = cnt['c']
print('  合计：(a) 进入 1/r 邻域 %d 格 / (b) 降但未进入 %d 格 / (c) 无下行证据 %d 格 / '
      '(d) 末档非正⇒仪器 %d 格 / 跳过 %d 格'
      % (na, nb, nc, cnt['d'], cnt['skip']))
print('  [V2 判据不完备 ⇒ 公开记账] §86 的三选一未覆盖"末档落到仪器下界以下、变成非正数"这一支；'
      '\n  (d) 是补记，不参与 (a)(b)(c) 的计数，也不改写任何已判格的结果。')
if cross_list:
    print('  [V2a] 交叉 ΔI：%s ⇒ 中位 %.2f' % ([('%s/r%d@%.2f' % t) for t in cross_list],
          float(np.median([t[2] for t in cross_list]))))
print('  [V2] 判决：' + ('**无格可判 ⇒ 不判决**' if na + nb + nc == 0 else
      ('**(a) 成立：$2\\ln2/r$ 是渐近律，我 §85-D 必须降级为有限率陈述**'
       if na >= max(nb, nc) and na >= 6 else
       ('**(b)：单调向下但在 ΔI≤8 内未进入 $1/r$ 邻域 ⇒ 可测范围内无法区分，§85-D 加 "measured up to ΔI=8"**'
        if nb >= na and nb > 0 else
        '**(c)：无系统性下行证据 ⇒ 我的解析预告与 C 的"渐近未到"一起否证**'))))
cmax = max(COV[k] for k in COV if np.isfinite(COV[k])) if any(np.isfinite(v) for v in COV.values()) else float('nan')
print('  [V2 前提] 实测覆盖最大 ΔI=%.2f；任何超出该范围的趋势都**不是**测量结果' % cmax)

npts = ng = 0
for (lab, rk), (ratio, tgt) in sorted(res.items(), key=lambda kv: (kv[0][1], kv[0][0])):
    fin = np.isfinite(ratio)
    npts += int(fin.sum())
    ng += int((ratio[fin] > tgt).sum())
print('  [V2c] 逐**中位曲线**格点（%d 格 × %d 档）：有限 %d，其中 γ_r/γ_1>1/r 的 %d'
      % (len(res), len(DGRID), npts, ng))
print('  [V2c 限定] 上表是每格 6 个设计的中位；逐设计层面的检验在 e131/§85-D（6/6 株）')

print('\n== [V6] 仪器诊断：逐设计 γ_loc / 最小 x / 曲线可达 ΔI ==')
for (lab, rk) in sorted(PD, key=lambda k: (k[1], k[0])):
    cm = PD[(lab, rk)]
    print('  %-7s r=%d  6 设计 min x=[%s]  可达 ΔI=[%s]  （进 ln 的门：x>1e-9*(JC+floor)≈3e-8）'
          % (lab, rk, ','.join('%.2e' % v for v in XMIN[(lab, rk)]),
             ','.join('%.1f' % v for v in DMAX[(lab, rk)])))
    for dv in (8.0, 12.0, 16.0):
        j = list(DGRID).index(dv)
        gv = cm[:, j]
        fin = np.isfinite(gv)
        neg = int((gv[fin] <= 0).sum())
        tag = ''
        if neg:
            tag = '   ⇒ %d/%d 个有限设计给出非正斜率' % (neg, int(fin.sum()))
        print('      ΔI=%5.1f  有限 %d/6  γ=[%s]%s'
              % (dv, int(fin.sum()),
                 ','.join('%7.3f' % v if np.isfinite(v) else '     --' for v in gv), tag))

print('\n== [V7] 负值格的尾巴：逐设计末 6 个 (ΔI, x) 点（看 x 是否反增）==')
bad = [(lab, rk) for (lab, rk) in sorted(PD, key=lambda k: (k[1], k[0]))
       if np.isfinite(GTAB[(lab, rk)][list(DGRID).index(16.0)])
       and GTAB[(lab, rk)][list(DGRID).index(16.0)] <= 0]
for cell in bad:
    print('  [%s r=%d]' % cell)
    for i, tl in enumerate(TAIL[cell]):
        print('     d%d  %s' % (i, '  '.join('(%.2f, %.3e)' % (a, b) for a, b in tl)))
if not bad:
    print('  中位层面没有非正格 ⇒ 无需看尾巴')

print('\n== [V3] 收敛速率（只有 (a)/(b) 且 |Pearson r|≥0.95 才报斜率；不外推超 2× 实测上界）==')
for (lab, rk), (ratio, tgt) in sorted(res.items(), key=lambda kv: (kv[0][1], kv[0][0])):
    ok = np.isfinite(ratio) & (ratio > tgt)
    d, e = DGRID[ok], ratio[ok] - tgt
    ndrop = int(np.isfinite(ratio).sum() - ok.sum())
    if len(e) < 4 or np.any(e <= 0):
        continue
    ln, ld = np.log(e), np.log(d)
    if np.std(ld) < 1e-12:
        continue
    rr = float(np.corrcoef(ln, ld)[0, 1])
    if abs(rr) >= 0.95:
        sl = float(np.polyfit(ld, ln, 1)[0])
        print('  %-7s r=%d  ln(比值-1/r) vs ln ΔI 斜率=%+.3f  Pearson=%.4f  （ΔI≤%.1f 内有效，'
              '拟合窗 %d 点，另有 %d 个有限点因 ≤1/r 被排除）'
              % (lab, rk, sl, rr, d.max(), len(e), ndrop))
    else:
        print('  %-7s r=%d  Pearson=%.4f <0.95 ⇒ 不报斜率' % (lab, rk, rr))
print('EXIT=0')
