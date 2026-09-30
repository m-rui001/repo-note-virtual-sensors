# -*- coding: utf-8 -*-
"""生成 p0/e137c_gamma_dI16_fixed.py = e137b 的四处改动，其余口径（网格、差分、地板、[V0] 自检）一律不动：

1) 判决支 (a) 按**预注册原文**收紧为双侧带 $[1/r,\,1/r+0.05]$。e137b:227 写的是
   `v[-1] <= tgt + 0.05`（只查上沿），于是末档为 $-29.7$ 这种纯噪声也被印成 "(a) 进入 1/r 邻域"
   ⇒ 判决行与自己那张表矛盾，代码账。
2) 新增支 (d)：末档 $\le0$ ⇒ 判为仪器下界、不参与判决。这一支在 §86 的三选一里**没有预注册**，
   按 [W4] 的规矩把"判据本身不完备"公开记账，但不改写任何已判格的结果。
3) [V6] 逐设计诊断：印每格 6 个设计在 $\Delta I\in\{8,12,16\}$ 的 $\gamma_{\rm loc}$、
   各设计进入 $\ln$ 的最小 $x=D-JC-D_{\rm floor}$、以及该设计曲线可达的最大 $\Delta I$。
4) [V7] 尾巴原始对：负值格逐设计印末 6 个 $(\Delta I, x)$ 点，直接看 $x$ 在尾部是否**反增**。
   （[V3] 那里另加一行"丢弃了多少非正点"，让拟合窗口的取舍可见。）
"""
SRC = 'p0/e137b_gamma_dI16.py'
DST = 'p0/e137c_gamma_dI16_fixed.py'
s = open(SRC, encoding='utf-8').read()


def rep(a, b, n=1):
    global s
    assert s.count(a) == n, '锚点 %r 命中 %d 次（应为 %d）' % (a[:48], s.count(a), n)
    s = s.replace(a, b)


# ---------- 3)+4) 逐设计记录 ----------
rep("""            di, lx = np.array(di), np.array(lx)
            if len(di) < 40:
                continue
            cols.append(local_gamma(di, lx))""",
    """            di, lx = np.array(di), np.array(lx)
            if len(di) < 40:
                continue
            xmin_l.append(float(np.exp(lx).min()))
            dmax_l.append(float(di.max()))
            tail_l.append(list(zip(di[-6:], np.exp(lx[-6:]))))
            cols.append(local_gamma(di, lx))""")
rep("""        cols = []
        for Z in pool:""",
    """        cols = []
        xmin_l, dmax_l, tail_l = [], [], []
        for Z in pool:""")
rep("""        GTAB[(lab, rk)] = med""",
    """        PD[(lab, rk)] = cm
        XMIN[(lab, rk)] = np.array(xmin_l)
        DMAX[(lab, rk)] = np.array(dmax_l)
        TAIL[(lab, rk)] = tail_l
        GTAB[(lab, rk)] = med""")
rep("""GTAB = {}   # (lab, rk) -> array over DGRID (median over designs)""",
    """GTAB = {}   # (lab, rk) -> array over DGRID (median over designs)
PD = {}     # (lab, rk) -> designs x DGRID 的逐设计 γ_loc
XMIN = {}   # (lab, rk) -> 每个设计进入 ln 的最小 x = D - JC - D_floor
DMAX = {}   # (lab, rk) -> 每个设计曲线可达的最大 ΔI
TAIL = {}   # (lab, rk) -> 每个设计末 6 个 (ΔI, x) 点""")

# ---------- 1)+2) 判决支 ----------
rep("""cnt = {'a': 0, 'b': 0, 'c': 0, 'skip': 0}""",
    """cnt = {'a': 0, 'b': 0, 'c': 0, 'd': 0, 'skip': 0}""")
rep("""    mono = bool(np.all(np.diff(v) <= 1e-9))
    if v[-1] <= tgt + 0.05 and v[-1] < v[0]:""",
    """    mono = bool(np.all(np.diff(v) <= 1e-9))
    if v[-1] <= 0.0:                      # (d) 末档非正 ⇒ 仪器下界，不判
        cnt['d'] += 1
        print('  %-7s r=%d 末档(ΔI=%.1f)比值=%+.3f ≤0 ⇒ (d) 仪器下界，**不参与判决**'
              % (lab, rk, d[-1], v[-1]))
        continue
    if tgt - 1e-12 <= v[-1] <= tgt + 0.05 and v[-1] < v[0]:""")
rep("""    cnt['a' if v_ == 'a' else ('b' if v_ == 'b' else 'c')] += 1""",
    """    cnt[v_] += 1""")
rep("""    print('  %-7s r=%d 全程单调降=%-5s 末档(ΔI=%.1f)比值=%.3f 对 1/%d 的余隙=%+.3f ⇒ (%s)'
          % (lab, rk, mono, d[-1], v[-1], rk, tail_gap, v_))""",
    """    print('  %-7s r=%d 全程单调降=%-5s 末档(ΔI=%.1f)比值=%.3f 对 1/%d 的余隙=%+.3f ⇒ (%s)'
          % (lab, rk, mono, d[-1], v[-1], rk, tail_gap, v_))
    if v_ == 'b':
        print('       （末档 %.3f 在带沿 %.3f 之上 ⇒ 按预注册原文不是 (a)）' % (v[-1], tgt + 0.05))""")
rep("""print('  合计：(a) 进入 1/r 邻域 %d 格 / (b) 降但未进入 %d 格 / (c) 无下行证据 %d 格 / 跳过 %d 格'
      % (na, nb, nc, cnt['skip']))""",
    """print('  合计：(a) 进入 1/r 邻域 %d 格 / (b) 降但未进入 %d 格 / (c) 无下行证据 %d 格 / '
      '(d) 末档非正⇒仪器 %d 格 / 跳过 %d 格'
      % (na, nb, nc, cnt['d'], cnt['skip']))
print('  [V2 判据不完备 ⇒ 公开记账] §86 的三选一未覆盖"末档落到仪器下界以下、变成非正数"这一支；'
      '\\n  (d) 是补记，不参与 (a)(b)(c) 的计数，也不改写任何已判格的结果。')""")

# ---------- [V3] 的取舍可见化 ----------
rep("""    d, e = DGRID[ok], ratio[ok] - tgt
    if len(e) < 4 or np.any(e <= 0):
        continue""",
    """    d, e = DGRID[ok], ratio[ok] - tgt
    ndrop = int(np.isfinite(ratio).sum() - ok.sum())
    if len(e) < 4 or np.any(e <= 0):
        continue""")
rep("""        print('  %-7s r=%d  ln(比值-1/r) vs ln ΔI 斜率=%+.3f  Pearson=%.4f  （ΔI≤%.1f 内有效）'
              % (lab, rk, sl, rr, d.max()))""",
    """        print('  %-7s r=%d  ln(比值-1/r) vs ln ΔI 斜率=%+.3f  Pearson=%.4f  （ΔI≤%.1f 内有效，'
              '拟合窗 %d 点，另有 %d 个有限点因 ≤1/r 被排除）'
              % (lab, rk, sl, rr, d.max(), len(e), ndrop))""")

# ---------- [V6] [V7] ----------
rep("""print('\\n== [V3] 收敛速率""",
    """print('\\n== [V6] 仪器诊断：逐设计 γ_loc / 最小 x / 曲线可达 ΔI ==')
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

print('\\n== [V7] 负值格的尾巴：逐设计末 6 个 (ΔI, x) 点（看 x 是否反增）==')
bad = [(lab, rk) for (lab, rk) in sorted(PD, key=lambda k: (k[1], k[0]))
       if np.isfinite(GTAB[(lab, rk)][list(DGRID).index(16.0)])
       and GTAB[(lab, rk)][list(DGRID).index(16.0)] <= 0]
for cell in bad:
    print('  [%s r=%d]' % cell)
    for i, tl in enumerate(TAIL[cell]):
        print('     d%d  %s' % (i, '  '.join('(%.2f, %.3e)' % (a, b) for a, b in tl)))
if not bad:
    print('  中位层面没有非正格 ⇒ 无需看尾巴')

print('\\n== [V3] 收敛速率""")

open(DST, 'w', encoding='utf-8', newline='\n').write(s)
b = open(DST, 'rb').read()
assert b.count(b'\r') == 0
print('WROTE %s (%d bytes, CR=%d)' % (DST, len(b), b.count(b'\r')))
