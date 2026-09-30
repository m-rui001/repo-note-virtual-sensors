# -*- coding: utf-8 -*-
r"""e162b（合并版）：两件都在判"C 的秩报告"，用的是**证书逻辑**而不是求解器互信。

证书逻辑（本脚本的全部要点，跑前写死）：
  我的可行集（秩 $\le r$ 的设计）$\\subseteq$ C 的可行集（全体 $S\\succeq0$），而 C 的 $I(S^*)$ 来自凸 SDP（CLARABEL）$\\Rightarrow$
  **若我的秩-$r$ 搜索达到 C 的 SDP 值（差 $\\le\\epsilon$），就证明全局最优被某个秩 $\\le r$ 的点取到** $-$ $-$
  这是"存在秩-$r$ 全局最优解"的**证书**，不依赖 C 的数值秩诊断；反之若我的秩-$r$ 严格高于 C 的值，
  只能说明我的搜索欠优化（e162 已实测到：秩$\\le$3 的 min 竟比秩$\\le$2 差 $6.9\\times10^{-3}$ $\\Rightarrow$ 包含关系违反 $\\Rightarrow$ 我方参照作废）。

三块：
 [P] 秩-1/秩-2 的**可达代价下界**（$s=10^7$）$\\Rightarrow$ 每格的可达地板（板账 #43 的口径，不许把不可达当结论）。
 [A2] $D{=}34.41$：C 报 $k{=}3$、$I(S^*)=4.893329$。我方秩-2 已读到 $4.893399$（差 $+7.0\\times10^{-5}$）$\\Rightarrow$  Intensify 秩-2：
      若 $\le 4.893329+10^{-5}$ $\\Rightarrow$ **存在秩-2 全局最优** $\\Rightarrow$ C 的 $k{=}3$ 那一列要么报错、要么 $\\mathrm{Opt}(D)$ 含秩-2 成员
      （正是 C 补记十二 73-B 说的"面非平凡"的具体实例）；若停在 $>10^{-4}$ 之上 $\\Rightarrow$ 只说明我仍欠优化，不判。
 [B] 强秩-1 重测 e152 全格（锚点两格已由 e162 修正为 $7.96°/6.02°$）：rand-1 四格 $d I$ 是否跨过 $10^{-2}$ bit。

预注册判据：
 [Q0] 只用 $|D-D_{\\rm tgt}|\\le10^{-7}$ 的点；逐行打印有效读数条数与起点一致性（#46）。
 [Q1] [A2] 的结论只允许三选一：证书成立（差 $\le10^{-5}$）/ 欠优化（$10^{-5}\sim10^{-3}$，不判）/ 不可达（无有效读数）。
 [Q2] [B]：修正后 $d I\ge10^{-2}$ $\\Rightarrow$ 该格改判；$<10^{-2}$ $\\Rightarrow$ 原结论不变、只更数值。
      若某格修正后的 $I_{\rm free}$ **高于** e152 的读数 $\\Rightarrow$ 不可复现（同可行集、更强优化器不可能更差）$\\Rightarrow$ 按 e162 [K4] 处理，先查实现。
"""
import sys
import time
import numpy as np
sys.stdout.reconfigure(encoding='utf-8')

out = []
T0 = time.time()


def p(*a):
    s = ' '.join(str(x) for x in a)
    out.append(s)
    print(s)
    sys.stdout.flush()


P = {'__name__': 'p'}
exec(compile(open('p0/e162_ray_closure.py', encoding='utf-8').read().split('CRAB = {')[0],
             'p0/e162_ray_closure.py[prefix]', 'exec'), P)
AUD = P['_ns']                                    # exp_c_audit 前言：A, B, W 在这里


def set_plant(A, B, W):
    Pc, K, TH, JC = P['ctrl_full'](A, B, W, np.eye(4), np.eye(4))
    P['A'], P['B'], P['W'], P['TH'], P['JC'] = A, B, W, TH, JC
    wth, Uth = np.linalg.eigh(TH)
    P['Uth'] = Uth[:, np.argsort(-wth)]
    return TH, JC, P['Uth']


def rand_plant(seed):
    rgn = np.random.default_rng(seed)
    A = np.linalg.qr(rgn.normal(0, 1, (4, 4)))[0] @ np.diag([1.6, 1.25, 0.9, 0.7]) @ np.linalg.qr(rgn.normal(0, 1, (4, 4)))[0].T
    B = rgn.normal(0, 1, (4, 4))
    M = rgn.normal(0.0, 1.0, (4, 4))
    return A, B, P['sym'](M @ M.T + 0.6 * np.eye(4))


def floor_D(Uth, r, k=120, hi=7.0):
    """秩-$r$ 的可达代价下界：$s=10^{hi}$ 下 $\\min$ 过 $k$ 个标架的 $D$。"""
    best_v = np.inf
    frames = [Uth[:, :r]] + [np.linalg.qr(np.random.default_rng(1000 + q).normal(0, 1, (4, r)))[0] for q in range(k)]
    for Z in frames:
        Zn = Z / np.linalg.norm(Z)
        _, D = P['make_cur'](Zn, Z.shape[1])(hi)
        if np.isfinite(D) and D < best_v:
            best_v = D
    return best_v


SDP = {('anchor', 34.41): 4.893329, ('anchor', 56.66): 2.053565, ('anchor', 80.00): 1.603105}
E152 = {('anchor', 56.66): (2.063303, 8.51e-2, 5.13), ('anchor', 80.00): (1.608043, 4.50e-2, 4.17),
        ('rand-1', 32.31): (1.743144, 1.95e-2, 3.59), ('rand-1', 34.41): (1.599991, 8.71e-3, 4.11),
        ('rand-1', 56.66): (1.097247, 8.93e-3, 2.68), ('rand-1', 80.00): (0.961934, 2.40e-2, 6.28)}

p('== [P] 锚点各秩的可达代价下界（$s=10^7$，120 标架 + $\\Theta$ 前 $r$ 列）==')
A, B, W = AUD['A'].copy(), AUD['B'].copy(), AUD['W'].copy()
TH, JC, Uth = set_plant(A, B, W)
FL = {}
for r in (1, 2, 3):
    FL[r] = floor_D(Uth, r)
    p('  秩-%d 地板 $\\approx$ %.4f（高出 $J_C$ %.4f）%s' % (r, FL[r], FL[r] - JC,
      ' $\\Rightarrow$ $D\\in\\{32.31,34.41\\}$ 对秩-1 不可达 $\\Rightarrow$ 空判据' if r == 1 and FL[r] > 34.41 else ''))

p('')
p('== [A2] $D{=}34.41$ 的秩-2 证书：C 的 SDP $I(S^*)=%.6f$（C 报 $k{=}3$）==' % SDP[('anchor', 34.41)])
t = 34.41
rng = np.random.default_rng(31337)
starts2 = [Uth[:, :2]] + [np.linalg.qr(rng.normal(0, 1, (4, 2)))[0] for _ in range(24)]
res2 = P['best'](t, starts2, coarse=700, fine=3000, ntop=5)
if not res2:
    p('  有效 0 条 $\\Rightarrow$ [Q1] 不可达/无读数，不判')
else:
    I2, rb, G2 = res2[0]
    gap = I2 - SDP[('anchor', t)]
    pv = np.degrees(np.arccos(np.clip(np.linalg.svd(Uth[:, :2].T @ G2, compute_uv=False), -1, 1)))
    p('  秩-2 $\\min I=%.6f$（有效 %d 条，最差残差 %.1e）vs SDP %.6f $\\Rightarrow$ 差 %+.2e bit' % (I2, len(res2), rb, SDP[('anchor', t)], gap))
    p('  [Q1] 结论：%s' % ('**证书成立 $\\Rightarrow$ 存在秩-2 全局最优 $\\Rightarrow$ C 的 $k{=}3$ 列须重判（或 $\\mathrm{Opt}(D)$ 含秩-2 成员）**'
                            if gap <= 1e-5 else ('仍欠优化（$10^{-5}<\\cdot\\le10^{-3}$）$\\Rightarrow$ 不判' if gap <= 1e-3 else '差距 $>10^{-3}$ $\\Rightarrow$ 我方秩-2 搜索不可用，只报数')))
    p('  最优 2-平面与 $\\Theta$ 前 2 平面的主角度 %s deg' % np.array2string(pv, precision=2))

p('')
p('== [B] 强秩-1 重测 e152 全格 ==')
p('%-8s %7s %10s %10s %7s %8s %8s %6s' % ('株', 'D', 'e152 dI', '本次 dI', '比值', 'e152角', '新角', '有效'))
RES = []
for nm in ('anchor', 'rand-1'):
    if nm == 'anchor':
        A2, B2, W2 = AUD['A'].copy(), AUD['B'].copy(), AUD['W'].copy()
    else:
        A2, B2, W2 = rand_plant(1)
    TH, JC, Uth = set_plant(A2, B2, W2)
    for t in (32.31, 34.41, 56.66, 80.00):
        Ia, Da, _ = P['read_pt'](Uth[:, [0]], t)
        if not np.isfinite(Ia):
            flr = floor_D(Uth, 1)
            p('%-8s %7.2f  Theta 首轴不可达（本株秩-1 地板 %.2f > 该 D）=> 空判据' % (nm, t, flr))
            continue
        st = [Uth[:, [i]] for i in range(4)] + [np.linalg.qr(rng.normal(0, 1, (4, 1)))[0] for _ in range(16)]
        res = P['best'](t, st, coarse=300, fine=1200, ntop=3)
        if not res:
            p('%-8s %7.2f  有效 0 条 => 不判' % (nm, t))
            continue
        If, rb, G = res[0]
        u = G[:, 0] / np.linalg.norm(G[:, 0])
        ang = np.degrees(np.arccos(min(1.0, abs(float(u @ Uth[:, 0])))))
        dI = Ia - If
        old = E152.get((nm, t))
        flag = ''
        if old and If > old[0] + 1e-5:
            flag = ' [K4!] 新搜索比 e152 差 => 不可复现'
        p('%-8s %7.2f %10.2e %10.2e %7s %8s %8.2f %6d%s'
          % (nm, t, (old[1] if old else np.nan), dI, ('%.2f' % (dI / old[1])) if old else 'n/a',
             ('%.2f' % old[2]) if old else 'n/a', ang, len(res), flag))
        RES.append((nm, t, dI, ang, old))

p('')
p('== [Q2] 判决 ==')
for nm, t, dI, ang, old in RES:
    tag = '跨过 $10^{-2}$ $\\Rightarrow$ 该格改判入"非对角"' if dI >= 1e-2 else '仍 $<10^{-2}$ $\\Rightarrow$ 原结论不变，只更数值'
    p('  %-8s D=%6.2f：dI %s -> %.2e（%s）；离开角 %s -> %.2f°' % (nm, t, ('%.2e' % old[1]) if old else 'n/a', dI, tag, ('%.2f' % old[2]) if old else 'n/a', ang))
p('  锚点两格已由 e162 独立修正（$8.51e{-}2\\to9.49e{-}2$、$4.50e{-}2\\to5.00e{-}2$ bit），本次只补 rand 株。')
p('')
p('用时 %.0f s | DARE %d 次 | 有效读数 %d | 回退 %d | 异常 %s' % (time.time() - T0, P['NDB'][0], P['NVALID'][0], P['NFALLBACK'][0], P['EXC'] or '无'))
open('p0/e162b_out.txt', 'w', encoding='utf-8').write('\n'.join(out) + '\n')
