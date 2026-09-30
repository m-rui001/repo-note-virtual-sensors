# -*- coding: utf-8 -*-
r"""e162c：用 **e152 自己的 rand 株定义**重测秩-1（e162b 的 [B] 段作废的根因）。

作废说明（本车道自己的账，#43 家族）：e162b 里我把 rand 株写成了 e159 的
$A=Q\\,\\mathrm{diag}(1.6,1.25,0.9,0.7)\\,Q^{\\mathsf T}$ 构造，而 e152 用的是
$A\\sim N(0,1.05)$、$B\\sim N(0,1.15)$、$W=(MM^{\\mathsf T}+0.6I)$（seed=1）。
$\\Rightarrow$ e162b 的 rand-1 行（"$d I:2.40e{-}2\\to9.42e{-}1$、$19.15°$"、"秩-1 地板 39.41"）**是另一株的数**，不得引用。
本脚本只换株的定义，机器（割线命中 $J=D$、$|D-D_{\\rm tgt}|\\le10^{-7}$、Powell 粗 300/细 1200）不动。

判据（沿用 e162/e162b，跑前写死）：
 [R0] 逐格打印有效读数条数与 $D$ 残差；无读数不判（#46）。锚点两格已由 e162 六位锁死，本轮只补 rand 株。
 [R1] 若修正后 $I_{\rm free}$ **高于** e152 同格读数 $\\Rightarrow$ 不可复现（同可行集、更强优化器不可能更差）$\\Rightarrow$ 判实现错，不判物理。
 [R2] 修正后 $d I\\ge10^{-2}$ bit $\\Rightarrow$ 该格进入 C 的"非对角"判读带；仍 $<10^{-2}$ $\\Rightarrow$ 原结论不变、只更数值。
 [R3] 顺带打印秩-1 的可达代价下界（$s=10^7$、200 方向）$\\Rightarrow$ "不可达"必须有地板佐证（#43），
      并检验是否等于 e153 口径的带权地板 $J_C+\\min_{\\mathrm{rank}Z=1}\\mathrm{tr}(\\Theta P_p)$。
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
exec(compile(open('p0/e162_ray_closure.py', encoding='utf-8').read().split("p('== 锚点")[0],
             'p0/e162_ray_closure.py[prefix]', 'exec'), P)
AUD = P['_ns']


def sym(X):
    return 0.5 * (X + X.T)


def e152_rand(seed):
    """逐字照 e152.make_plant('rand4', seed)。"""
    rgn = np.random.default_rng(seed)
    m = 4
    A = rgn.normal(0.0, 1.05, (m, m))
    B = rgn.normal(0.0, 1.15, (m, m))
    M = rgn.normal(0.0, 1.0, (m, m))
    return A, B, sym(M @ M.T + 0.6 * np.eye(m))


E152 = {('rand-1', 32.31): (1.743144, 1.95e-2, 3.59), ('rand-1', 34.41): (1.599991, 8.71e-3, 4.11),
        ('rand-1', 56.66): (1.097247, 8.93e-3, 2.68), ('rand-1', 80.00): (0.961934, 2.40e-2, 6.28)}

p('== [R3] 秩-1 可达代价下界（e152 的 rand-1 株，$s=10^7$，200 方向）==')
rng = np.random.default_rng(9090)
for nm, get in [('anchor', lambda: (AUD['A'].copy(), AUD['B'].copy(), AUD['W'].copy())),
                ('rand-1', lambda: e152_rand(1))]:
    A, B, W = get()
    Pc, K, TH, JC = P['ctrl_full'](A, B, W, np.eye(4), np.eye(4))
    wth, U = np.linalg.eigh(TH)
    Uth = U[:, np.argsort(-wth)]
    P['A'], P['B'], P['W'], P['TH'], P['JC'], P['Uth'] = A, B, W, TH, JC, Uth
    fl = np.inf
    dirs = [Uth[:, [i]] for i in range(4)] + [np.linalg.qr(rng.normal(0, 1, (4, 1)))[0] for _ in range(200)]
    for Zr in dirs:
        Zn = Zr / np.linalg.norm(Zr)
        _, D = P['make_cur'](Zn, 1)(7.0)
        if np.isfinite(D) and D < fl:
            fl = D
    p('  %-8s $J_C=%.4f$  秩-1 地板 $\\min D=%.4f$（高出 $J_C$ %.4f）' % (nm, JC, fl, fl - JC))
    p('')
    p('  --- %s：强秩-1 重测 ---' % nm)
    p('  %7s %10s %10s %8s %8s %8s %6s' % ('D', 'e152 dI', '本次 dI', 'e152角', '新角', '残差', '有效'))
    for t in (32.31, 34.41, 56.66, 80.00):
        P['GUESS'][0] = 0.0
        Ia, Da, _ = P['read_pt'](Uth[:, [0]], t)
        st = [Uth[:, [i]] for i in range(4)] + [np.linalg.qr(rng.normal(0, 1, (4, 1)))[0] for _ in range(16)]
        res = P['best'](t, st, coarse=300, fine=1200, ntop=3)
        if not res:
            p('  %7.2f  有效 0 条（地板 %.2f）=> 不判' % (t, fl))
            continue
        If, rb, G = res[0]
        u = G[:, 0] / np.linalg.norm(G[:, 0])
        ang = np.degrees(np.arccos(min(1.0, abs(float(u @ Uth[:, 0])))))
        old = E152.get((nm, t))
        if not np.isfinite(Ia):
            p('  %7.2f  Theta 首轴不可达（残差门内无读数）$\\Rightarrow$ 该格 $d I$ 无定义，只报自由端 $I=%.6f$' % (t, If))
            continue
        dI = Ia - If
        flag = ' [R1!] 新 $I_{\\rm free}$ 高于 e152 $\\Rightarrow$ 不可复现' if (old and If > old[0] + 1e-5) else ''
        p('  %7.2f %10s %10.2e %8s %8.2f %8.1e %6d%s'
          % (t, ('%.2e' % old[1]) if old else 'n/a', dI, ('%.2f' % old[2]) if old else 'n/a', ang, rb, len(res), flag))
        if old:
            p('          [R2] $d I$ %.2e $\\to$ %.2e（%.1f$\\times$）$\\Rightarrow$ %s'
              % (old[1], dI, dI / old[1], '跨过 $10^{-2}$ $\\Rightarrow$ 该格改判' if dI >= 1e-2 else '仍 $<10^{-2}$ $\\Rightarrow$ 原结论不变'))
p('')
p('用时 %.0f s | DARE %d | 有效读数 %d | 回退 %d | 异常 %s' % (time.time() - T0, P['NDB'][0], P['NVALID'][0], P['NFALLBACK'][0], P['EXC'] or '无'))
open('p0/e162c_out.txt', 'w', encoding='utf-8').write('\n'.join(out) + '\n')
