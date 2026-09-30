# -*- coding: utf-8 -*-
r"""实验 148：把 §96 要用到的**全部数字**先算出来落盘，不许口算。

对齐 1810.00298（IEEE JSTSP 12(5):841--856, 2018）全文实读的三处：
 [G1] Thm 5 / Cor 2 式 (44)(45)：标量量化（SDUSQ）的上界间隙 $=\frac r2\log_2(\pi e/6)+1$ bits/vector，
      文中把 $\frac12\log_2(\pi e/6)$ 印成 "0.254" ⇒ 核它到底等于多少，"0.254" 是舍入还是笔误。
 [G2] Remark 6 式 (49)(51)：向量格量化的间隙 $=\frac r2\log_2(2\pi e\,G_r)+1$，$G_p\to1/(2\pi e)$ 时**按维平均**归零；
      用文中 Example 3 的 $G_4=0.076603$（D4 格）算实际值，并与 [G1] 的标量系数比。
 [G3] Remark 1 式 (15)：$R^{na}_{GM}(D)\ge\sum_{|\mu_{A,i}|>1}\log|\mu_{A,i}|$（文中归给 ref [37]
      = Tatikonda & Mitter, TAC 49(7):1056--1068, 2004 ⇒ **不是**我们 bibitem 里那篇 TAC 49(9) 1549--1561）。
      在锚株上算这条下界，与我车道正文的 $R_{\rm exp}$ 逐位对比。
 [G4] 汇率换算：把 [G1] 的"加性比特间隙"翻译进我们自己的定价器坐标（$\gamma_r=2\ln2/r$，
      地板上方超量 $x\propto e^{-\gamma\Delta I}$）⇒ 同一间隙在代价货币里是地板上方超量的
      放大因子 $4^{0.254+1/r}$。只用代数，不碰任何设计。
"""
import sys
import numpy as np
from scipy.linalg import eigh
sys.stdout.reconfigure(encoding='utf-8')

src = open('p0/exp_c_audit.py', encoding='utf-8').read().split("print(r'== E53")[0]
_ns = {'__name__': 'p'}
exec(compile(src, 'p0/exp_c_audit.py[preamble]', 'exec'), _ns)
A = _ns['A']

ln2 = np.log(2.0)
print('== [G1] 标量量化间隙：文中印 "0.254" 的那个系数 ==')
c_scalar = 0.5 * np.log2(np.pi * np.e / 6.0)
print('  (1/2)·log2(πe/6) = %.6f  ⇒ 文中 "0.254" 是 %.4f 的截断（不是四舍五入到 3 位）' % (c_scalar, c_scalar))
for r in (1, 2, 3, 4):
    print('  r=%d  间隙 = %.4f·r + 1 = %.4f bits/vector' % (r, c_scalar, c_scalar * r + 1.0))

print('\n== [G2] 向量/格量化间隙：½log2(2πe·G_r) ==')
G4 = 0.076603
c_lat4 = 0.5 * np.log2(2.0 * np.pi * np.e * G4)
print('  D4 (文中 G_4=%.6f)：½log2(2πe·G_4) = %.6f bits/dim  vs 标量 %.6f ⇒ 省 %.1f%%'
      % (G4, c_lat4, c_scalar, 100 * (1 - c_lat4 / c_scalar)))
print('  极限：G_p → 1/(2πe) ⇒ ½log2(2πe·G_p) → %.6f（式 (51) 说的就是这个消失）'
      % (0.5 * np.log2(2.0 * np.pi * np.e * (1.0 / (2.0 * np.pi * np.e)))))

print('\n== [G3] Remark 1 式 (15) 的下界 vs 我车道的 R_exp（锚株） ==')
mod = np.abs(np.linalg.eigvals(A))
uns = mod[mod > 1.0]
bound_bits = float(np.log2(uns).sum())
print('  锚株 A 的 |λ|：%s' % np.array2string(np.sort(mod)[::-1], precision=6))
print('  Σ_{|λ|>1} log2|λ| = %.6f bits/vector   （我正文 R_exp = 1.168539 ⇒ 相对差 %.2e）'
      % (bound_bits, abs(bound_bits - 1.168539) / 1.168539))
ge = mod[mod >= 1.0]
n1 = int(np.sum(np.abs(mod - 1.0) < 1e-9))
print('  两种写法对比：|λ|>1 取 %d 个模式、|λ|≥1 取 %d 个，二者之差 = log2(Π_{≥1}/Π_{>1}) = %.3e；'
      '恰好 |λ|=1 的模式有 %d 个 ⇒ %s'
      % (len(uns), len(ge), float(np.log2(np.prod(ge) / np.prod(uns))), n1,
         '本株两种写法逐位相同' if n1 == 0 else '本株两种写法不同，必须点名用哪一种'))
print('  该文 ref [37] = Tatikonda & Mitter, TAC 49(7):1056--1068 (2004)；'
      '我们 bibitem tatikonda2004stochastic = Tatikonda–Sahai–Mitter, TAC 49(9):1549--1561 (2004) ⇒ 两篇不同')

print('\n== [G4] 同一间隙换算到代价货币（γ_r = 2ln2/r，x ∝ e^{−γΔI}） ==')
for r in (1, 2, 3, 4):
    dI = c_scalar * r + 1.0
    gam = 2.0 * ln2 / r
    print('  r=%d  ΔI=%.4f bits ⇒ 地板上方超量放大 e^{γΔI} = %.3f×   (= 4^{0.254+1/r})'
          % (r, dI, np.exp(gam * dI)))
print('  闭合检验：e^{(2ln2/r)(c r+1)} = 4^{c+1/r}，c=%.6f ⇒ r=1..4 的因子 %s'
      % (c_scalar, ['%.3f' % (4.0 ** (c_scalar + 1.0 / r)) for r in (1, 2, 3, 4)]))

print('\n== 附带：0.254 的身份核对（文中式 (38)：µ_{Σv,i}=Δ_i²/12 ⇒ 均匀量化器的二阶矩是 1/12） ==')
g_unif = 1.0 / 12.0
print('  ½log2(2πe·(1/12)) = %.6f  vs  ½log2(πe/6) = %.6f  ⇒ 差 %.3e（同一个量：0.254 就是均匀量化器的'
      '每维熵幂/空间填充损失）' % (0.5 * np.log2(2 * np.pi * np.e * g_unif), c_scalar,
                                  abs(0.5 * np.log2(2 * np.pi * np.e * g_unif) - c_scalar)))
print('  D4：½log2(2πe·G_4) = %.6f ⇒ 每维比标量少 %.6f bits（文中 Fig.7 的实测间隙随 r 变化正来自这一项）'
      % (c_lat4, c_scalar - c_lat4))
print('EXIT=0')
