# -*- coding: utf-8 -*-
# e183：Stavrou--Skoglund TAC-2022（read/ 全文，papers/txt/tac2022_stavrou_rwf.txt）的 Prop.1 结构条件
# 是否覆盖本车道的锚点植物。判决只许落在**可算残差**上：命题 1 的三种情形、正规性、A^2 的实性。
# 纪律：#46 每个判据打印数值与残差；#53 不许用上界去否证上界律；这里不是判决科学结论，是**适用域**。
import io
import sys
import numpy as np

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')
OUT = []


def p(*a):
    s = ' '.join(str(x) for x in a)
    OUT.append(s)
    print(s, flush=True)


P = {'__name__': 'p'}
exec(compile(open('p0/e162_ray_closure.py', encoding='utf-8').read().split('CRAB = {')[0],
             'prefix', 'exec'), P)
AUD = P['_ns']
A, W = AUD['A'].copy(), AUD['W'].copy()
TH = AUD['TH']
n = A.shape[1]
sym = lambda X: 0.5 * (X + X.T)
nrm = lambda X: float(np.max(np.abs(X)))

p('== [V0] 锚点数据（从 p0/e162_ray_closure.py 的 _ns 取，不另起数）==')
p('  A 维度 %d；eig(A)=%s' % (n, np.array2string(np.linalg.eigvals(A), precision=4)))
p('  W=%s' % np.array2string(W, precision=4))
p('  eig(W)=%s' % np.array2string(np.linalg.eigvalsh(sym(W)), precision=4))
p('  eig(Theta)=%s' % np.array2string(np.linalg.eigvalsh(sym(TH))[::-1], precision=4))

p('')
p('== [V1] TAC-2022 Proposition 1 三情形逐一带残差（他们的式 (14) 只在情形成立时才是那个凸规划）==')
alpha = np.trace(A) / n
r1 = nrm(A - alpha * np.eye(n))
p('  情形 1  $A=\\alpha I_p$：$\\alpha=%.4f$，残差 $\\max|A-\\alpha I|=%.4e$ ⇒ %s'
  % (alpha, r1, '成立' if r1 < 1e-9 else '不成立'))
r2a = nrm(A - A.T)
sw = np.trace(W) / n
r2b = nrm(sym(W) - sw * np.eye(n))
p('  情形 2  $A=A^{\\top}$ 且 $W=\\sigma_w^2 I_p$：$\\max|A-A^{\\top}|=%.4e$，'
  '$\\sigma_w^2=%.4f$，$\\max|W-\\sigma_w^2 I|=%.4e$ ⇒ %s'
  % (r2a, sw, r2b, '成立' if (r2a < 1e-9 and r2b < 1e-9) else '不成立'))
r3 = nrm(A - W)
p('  情形 3  $A=W\\succeq0$：$\\max|A-W|=%.4e$ ⇒ %s'
  % (r3, '成立' if r3 < 1e-9 else '不成立'))
cov = r1 < 1e-9 or (r2a < 1e-9 and r2b < 1e-9) or r3 < 1e-9
p('  三情形合计命中 **%d／3** ⇒ 锚点植物 %s他们的结构类'
  % (int(r1 < 1e-9) + int(r2a < 1e-9 and r2b < 1e-9) + int(r3 < 1e-9),
     '落在' if cov else '**不在**'))

p('')
p('== [V2] 他们论证真正依赖的两件事：成对交换 与 $A^2$ 的实特征值 ==')
p('  $\\max|AW-WA|=%.4e$（情形 1/3 的必要前提；非零即不能同时正交对角化）'
  % nrm(A @ W - W @ A))
p('  $\\max|AA^{\\top}-A^{\\top}A|=%.4e$ ⇒ $A$ %s（正规才谈得上式 (19) 的谱分解口径）'
  % (nrm(A @ A.T - A.T @ A), '正规' if nrm(A @ A.T - A.T @ A) < 1e-9 else '**非正规**'))
eigA2 = np.linalg.eigvals(A @ A)
p('  eig$(A^2)=%s$；虚部最大 $%.3e$ ⇒ 式 (14) 的 $\\mu_{A^2,i}$ 在本锚点 %s实数'
  % (np.array2string(eigA2, precision=4), float(np.max(np.abs(eigA2.imag))),
     '全是' if float(np.max(np.abs(eigA2.imag))) < 1e-9 else '**不是**'))
p('  他们式 (14) 要求 $\\mu_{\\Lambda,i}=\\mu_{A^2,i}\\mu_{\\Delta,i}+\\mu_{\\Sigma_W,i}$ '
   '是一组**实**数并对 $i$ 排序求和；上面两条残差 ⇒ 该式在本锚点上不是良定义的实规划。')

p('')
p('== [V3] 与本车道式 (8)（eq:wf）的结构对照（只列形式，不判优劣）==')
p('  TAC-2022 式 (14)：目标 $\\tfrac12\\sum_i\\log(\\mu_{\\Lambda,i}/\\mu_{\\Delta,i})$，'
  '约束 $\\sum_i\\mu_{\\Delta,i}\\le D$，而 $\\mu_{\\Lambda}$ **依赖设计** $\\Delta$（自洽）。')
p('  本车道 eq:wf：目标 $\\sum_i\\min(\\gamma_i,\\nu)$ 即 $\\operatorname{tr}(\\Gamma R)$ **线性**，'
  '约束 $\\det R\\ge2^{-2I}$（率），$\\gamma_i$ 是**冻结先验 $X$ 下**对称阵 $\\Gamma$ 的特征值，不随 $R$ 变。')
p('  两者不是同一个规划：变量口径（失真 vs 收缩比）、目标（对数比 vs 线性）、约束（trace vs log-det）三处都不同；'
  '相同点只有"阈值形如 $\\min(\\cdot,\\nu)$ 的反向水填"和"用 Hadamard/AM--GM 把行列式换成逐维"这一步。')
p('  对角化依据也不同：他们靠 Prop.1 的成对交换（附录 A 定理 3）；'
  '本车道靠 $\\Gamma$ 本身对称（$\\Gamma=G^{-1/2}(Z^{\\top}X\\Theta XZ)G^{-1/2}\\succeq0$），'
  '与 $A$ 是否对称/正规无关。')

p('')
p('== [V4] 条数与退出（#46）==')
p('  结构情形命中 %d/3；交换残差 %.3e；正规残差 %.3e；eig($A^2$) 虚部 %.3e'
  % (int(r1 < 1e-9) + int(r2a < 1e-9 and r2b < 1e-9) + int(r3 < 1e-9),
     nrm(A @ W - W @ A), nrm(A @ A.T - A.T @ A),
     float(np.max(np.abs(eigA2.imag)))))
p('  ⇒ 可写进正文的**最窄**一句：本锚点不满足他们的结构条件，其定理不为本对象定价；'
  '但"反向水填"这一步的形状与 AM--GM/Hadamard 这一步**不是我们的发现**（已在 §VII/§X 口径里降级）。')
open('p0/e183_out.txt', 'w', encoding='utf-8').write('\r\n'.join(OUT) + '\r\n')
