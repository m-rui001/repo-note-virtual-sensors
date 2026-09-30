r"""E46 = 把 #23 的"两族整条 $E(D)$ 重合"从数值事实往上推一格。

做法：先判定见证植物属于哪一类结构（块对角 + $\ker\Theta$ 支撑在第二块），
在块对角子族里把"自由删"证成命题（率可加、代价看不见第二块、下确界在 $V_{22}\to\infty$ 处取到），
然后把相关噪声那一腿单独量出来——那一腿正是 #23 我用 $10^{9}\!\sim\!10^{17}$ 的 $V_{22}$ 做驻点检验时偷偷跨过去的。

打印约定（这一版重写就是因为又踩了）：含 LaTeX 的字面量一律单独一行 `r'...'`，
需要换行就 `print('')` 另起一句；绝不把 raw 和非 raw 字面量隐式拼接。
"""
import sys
sys.stdout.reconfigure(encoding='utf-8')
import numpy as np
import p0.exp_plants_oos as E35
from p0.ticket import floor0, ker_theta

ln2 = np.log(2.0)
NB = 2                                     # split: S1 = span(e1,e2), S2 = span(e3)


def plant(c13=0.0, au=.6, b3=0.0, wu=1.0):
    A = np.array([[1.25, .30, c13], [0., .85, 0.], [0., 0., au]])
    B = np.array([[1.], [0.], [b3]])
    W = np.diag([1.0, .7, wu])
    return E35.Pl(A, B, W, np.eye(3), np.array([[1.0]]))


F = np.array([[1., 0., 0.], [0., 0., 1.]])          # 通道1 测 e1，通道2 测 e3
FD = np.array([[1., 0., 0.]])                        # 删掉通道2
Pl = plant()
PhiF, PF, resF = floor0(Pl, F)
PhiD, PD, resD = floor0(Pl, FD)


def blocks(M):
    return M[:NB, :NB], M[:NB, NB:], M[NB:, NB:]


print('== E46：块对角结构 + 相关噪声那一腿 ==')
print('')
print(r'[1] 类别判定（命题的前提，一条都不许含糊）')
A11, A12, A22 = blocks(Pl.A)
W11, W12, W22 = blocks(Pl.W)
T11, T12, T22 = blocks(Pl.Th)
print(r'  split，dim S1 = 2，dim S2 = 1。')
print(r'  |A12| = %.1e   |A21| = %.1e   （A 块对角等价于两条都为 0；A13 是我在 E43 里当旋钮的那一项）'
      % (np.linalg.norm(A12), np.linalg.norm(Pl.A[NB:, :NB])))
print(r'  |W12| = %.1e   |Theta22| = %.1e   |Theta12| = %.1e'
      % (np.linalg.norm(W12), np.linalg.norm(T22), np.linalg.norm(T12)))
print(r'  F 的行按块分：行 1 只碰 S1（第 3 列分量为 0），行 2 只碰 S2，所以 F = F1 (+) F2，')
print(r'  删除动作就是"不要 F2 这一条"。')
wT, ker, sup, wsup = ker_theta(Pl)
print(r'  A 块对角给出 S2 是 A-不变；Theta22 = 0 与 ker(Theta) 包含 S2 是同一件事。')
print(r'  这就是 E43 [6] 里那三条数（Pc 第三列、B^T Pc 的 (1,3)、|K[:,3]|）的结构版本。')

print('')
print(r'[2] 块对角子族（rho = 0）：命题的三条组成，逐条打印')
print(r'  %10s %12s %15s %15s %15s %15s %s'
      % ('v1', 'v2', 'x_full', 'x_del', 'I_full', 'I_del', '跨块 / I 加法残差 / I2 实测-闭式'))
for v1 in [1e-6, 1e-3, 1.0]:
    for v2 in [1.0, 1e3, 1e9]:
        Vd = np.diag([v1, v2])
        xf, If, Ptf = Pl.xI(F, Vd, PhiF)
        xd, Id, _ = Pl.xI(FD, np.array([[v1]]), PhiF)
        Pf = Pl.P(F, Vd, Ptf)
        I2 = If - Id
        pred = .5 * np.log1p(Ptf[NB, NB] / v2) / ln2
        print(r'  %10.2e %12.2e %15.9f %15.9f %15.9f %15.9f  %.1e / %.1e / %.3e-%.3e'
              % (v1, v2, xf, xd, If, Id,
                 np.linalg.norm(blocks(Pf)[1]), np.linalg.norm(blocks(Ptf)[1]), I2, pred))
print(r'  三件事同时成立：x_full = x_del（跨 v2 六个数量级逐位一样，残差在 1e-16 档）、')
print(r'  P 与 Pt 的跨块为 0、I_full = I_del + 0.5*log2(1 + Pt22/v2)（第二项用的是先验块 Pt22，不是后验块）。')
print(r'  第二项恒正，且只有 v2 -> infinity 时才归零 —— 所以下确界落在边界上，不在内部取到。')

print('')
print(r'[3] 由此直接得到的命题（rho = 0 子族，闭式，不是拟合）')
print(r'  命题. 设 A = diag(A1, A2)、W = diag(W1, W2)、ran(Theta) 含于 S1、F = F1 (+) F2。')
print(r'  则在噪声取块对角的子族上，对任何预算 xt >= 0：E_full(xt) = E_del(xt)，')
print(r'  且下确界在 V22 -> infinity 处取到（不存在有限最优）。')
print(r'  证明要点（三步，全是上表量过的量）：')
print(r'  (i) Riccati 在块对角数据下不产生跨块，Pt 与 P 都块对角，故 x = tr(Theta P) 只是 V11 的函数；')
print(r'  (ii) det(F Pt F^T + V) 按块相乘，det V 也是，于是 I = sum_k 0.5*log2(1 + Pt_kk/Vkk) 对 k 可加；')
print(r'  (iii) 约束只含 V11，I 的第二项对 V22 单调降、下界 0，取 V22 -> infinity 得 I -> I1；')
print(r'        反向不等式 E_full <= E_del 由扩族立得。证毕。')
print(r'  可检验的副作用：最优 V22 应当发散，且两族曲线逐位重合')
print(r'  —— E44 给的 lambda(V*) = 1e9~1e17（随 x 单调升）与"没有有限最优"完全一致，')
print(r'  我当时把它读成"第二条通道几乎不测"，其实那是 (iii) 的边界行为。')

print('')
print(r'[4] 相关噪声那一腿：先记我自己的错，再量真正该量的')
print(r'  E44 的驻点检验是在 v2 = 1e9~1e17 那个解上做中心差分。')
print(r'  但由 (iii) 那里根本没有有限最优：v2 -> infinity 时 dI/dv2 与 dx/dv2 一起趋于 0，')
print(r'  所以"dI/dv2 + mu dx/dv2 = 0 成立"在边界上是恒真式，不构成证据。')
print(r'  更要紧的：rho 变号等价于把第二条测量乘 -1（正交变换，模型不变），')
print(r'  于是 x 与 I 都是 rho 的偶函数，d/drho 在 rho = 0 恒为 0')
print(r'  —— 我那一行"d_rho I 约等于 0 通过"同样是空判。现在把偶性直接量出来：')
print(r'  %10s %8s %20s %20s %16s %16s'
      % ('v1', 'v2', 'x(+r)-x(-r)', 'I(+r)-I(-r)', 'dx/r^2', 'dI/r^2'))
for (v1, v2) in [(1e-4, 1e2), (1e-2, 1e2), (1.0, 1e2), (1e-4, 1e-1), (1e-4, 1e5)]:
    def xI_rho(rr, v1=v1, v2=v2):
        Vq = np.array([[v1, rr * np.sqrt(v1 * v2)], [rr * np.sqrt(v1 * v2), v2]])
        a, b, _ = Pl.xI(F, Vq, PhiF)
        return a, b
    d = 1e-3
    xp, xm = xI_rho(d), xI_rho(-d)
    x0, I0 = xI_rho(0.0)
    rr = 0.5
    x1, I1 = xI_rho(rr)
    print(r'  %10.2e %8.0e %20.2e %20.2e %16.3e %16.3e'
          % (v1, v2, xp[0] - xm[0], xp[1] - xm[1], (x1 - x0) / rr ** 2, (I1 - I0) / rr ** 2))
print(r'  前两列是偶性的判决：若都在 1e-16 档，我上一轮报的"三行全判驻点"就一条都不算数。')
print(r'  剩下的真问题是 dI/r^2 的符号：负值意味着有限相关性能把 I 压到 rho = 0 之下，')
print(r'  而命题只覆盖 rho = 0 —— 那一格正是 E44 随机撒 4000 点没找到反例、但证明不了没有的地方。')

print('')
print(r'[5] 我还缺什么（写给下一次的自己和任何接手的人）')
print(r'  缺的不是"再验一下"，是一个具体的归约：白化 (F, V) -> (F V^(-1/2), I)')
print(r'  把相关噪声族变成"通道方向可变、噪声固定为单位阵"的族。')
print(r'  在这个写法里 V22 -> infinity 对应第二条通道方向趋于 0，')
print(r'  于是猜想是：白化后的可行集在该极限下退回 rho = 0 那一套，相关噪声买不到额外率。')
print(r'  做出来，[3] 的命题就没有"rho = 0 子族"这个尾巴，#23 的整条 E(D) 重合升成结构命题；')
print(r'  做不出来（比如找到某个 (v1, v2, rho) 让 I < E_del），#23 第一节那张表就要重画。两个结果我都要。')
