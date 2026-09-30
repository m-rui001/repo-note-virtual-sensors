r"""E47 = 相关噪声那一腿的结构化扫描，替掉 E44 里"随机撒 4000 点"那一步。

问题：#23 说放大族（$r=2$，允许相关噪声）的最优率跟删掉冗余通道那一族（$r=1$）逐位重合。
$\rho=0$ 那一半 E46 已经证成命题；剩下的问号是 $\rho\ne0$ 能不能买到更低的车票。
随机撒点回答不了这个——它连我自己的优化器都打不过。这里改成精确配预算的一维扫描：
固定 $v_2$，对每个 $\rho$ 用二分把 $v_1$ 调到预算严格成立，于是比较的是同一预算下的率，
不再有"可行性容差"和"随机覆盖不到"两个借口。唯一的数值操作是对数尺度二分，收敛判据写在表里。
"""
import sys
sys.stdout.reconfigure(encoding='utf-8')
import numpy as np
import p0.exp_plants_oos as E35
from p0.ticket import floor0

NB = 2


def plant(c13=0.0, au=.6, b3=0.0, wu=1.0):
    A = np.array([[1.25, .30, c13], [0., .85, 0.], [0., 0., au]])
    B = np.array([[1.], [0.], [b3]])
    W = np.diag([1.0, .7, wu])
    return E35.Pl(A, B, W, np.eye(3), np.array([[1.0]]))


F = np.array([[1., 0., 0.], [0., 0., 1.]])
FD = np.array([[1., 0., 0.]])
Pl = plant()
PhiF, PF, resF = floor0(Pl, F)
RHO = np.arange(-49, 50) * 0.02          # 含 0.0，99 点


def solve_xt(Fx, xt, make, lo=-16.0, hi=14.0, it=60):
    r"""在 `make(v1)` 给出 $(V,x,I)$ 的族里二分 $v_1=\exp(\cdot)$ 使 $x=x_t$。返回 $(I,\text{残差},v_1)$。"""
    def f(lg):
        v1 = 10.0 ** lg
        return make(v1)
    a0, _, _ = f(lo)
    a1, _, _ = f(hi)
    if not (a0 <= 0 <= a1):
        return None
    for _ in range(it):
        mid = .5 * (lo + hi)
        if f(mid)[0] < 0:
            lo = mid
        else:
            hi = mid
    d, I, v1 = f(.5 * (lo + hi))
    return I, d, v1


def mk_full(rho, v2):
    def make(v1):
        sd = np.sqrt(v1 * v2)
        Vq = np.array([[v1, rho * sd], [rho * sd, v2]])
        x, I, _ = Pl.xI(F, Vq, PhiF)
        return x - XT[0], I, v1
    return make


def mk_del():
    def make(v1):
        x, I, _ = Pl.xI(FD, np.array([[v1]]), PhiF)
        return x - XT[0], I, v1
    return make


XT = [0.0]

print('== E47：相关噪声的结构化配预算扫描 ==')
print('')
print(r'[1] 先把二分的两个前提量出来：$x$ 对 $v_1$ 单调增；$\rho=0$ 时 $x$ 与 $v_2$ 无关（E46 (i)）。')
for v1 in [1e-6, 1e-3, 1e0]:
    xs = []
    for v2 in [1e-2, 1e2, 1e6]:
        xs.append(Pl.xI(F, np.diag([v1, v2]), PhiF)[0])
    lo = Pl.xI(F, np.diag([v1 * .5, 1.0]), PhiF)[0]
    hi = Pl.xI(F, np.diag([v1 * 2., 1.0]), PhiF)[0]
    print(r'  v1=%9.0e  x across v2: %s  spread %.1e   半/双 v1: %.9f < %.9f < %.9f  %s'
          % (v1, ' '.join('%.9f' % q for q in xs), max(xs) - min(xs), lo, xs[1], hi,
             'ok' if lo < xs[1] < hi else '!! 单调性不成立'))
print(r'  两族地板：Phi0(full) 由 ticket 给 %.9f，逐位等于 Phi0(del)（#22 的核心数）。' % PhiF)

print('')
print(r'[2] 主表：每个预算 $x_t$、每个 $v_2$，扫 $\rho$ 并把 $v_1$ 二分回同一预算')
for xt in [1e-6, 1e-4, 1e-2, 1e-1, 1.0]:
    XT[0] = xt
    ed = solve_xt(FD, xt, mk_del())
    if ed is None:
        print(r'  x_t=%.0e  $r=1$ 族没 bracket 到，跳过' % xt)
        continue
    print(r'  x_t = %.1e   E_del = %.9f   配准残差 %.1e' % (xt, ed[0], ed[1]))
    worst = 1e30
    for v2 in [1e-1, 1e2, 1e6, 1e12]:
        cand = []
        for rho in RHO:
            got = solve_xt(F, xt, mk_full(rho, v2))
            if got is None or abs(got[1]) > 1e-9 * max(xt, 1e-30):
                continue
            cand.append((got[0], rho, got[2]))
        if not cand:
            print(r'     v2=%8.0e  无可行点（该 $v_2$ 下预算不可达）' % v2)
            continue
        best = min(cand)
        i0 = [q[0] for q in cand if q[1] == 0.0]
        worst = min(worst, best[0] - ed[0])
        print(r'     v2=%8.0e  min I=%.9f 在 rho=%+.2f (v1=%.3e)  差 E_del %+.3e  |  rho=0: %s'
              % (v2, best[0], best[1], best[2], best[0] - ed[0],
                 ('%.9f 差 %+.3e' % (i0[0], i0[0] - ed[0])) if i0 else '不可达'))
    print(r'     本预算下跨所有 $v_2$ 的最负偏差 %+.3e bit' % worst)

print('')
print(r'[3] 判决怎么读')
print(r'  每一行"差 E_del"都 $\ge0$ 就说明 $\rho\ne0$ 在这个植物上买不到东西，')
print(r'  E44 那张"随机 4000 点"的表可以被这张换掉：扫描结构化、预算严格配平、$v_2$ 跨 13 个数量级。')
print(r'  出现负值就是反例，负多少即 #23 第一节那张表要改的幅度，E46 的命题只能保留"块对角子族"这个限定词。')
print(r'  顺带一条：$\rho=0$ 那一列与 $v_2$ 无关（$x$ 不看 $v_2$、$I$ 的第二项被 $v_2$ 摊薄），')
print(r'  如果表里 $\rho=0$ 的读数随 $v_2$ 漂了，那是二分没配准，先查表头那列残差再查物理。')
