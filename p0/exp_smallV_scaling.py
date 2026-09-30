"""E27：直接测 $x(\gamma)=\operatorname{tr}(\Theta[P_\infty(\gamma^{-1}I)-P_0])$ 的对数标度指数。

Prop E 的步(b）断言 $P_\infty(V)$ 在 $V=0$ 沿对角方向可微，于是
    $x=\sum_j c_j V_{jj}+O(|V|^2)$，各向同性支上 $x\propto 1/\gamma$。
这条断言有真实的反面文献：奇异摄动的小测量噪声展开里常用 $\varepsilon^2=R$ 作小量
（克里希纳-森古普塔 1980；哈利勒的降阶观测器线），那种标度下 $x\propto 1/\sqrt\gamma$，
拉格朗日给 $\gamma_j\propto1/x^2$，系数变成 $r$ 而非 $r/2$——差整整一倍。

所以这不是修辞问题，是可测的：拟合 $\log x$ vs $\log\gamma$ 的斜率，$-1$ 支持步(b)，$-0.5$ 支持 $\sqrt{\ }$ 支。
顺带用标量例子里的 ARE 闭式解做同尺度的对照（标量情形 $P$ 有解析式，指数一测就知道）。
"""
import sys
import numpy as np
sys.stdout.reconfigure(encoding='utf-8')
from p0.p0_replicate_letter import A, W, n, sym, ln2, ctrl, rate_cost
from p0.exp_authoritative import noiseless
from p0.exp_nearfloor_law import Fs, floor_of, Th, Pc

UNC = float(np.trace(W @ Pc))
lgs = np.linspace(2.0, 13.0, 45)          # gamma = e^lg，跨 11 个数量级


def slope(logg, logx):
    M = np.vstack([logg, np.ones_like(logg)]).T
    a, b = np.linalg.lstsq(M, logx, rcond=None)[0]
    r = logx - (a * logg + b)
    return a, np.max(np.abs(r))


if __name__ == '__main__':
    print('D_min^unc = %.4f；步(b) 预言 x ∝ 1/γ（指数 −1），竞争支 x ∝ 1/√γ（指数 −0.5）\n')

    for tag, F in Fs.items():
        P0, res, it = noiseless(F, A, W, iters=80000, tol=1e-15)
        phi = float(np.trace(Th @ P0))
        r = F.shape[0]
        xs = []
        for lg in lgs:
            I, J, _, _ = rate_cost(F, np.exp(lg) * np.eye(r), Th, Pc)
            xs.append(J - UNC - phi)
        xs = np.array(xs)
        ok = xs > 1e-9
        lg, lx = lgs[ok], np.log(xs[ok])
        a_full, r_full = slope(lg, lx)
        # 最深的三分之一（真正接近 V=0 的渐近段）
        k = max(3, len(lg) // 3)
        a_deep, r_deep = slope(lg[-k:], lx[-k:])
        print('[%s]  rank=%d  Φ=%.4f  D_min(K_F)=%.4f  可用 %d 点 (γ 到 e^%.1f)'
              % (tag, r, phi, UNC + phi, ok.sum(), lgs[ok][-1]))
        print('   全窗指数 %.4f (残差 %.3f)   深端 1/3 指数 %.4f (残差 %.3f)   '
              '距 −1 偏 %+.1f%%，距 −0.5 偏 %+.1f%%'
              % (a_full, r_full, a_deep, r_deep, 100 * (a_deep + 1), 100 * (a_deep + .5)))
        print('   x 量级: %.3e → %.3e' % (xs[ok][0], xs[ok][-1]))

    # 标量 ARE 闭式对照：P = a^2 P + w - a^2 P^2/(P+v)
    print('\n标量例子里 $P$ 有闭式解，同一检验不用数值 Riccati：')
    a2, w = 1.3127 ** 2, 1.0
    P0 = 0.0                                   # F=1 时 V=0 的极限解（完美观测，P_-=0）
    print('%8s %12s %12s %12s' % ('gamma', 'x(数值解)', 'x/(1/γ)', 'x/(1/√γ)'))
    for lg in (4., 8., 12., 15.):
        v = np.exp(-lg)
        # H=I：P = Σv/(Σ+v)，Σ=a²P+w  ⟹  a²P² + [w+v(1-a²)]P − wv = 0
        s = np.roots([a2, w + v * (1 - a2), -w * v])
        P = max([s_.real for s_ in s if abs(s_.imag) < 1e-9])
        x = P - P0
        print('%8.0e %12.5e %12.4f %12.4f' % (np.exp(lg), x, x * np.exp(lg), x * np.exp(.5 * lg)))
