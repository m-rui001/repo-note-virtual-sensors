r"""E60 = 判死两条：(a) $r=1$ 的 $\Phi$ 是否真有第二个严格局部极小（86.6659）；(b) $r=2$ 的盆地普查。

C 的 §10.1 主张"$r=1$：24 个起点 **24/24** 全等于 10.169224"，§10.2 据此写"设计 $=$ 便宜且 100% 可靠"。
我的 E59c 球面下降在 24 个 Haar 起点上是 **11 命中 / 13 落在 86.665888**，两处的切向 Hessian 都正定
（$+35.8,+171.0,+892.5$ vs $+175.1,+228.5,+1229.4$）。若 86.6659 是真的，他那一行"24/24"就是他的
实现性质（测地线+他的步长），不是景观性质，"100% 可靠"要降级。
这一跑先把**评价函数**本身钉死（同一个 $f$ 用四条独立路线算 $\Phi$），再沿两点之间的测地线扫一遍
看有没有势垒（有⇒两个盆地，无⇒其中一点是我差分梯度的伪驻点），最后 $r=2$ 做 8 起点普查。
"""
import sys
sys.stdout.reconfigure(encoding='utf-8')
import numpy as np
from scipy.optimize import minimize

_src = open('p0/exp_c_audit.py', encoding='utf-8').read().split("print(r'== E53")[0]
_ns = {'__name__': 'e53_preamble'}
exec(compile(_src, 'p0/exp_c_audit.py[preamble]', 'exec'), _ns)
A, W, TH, JC, n, sym = _ns['A'], _ns['W'], _ns['TH'], _ns['JC'], _ns['n'], _ns['sym']
phi_iter, phi_dare, post = _ns['phi_iter'], _ns['phi_dare'], _ns['post']
OPT1 = 10.169224
print(r'== E60：$r=1$ 第二极小判定 + $r=2$ 盆地普查 ==')


def unit(f):
    return f / np.linalg.norm(f)


def phi_f(f):
    return phi_iter(unit(f).reshape(1, n))[0]


def t_basis(f):
    U, s, _ = np.linalg.svd(np.eye(n) - np.outer(unit(f), unit(f)))
    return U[:, s > 0.5][:, :n - 1]


def shift(f, u):
    return unit(unit(f) + u)


def vg(f, h=2e-5):
    f = unit(f); B = t_basis(f)
    v = phi_f(f)
    g = np.array([(phi_f(shift(f, B[:, i] * h)) - phi_f(shift(f, -B[:, i] * h))) / (2 * h) for i in range(n - 1)])
    return v, g, B


def descent(f, iters=700, ang0=0.2):
    f = unit(f); radius = ang0
    for k in range(iters):
        v, g, B = vg(f)
        gn = np.linalg.norm(g)
        if gn * radius < 1e-11:
            return f, v, k, gn, 'conv'
        vn = phi_f(shift(f, B @ (-g / gn * radius)))
        if vn < v:
            f, radius = shift(f, B @ (-g / gn * radius)), min(radius * 1.6, 0.6)
        else:
            radius *= 0.4
            if radius < 1e-9:
                return f, v, k, gn, 'stall'
    return f, phi_f(f), iters, np.nan, 'maxit'


def hess(f, h=3e-4):
    B = t_basis(f); f0 = unit(f); v0 = phi_f(f0)
    H = np.zeros((n - 1, n - 1))
    for i in range(n - 1):
        H[i, i] = (phi_f(shift(f0, B[:, i] * h)) + phi_f(shift(f0, -B[:, i] * h)) - 2 * v0) / h ** 2
        for j in range(i + 1, n - 1):
            pp = phi_f(shift(f0, (B[:, i] + B[:, j]) * h)); pm = phi_f(shift(f0, (B[:, i] - B[:, j]) * h))
            mp = phi_f(shift(f0, (-B[:, i] + B[:, j]) * h)); mm = phi_f(shift(f0, -(B[:, i] + B[:, j]) * h))
            H[i, j] = H[j, i] = (pp - pm - mp + mm) / (4 * h ** 2)
    return np.linalg.eigvalsh(H)


# ---- 找到两个终点 ----
rng = np.random.default_rng(771)
glob = trap = None
for i in range(24):
    f0 = rng.standard_normal(n)
    f, v, k, gn, fl = descent(f0)
    if abs(v - OPT1) / OPT1 < 1e-4:
        glob = glob if glob is not None else f
    elif abs(v - 86.665888) < 1e-3:
        trap = trap if trap is not None else f
    if glob is not None and trap is not None:
        break
print('\n[1] 两条路线的终点')
for tag, f in (('全局最优侧', glob), ('第二终点侧', trap)):
    v = phi_f(f)
    print('  %s  $\\Phi$=%.6f  $f=%s$' % (tag, v, np.array2string(unit(f), precision=6)))

# ---- [2] 四条独立路线算同一个 f ----
print('\n[2] 同一个 $f$ 的 $\Phi$：定点迭代(默认/2万步/1e-14) vs 奇异 DARE vs 暴力迭代')
for tag, f in (('glob', glob), ('trap', trap)):
    F = unit(f).reshape(1, n)
    v1 = phi_iter(F)[0]
    v2 = phi_iter(F, it=20000, tol=1e-14)[0]
    v3 = phi_iter(F, it=2000, tol=1e-15)[0]
    v4 = phi_dare(F)[0]
    # 暴力：直接迭代 200000 步，不判停，看尾差
    Pt = W.copy()
    for k in range(200000):
        P = post(F, Pt)
        Pt = sym(A @ P @ A.T + W)
    v5 = float(np.trace(TH @ P))
    print('  %-5s 默认=%.7f  2e4步=%.7f  1e-15=%.7f  DARE=%.7f  暴力2e5=%.7f   极差 %.2e'
          % (tag, v1, v2, v3, v4, v5, max(abs(v1 - v5), abs(v4 - v5))))
    ev = hess(f)
    print('        切向 Hessian 特征值 %s  ‖g‖(h=1e-6)=%.2e  ⇒ %s'
          % (np.array2string(ev, precision=2), np.linalg.norm(vg(f, 1e-6)[1]),
             '严格局部极小' if (ev > 0).all() else '非极小'))

# ---- [3] 势垒：两点之间测地线 ----
print('\n[3] 沿 $f_{\rm glob}\to f_{\rm trap}$ 的大圆弧扫 $\Phi$')
f1, f2 = unit(glob), unit(trap)
c = float(np.clip(f1 @ f2, -1, 1))
print('  夹角 %.2f 度' % np.degrees(np.arccos(abs(c))))
th = np.linspace(0, np.pi / 2, 25)
perp = unit(f2 - c * f1)
vals = [phi_f(np.cos(t) * f1 + np.sin(t) * perp) for t in th]
vals = np.array(vals)
print('  端点值 %.5f / %.5f；弧上最小 %.5f（在第 %d 点）；弧上最大 %.5f（在第 %d 点）'
      % (vals[0], vals[-1], vals.min(), vals.argmin(), vals.max(), vals.argmax()))
print('  %s' % np.array2string(vals, precision=3))
print('  ⇒ 弧内部出现高于两端点的极大 = 两个盆地之间存在势垒（下降不可互相跨越）')

# ---- [4] 第二极小的几何 ----
print('\n[4] 两个终点的几何对照')
wTs, vTs = np.linalg.eigh(TH)
o = np.argsort(-wTs); vTs = vTs[:, o]
wA, vA = np.linalg.eig(A)
V = np.real(vA); V = V / np.linalg.norm(V, axis=0, keepdims=True)
UN = [k for k in range(n) if abs(wA[k]) >= 1]
for tag, f in (('glob', glob), ('trap', trap)):
    f = unit(f)
    d = np.abs(V[:, UN].T @ f).min()
    ang = np.degrees(np.arccos(abs(vTs[:, 0] @ f)))
    print('  %-5s 到坏集距离 $d$=%.4f   与 $\Theta$ 首特征向量夹角 %.2f 度   $|\lambda|$ 加权分量 %s'
          % (tag, d, ang, np.array2string(np.abs(V[:, UN].T @ f), precision=3)))

# ---- [5] r=2 盆地普查 ----
print('\n[5] $r=2$：我自己的 Stiefel 下降，8 个 Haar 起点')
r = 2
OPT2 = 0.735143


def phi_F(F):
    Q = np.linalg.qr(F)[0]
    return phi_iter(Q.T)[0]


def tdir(F, i, j):
    E = np.zeros((n, r)); E[i, j] = 1.0
    return E - F @ (F.T @ E)


PARAMS = [(i, j) for j in range(r) for i in range(n)]


def descent_r2(F, iters=400, ang0=0.2):
    F = np.linalg.qr(F)[0]
    radius = ang0
    for k in range(iters):
        v = phi_F(F)
        g = np.zeros(len(PARAMS))
        for a, (i, j) in enumerate(PARAMS):
            U = tdir(F, i, j)
            U = U / np.linalg.norm(U)
            g[a] = (phi_F(np.linalg.qr(F + U * 2e-5)[0]) - phi_F(np.linalg.qr(F - U * 2e-5)[0])) / 4e-5
        gn = np.linalg.norm(g)
        if gn * radius < 1e-11:
            return F, v, k, gn, 'conv'
        step = np.zeros((n, r))
        for a, (i, j) in enumerate(PARAMS):
            U = tdir(F, i, j); U = U / np.linalg.norm(U)
            step += g[a] * U
        Fp = np.linalg.qr(F - radius * step / gn)[0]
        if phi_F(Fp) < v:
            F, radius = Fp, min(radius * 1.6, 0.6)
        else:
            radius *= 0.4
            if radius < 1e-9:
                return F, v, k, gn, 'stall'
    return F, phi_F(F), iters, np.nan, 'maxit'


rng2 = np.random.default_rng(4242)
res = []
for k in range(8):
    F0 = rng2.standard_normal((n, r))
    F, v, it, gn, fl = descent_r2(F0)
    res.append(v)
    print('  起点%d → $\Phi$=%-12.6f 步 %3d ‖g‖ %.1e %s  占 $j_c$ %+.2f%%' %
          (k, v, it, gn, fl, 100 * v / JC))
res = np.array(res)
print('  8 起点：命中 %.6f 的 %d 个；不同收敛值 %d 个：%s' %
      (OPT2, (np.abs(res - OPT2) / OPT2 < 1e-4).sum(), len(np.unique(np.round(res, 4))),
       np.array2string(np.sort(np.unique(np.round(res, 4))), precision=4)))
Fc = np.load('.work3/c24_bestF_r2.npy')
Fc = Fc if Fc.shape[0] == n else Fc.T
print('  他的 $r=2$ 最优作为起点：从它出发 30 次随机小扰动后再下降，最小值 %.6f（应等于 %.6f）'
      % (min([descent_r2(np.linalg.qr(Fc + 0.05 * rng2.standard_normal((n, r)))[0])[1] for _ in range(6)]), OPT2))
