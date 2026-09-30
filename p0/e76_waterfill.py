r"""E76 = 把 E75 的 AM-GM 换成**精确的逆注水**，并用它直接造传感器。

E75 的读数（p0/e75_out.txt）：
 [0] 墙复现 OK；[1] 界在 8 个实测前沿点全部不越界 OK；[3] detGamma 单调性无反例；
 [2] 蓝平面 I=3 的 (*) 界 = 40.37445 < 43.711 => 验收线**没有**被闭式界认证。
但同一份输出里有一条结构性线索：冻结 X 版界（61.02/53.20/50.08/47.75/46.94）离实测
（62.90/54.29/50.60/47.94/47.10）只差 0.5~1.9，而墙版 (*) 差 3~14。
=> 瓶颈不是不等式，是 X -> X_w 这一步松弛；AM-GM 本身还有可提的精度。

== 精确化：min tr(Gamma R) s.t. det R >= delta, 0 <= R <= I ==
设 gamma_1..gamma_r 是 Gamma(X) 的特征值，则逐特征方向可分：
    min sum gamma_i r_i  s.t.  prod r_i >= delta, 0 <= r_i <= 1
    =>  r_i = min(1, nu/gamma_i),   nu^m = delta * prod_{gamma_i>nu} gamma_i
（m = 未截断的方向数）。这就是**逆注水**。关键的一点：截断 r_i=1 等价于 T_i=0，即
"该方向不传感" —— 于是 E63f/E69 实测到的前沿胜者是 rank-1（rho=0，关掉一条通道）
不再是数值巧合，而是注水解的结构性预言：率预算越小，nu 越大，被截断的小 gamma 方向越多。
AM-GM 是它的下界（等号当且仅当 R 在 Gamma 的本征基下各向同性，即无截断），所以
墙版 (*) 只会更弱，这是可以量化比较的。

== 逆注水给出的闭式传感器 ==
K := G^{1/2} M G^{1/2}，则 T = K(I+K)^{-1}、R = (I+K)^{-1}（对奇异 M 由连续性仍成立），
所以注水解直接给出 M = G^{-1/2} K G^{-1/2}，K = V diag((gamma_i-nu)^+/nu) V^t，
V = Gamma 的本征向量。给定 X 就有一个传感器；X 又是传感器的定点，所以套一层自洽迭代，
再用 nu 二分命中目标率。**若这条闭式路线能重现 E63f 穷尽网格/Powell 拿到的前沿值，
就等于给出"前沿设计 = 逆注水"的正面定理**（这是我想要的定理级结论，不是微雕）。

判据（先写下，跑完不许改）：
 [A] 代数恒等式自检：对 E75 的实测红平面设计，核
     D = j_c + tr(Theta X) - tr(Gamma) + tr(Gamma R) 与 det R = 2^{-2I}。
     差 > 1e-6 说明我这套 T-约化是错的，后面全废。
 [B] 单调性压力测试：det Gamma(X) 在 X = X_w + P (P 随机 PSD, 跨三个谱宽尺度) 上是否
     总 >= det Gamma(X_w)。出现反例 => (*) 的墙版只能降级为"冻结 X 版"。
 [C] 界精度对表：同一 X 上比较 AM-GM 值与注水精确值，再与实测 D 比。
     注水界必须 <= 实测 D（越界=我的注水解算错）；若差距普遍 < 1，则"冻结 X 松弛基本紧"。
 [D] 硬门（最重要）：自洽注水传感器在 (red I=2,2.5,3,4,4.90)、(blue I=2,3,4)、
     (free I=2,3,4) 上的 D 值。它是**可达上界**，所以必须 >= 该点已知最好前沿值
     （red 62.8975/54.2924/50.6023/47.9396/47.1006, blue 59.5611/45.4537/40.8072,
      free 58.1931/42.0427/36.5386）。低于已知值 = 我的 (D,I) 评估口径错，不是新纪录。
     与已知值之差 < 1% => 前沿被闭式解释，[E] 的排名预言才有意义。
 [E] 秩亏损预言：逐点报被截断方向数，与实测形状对照（红平面 rho=0 出现在 I=2,2.5,3,
     rho=0.01@I=4, rho=0.1@I=4.90，即率越大开的方向越多）。预言错就明说。
"""
import sys, time
sys.stdout.reconfigure(encoding='utf-8')
import numpy as np
from scipy.linalg import schur

np.set_printoptions(precision=4, suppress=True, linewidth=170)
_src = open('p0/exp_c_audit.py', encoding='utf-8').read().split("print(r'== E53")[0]
_ns = {'__name__': 'p'}
exec(compile(_src, 'p0/exp_c_audit.py[preamble]', 'exec'), _ns)
A, W, TH, JC, n, sym = _ns['A'], _ns['W'], _ns['TH'], _ns['JC'], _ns['n'], _ns['sym']
ln2 = np.log(2.0)
t0 = time.time()
print('== E76：逆注水前沿下界 + 闭式传感器 ==')


def stat(S, maxit=30000):
    """设计的定点 (X_prior, P_post)；不收敛返回 None。与 E75 同一实现。"""
    Pt = W.copy()
    for _ in range(maxit):
        Pm = sym(np.linalg.inv(np.linalg.inv(Pt) + S))
        Pn = sym(A @ Pm @ A.T + W)
        if not np.all(np.isfinite(Pn)) or np.max(np.abs(Pn)) > 1e14:
            return None
        if np.max(np.abs(Pn - Pt)) < 1e-13 * max(1.0, np.max(np.abs(Pn))):
            return Pt, Pm
        Pt = Pn
    return None


def wall_prior(Z, maxit=60000):
    X = W.copy()
    for _ in range(maxit):
        G = sym(Z.T @ X @ Z)
        P = sym(X - X @ Z @ np.linalg.solve(G, Z.T @ X))
        Xn = sym(A @ P @ A.T + W)
        if not np.all(np.isfinite(Xn)):
            return None
        if np.max(np.abs(Xn - X)) < 1e-13 * max(1.0, np.max(np.abs(Xn))):
            return Xn
        X = Xn
    return X


def gpow(G, p):
    """G^{p}，p = +1/2 用于 K = G^{1/2} M G^{1/2}，p = -1/2 用于 M = G^{-1/2} K G^{-1/2}。"""
    w, V = np.linalg.eigh(G)
    return (V * np.clip(w, 1e-300, None) ** p) @ V.T


def gamma_of(Z, X):
    """Gamma = G^{-1/2}(Z^t X Theta X Z)G^{-1/2}，连同 G、G^{-1/2}。"""
    G = sym(Z.T @ X @ Z)
    Lam = sym(Z.T @ X @ TH @ X @ Z)
    Gi2 = gpow(G, -0.5)
    return sym(Gi2 @ Lam @ Gi2), G, Gi2


def wf_R(gam, nu):
    """逆注水：r_i = min(1, nu/gamma_i)。返回 (R 的谱, 截断数)。"""
    r = np.minimum(1.0, nu / np.clip(gam, 1e-300, None))
    return r, int(np.sum(gam <= nu))


def DI(S):
    sp = stat(S)
    if sp is None:
        return None
    X, Pm = sp
    D = JC + float(np.trace(TH @ Pm))
    I = 0.5 * (np.linalg.slogdet(X)[1] - np.linalg.slogdet(Pm)[1]) / ln2
    return D, I, X, Pm


Ts, U = schur(A, output='real', sort=lambda a: abs(a) < 1.0)[:2]
PLANES = {'free': U[:, :], 'blue': U[:, n - 3:], 'red': U[:, n - 2:]}
Z2, Z3 = PLANES['red'], PLANES['blue']
RHO_TH = {2.0: (0.0, 138.0), 2.5: (0.0, 132.0), 3.0: (0.0, 130.0), 4.0: (0.01, 128.0),
          4.90: (0.1, 128.0)}
KNOWN = {'red': {2.0: 62.8975, 2.5: 54.2924, 3.0: 50.6023, 4.0: 47.9396, 4.90: 47.1006},
         'blue': {2.0: 59.5611, 3.0: 45.4537, 4.0: 40.8072},
         'free': {2.0: 58.1931, 3.0: 42.0427, 4.0: 36.5386}}
WALLS = {'red': 46.12314, 'blue': 36.60320, 'free': 31.48334}

Xw = {k: wall_prior(Z) for k, Z in PLANES.items()}


def scale_to(Z, S0, It, steps=44):
    """沿射线 10^m * S0 二分命中 I(S)=It，返回 (m, D, I)。E75 的口径。"""
    lo, hi = -3.0, 13.0
    for _ in range(steps):
        m = 0.5 * (lo + hi)
        a = DI((10.0 ** m) * S0)
        if a is None:
            hi = m
            continue
        if a[1] >= It:
            hi = m
        else:
            lo = m
    a = DI((10.0 ** hi) * S0)
    return hi, a


print('\n[A] 代数恒等式自检（E75 的实测红平面设计）')
ok = True
for It in sorted(KNOWN['red']):
    rho, th_deg = RHO_TH[It]
    th = np.radians(th_deg)
    Rv = np.array([[np.cos(th), -np.sin(th)], [np.sin(th), np.cos(th)]])
    M0 = sym(Rv @ np.diag([1.0, rho]) @ Rv.T)
    m, a = scale_to(Z2, Z2 @ M0 @ Z2.T, It)
    if a is None:
        print('  I=%.2f 定点不收敛' % It); ok = False; continue
    D, Iv, X, Pm = a
    gam, G, Gi2 = gamma_of(Z2, X)
    # R = (I+K)^{-1}, K = G^{1/2} M G^{1/2}; 用 S 的尺度还原 M
    M = (10.0 ** m) * M0
    K = sym(gpow(G, 0.5) @ M @ gpow(G, 0.5))
    R = sym(np.linalg.inv(np.eye(2) + K))
    trGR = float(np.trace(gam @ R))
    lhs = D
    rhs = JC + float(np.trace(TH @ X)) - float(np.trace(gam)) + trGR
    dR = float(np.linalg.det(R))
    e1 = abs(lhs - rhs)
    e2 = abs(dR - 2.0 ** (-2.0 * Iv))
    print('  I=%4.2f: D=%9.4f | j_c+tr(TH X)-tr(Gamma)+tr(Gamma R)=%9.4f 差 %.2e | detR=%.6e vs 2^-2I=%.6e 差 %.2e'
          % (It, lhs, rhs, e1, dR, 2.0 ** (-2.0 * Iv), e2))
    ok &= e1 < 1e-6 and e2 < 1e-6
print('  => %s' % ('T-约化恒等式成立' if ok else '**恒等式不成立，E76 全废**'))

print('\n[B] det(Gamma) 单调性压力测试（X = X_w + P，P 随机 PSD，三种谱宽）')
rng = np.random.default_rng(20260929)
for tag in ('red', 'blue', 'free'):
    Z = PLANES[tag]
    X_ = Xw[tag]
    r0 = float(np.linalg.det(gamma_of(Z, X_)[0]))
    nbad = 0
    ntot = 0
    worst = 0.0
    for scale in (1e-3, 1.0, 1e3):
        for _ in range(60):
            B = rng.normal(size=(n, n)) * scale
            P = sym(B @ B.T)
            r1 = float(np.linalg.det(gamma_of(Z, sym(X_ + P))[0]))
            ntot += 1
            if r1 < r0:
                nbad += 1
                worst = max(worst, r0 - r1)
    print('  %-5s detGamma(X_w)=%.4e | 随机 PSD 扰动 %d 次，出现 detGamma 变小的 %d 次，最大违反 %.3e'
          % (tag, r0, ntot, nbad, worst))

print('\n[C] 同一 X 上：AM-GM vs 注水精确值 vs 实测（红平面墙先验与前沿先验）')
for tag in ('red', 'blue', 'free'):
    r = PLANES[tag].shape[1]
    gam, G, Gi2 = gamma_of(PLANES[tag], Xw[tag])
    g = np.sort(np.linalg.eigvalsh(gam))[::-1]
    print('  %s lam(Gamma(X_w)) = %s  detGamma=%.4e' % (tag, np.array2string(g, precision=4),
                                                        float(np.linalg.det(gam))))
    for It in (2.0, 3.0, 4.0):
        delta = 2.0 ** (-2.0 * It)
        # 注水：nu 由 prod min(1,nu/g_i) = delta 决定，单调，二分
        lo, hi = 1e-12 * g[0], g[0] * 1e6
        for _ in range(200):
            nu = 0.5 * (lo + hi)
            rr, nc = wf_R(g, nu)
            # 截断方向贡献 det=1，未截断贡献 nu/g
            m = r - nc
            prod = (nu ** m) / np.prod(g[:len(g) - nc]) if m > 0 else 1.0
            # g 降序，未截断的是最大的 m 个
            if prod < delta:
                lo = nu
            else:
                hi = nu
        nu = 0.5 * (lo + hi)
        rr, nc = wf_R(g, nu)
        wf = float(np.dot(g, rr))
        am = r * float(np.linalg.det(gam)) ** (1.0 / r) * delta ** (1.0 / r)
        kn = KNOWN[tag].get(It)
        base = WALLS[tag]
        print('    I=%.2f: 墙+AMGM=%9.4f  墙+注水=%9.4f  (注水-AMGM=%+.4f, 截断 %d/%d 方向)%s'
              % (It, base + am, base + wf, wf - am, nc, r,
                 '' if kn is None else '  | 实测可达 %9.4f' % kn))

print('\n[D] 自洽注水传感器（闭式造设计）对已知前沿')
gate = True
for tag in ('red', 'blue', 'free'):
    Z = PLANES[tag]
    r = Z.shape[1]
    for It in sorted(KNOWN[tag]):
        # 内层自洽：X -> Gamma(X) -> 注水 M(nu) -> 定点 X
        def build(nu, X0):
            X = X0
            S = None
            for _ in range(30):
                gam, G, Gi2 = gamma_of(Z, X)
                gg = np.sort(np.linalg.eigvalsh(gam))[::-1]
                w, V = np.linalg.eigh(sym(gam))
                kk = np.maximum(0.0, w - nu) / nu
                K = sym(V @ np.diag(kk) @ V.T)
                M = sym(Gi2 @ K @ Gi2)
                S = sym(Z @ M @ Z.T)
                a = DI(S)
                if a is None:
                    return None
                Xn = a[2]
                if np.max(np.abs(Xn - X)) < 1e-10 * max(1.0, np.max(np.abs(Xn))):
                    X = Xn
                    break
                X = Xn
            return (DI(S), nu, int(np.sum(np.linalg.eigvalsh(gam) <= nu)))
        # 外层二分 nu 命中 I：nu 越大 => 截断越多 => I 越小
        gam0, _, _ = gamma_of(Z, Xw[tag])
        gmax = float(np.max(np.linalg.eigvalsh(gam0)))
        lo, hi = 1e-6 * gmax, 50.0 * gmax
        res = None
        for _ in range(34):
            nu = 0.5 * (lo + hi)
            b = build(nu, Xw[tag])
            if b is None:
                hi = nu
                continue
            a, nu_, nc = b
            if a[1] > It:
                lo = nu
            else:
                hi = nu
        b = build(0.5 * (lo + hi), Xw[tag])
        if b is None:
            print('  %-5s I=%.2f: 自洽迭代失败' % (tag, It))
            gate = False
            continue
        a, nu_, nc = b
        D, Iv = a[0], a[1]
        kn = KNOWN[tag][It]
        rel = (D - kn) / kn
        bad = D < kn - 1e-9
        print('  %-5s I=%.2f: 注水 D=%9.4f I=%7.4f | 已知可达 %9.4f | 差 %+.4f (%+.2f%%) 截断 %d/%d %s'
              % (tag, It, D, Iv, kn, D - kn, 100 * rel, nc, r,
                 '**低于已知：口径错**' if bad else 'OK'))
        gate &= not bad
print('  => 硬门 %s' % ('通过：闭式注水没有伪造新纪录' if gate else '**未通过**'))
print('\n[E] 秩亏损预言 vs 实测形状（红平面 rho）')
print('  实测（E63f/E69）：I=2,2.5,3 时 rho=0（关一条通道）；I=4 rho=0.01；I=4.90 rho=0.1')
print('  耗时 %.0f s' % (time.time() - t0))
