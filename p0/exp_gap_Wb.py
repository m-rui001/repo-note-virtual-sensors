r"""E38b：把 E38 的两条推论钉死。

[4] **成对不变的解析恒等式**。$P(\varepsilon W,V)=\varepsilon P(W,V/\varepsilon)$ 微分一次给
$$\mathrm dP_{\varepsilon W}(V)[Z]=\mathrm dP_{W}(V/\varepsilon)[Z]$$
    即 $V\mapsto P$ 的 Fréchet 导数**完全不变**，只要把自变量一起缩放。这比 E37/E36 的两点外推强：
    它用闭式算子（E37 的 `dP_op_p`）在有限 $V$ 上验，不涉 extrapolation，应当打到机器精度。
    同时它给出 $C$ 的零次齐次性是**恒等式**而非渐近断言。

[5] **台阶下移的可反驳预言**。E38 测到 $\mathrm{gap}_A\propto\varepsilon$、$\lambda_i(C)$ 不动，
    于是锚点那一格 $\mathrm{gap}_{\{1\}}=1.084$ 应当在 $\varepsilon=10^{-2}$ 处变成 $1.084\times10^{-2}$。
    换成观测量：在 $x\in[5\times10^{-3},2\times10^{-2}]$ 这段窗口，包络斜率
    $\mathrm dI/\mathrm d\log_2(1/x)$ 在 $\varepsilon=1$ 还在过渡（$\approx0.84$），
    在 $\varepsilon=10^{-2}$ 应当已经坐满 $\operatorname{rank}C/2=1$。
    这条如果反了，齐次性就只能停在代数层面、不算物理结论。
"""
import sys
import numpy as np
sys.stdout.reconfigure(encoding='utf-8')
import p0.exp_plants_oos as E35
import p0.exp_ladder as X37
import p0.exp_C_exact as X36

sym = X36.sym
ln2 = np.log(2.0)


def rebuild(Pl, Wnew):
    return E35.Pl(Pl.A, Pl.B, Wnew, Pl.Q, Pl.R)


base = None
for nm, pl, F in E35.plants():
    if nm.find('锚点') >= 0:
        base = (nm.replace('$', ''), pl, F)
        break
nm0, Pl0, F0 = base
r = F0.shape[0]
print('== E38b ==  算例：%s  $n=%d$ $r=%d$' % (nm0, Pl0.nn, r))

# ---------- [4] dP 的成对不变量 ----------
print('\n[4] $\mathrm dP_{\\varepsilon W}(V)=\mathrm dP_{W}(V/\\varepsilon)$（闭式算子，无外推）')
ts = (1e-3, 1e-5, 1e-7)
for eps in (1e-2, 1e2):
    Pe = rebuild(Pl0, eps * Pl0.W)
    rows = []
    for t in ts:
        V = t * np.eye(r)
        Ja, _p, _b, dev_a = X37.dP_op_p(Pe, F0, V)
        Jb, _p2, _b2, dev_b = X37.dP_op_p(Pl0, F0, V / eps)
        rows.append(np.linalg.norm(Ja - Jb) / np.linalg.norm(Jb))
    print('    $\\varepsilon=%-5.0e$ 相对偏差 %s   （算子残差 dev：%s / %s）'
          % (eps, ' '.join('%.2e' % q for q in rows),
             ' '.join('%.1e' % q for q in (X37.dP_op_p(Pe, F0, t * np.eye(r))[3] for t in ts[:1])),
             ' '.join('%.1e' % q for q in (X37.dP_op_p(Pl0, F0, t / eps * np.eye(r))[3] for t in ts[:1]))))

# ---------- [5] 台阶下移 ----------
print('\n[5] 包络斜率窗口 $x\in[5\times10^{-3},2\times10^{-2}]$：预言 $\\varepsilon=10^{-2}$ 已到 $\\operatorname{rank}C/2=1$')


def env_slope(Pl, eps, X=(2e-2, 1e-2, 5e-3)):
    fx = X37.floor_exact(Pl, F0, None)
    phi0 = fx['闭式']
    C = X37.C_of(Pl, F0)
    fams, e, U, rho = X37.fams_of(Pl, F0, C, phi0)
    ts = tuple(np.geomspace(3e-2, 3e-8, 60))
    curves = []
    for f in fams:
        rows = X37.curve(Pl, f, ts, phi0)
        if len(rows) < 8:
            continue
        lx = np.array([np.log2(1.0 / q[1]) for q in rows])
        I = np.array([q[2] for q in rows])
        ok = np.isfinite(lx) & np.isfinite(I)
        curves.append((lx[ok], I[ok], f['name'], f['gap']))
    out = []
    for xq in X:
        uq = np.log2(1.0 / xq)
        cand = []
        for lx, I, nm, gp in curves:
            if lx.min() <= uq <= lx.max():
                cand.append((np.interp(uq, lx, I), nm))
        best = min(cand, key=lambda q: q[0]) if cand else (np.nan, '-')
        out.append((xq, uq, best[0], best[1]))
    sl = [(out[j + 1][2] - out[j][2]) / (out[j + 1][1] - out[j][1]) for j in range(len(out) - 1)]
    print('  $\\varepsilon=%-6.0e$ 地板 $%.4e$  各 $\\mathrm{gap}$：%s'
          % (eps, phi0, ' '.join('%s:%.3e' % (f['name'], f['gap']) for f in fams if len(f['A']) == 1)))
    print('      ' + '   '.join('$x=%.0e$:%.4f($A$%s)' % (o[0], o[2], o[3]) for o in out)
          + '   斜率 ' + ' '.join('%.4f' % q for q in sl))


env_slope(Pl0, 1.0)
env_slope(rebuild(Pl0, 1e-2 * Pl0.W), 1e-2)
env_slope(rebuild(Pl0, 1e-3 * Pl0.W), 1e-3)
