r"""E59c = $r=1$ 的 $\Phi$ 在 $\mathbb{P}^3$ 上的驻点普查（回答 C §10.1 的"24/24 命中同一个值"）。

E59b 暴露两件事：
  (1) 我的下降**不是**停在驻点：停滞处 $\|g\|=15\sim 4.5\times10^3$，是我的线搜索退出条件写坏了
      （$\eta=\text{const}$ 的绝对阈值配上 step$=L/(1+\|g\|)$，在 $\|g\|$ 大时预测下降掉到 $10^{-12}$ 以下）。
      所以 E59b 那 0/24 **不能**当反驳用，只能当我自己的缺陷登出。
  (2) 但 Powell 从三个停滞点收敛到 **86.6659**，与全局最优 $10.169224$ 之间隔着一个真正的结构。
这一跑：可信域式球面下降（固定切向步长=角度，失败收缩、成功放大，停止条件用预测下降量），
24 个 Haar 起点 + 8 个"作弊"起点（$\Theta$ 首向量、标准基、他给的三个最优 $F$ 的扰动），
对每个收敛点做 $3\times3$ 切向 Hessian 中心差分 $\to$ 判极小/鞍/极大。
"""
import sys
sys.stdout.reconfigure(encoding='utf-8')
import numpy as np
import itertools
from scipy.linalg import schur

_src = open('p0/exp_c_audit.py', encoding='utf-8').read().split("print(r'== E53")[0]
_ns = {'__name__': 'e53_preamble'}
exec(compile(_src, 'p0/exp_c_audit.py[preamble]', 'exec'), _ns)
A, W, TH, JC, n, sym = _ns['A'], _ns['W'], _ns['TH'], _ns['JC'], _ns['n'], _ns['sym']
phi_iter = _ns['phi_iter']
OPT1 = 10.169224
print(r'== E59c：$r=1$ 驻点普查 ==')


def unit(f):
    return f / np.linalg.norm(f)


def phi_f(f):
    return phi_iter(unit(f).reshape(1, n))[0]


def t_basis(f):
    U, s, _ = np.linalg.svd(np.eye(n) - np.outer(unit(f), unit(f)))
    return U[:, s > 0.5][:, :n - 1]


def shift(f, u):
    return unit(unit(f) + u)


def val_grad(f, h=2e-5):
    f = unit(f)
    B = t_basis(f)
    v = phi_f(f)
    g = np.array([(phi_f(shift(f, B[:, i] * h)) - phi_f(shift(f, -B[:, i] * h))) / (2 * h)
                  for i in range(n - 1)])
    return v, g, B


def descent(f, iters=600, ang0=0.2):
    """可信域式：每步沿 $-g$ 走"角度 radius"，成功放大 1.6、失败收缩 0.4。"""
    f = unit(f)
    radius = ang0
    for k in range(iters):
        v, g, B = val_grad(f)
        gn = np.linalg.norm(g)
        if gn * radius < 1e-11:                     # 预测下降量小于评价噪声
            return f, v, k, gn, 'conv'
        u = -g / gn * radius
        vn = phi_f(shift(f, B @ u))
        if vn < v:                                  # Armijo 用纯下降比较
            f, radius = unit(unit(f) + B @ u), min(radius * 1.6, 0.6)
        else:
            radius *= 0.4
            if radius < 1e-9:
                return f, v, k, gn, 'stall'
    return f, phi_f(f), iters, np.nan, 'maxit'


def hess(f, h=3e-4):
    B = t_basis(f)
    H = np.zeros((n - 1, n - 1))
    f0 = unit(f)
    v0 = phi_f(f0)
    for i in range(n - 1):
        vi = phi_f(shift(f0, B[:, i] * h))
        for j in range(i, n - 1):
            if i == j:
                vj = phi_f(shift(f0, -B[:, i] * h))
                H[i, i] = (vi + vj - 2 * v0) / h ** 2
            else:
                pp = phi_f(shift(f0, (B[:, i] + B[:, j]) * h))
                pm = phi_f(shift(f0, (B[:, i] - B[:, j]) * h))
                mp = phi_f(shift(f0, (-B[:, i] + B[:, j]) * h))
                mm = phi_f(shift(f0, -(B[:, i] + B[:, j]) * h))
                H[i, j] = H[j, i] = (pp - pm - mp + mm) / (4 * h ** 2)
    return H


# ---- [1] 24 个 Haar 起点 ----
rng = np.random.default_rng(771)
starts = [('haar%d' % i, rng.standard_normal(n)) for i in range(24)]
wT, vT = np.linalg.eigh(TH)
vT = vT[:, np.argsort(-wT[np.argsort(-wT)])]
wTs, vTs = np.linalg.eigh(TH)
o = np.argsort(-wTs)
vTs = vTs[:, o]
starts += [('theta_top', vTs[:, 0]), ('std_e1', np.eye(n)[:, 0]), ('std_e4', np.eye(n)[:, 3]),
           ('W_top', np.linalg.eigh(W)[1][:, -1]),
           ('unstab1', np.real(np.linalg.eig(A)[1][:, np.argmax(np.abs(np.linalg.eig(A)[0]))]))]
for i in (2, 3, 4, 5):
    starts += [('haarx%d' % i, rng.standard_normal(n))]
F1 = np.load('.work3/c24_bestF_r1.npy')
F1 = F1 if F1.shape[0] == n else F1.T
starts += [('C_optF1', F1[:, 0]), ('C_optF1+30deg', unit(F1[:, 0] + 0.5 * rng.standard_normal(n)))]

print('\n[1] %-12s %-11s %-10s %-10s %-11s %s' %
      ('起点', '收敛Φ', '步', '终点‖g‖', '旗标', 'Hessian 特征值 / 类型'))
tab = []
for tag, f0 in starts:
    f, v, k, gn, fl = descent(f0)
    ev = np.linalg.eigvalsh(hess(f))
    typ = '极小' if (ev > 0).all() else ('鞍' if (ev > 0).any() else ('极大' if (ev < 0).all() else '退化'))
    tab.append((tag, v, k, gn, fl, ev, typ, f))
    print('  %-12s %-11.5f %-10d %-10.2e %-11s [%s] %s' %
          (tag, v, k, gn, fl, ' '.join('%+.1f' % x for x in ev), typ))

vals = np.array([t[1] for t in tab])
nh = (np.abs(vals - OPT1) / OPT1 < 1e-4).sum()
print('\n  全部 %d 个起点：落到全局最优 %.6f 的 %d 个；不同的收敛值 %d 个' %
      (len(tab), OPT1, nh, len(np.unique(np.round(vals, 4)))))
print('  收敛值列表：%s' % np.array2string(np.sort(np.round(vals, 4)), precision=4))
mins = [t for t in tab if t[6] == '极小' and np.abs(t[1] - OPT1) / OPT1 > 1e-4]
print('  **非全局的局部极小**：%s' % (['%.5f' % t[1] for t in mins] or '无'))

# ---- [2] 86.6659 那一格到底是什么 ----
print('\n[2] E59b 的 86.6659 复现与定性')
cand = [t for t in tab if abs(t[1] - 86.6659) < 2.0]
if cand:
    for t in cand:
        print('  起点 %s → %.6f  ‖g‖=%.2e  类型 %s  特征值 %s' % (t[0], t[1], t[3], t[6], np.array2string(t[5], precision=1)))
else:
    print('  本轮 32 个起点没有一个收敛到 86.6659 附近 ⇒ 那是 E59b 里 Powell 在坏线搜索终点上的提前终止')
    near = sorted(tab, key=lambda t: abs(t[1] - 86.6659))[:3]
    for t in near:
        print('    最近的三个收敛值：%s → %.5f（%s）' % (t[0], t[1], t[6]))

# ---- [3] Schur 14.6398 的对象判定 ----
print('\n[3] $\Phi=14.6398$ 是哪一个两维子空间')
wA, vA = np.linalg.eig(A)
V = np.real(vA)
V = V / np.linalg.norm(V, axis=0, keepdims=True)
print('  $\operatorname{eig}(A)=%s$' % np.array2string(wA, precision=4))
UN = [k for k in range(n) if abs(wA[k]) >= 1]
print('  不稳下标 %s（模 %s）；实 Schur 与它张成的不变子空间相同，故"不稳模"只有一种读法' %
      (UN, np.array2string(np.abs(wA[UN]), precision=4)))
best = []
for cols in itertools.combinations(range(n), 2):
    Q = np.linalg.qr(V[:, list(cols)])[0]
    v = phi_iter(Q.T)[0]
    best.append((abs(v - 14.6398), 'evec%s' % str(cols), v))
Tm, Z = schur(A, output='real')
for cols in itertools.combinations(range(n), 2):
    Q = np.linalg.qr(Z[:, list(cols)])[0]
    v = phi_iter(Q.T)[0]
    best.append((abs(v - 14.6398), 'schurZ%s' % str(cols), v))
for cols in itertools.combinations(range(n), 2):
    Q = np.linalg.qr(np.eye(n)[:, list(cols)])[0]
    v = phi_iter(Q.T)[0]
    best.append((abs(v - 14.6398), 'std%s' % str(cols), v))
best.sort()
for e, tag, v in best[:8]:
    print('  %-12s $\Phi$=%10.4f   |差|=%.4f%s' % (tag, v, e, '   <== 命中' if e < 0.05 else ''))
print('  判读：若 15 个自然候选都对不上，他那一格就需要把 $F$ 矩阵本身印出来（否则 46.5% 这个口径锚点不可复核）。')

# ---- [4] 495× 出处 ----
print('\n[4] 他 §10.4 的 $495\times$')
print('  $78.6010/0.158701=%.1f$（r=2 的"A 不稳模"÷ r=3 的最优 ⇒ 串列）；'
      'r=3 用不稳两模+补一维：%s' % (78.6010 / 0.158701))
Q3 = np.linalg.qr(np.hstack([V[:, UN], rng.standard_normal((n, 1))]))[0]
print('  $\Phi(r=3,\,\text{span}(E_u)+\text{随机补维})=%.4f$ = 最优的 %.1f 倍' %
      (phi_iter(Q3.T)[0], phi_iter(Q3.T)[0] / 0.158701))
