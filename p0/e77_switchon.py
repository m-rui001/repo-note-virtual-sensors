r"""E77 = 注水转置率（switch-on rate）的闭式表，给 note 用的数字。

E76 的结论已经摆在那儿：冻结 Gamma 上的逆注水给出 r_i = min(1, nu/gamma_i)，
方向 i 被传感 <=> gamma_i > nu。令前 k 个方向开（gamma 降序），det R = delta 给出
nu^k = delta * prod_{i<=k} gamma_i；第 k+1 个方向在 nu = gamma_{k+1} 时打开，于是
    I*(k) = 1/2 log2( prod_{i<=k} gamma_i / gamma_{k+1}^k ).
这条 I* 是**可证伪的预言**：E63f/E69 实测的前沿胜者秩应当是 k = #{gamma_i > nu(I)}。

判据（先写下）：
 [1] 逐平面报 gamma(X_w) 与全部 I*(k)，并与实测胜者秩对照：
     red  rho=0 @ I=2,2.5,3；rho=0.01 @ 4；rho=0.1 @ 4.90  => 秩 1,1,1,2,2
     free rank(>=1e-6)=1 @ I=2；"=3" 但形状 [0.927 0.073 0 0] @ I=3、[0.906 0.094 0 0] @ 4
          => 有效秩（按 1% 阈值）1 @ 2、2 @ 3、2 @ 4。
     命中/错位都要报，不许只报命中的。
 [2] 与 E76[D] 的代价差距放一起：错位的那一格（free I=2）是不是代价差最大的一格。
     若"预言错 <=> 代价差大"同向，则这套冻结谱的说法是自洽的；否则我说它"解释了前沿"就是over-claim。
"""
import sys
sys.stdout.reconfigure(encoding='utf-8')
import numpy as np
from scipy.linalg import schur

_src = open('p0/exp_c_audit.py', encoding='utf-8').read().split("print(r'== E53")[0]
_ns = {'__name__': 'p'}
exec(compile(_src, 'p0/exp_c_audit.py[preamble]', 'exec'), _ns)
A, W, TH, JC, n, sym = _ns['A'], _ns['W'], _ns['TH'], _ns['JC'], _ns['n'], _ns['sym']
Ts, U = schur(A, output='real', sort=lambda a: abs(a) < 1.0)[:2]
PLANES = {'free': U[:, :], 'blue': U[:, n - 3:], 'red': U[:, n - 2:]}


def wall_prior(Z, maxit=60000):
    X = W.copy()
    for _ in range(maxit):
        G = sym(Z.T @ X @ Z)
        P = sym(X - X @ Z @ np.linalg.solve(G, Z.T @ X))
        Xn = sym(A @ P @ A.T + W)
        if np.max(np.abs(Xn - X)) < 1e-13 * max(1.0, np.max(np.abs(Xn))):
            return Xn
        X = Xn
    return X


def gam_of(Z, X):
    G = sym(Z.T @ X @ Z)
    w, V = np.linalg.eigh(G)
    Gi2 = (V / np.sqrt(np.clip(w, 1e-300, None))) @ V.T
    return sym(Gi2 @ sym(Z.T @ X @ TH @ X @ Z) @ Gi2)


MEAS = {'red': [(2.0, 1), (2.5, 1), (3.0, 1), (4.0, 2), (4.90, 2)],
        'free': [(2.0, 1), (3.0, 2), (4.0, 2)]}
# free 的实测秩有两个口径：E63f 打印的 rank(>=1e-6)=1,3,3，以及按 1% 权重的有效秩 1,2,2
# （形状 [1 0 0 0] / [.927 .073 0 0] / [.906 .094 0 0]，第三个方向权重 <0.001）。
STRICT = {2.0: 1, 3.0: 3, 4.0: 3}
GAP = {('red', 2.0): 0.62, ('red', 2.5): 0.86, ('red', 3.0): 0.96, ('red', 4.0): 0.86,
       ('red', 4.9): 0.32, ('free', 2.0): 3.90, ('free', 3.0): 1.11, ('free', 4.0): 0.47,
       ('blue', 2.0): 1.36, ('blue', 3.0): 0.33, ('blue', 4.0): 0.16}

print('== E77：逆注水的开关率闭式表 ==')
print('\n[1] gamma(X_w)、闭式开关率 I*(k)、实测胜者秩')
pred = {}
for tag, Z in PLANES.items():
    X_ = wall_prior(Z)
    g = np.sort(np.linalg.eigvalsh(gam_of(Z, X_)))[::-1]
    r = len(g)
    print('  %-5s r=%d  gamma(X_w) = %s' % (tag, r,
          ' '.join('%.4f' % x for x in g)))
    Is = []
    for k in range(1, r):
        Ik = 0.5 * np.log2(np.prod(g[:k]) / g[k] ** k)
        Is.append(Ik)
        print('        k=%d -> k+1:  I* = %7.4f   (秩 %d 传感在 I > %.4f 才划算)'
              % (k, Ik, k + 1, Ik))
    pred[tag] = Is
    for It, mr in MEAS.get(tag, []):
        pk = 1 + sum(1 for x in Is if It > x)
        extra = '' if tag != 'free' else '（严格口径 rank=%d，第三方向权重<0.1%%）' % STRICT[It]
        print('     I=%4.2f: 预言秩 %d | 实测有效秩 %d | %s | 代价差 %+.2f%% %s'
              % (It, pk, mr, '命中' if pk == mr else '**错位**', GAP[(tag, It)], extra))
print('  blue 的实测胜者秩我这轮没有可信读数（E63f 只给了 Powell 的 D 值，没存形状），'
      '所以 blue 只报预言、不报对照：')
print('        blue 预言：%s' % ', '.join('I>%.4f 时秩>=%d' % (x, k + 2) for k, x in enumerate(pred['blue'])))

print('\n[2] 错位与代价差是否同向')
rows = []
for tag in ('red', 'free'):
    for It, mr in MEAS[tag]:
        pk = 1 + sum(1 for x in pred[tag] if It > x)
        rows.append((tag, It, pk == mr, GAP[(tag, It)]))
bad = [r for r in rows if not r[2]]
print('  错位格：%s' % (', '.join('%s I=%.2f (差 %+.2f%%)' % (t, i, g) for t, i, _, g in bad) or '无'))
print('  命中格最大代价差 %+.2f%% | 错位格 %s' % (
    max(g for r in rows if r[2] for g in [r[3]]),
    '，'.join('%+.2f%%' % r[3] for r in bad) or '—'))
print('  => 若错位格恰是代价差最大的格，则"冻结谱解释前沿"是自洽的；否则要降级。')
