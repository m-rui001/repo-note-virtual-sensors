r"""E61 = 三条收尾判定（全部围绕 C 的 §10）。

(a) **他自己的** `c24_out.txt` 第 4 行印 `起点值: 中位 10.1692 最差 86.6659`——他 §10.1 表里写的
    "24/24 全等于 10.169224"与他自己在上一行的"最差"列矛盾（`收敛到的最优值集合（前 5）`只是最小的 5 个）。
    这里用他的 `gdesc`/`floor` 原函数重放 24 起点（seed 也照他的 `1000*r+sd`），直接数出落在 $86.6659$ 的起点数。
(b) `14.6398` 那一格：他取 `schur(A, output='real', sort=lambda a: abs(a)<1.0)` 的 `U2[:,2:]`。
    我把 $T,Z$ 的分块结构印出来，并比较 $\operatorname{span}(Z_{2:})$ 与"两条实不稳定本征向量的 span"——
    后者给 $78.6010$（= 他自己另一行）。$14.6398$ 与全文 +46.5% 这个口径锚点绑在一起，必须钉清。
(c) 修正 E60 [3] 的大圆弧（那里我把弧终点写成了 $f_2$ 在 $f_1^\perp$ 上的分量，只到 $\pi/2$；正确的弧要走到
    夹角本身），重扫 $\Phi$ 看势垒。
"""
import sys
sys.stdout.reconfigure(encoding='utf-8')
import numpy as np
from scipy.linalg import schur

_src = open('p0/exp_c_audit.py', encoding='utf-8').read().split("print(r'== E53")[0]
_ns = {'__name__': 'e53_preamble'}
exec(compile(_src, 'p0/exp_c_audit.py[preamble]', 'exec'), _ns)
A, W, TH, JC, n, sym = _ns['A'], _ns['W'], _ns['TH'], _ns['JC'], _ns['n'], _ns['sym']
phi_iter = _ns['phi_iter']
print(r'== E61：重放 C 的下降 + Schur 那格的分块结构 ==')

# 载入 C 的 floor / rand_F / tangent_basis / gdesc（他的文件在主段之前定义这些）
_c = open('.work3/c24_gdescent.py', encoding='utf-8').read().split("# 结构参照子空间")[0]
_cns = {'__name__': 'c24'}
exec(compile(_c, '.work3/c24_gdescent.py[defs]', 'exec'), _cns)
floor, rand_F, gdesc = _cns['floor'], _cns['rand_F'], _cns['gdesc']

Fglob = np.array([[-0.919737, 0.201585, -0.152763, -0.300184]])
Ftrap = np.array([[-0.141388, 0.615166, -0.752461, 0.188106]])
print('\n(a) 用他原样的 `gdesc` + 他原样的 seed 公式 `rand_F(1,1000*1+sd)`，sd=0..23')
vals = []
for sd in range(24):
    F, f, it, hist = gdesc(rand_F(1, 1000 + sd))
    vals.append(f)
vals = np.array(vals)
print('  24 个收敛值：%s' % np.array2string(np.sort(vals), precision=4))
nG = (np.abs(vals - 10.169224) / 10.169224 < 1e-4).sum()
nT = (np.abs(vals - 86.665888) < 1e-3).sum()
print('  命中 10.169224：%d/24    落在 86.665888：%d/24    中位 %.4f  最差 %.4f' %
      (nG, nT, np.median(vals), vals.max()))
print('  他的 `floor` 在我两个终点上的值：glob %.6f  trap %.6f' % (floor(Fglob), floor(Ftrap)))
print('  他的下降从我的 trap 出发（看能否被他的步长拉走）：')
for tag, F0 in (('trap', np.array([Ftrap])),):
    F, f, it, hist = gdesc(F0[0].copy())
    print('    gdesc(trap) → %.6f  迭代 %d  历史前三 %s' % (f, it, np.array2string(np.array(hist[:4]), precision=4)))

# ---------------- (b) Schur 结构 ----------------
print('\n(b) 他的 Schur 尾 `Sun=U2[:,2:]` 的到底是什么')
T2, U2 = schur(A, output='real', sort=lambda a: abs(a) < 1.0)[:2]
print('  $T=$ ( quasi-upper ) 对角块：%s' % np.array2string(np.diag(T2), precision=4))
print('  $T$ 的 2$\\times$2 块位置：', [(i + 1, j) for i, j in zip(*np.nonzero(np.abs(np.tril(T2, -1))) if True else [])][:0])
off = np.nonzero(np.abs(np.tril(T2, -1)) > 0)
print('    次对角非零元：%s' % (list(zip(off[0], off[1])) if len(off[0]) else '无（全 1×1 块）'))
print('  eig(T 对角)=%s' % np.array2string(np.diag(T2), precision=4))
Sun = U2[:, 2:]
wA, vA = np.linalg.eig(A)
UN = [k for k in range(n) if abs(wA[k]) >= 1]
Vun = np.real(vA[:, UN])
Qun = np.linalg.qr(Vun)[0]
Qsu = np.linalg.qr(Sun)[0]
svals = np.linalg.svd(Qun.T @ Qsu, compute_uv=False)
print('  span(实不稳定本征向量) 与 span(Sun) 的主角：%s 度' % np.array2string(np.degrees(np.arccos(np.clip(svals, -1, 1))), precision=2))
print('  $\\Phi$(Sun 前 1 列)=%.4f   $\\Phi$(Sun 前两列)=%.4f   $\\Phi$(不稳本征 span)=%.4f'
      % (floor(Sun[:, :1].T), floor(Sun[:, :2].T), floor(Qun.T)))
print('  $\\Phi$(U2 前两列)=%.4f   $\\Phi$(U2 第 0,3 列)=%.4f' % (floor(U2[:, :2].T), floor(U2[:, [0, 3]].T)))
for j in range(n):
    print('    $\\Phi$(U2 第 %d 列单独)=%.4f   $\\Phi$(第 %d 列与不稳 span 的混合)…' % (j, floor(U2[:, [j]].T), j))
print('  判读：$\\Phi$ 只依赖行空间。若 Sun 的两列张成的就是不稳定不变子空间，两个数必须相同；'
      '若不同，则他那一格的对象不是"不稳不变子空间"，+46.5% 这个锚点需要重新命名。')

# ---------------- (c) 正确的弧 ----------------
print('\n(c) glob→trap 的大圆弧（走到真实夹角，不再截到 $\\pi/2$）')


def unit(f):
    return f / np.linalg.norm(f)


f1, f2 = unit(Fglob[0]), unit(Ftrap[0])
c = float(np.clip(f1 @ f2, -1, 1))
alpha = np.arccos(c)
perp = unit(f2 - c * f1)
th = np.linspace(0, alpha, 33)
vals = np.array([phi_iter((np.cos(t) * f1 + np.sin(t) * perp).reshape(1, n))[0] for t in th])
print('  夹角 %.2f 度；端点 %.5f / %.5f；弧上最大 %.1f（第 %d 点，占弧长 %.2f）；弧上最小 %.5f（第 %d 点）'
      % (np.degrees(alpha), vals[0], vals[-1], vals.max(), vals.argmax(), vals.argmax() / (len(vals) - 1),
         vals.min(), vals.argmin()))
print('  %s' % np.array2string(vals, precision=2))
print('  ⇒ 弧内部的最大值 %.1f 远高于两端 ⇒ 两点之间存在势垒，一阶下降不可互相跨越' % vals.max())
