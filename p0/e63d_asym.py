r"""E63d = 射线口径下成本逼近"墙"有多慢：判掉"图上左端点=墙"这类读图推断。
结论用一句话：要把 D 读到离墙 0.01 以内需要 I>=15 bit，而原文 Fig.1 的 I 轴只画到 5。
跑法：python -X utf8 p0/e63d_asym.py > p0/e63d_out.txt 2>&1
"""
import sys; sys.stdout.reconfigure(encoding='utf-8')
import numpy as np
from scipy.linalg import schur, solve_discrete_are as sdare
_src = open('p0/exp_c_audit.py', encoding='utf-8').read().split("print(r'== E53")[0]
_ns = {'__name__': 'p'}
exec(compile(_src, 'x', 'exec'), _ns)
A, W, TH, JC, n, sym = _ns['A'], _ns['W'], _ns['TH'], _ns['JC'], _ns['n'], _ns['sym']
phi = _ns['phi_iter']
Ts, U = schur(A, output='real', sort=lambda a: abs(a) < 1.0)[:2]
F = U[:, 2:].T; r = 2


def IJ(V):
    Pt = sym(sdare(A.T, F.T, sym(W), sym(V)))
    FPt = F @ Pt
    P = sym(Pt - FPt.T @ np.linalg.solve(FPt @ F.T + V, FPt))
    s, ld = np.linalg.slogdet(np.eye(r) + np.linalg.solve(V.T, (F @ Pt @ F.T)).T)
    return 0.5 * ld / np.log(2), JC + float(np.trace(TH @ P))


def ray(I_t):
    lo, hi = -8.0, 20.0
    for _ in range(70):
        m = 0.5 * (lo + hi)
        a = IJ(np.eye(r) / 10.0 ** m)
        if a is None or a[0] < I_t:
            lo = m
        else:
            hi = m
    return IJ(np.eye(r) / 10.0 ** hi)[1]


wall = JC + phi(F)[0]
print('E63d: 原文 F_2（Schur 尾 2 列）射线口径逼近墙的过程   墙 = %.4f' % wall)
for I in [4.926, 5.0, 6, 8, 10, 15, 20, 40, 200]:
    d = ray(I)
    print('   I=%8.3g bit -> D = %.4f  (超墙 %+.4f, 相对 %+.2f%%)' % (I, d, d - wall, 100 * (d - wall) / wall))
print('要读到离墙 0.01 以内需 I >= 15 bit；Fig.1 的 I 轴上限 5 ⇒ 图上的左端点不可能读出一个不可行墙。')
