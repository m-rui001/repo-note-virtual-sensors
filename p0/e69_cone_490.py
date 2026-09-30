r"""E69 = 把 Table I 的 cone 列补满：E63f[A] 的穷尽网格在 I=4.90 上没跑过。

判据先写下：
 [0] 先原样复现 E63f[A] 的四个已知锥值（2.0/2.5/3.0/4.0）。任一格差 >5e-4 就说明这个脚本
     和那条证据链不是同一口径，4.90 的数作废，不进口径。
 [1] 复现通过后，同一网格多跑 I=4.90 一格，并用信息型定点 JI 独立复核胜者设计
     （JIfast 走 sdare + 闭环谱半径判据，是 debt #6 里我公开过的旧 guard；
      复核是为了确认胜者不是那条 guard 的漏检/误检点）。
 [2] 只报符号与量级，不加小数位：4.90 这一格只用来判定"锥值在图线上方还是下方"。
"""
import sys, time
sys.stdout.reconfigure(encoding='utf-8')
import numpy as np
from scipy.linalg import schur, solve_discrete_are as sdare

np.set_printoptions(precision=4, suppress=True, linewidth=170)
_src = open('p0/exp_c_audit.py', encoding='utf-8').read().split("print(r'== E53")[0]
_ns = {'__name__': 'p'}
exec(compile(_src, 'p0/exp_c_audit.py[preamble]', 'exec'), _ns)
A, W, TH, JC, n, sym = _ns['A'], _ns['W'], _ns['TH'], _ns['JC'], _ns['n'], _ns['sym']
In_ = np.eye(n)
EPSR = 1e-11
print('== E69：Table I 锥列补格（I=4.90）==')


def JI(S, tol=1e-12, maxit=20000):
    """信息型定点，对任意半正定 S 安全；不收敛=inf。与 E63f 同一实现。"""
    Pt = W.copy()
    for _ in range(maxit):
        Pm = sym(np.linalg.inv(np.linalg.inv(Pt) + S))
        Pn = sym(A @ Pm @ A.T + W)
        if not np.all(np.isfinite(Pn)) or np.max(np.abs(Pn)) > 1e14:
            return np.inf, np.inf
        if np.max(np.abs(Pn - Pt)) < tol * max(1.0, np.max(np.abs(Pn))):
            Pt = Pn; break
        Pt = Pn
    else:
        return np.inf, np.inf
    return (float(np.trace(TH @ Pm)) + JC,
            float(0.5 * (np.linalg.slogdet(Pt)[1] - np.linalg.slogdet(Pm)[1]) / np.log(2)))


def JIfast(S):
    """S 的秩 k 因式 -> F(k x n), V=I_k -> sdare 拿 P~ -> (J,I)；闭环不稳定则 inf。与 E63f 逐字相同。"""
    ev, EV = np.linalg.eigh(sym(S))
    keep = ev > max(EPSR * ev.max(), 1e-300)
    if not keep.any():
        return np.inf, np.inf
    F = (EV[:, keep] * np.sqrt(ev[keep]))[:, :].T
    k = F.shape[0]
    try:
        Pt = sym(sdare(A.T, F.T, sym(W), np.eye(k)))
    except Exception:
        return np.inf, np.inf
    if not np.all(np.isfinite(Pt)):
        return np.inf, np.inf
    FPt = F @ Pt
    try:
        P = sym(Pt - FPt.T @ np.linalg.solve(FPt @ F.T + np.eye(k), FPt))
    except Exception:
        return np.inf, np.inf
    cl = A - A @ Pt @ F.T @ np.linalg.solve(FPt @ F.T + np.eye(k), F)
    if np.max(np.abs(np.linalg.eigvals(cl))) >= 1.0 - 1e-9:
        return np.inf, np.inf
    ld = np.linalg.slogdet(Pt)[1] - np.linalg.slogdet(P)[1]
    if ld <= 0:
        return np.inf, np.inf
    return JC + float(np.trace(TH @ P)), 0.5 * ld / np.log(2)


Ts, U = schur(A, output='real', sort=lambda a: abs(a) < 1.0)[:2]
Z2 = U[:, n - 2:]


def J_hit(S1, It, lo=-3.0, hi=13.0, steps=26):
    a = JIfast(10.0 ** hi * S1)
    if not np.isfinite(a[1]) or a[1] < It:
        return np.nan
    for _ in range(steps):
        m = 0.5 * (lo + hi)
        b = JIfast(10.0 ** m * S1)
        if np.isfinite(b[1]) and b[1] >= It:
            hi = m
        else:
            lo = m
    return JIfast(10.0 ** hi * S1)[0]


RHO = [0.0, 1e-4, 1e-2, 1e-1, 0.5, 1.0]
TH_ = np.linspace(0, np.pi, 91)[:-1]
IS = (2.0, 2.5, 3.0, 4.0, 4.90)
E63F = {2.0: 62.8975, 2.5: 54.2924, 3.0: 50.6023, 4.0: 47.9396}

best = {It: (np.inf, None) for It in IS}
t0 = time.time()
for rho in RHO:
    for th in TH_:
        Rv = np.array([[np.cos(th), -np.sin(th)], [np.sin(th), np.cos(th)]])
        S1 = sym(Z2 @ sym(Rv @ np.diag([1.0, rho]) @ Rv.T) @ Z2.T)
        for It in IS:
            j = J_hit(S1, It)
            if np.isfinite(j) and j < best[It][0]:
                best[It] = (j, (rho, th))

print('\n[0] 复现 E63f[A] 的四个已知锥值')
ok = True
for It in (2.0, 2.5, 3.0, 4.0):
    j = best[It][0]
    d = j - E63F[It]
    print('  I=%.2f: E69 %9.4f | E63f %9.4f | 差 %+.5f %s'
          % (It, j, E63F[It], d, 'OK' if abs(d) < 5e-4 else '**不一致**'))
    ok &= abs(d) < 5e-4
print('  => %s（耗时 %.0f s）' % ('口径一致，[1] 可信' if ok else '口径不一致，4.90 不进口径', time.time() - t0))

print('\n[1] 新增的一格 I=4.90')
j, arg = best[4.90]
rho, th = arg
Rv = np.array([[np.cos(th), -np.sin(th)], [np.sin(th), np.cos(th)]])
S1 = sym(Z2 @ sym(Rv @ np.diag([1.0, rho]) @ Rv.T) @ Z2.T)
print('  胜者设计: rho=%.3g, theta=%.1f deg（rho=0 即关掉一条通道，rank-1 边界）'
      % (rho, np.degrees(th)))

# 用信息型定点 JI 独立复核：重新二分尺度命中 I=4.90，比较两条路线的 J
lo, hi = -3.0, 13.0
for _ in range(40):
    m = 0.5 * (lo + hi)
    a = JI((10.0 ** m) * S1)
    if np.isfinite(a[1]) and a[1] >= 4.90:
        hi = m
    else:
        lo = m
Jv, Iv = JI((10.0 ** hi) * S1)
print('  JI   复核: J=%9.4f I=%7.4f' % (Jv, Iv))
Jf, If_ = JIfast((10.0 ** hi) * S1)
print('  JIfast 同点: J=%9.4f I=%7.4f | dJ=%+.5f dI=%+.5f' % (Jf, If_, Jv - Jf, Iv - If_))
print('  穷尽网格值 %9.4f | JI 复核值 %9.4f | 差 %+.5f' % (j, Jv, Jv - j))

print('\n[2] 与规范化读图值（repro/readoff.py）对表')
DIG = {2.00: (62.902, 0.021), 2.50: (54.311, 0.014), 3.00: (50.651, 0.003),
       4.00: (47.978, 0.007), 4.90: (47.040, 0.006)}
for It in IS:
    Dd, sp = DIG[It]
    d = best[It][0] - Dd
    print('  I=%.2f: 锥 %9.4f | 图 %8.3f+-%.3f | 锥-图 %+.4f | 读图半宽 %.3f | %s'
          % (It, best[It][0], Dd, sp, d, sp,
             '符号可读' if abs(d) > sp else '符号不可读'))
