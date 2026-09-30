r"""E67 = 冻结前的最后一次定点兜底：只做两件事，都是"已经写进板子的标题句"的验证，不加小数。

[A] §15.4 第 3 条：C 的对偶分支 S=c*Theta，他 c33 的 10^0.5 网格线性插值 vs 二分命中（7 个率，~150 次迭代）
[B] §15.3 标题句"等地板不等率"：120 个随机方向，全程只用信息型定点（无 sdare、无我的坏判据），
    地板=phi_iter，率=J_hit(I=2/3)。这是把 E64 那条结论从"被我公开标记为不可靠的路线"
    搬到"两条路线都对到的路线"上。不做分位数、不做更细容差、不再加密。
"""
import sys, time
sys.stdout.reconfigure(encoding='utf-8')
import numpy as np

np.set_printoptions(precision=4, suppress=True, linewidth=170)
_src = open('p0/exp_c_audit.py', encoding='utf-8').read().split("print(r'== E53")[0]
_ns = {'__name__': 'p'}
exec(compile(_src, 'p0/exp_c_audit.py[preamble]', 'exec'), _ns)
A, W, TH, JC, n, sym = _ns['A'], _ns['W'], _ns['TH'], _ns['JC'], _ns['n'], _ns['sym']
phi_iter = _ns['phi_iter']
print('== E67 ==')


def JI(S, tol=1e-11, maxit=6000):
    Pt = W.copy(); Pm = None
    try:
        for _ in range(maxit):
            Pm = sym(np.linalg.solve(np.linalg.inv(Pt) + S, np.eye(n)))
            Pn = sym(A @ Pm @ A.T + W)
            if not np.all(np.isfinite(Pn)) or np.max(np.abs(Pn)) > 1e14:
                return np.inf, np.inf
            if np.max(np.abs(Pn - Pt)) < tol * max(1.0, np.max(np.abs(Pn))):
                Pt = Pn; break
            Pt = Pn
        else:
            return np.inf, np.inf
    except np.linalg.LinAlgError:
        return np.inf, np.inf
    return (float(np.trace(TH @ Pm)) + JC,
            float(0.5 * (np.linalg.slogdet(Pt)[1] - np.linalg.slogdet(Pm)[1]) / np.log(2)))


def J_hit(S1, It, lo=-4.0, hi=13.0, steps=20):
    if not np.isfinite(JI(10.0 ** hi * S1)[0]):
        return np.nan
    for _ in range(steps):
        m = 0.5 * (lo + hi)
        j, i_ = JI(10.0 ** m * S1)
        if np.isfinite(i_) and i_ >= It:
            hi = m
        else:
            lo = m
    return JI(10.0 ** hi * S1)[0]


print('\n[A] 对偶分支：插值 vs 二分')
his = {2.0: 70.9922, 2.5: 54.0421, 3.0: 45.3757, 3.5: 40.9311, 4.0: 37.5685, 4.5: 36.2756, 5.0: 34.9827}
t0 = time.time()
for It, jh in his.items():
    jb = J_hit(TH / np.trace(TH), It)
    print('  I=%.2f: c33 插值 %9.4f | 二分 %9.4f | 差 %+8.4f（插值偏高 %+.2f%%）'
          % (It, jh, jb, jb - jh, 100 * (jh - jb) / jb))
print('  (%.0f s)' % (time.time() - t0))

print('\n[B] 等地板配对（纯定点，120 方向）')
rng = np.random.default_rng(5)
Zs = rng.normal(0, 1, (120, n)); Zs /= np.linalg.norm(Zs, axis=1, keepdims=True)
t0 = time.time()
FL, J2, J3 = [], [], []
for z in Zs:
    f, k, s, tag = phi_iter(z[None, :])
    if tag != 'conv' or not np.isfinite(f):
        continue
    a, b = J_hit(np.outer(z, z), 2.0), J_hit(np.outer(z, z), 3.0)
    if not (np.isfinite(a) and np.isfinite(b)):
        continue
    FL.append(f); J2.append(a); J3.append(b)
FL = np.array(FL); J2 = np.array(J2); J3 = np.array(J3)
print('  可用 %d/120（%.0f s），地板 [%.3f, %.3f]，J(I=2) [%.3f, %.3f]'
      % (len(FL), time.time() - t0, FL.min(), FL.max(), J2.min(), J2.max()))
for tol in (0.005, 0.001):
    pairs = []
    for i in range(len(FL)):
        d = np.abs(FL - FL[i]) / np.maximum(FL, FL[i])
        for k in np.where((d < tol) & (np.arange(len(FL)) > i))[0]:
            pairs.append((abs(J2[i] - J2[k]) / min(J2[i], J2[k]), i, k, d[k]))
    if not pairs:
        print('  地板相对差 <%.4f: 无配对' % tol); continue
    pairs.sort(reverse=True)
    sp, i, k, dd = pairs[0]
    print('  地板相对差 <%.4f: %d 对；率(I=2)相对极差 中位 %.1f%% / 最大 %.1f%%'
          % (tol, len(pairs), 100 * np.median([p[0] for p in pairs]), 100 * sp))
    print('    最大一对：地板 %.3f vs %.3f（差 %.4f%%）→ J(I=2) %.3f vs %.3f；同配对 J(I=3) %.3f vs %.3f'
          % (FL[i], FL[k], 100 * dd, J2[i], J2[k], J3[i], J3[k]))
for nm, JV in (('I=2', J2), ('I=3', J3)):
    cross = 0; ex = None; worst = 0
    for i in range(len(FL)):
        w = np.where((FL > FL[i]) & (JV < JV[i]))[0]
        cross += len(w)
        if len(w):
            k = w[np.argmax(FL[w] - FL[i])]
            if FL[k] - FL[i] > worst:
                worst, ex = FL[k] - FL[i], (i, k)
    if ex is None:
        print('  %s: 无交叉对' % nm)
    else:
        print('  %s: 交叉对 %d 个；最极端：地板高 %.3f（%.1f%%）却省率 %.3f -> %.3f（%.1f%%）'
              % (nm, cross, worst, 100 * worst / FL[ex[0]], JV[ex[0]], JV[ex[1]],
                 100 * (JV[ex[1]] - JV[ex[0]]) / JV[ex[0]]))
q = np.percentile(FL, np.arange(11) * 10)
print('  分桶（地板十等分，只报桶内 J(I=2) 极差的相对量级）：')
for b in range(10):
    m = (FL >= q[b]) & (FL <= q[b + 1])
    if m.sum() >= 3:
        print('    桶 %d 地板[%.2f,%.2f] n=%d: J2 极差 %.1f%%' % (b, q[b], q[b + 1], m.sum(),
                                                                  100 * (J2[m].max() - J2[m].min()) / J2[m].min()))
print('== E67 done ==')
