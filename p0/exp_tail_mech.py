r"""E56c = 尾的机制：$r=1$ 时 $\Phi$ 到底是"到坏集的距离"的什么函数？

E56b 量到 $r=1$ 的 $\Pr(\Phi>t)\sim t^{-0.51}$、$r=2$ 的 $\sim t^{-0.98}$，而我自己的单极点
$\Phi\sim\kappa/(\varepsilon^*-\varepsilon)$（E55，扰动 $A$）换成"随机方向"就该给出：坏集
$\{f:f\perp v_u\}$ 在 $S^3$ 上是**余维 1**（$\ker F$ 含一条不稳定本征向量即破坏可检测性——我先前写余维 3
是错的，那要求整个超平面不变），$\Pr(d<s)\sim s$，于是尾指数应为 $-1/q$，$q$ 是 $\Phi\sim d^{-q}$ 的幂。
$-0.51$ 反解出 $q\approx1.95$，$-0.98$ 反解出 $q\approx1.02$。**同一株植物上 $r=1$ 像双极点、$r=2$ 像单极点。**
这一跑直接量 $q$，不靠反解。
"""
import sys
sys.stdout.reconfigure(encoding='utf-8')
import numpy as np

_src = open('p0/exp_c_audit.py', encoding='utf-8').read().split("print(r'== E53")[0]
_ns = {'__name__': 'e53_preamble'}
exec(compile(_src, 'p0/exp_c_audit.py[preamble]', 'exec'), _ns)
A, W, TH, JC, n, sym = _ns['A'], _ns['W'], _ns['TH'], _ns['JC'], _ns['n'], _ns['sym']
phi_iter = _ns['phi_iter']
bsym = lambda M: 0.5 * (M + M.transpose(0, 2, 1))
CAP = 1e14
print(r'== E56c：$r=1$ 的 $\Phi$ 对"到坏集距离"的幂 ==')


def _run(Fs, it):
    N = Fs.shape[0]
    Pt = np.repeat(W[None], N, 0)
    for k in range(it):
        FPtFT = np.einsum('nia,nij,njb->nab', Fs, Pt, Fs)
        PFt = np.einsum('nij,njr->nir', Pt, Fs)
        X = np.linalg.solve(FPtFT, np.transpose(PFt, (0, 2, 1)))
        P = bsym(Pt - np.einsum('nir,nrj->nij', PFt, X))
        Pt = bsym(bsym(A @ P @ A.T) + W)
    return np.einsum('ij,nji->n', TH, P)


def batch_phi(Fs, it=600, chunk=1000):
    out = []
    for s in range(0, Fs.shape[0], chunk):
        try:
            out.append(_run(Fs[s:s + chunk], it))
        except np.linalg.LinAlgError:
            out.append(np.full(min(chunk, Fs.shape[0] - s), np.inf))
    return np.where(np.isnan(np.concatenate(out)), np.inf, np.concatenate(out))


# 坏集 = 与某条不稳定右本征向量正交的超平面
w_ev, V_ev = np.linalg.eig(A)
VU = []
for k, lam in enumerate(w_ev):
    if abs(lam) >= 1.0 and abs(lam.imag) < 1e-12:
        v = np.real(V_ev[:, k])
        VU.append(v / np.linalg.norm(v))
        print('  不稳定实模态 $\\lambda=%.4f$，右本征向量 $v=%s$' % (lam.real, np.array2string(v, precision=4)))
VU = np.stack(VU)                      # (m,n)  m=2
_cpx = [abs(l) for l in w_ev if abs(l) >= 1.0 and l.imag != 0]
if _cpx:
    print('  复对（模 %s）$\\ge1$，不进实坏集' % _cpx)
else:
    print('  无模 $\\ge1$ 的复对本征值：坏集仅由两条实不稳定模态张成（$m=2$ 个超平面）')

rng = np.random.default_rng(20260929)
G = rng.standard_normal((20000, n))
Fm = G / np.linalg.norm(G, axis=1, keepdims=True)
d = np.abs(Fm @ VU.T).min(axis=1)                 # 到最近坏超平面的距离 = |cos|
d1, d2 = np.abs(Fm @ VU.T).T
print('\n  [A] $\\Pr(d<s)$ 的幂（余维 1 ⇒ 斜率 1）')
ss = np.geomspace(0.002, 0.08, 8)
pr = np.array([(d < s).mean() / s for s in ss])   # 除 $s$ 后应为常数
print('    s      ' + ' '.join('%9.4f' % s for s in ss))
print('    ratio  ' + ' '.join('%9.4f' % p for p in pr))
print('    ⇒ 斜率 $\\log\\Pr/\\log s$ = %.3f（相邻点）；余维 1 给 1.000，余维 2 给 2.000' %
      np.polyfit(np.log(ss), np.log([(d < s).mean() for s in ss]), 1)[0])

sel = np.where(d < 0.35)[0]                       # 小距离子样本，标量路线逐个算（批量 solve 在这些方向会奇异）
ph = np.full(Fm.shape[0], np.inf)
flag = np.zeros(Fm.shape[0], int)
for i in sel:
    out = phi_iter(Fm[i][None, :], it=900, cap=1e16)
    ph[i] = out[0]
    flag[i] = {'conv': 1, 'maxit': 2, 'cap': 3}.get(out[3], 0)
ok = np.isfinite(ph) & (d > 0)
print('\n  子样本 %d（$d<0.35$）：conv %d  maxit %d  cap(视为发散) %d' %
      (len(sel), (flag[sel] == 1).sum(), (flag[sel] == 2).sum(), (flag[sel] == 3).sum()))
ld, lp = np.log(d[ok]), np.log(np.maximum(ph[ok], 1e-12))
print('\n  [B] $\\log\\Phi$ 对 $\\log d$ 的幂（分小距离窗口）')
for cut in (0.4, 0.2, 0.1, 0.05):
    m = ld < np.log(cut)
    b = np.polyfit(ld[m], lp[m], 1)
    print('    $d<%.2f$  点数 %5d  斜率 $q$ = %+.3f  截距 %.3f  相关 %.3f' %
          (cut, m.sum(), b[0], b[1], np.corrcoef(ld[m], lp[m])[0, 1]))

print('\n  [C] 尾指数的闭环检查：$\\Pr(\\Phi>t)$ 的斜率应 $=-$(A 的斜率)$/q$')
TS = np.geomspace(1.0, 1e6, 10)
prr = np.array([(ph[ok] > t).mean() for t in TS])
m = (prr > 0.002) & (prr < 0.9)
sl = np.polyfit(np.log(TS[m]), np.log(prr[m]), 1)[0]
q = -np.polyfit(ld, lp, 1)[0]
print('    实测 $\\Pr(\\Phi>t)$ 斜率 %.3f；由 $-1/q$ 预言 %.3f（$q$ 取全窗口 %.3f）' % (sl, -1.0 / q, q))
print('    只用 $d<0.1$ 的小距离窗口：$q$ = %.3f → 预言 %.3f' %
      tuple([-np.polyfit(ld[ld < np.log(.1)], lp[ld < np.log(.1)], 1)[0]] * 2))

print('\n  [D] 两模态分解：$\\Phi\\approx A/d_1^{-2}\\,$项 $+\\,B/d_2^{-2}\\,$项 $+\\,C$？')
for cut in (0.02, 0.05, 0.15):
    msk = ok & (d < cut)
    y = ph[msk]
    Xd = np.stack([d1[msk] ** -2, d2[msk] ** -2, np.ones(msk.sum())], 1)
    coef, *_ = np.linalg.lstsq(Xd, y, rcond=None)
    rms = lambda pred: np.sqrt(((pred - y) ** 2).mean()) / y.mean()
    a1 = np.polyfit(d1[msk] ** -2, y, 1)
    a2 = np.polyfit(d[msk] ** -2, y, 1)
    print('    $d<%.2f$ 样本 %4d：两路 $A=%.3g\\ B=%.3g\\ C=%.3g$  RMS %.3f | 单路 $d_1^{-2}$ RMS %.3f | '
          '$\\min$ 路 $d^{-2}$ RMS %.3f'
          % (cut, msk.sum(), *coef, rms(Xd @ coef), rms(np.polyval(a1, d1[msk] ** -2)),
             rms(np.polyval(a2, d[msk] ** -2))))

print('\n  判读')
print("""    (1) [A] 斜率若为 1.0，坏集余维 1 坐实（我先前"余维 3"作废）。
    (2) [B]/[C] 若 $q\\approx2$ 且 $-1/q$ 与实测尾斜率同号同量级，则 $r=1$ 的重尾**完全**由
        "双极点 $\\times$ 余维 1"给出，$\\mathbb E\\Phi=\\infty$ 是一句话可证的推论。
    (3) [D] 若两模态各自 $d_i^{-2}$ 的线性组合就够（RMS 小），则不存在跨模态耦合项，
        $\\Phi$ 的近坏集主部 = $\\sum_i c_i\\,d_i^{-2}$；若 RMS 大，则还有角度依赖，需保留 $q$ 的窗口性。""")
