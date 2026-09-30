r"""E56 = 对 C §7.1"崎岖度"表的独立复算（同一株 $n=4$ 植物，不同实现），并测他新给的两个百分位读数。

C 在 §7.1 交回一张分位数表（0/0.1/1/10/50%），并给出两条**可推翻**的读数：
 (i) 原文的 Schur 尾 2 通道（$\Phi=14.6398$，$+46.5\%$）≈ 随机分布的**第 11 百分位**；Schur 尾 3 通道
     （$5.1199$，$+16.3\%$）≈ **第 40 百分位**；
 (ii) $r=2$ 的随机中位是 $+236\%$、$r=1$ 是 $+1131\%$，跨 5 个数量级。
我只用 200 个样本（E53[3]）撞过 min/median，分位数分辨率不够。这一跑把 $N$ 提到 4000，直接检验 (i)。

按 §7.3（C 自己补的方法学教训）批量迭代**每步对称化**，并且先在 8 个方向上与标量路线逐位核对再放量。
"""
import sys
sys.stdout.reconfigure(encoding='utf-8')
import numpy as np
import time

_src = open('p0/exp_c_audit.py', encoding='utf-8').read().split("print(r'== E53")[0]
_ns = {'__name__': 'e53_preamble'}
exec(compile(_src, 'p0/exp_c_audit.py[preamble]', 'exec'), _ns)
A, B, W, TH, JC, n, sym = _ns['A'], _ns['B'], _ns['W'], _ns['TH'], _ns['JC'], _ns['n'], _ns['sym']
post_scalar, phi_iter = _ns['post'], _ns['phi_iter']
phi_dare = _ns['phi_dare']
bsym = lambda M: 0.5 * (M + M.transpose(0, 2, 1))

CAP = 1e14


def batch_phi(Fs, it=1500):
    """Fs: (N,r,n) 正交基 → (Phi, 末步号)。每步 sym（C §7.3 的教训）。"""
    N = Fs.shape[0]
    Pt = np.repeat(W[None], N, 0)
    phi = np.full(N, np.inf)
    with np.errstate(over='ignore', invalid='ignore'):
      for k in range(it):
          FPtFT = np.einsum('nia,nij,njb->nab', Fs, Pt, Fs)
          PFt = np.einsum('nij,njr->nir', Pt, Fs)                       # (N,n,r)
          X = np.linalg.solve(FPtFT, np.transpose(PFt, (0, 2, 1)))      # (N,r,n)
          P = Pt - np.einsum('nir,nrj->nij', PFt, X)
          P = bsym(P)
          phi = np.einsum('ij,nji->n', TH, P)
          Pt = bsym(bsym(A @ P @ A.T) + W)
    return np.where(np.isnan(phi), np.inf, phi), k


def haar(N, r, seed):
    g = np.random.default_rng(seed).standard_normal((N, n, r))
    return np.stack([np.linalg.qr(g[i])[0] for i in range(N)])


print('== E56：Grassmannian 上 Φ 的分布（独立实现，检验 C §7.1 的百分位读数）==')
print('  $j_c=\\operatorname{tr}(WP_c)=%.4f$   Schur 尾：$r$=2 → 14.6398，$r$=3 → 5.1199' % JC)

# ---------------------------------------------------------------- 0) 批量 vs 标量逐位核对
print('\n  [0] 批量（每步 sym）vs 两条标量路线：同一批 8 个方向')
F8 = haar(8, 2, 7)
pb, _ = batch_phi(F8, it=1500)
print('      %-4s %13s %13s %13s %10s %10s' % ('#', '批量', '标量定点', '标量 DARE', '批/定', '批/DARE'))
for i in range(8):
    ps, k, s, fl = phi_iter(F8[i].T)
    pd = phi_dare(F8[i].T)[0]
    print('    %-4d %13.6f %13.6f %13.6f %10.2e %10.2e'
          % (i, pb[i], ps, pd, abs(pb[i] - ps) / max(abs(ps), 1e-30),
             abs(pb[i] - pd) / max(abs(pd), 1e-30)))

# ---------------------------------------------------------------- 1) 分位数表
print('\n  [1] $N=4000$ 个 Haar 子空间的分位数（相对 $j_c$ 的 %），对照 C §7.1（$3\\times10^5$ 点）')
CS = {1: [32.37, 38.05, 68.03, 303.6, 1131.0], 2: [2.40, 5.46, 13.37, 38.74, 236.4],
      3: [0.52, 0.93, 1.56, 3.76, 31.27]}
SCHUR = {2: 14.639795, 3: 5.119861}
QS = [0, 0.1, 1, 10, 50]
RES = {}
for r in (1, 2, 3):
    t0 = time.time()
    ph, nk = batch_phi(haar(4000, r, 11 + r), it=1500)
    ph = np.sort(ph)
    fin = ph[ph < CAP]
    RES[r] = ph
    mine = [float(np.percentile(fin, q) / JC * 100 - 100) for q in QS]
    print('    r=%d  有限 %d/4000  %.1fs' % (r, len(fin), time.time() - t0))
    print('      %-8s %10s %10s %10s %10s %10s' % ('分位', 'min', '0.1%', '1%', '10%', '50%'))
    print('      %-8s %10.2f %10.2f %10.2f %10.2f %10.2f' % ('E56', *mine))
    print('      %-8s %10.2f %10.2f %10.2f %10.2f %10.2f' % ('C §7.1', *CS[r]))
    print('      比值     %10.3f %10.3f %10.3f %10.3f %10.3f' %
          tuple(m / c for m, c in zip(mine, CS[r])))

# ---------------------------------------------------------------- 2) Schur 尾落在第几百分位
print('\n  [2] 原文 Schur 尾规则的百分位位置（C 给的是 11% / 40%）')
for r in (2, 3):
    fin = RES[r][RES[r] < CAP]
    v = SCHUR[r]
    pct = 100.0 * (fin < v).sum() / len(fin)
    se = 100 * np.sqrt(max(pct / 100, 1e-9) * (1 - pct / 100) / len(fin))
    band = [100.0 * (fin < v * f).sum() / len(fin) for f in (0.9, 1.1)]
    print('    r=%d  $\\Phi$=%.4f → 第 %.1f 百分位  二项 SE %.2f 点  $\\pm10\\%%$ 带宽 %.1f–%.1f  C 报 %s%%'
          % (r, v, pct, se, band[0], band[1], 11 if r == 2 else 40))

# ---------------------------------------------------------------- 3) 长尾质量
print('\n  [3] 长尾：$\\Phi$ 超过阈值的随机比例（期望是否存在）')
for r in (1, 2, 3):
    fin = RES[r][RES[r] < CAP]
    print('    r=%d  ' % r + '   '.join('$\\Phi>%g$:%.1f%%' % (t, 100.0 * (fin > t).sum() / len(fin))
                                        for t in (100, 1000, 10000)) +
          '   最大 %.3e   均值 %.3e' % (fin.max(), fin.mean()))

print('\n  判读（写死，跑出来再核）')
print("""    (1) [2] 若我的百分位与 C 的 11%/40% 相差 >2 个二倍标准误，则 §7.1 的位置读数要重新标定；
        落在误差内 → 独立实现坐实，"同样两条通道选哪两条值 44 个百分点"可以进正文（附分布图）。
    (2) [1] 若我的分位数列系统性低于 C 的（同一 $q$），先看 [0] 是否逐位一致：一致则差值是样本量
        （我 4000 vs 他 300000）在 $q\\le1\\%$ 处的抽样误差，不是实现差异。
    (3) [3] 若 $r=1$ 的长尾质量集中在 $10^3$ 以上，"跨 5 个数量级"就是重尾（均值被尾部主导）的性质，
        任何"随机方向差不多"的写法都会被一句"$\\mathbb E\\Phi$ 被尾部主导"杀掉；对 §2.4 措辞是硬约束。""")
