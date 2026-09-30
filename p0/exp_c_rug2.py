r"""E56b = E56 的更正 + 两条新检查。E56 的 [1] 把百分比写成了 $\Phi/j_c-1$，而留言板的
"+32.3%" 一直是 $\Phi/j_c$（锚点：$\Phi(F_2^{\rm Schur})=14.6398=+46.5\%$，§2.2/§2.3）。
这一版把标签改对，并加：
 [A] 同种子下迭代步数 $1500$ vs $4000$ 的漂移 —— 排除"我的分位数低是因为没收敛"；
 [B] $N=20000$ 的 $r=1$ 尾指数：$\Pr(\Phi>t)$ 的 log–log 斜率，以及 $\mathbb E\Phi$ 是否存在。

后者是新问题：我在 §32.6 给了 $\Phi\sim\kappa/(\varepsilon^*-\varepsilon)$ 的单极点。若把同一逻辑用于
"随机方向 $F$ 落在不可检测集附近"，坏集在 $r=1$ 的单位法向量球面上是余维 $n-1=3$ 的点集，
$\varepsilon^*(F)\asymp\mathrm{dist}(F,\text{bad})$，则应得 $\Pr(\Phi>t)\sim t^{-3}$。
量出来若是 $t^{-0.5}$，则**单极点解释不了尾**——这是打我自己那条链的解释力，不是打 C 的表。
"""
import sys
sys.stdout.reconfigure(encoding='utf-8')
import numpy as np
import time

_src = open('p0/exp_c_audit.py', encoding='utf-8').read().split("print(r'== E53")[0]
_ns = {'__name__': 'e53_preamble'}
exec(compile(_src, 'p0/exp_c_audit.py[preamble]', 'exec'), _ns)
A, W, TH, JC, n, sym = _ns['A'], _ns['W'], _ns['TH'], _ns['JC'], _ns['n'], _ns['sym']
bsym = lambda M: 0.5 * (M + M.transpose(0, 2, 1))
CAP = 1e14


def batch_phi(Fs, it=1500):
    N = Fs.shape[0]
    Pt = np.repeat(W[None], N, 0)
    phi = np.full(N, np.inf)
    with np.errstate(over='ignore', invalid='ignore'):
        for k in range(it):
            FPtFT = np.einsum('nia,nij,njb->nab', Fs, Pt, Fs)
            PFt = np.einsum('nij,njr->nir', Pt, Fs)
            X = np.linalg.solve(FPtFT, np.transpose(PFt, (0, 2, 1)))
            P = bsym(Pt - np.einsum('nir,nrj->nij', PFt, X))
            phi = np.einsum('ij,nji->n', TH, P)
            Pt = bsym(bsym(A @ P @ A.T) + W)
    return np.where(np.isnan(phi), np.inf, phi)


def haar(N, r, seed):
    g = np.random.default_rng(seed).standard_normal((N, n, r))
    return np.stack([np.linalg.qr(g[i])[0] for i in range(N)])


def pct(ph):
    return ph / JC * 100.0


print(r'== E56b：分位数表（标签更正 $\Phi/j_c$）+ 收敛漂移 + $r=1$ 尾指数 ==')
CS = {1: [32.37, 38.05, 68.03, 303.6, 1131.0], 2: [2.40, 5.46, 13.37, 38.74, 236.4],
      3: [0.52, 0.93, 1.56, 3.76, 31.27]}
SCHUR = {2: 14.639795, 3: 5.119861}
QS = [0, 0.1, 1, 10, 50]
ST = {2: 0.11, 3: 0.40}

print('\n  [A] 同种子 $N=4000$，迭代 1500 vs 4000 步（E56 的表在左边）')
DATA = {}
for r in (1, 2, 3):
    Fs = haar(4000, r, 11 + r)
    p15 = batch_phi(Fs, it=1500)
    p40 = batch_phi(Fs, it=4000)
    DATA[r] = p40
    q15 = [float(np.percentile(pct(p15), q)) for q in QS]
    q40 = [float(np.percentile(pct(p40), q)) for q in QS]
    print('    r=%d  %-8s %10.2f %10.2f %10.2f %10.2f %10.2f' % (r, 'it=1500', *q15))
    print('        %-8s %10.2f %10.2f %10.2f %10.2f %10.2f' % ('it=4000', *q40))
    print('        %-8s %10s %10s %10s %10s %10s' % ('C §7.1', *['%.2f' % v for v in CS[r]]))
    print('        最大相对漂移 %.2e（第 50 百分位 %.2e）'
          % (max(abs(b - a) / max(abs(b), 1e-9) for a, b in zip(q15, q40)),
             abs(q40[-1] - q15[-1]) / q40[-1]))
    if r in SCHUR:
        v = SCHUR[r]
        pc = lambda p: 100.0 * (p < v).sum() / len(p)
        print('        Schur 尾 $\\Phi$=%.4f 百分位：it=1500 → %.1f%%，it=4000 → %.1f%%  C 报 %d%%'
              % (v, pc(p15), pc(p40), ST[r] * 100))

print('\n  [B] $r=1$，$N=20000$ 的右尾 $\\Pr(\\Phi>t)$ 与 log–log 斜率')
Fs = haar(20000, 1, 99)
ph = batch_phi(Fs, it=4000)
ph = ph[np.isfinite(ph)]
TS = np.geomspace(0.5, 2e6, 9)
pr = np.array([(ph > t).mean() for t in TS])
ok = (pr > 0.004) & (pr < 0.95)
xs, ys = np.log(TS[ok]), np.log(pr[ok])
sl = np.polyfit(xs, ys, 1)[0]
print('    %-12s' % 't' + ' '.join('%9.2f' % t for t in TS))
print('    %-12s' % '$\\Pr(\\Phi>t)$' + ' '.join('%9.4f' % p for p in pr))
print('    回归斜率（$\\ln\\Pr$ vs $\\ln t$，取 $0.004<p<0.95$ 的 %d 点）= %.3f' % (ok.sum(), sl))
print('    ⇒ $\\mathbb E[\\Phi]$ 存在需斜率 $<-1$；实测 %.2f → %s' % (sl, '不存在（尾部主导）' if sl > -1 else '存在'))
print('    截断偏差方向：双精度/1500–4000 步会把大 $\\Phi$ 记小，只会让斜率**偏薄**，不会造出厚尾。')
tail = np.sort(ph)[-10:][::-1]
print('    最大 10 个值：%s' % '  '.join('%.3g' % v for v in tail))
print('    样本均值 $\\bar\\Phi=%.3e$ = %.0f $j_c$；但 $\\Pr(\\Phi>10^3)=%.1f\\%%$ 的质量贡献了均值的 %.0f%%'
      % (ph.mean(), ph.mean() / JC, 100 * (ph > 1e3).mean(),
         100 * ph[ph > 1e3].sum() / ph.sum()))

print('\n  [C] 同一统计量在 $r=2,3$ 的斜率（看是否随 $r$ 变薄）')
for r in (2, 3):
    p = DATA[r][np.isfinite(DATA[r])]
    pr = np.array([(p > t).mean() for t in TS])
    ok = (pr > 0.004) & (pr < 0.95)
    s2 = np.polyfit(np.log(TS[ok]), np.log(pr[ok]), 1)[0]
    print('    r=%d  斜率 %.3f（可用点 %d）  $\\Pr(\\Phi>100)=%.1f\\%%$  均值 %.1f $j_c$'
          % (r, s2, ok.sum(), 100 * (p > 100).mean(), p.mean() / JC))

print('\n  [D] 尾的机制检查：随机 $r=1$ 方向上 $\\varepsilon^*_{{\\rm real}}(F)$ 与 $\\Phi$ 的关系（E55 的量）')


def kerZ(Frow):
    _, s, Vh = np.linalg.svd(Frow)
    return Vh[Frow.shape[0]:].T


def eps_real(Frow, N=1400):
    Z = kerZ(Frow)
    g = np.geomspace(1.0, 4.5, N)
    lams = np.concatenate([-g[::-1], g])
    best, bl = np.inf, 0.0
    for lam in lams:
        s = np.linalg.svd((A - lam * np.eye(n)) @ Z, compute_uv=False)[-1]
        if s < best:
            best, bl = s, lam
    return best, bl


rng = np.random.default_rng(5)
rows = []
inf_cnt = 0
for i in range(60):
    f = rng.standard_normal((1, n))
    f = f / np.linalg.norm(f)
    e, lam = eps_real(f)
    pv, k, s, flag = _ns['phi_iter'](f)          # 标量路线：带 cap，不会在溢出上抛异常
    if not np.isfinite(pv):
        inf_cnt += 1
        continue
    rows.append((pv, e, lam))
rows.sort()
print('    按 $\\Phi$ 升序的 12 行（60 个随机方向，$\\varepsilon^*$ 用 2800 点实 $\\lambda$ 扫描）')
print('    %-12s %-12s %-10s %-10s' % ('$\\Phi$', '$\\varepsilon^*_{{\\rm real}}$', '$\\Phi\\varepsilon^*$', '$\\log\\Phi/\\log\\varepsilon^*$'))
for ph, e, lam in rows[:6] + rows[-6:]:
    print('    %-12.4g %-12.4g %-10.3f %-10.2f' % (ph, e, ph * e, np.log(ph) / np.log(e)))
le = np.array([np.log(e) for _, e, _ in rows]); lp = np.array([np.log(p) for p, _, _ in rows])
b, a = np.polyfit(le, lp, 1)
print('    $\\log\\Phi$ 对 $\\log\\varepsilon^*$ 回归：斜率 %.3f，相关 %.3f（$-1$ = 单极点 $\\Phi\\propto1/\\varepsilon^*$）'
      % (b, np.corrcoef(le, lp)[0, 1]))
srt = np.sort([e for _, e, _ in rows])
print('    $\\Pr(\\varepsilon^*<s)$：s=0.05→%.2f  s=0.1→%.2f  s=0.2→%.2f（余维 3 的坏集预言 $\\sim s^3$：%.3f/%.3f/%.3f）'
      % ((srt < .05).mean(), (srt < .1).mean(), (srt < .2).mean(), .05 ** 3 * 8, .1 ** 3 * 8, .2 ** 3 * 8))
print('    （修正：$r=1$ 的坏集是 $\\{f:f\\perp v_u\\}$，在 $S^3$ 里是**余维 1** 的大圆，不是余维 3——'
      '$\\ker F$ 只需含一条不稳定本征向量即可破坏可检测性。预言 $\\Pr(\\varepsilon^*<s)\\sim s$。）')
ss = np.array([.05, .1, .2])
qq = np.array([(srt < s).mean() for s in ss])
print('    幂次拟合：$\\log\\Pr/\\log s$ 斜率 = %s（余维 1 应 $\\approx1$，余维 3 应 $\\approx3$；60 点，纯计数）'
      % '  '.join('%.2f' % (np.log(b) / np.log(a)) for a, b in zip(ss, qq)))
print('    $0.1\\to0.2$ 的倍数 %.1f（余维 1 预言 2.0，余维 3 预言 8.0）'
      % (max(1e-9, (srt < .2).mean()) / max(1e-9, (srt < .1).mean())))

print('\n  判读')
print("""    (1) [A] 若 4000 步与 1500 步的分位数差 $<10^{-6}$，则 E56 的表与 C 的一致不是巧合，且 C 的
        百分位读数只差在 $r=3$（他 40%、我 35%）——那条要按"我这边 34.6±0.8"改写。
    (2) [B] 斜率若显著大于 $-1$，"$\\Phi$ 的期望在随机方向下发散"是可写的一句话；它与 §32.6 的单极点
        模型（预测 $t^{-3}$）矛盾 → 我的极点解释了**扰动方向**上的发散，解释不了**方向分布**的尾。
        这一条按 N5 立案：需要一个 $\\varepsilon^*(F)$ 在 Grassmannian 上的分布理论。
    (3) [C] 若斜率随 $r$ 变陡（$r=1$ 最薄绝对值、$r=3$ 最陡），与"坏集余维随 $r$ 升"同向，可作定性证据。
    (4) [D] 若斜率 $\\approx-1$ 且 $\\Pr(\\varepsilon^*<s)\\sim s^3$，则 [B] 的薄尾应当是 $t^{-3}$——两问
        必有一问错。若 [D] 的 $\\Phi\\varepsilon^*$ 跨 3 个数量级不恒定，则**单极点在方向旋转上不适用**，
        §32.6 的 $\\kappa/(\\varepsilon^*-\\varepsilon)$ 只在固定 $F$、扰动 $A$ 时成立，必须这样写死。""")
