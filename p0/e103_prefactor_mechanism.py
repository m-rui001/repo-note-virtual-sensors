r"""E103 = e101b 的律（D−J_c = K(V)/(I−R_exp)，指数 1 普适、K 跨 6437×）的**机制与可计算前因子**。

三条腿（都预注册）：
 [M1] 标度 DARE 的极限谱：令 θ(s) := eig(s·Pm)。标量情形解 DARE 得 θ → |λ|²−1（不稳）、0（稳定），
      于是 I(0) = Σ ½log2(1+θ*_i) = Σ_{|λ|>1} log2|λ| = R_exp，**与 V 无关**。
      判据：对每个可检测 V，eig(sPm) 在 s=1e-5,1e-6 两级之差的每分量相对变化 ≤ 2%，
      且 Σ½log2(1+θ_i(s)) → R_exp（绝对差 ≤ 1e-5）。成 ⇒ 地板无关性有结构解释，不只是拟合。
 [M2] 前因子的乘积形式：a(V) := lim ΔI/s，b(V) := lim (D−J_c)·s ⇒ K(V) = a(V)·b(V)。
      判据：|K_direct − a·b|/(a·b) ≤ 10%（K_direct 取 ΔI=1e-4 锚点处的 (D−J_c)·ΔI）。
 [M3] b(V) 的身份：b = tr(Θ·σ*)，σ* := eig(s·Pp) 的极限矩阵。判据：tr(TH·sPp) 在 s=1e-4…1e-6
      收敛（相对变化 ≤ 3%）。若收敛，则"近地板代价前因子"是**给定 V 就能算**的量。

 [G1] 设计规则的可证伪检验（这条直接接 C 的 `.work3/paper_prep.md:57` row 7：
      "取 Θ=Kᵗ(R+BᵗP_cB)K 前 r 个特征向量，距最优 1.08–1.20×"）：
      在 rank r 的**全体**可检测子空间里，最小化 b(V) 的是不是 Θ 的本征方向族？
      两种读法都测（Θ 最小 r 个本征方向 / 最大 r 个），并与 60 个随机样本比。
      判据：若 Θ-规则（任一读法）的 b 落在随机样本的 5% 分位以下 ⇒ "地板侧规则"在近地板口径下**有机制依据**；
      否则当场记为否证，并且 C 那条规则只能按"在 I=3 口径的经验规则"卖，不许升格。
      （注意口径差别：C 的 1.08–1.20× 是对**已知最好值**量的，而 §65-C 已证那个"最好值"没有证书。）
"""
import numpy as np
from scipy.linalg import solve_discrete_are, schur

np.set_printoptions(precision=6, suppress=True, linewidth=170)
OUT = open('p0/e103_out.txt', 'w', encoding='utf-8')


def p(*a):
    OUT.write(' '.join(str(x) for x in a) + '\n')
    OUT.flush()


_src = open('p0/exp_c_audit.py', encoding='utf-8').read().split("print(r'== E53")[0]
_ns = {'__name__': 'p'}
exec(compile(_src, 'p0/exp_c_audit.py[preamble]', 'exec'), _ns)
A, W, TH, JC, n = _ns['A'], _ns['W'], _ns['TH'], _ns['JC'], _ns['n']
sym = _ns['sym']
ln2 = np.log(2.0)
ev = np.linalg.eigvals(A)
theta_star_target = np.sort(np.array([max(abs(l) ** 2 - 1.0, 0.0) for l in ev]))
R_exp = float(np.sum(np.log2(np.abs(ev[np.abs(ev) >= 1.0]))))


def parts(Z, s):
    r = Z.shape[1]
    C = np.sqrt(s) * Z.T
    Pm = sym(solve_discrete_are(A.T, C.T, W, np.eye(r)))
    Sm = C @ Pm @ C.T + np.eye(r)
    Lk = Pm @ C.T @ np.linalg.inv(Sm)
    Pp = sym(Pm - Lk @ Sm @ Lk.T)
    rho = float(np.abs(np.linalg.eigvals((np.eye(n) - Lk @ C) @ A)).max())
    D = JC + float(np.trace(TH @ Pp))
    I = 0.5 * (np.linalg.slogdet(Pm)[1] - np.linalg.slogdet(Pp)[1]) / ln2
    return rho, D, I, Pm, Pp


def vis(Z):
    w_, Vv = np.linalg.eig(A)
    out = []
    for j in range(n):
        if abs(ev[j]) < 1.0:
            continue
        vr = np.real(Vv[:, j])
        out.append(float(np.linalg.norm(Z.T @ vr) / max(np.linalg.norm(vr), 1e-300)))
    return np.array(out)


p('== (0) 目标极限谱 θ* = eig(sPm)|_{s→0} 的标量式预测 ==')
p(f' λ(A) = {np.sort_complex(ev)}')
p(f' 预测 θ*_i = max(|λ_i|²−1, 0) = {np.round(theta_star_target, 4)}')
p(f' 校验：Σ ½log2(1+θ*) = {np.sum(0.5*np.log2(1+theta_star_target)):.6f}  vs R_exp = {R_exp:.6f}'
  f'  差 {abs(np.sum(0.5*np.log2(1+theta_star_target)) - R_exp):.2e}')
p(' ⇒ 若 [M1] 成立，地板与 V 无关**不是巧合**：标度预测协方差把每个不稳定模态钉在 |λ|²−1，')
p('   稳定模态钉在 0（不贡献率），所以 I(0) 只由谱决定，与 range(S) 无关。')

S1, S2 = 1e-5, 1e-6
p('\n== (1) [M1/M2/M3] 逐设计 ==')
p(f" {'设计':<14}{'eig(sPm)@1e-5':<28}{'Δθ相对':>9}{'a=ΔI/s':>11}{'b=(D−Jc)s':>12}"
  f"{'tr(TH·sPp)漂移':>15}{'K_direct':>11}{'a·b':>11}{'偏离':>8}")
_, U = schur(A, output='real', sort=lambda a: abs(a) < 1.0)[:2]
rng = np.random.default_rng(20260931)
designs = [(U[:, n - k:], f'tail-{k}') for k in (4, 3, 2)]
for rk in (1, 2, 3):
    for c in range(4):
        Z, _ = np.linalg.qr(rng.standard_normal((n, rk)))
        designs.append((Z, f'r{rk}-{c:02d}'))
recs = []
for Z, nm in designs:
    t1 = np.sort(np.linalg.eigvalsh(sym(Z.T @ (S1 * parts(Z, S1)[3]) @ Z)))
    rho1, D1, I1, Pm1, Pp1 = parts(Z, S1)
    rho2, D2, I2, Pm2, Pp2 = parts(Z, S2)
    th1 = np.sort(np.linalg.eigvalsh(sym(Z.T @ (S1 * Pm1) @ Z)))
    th2 = np.sort(np.linalg.eigvalsh(sym(Z.T @ (S2 * Pm2) @ Z)))
    dth = np.max(np.abs(th2 - th1) / np.maximum(th1, 1e-12)) if len(th1) else np.nan
    a1, a2 = (I1 - R_exp) / S1, (I2 - R_exp) / S2
    b1, b2 = (D1 - JC) * S1, (D2 - JC) * S2
    drift = abs(b2 - b1) / b1
    # K_direct：在 ΔI≈1e-4 处取 (D−Jc)·ΔI
    Kd = np.nan
    for s in 10.0 ** -np.arange(2.0, 6.01, 0.2):
        rho, D, I = parts(Z, s)[:3]
        d = I - R_exp
        if rho < 1 and 3e-5 <= d <= 3e-4:
            Kd = (D - JC) * d
    K = a2 * b2
    recs.append(dict(nm=nm, r=Z.shape[1], th=th2, a=a2, b=b2, K=K, Kd=Kd, dth=dth, drift=drift,
                     visv=vis(Z), Z=Z))
    p(f' {nm:<14}{str(np.round(th1, 4)):<28}{dth*100:8.2f}%{a2:11.4g}{b2:12.5g}'
      f'{drift*100:14.2f}%{Kd:11.4g}{K:11.4g}{abs(K-Kd)/Kd*100 if np.isfinite(Kd) else np.nan:7.1f}%')

p('\n== (2) 判据 ==')
p(' [M1] Σ½log2(1+θ) 到 R_exp 的绝对差（s=1e-6）：')
for g in recs:
    s_ = 0.5 * np.sum(np.log2(1 + np.maximum(g['th'], 0)))
    p(f"   {g['nm']:<10} Σ½log2(1+θ)={s_:.6f}  差 {s_-R_exp:+.2e}  θ={np.round(g['th'],4)}")
p(f" [M3] b 的收敛漂移最大 {max(g['drift'] for g in recs)*100:.2f}%（阈值 3%）")
devs = [abs(g['K'] - g['Kd']) / g['Kd'] for g in recs if np.isfinite(g['Kd'])]
p(f" [M2] K_direct 与 a·b 的相对偏离：中位 {np.median(devs)*100:.1f}%，最大 {np.max(devs)*100:.1f}%"
  f'（阈值 10%）')

p('\n== (3) [G1] 近地板前因子 b(V) 的最小化：Θ 本征方向规则 vs 随机 60 样本 ==')
wT, VT = np.linalg.eigh(sym(TH))
p(f' λ(Θ) 升序 = {np.round(wT, 4)}')
for rk in (1, 2, 3):
    rnd = []
    for c in range(60):
        Z, _ = np.linalg.qr(rng.standard_normal((n, rk)))
        v = vis(Z)
        if np.any(v < 1e-6):
            continue
        try:
            rho, D, I = parts(Z, S2)[:3]
        except Exception:
            continue
        if rho >= 1:
            continue
        rnd.append(((D - JC) * S2, np.log((I - R_exp) / S2 * (D - JC) * S2), Z))
    bs = np.array([x[0] for x in rnd])
    rule_lo = VT[:, :rk]
    rule_hi = VT[:, n - rk:]
    outs = {}
    for lab, Zr in [('Θ最小r方向', rule_lo), ('Θ最大r方向', rule_hi)]:
        try:
            rho, D, I = parts(Zr, S2)[:3]
            outs[lab] = ((D - JC) * S2 if rho < 1 else np.nan, rho)
        except Exception:
            outs[lab] = (np.nan, np.nan)
    p(f'\n rank {rk}: 随机可用样本 {len(rnd)}，b 分布 min {bs.min():.5g} / 中位 {np.median(bs):.5g} '
      f'/ max {bs.max():.5g}')
    for lab, (bv, rho) in outs.items():
        if not np.isfinite(bv):
            p(f'   {lab}: ρ={rho:.4f} 不可检测 ⇒ 无 b')
            continue
        q = np.mean(bs <= bv) * 100
        p(f'   {lab}: b={bv:.5g}（ρ={rho:.4f}）落在随机样本的 {q:.1f}% 分位 → '
          f'{"支持：近地板前因子确实由 Θ 方向选择压住" if q <= 5 else "否证：该规则不占优，只能按经验规则卖"}')
    best = min(rnd, key=lambda x: x[0])
    p(f'   随机里最小的 b={best[0]:.5g}，其 vis={np.round(vis(best[2]), 4)}；'
      f'与 Θ-规则差距 {best[0]/np.nanmin([v[0] for v in outs.values() if np.isfinite(v[0])]):.3f}×')
    p(f'   b 的跨样本倍差 = {bs.max()/bs.min():.1f}×，而 ln(a) 的散布 = '
      f'{np.std([x[1] for x in rnd]):.3f}（对数）⇒ 判定 K 的散布主要来自哪条腿')
vv = np.array([np.min(g['visv']) for g in recs])
bb = np.array([np.log(g['b']) for g in recs])
if len(vv) > 3:
    p(f'\n corr(log b, log min-vis)（含坐标对照，n={len(vv)}）= '
      f'{np.corrcoef(bb, np.log(np.maximum(vv, 1e-12)))[0,1]:.3f}')
OUT.close()
print('written p0/e103_out.txt')
