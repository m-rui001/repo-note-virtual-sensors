r"""E101 = 把 §66-C 的"率地板与平面无关"从**坐标尾平面**推广到**任意（随机）子空间**，
并顺手量一次"地板附近的代价发散律"，作为动态舞台定理候选（板 §2.2 / §66-G-④）。

背景（已定）：
 · 2606.31396 Th.1（任意观测律、前提=均方可观测）给 I ≥ R_exp = Σ_{|λ|≥1} log2|λ| = 1.168539。
 · 本车道 §66-C 在 tail-4/3/2 三个**坐标**平面上测到 s→0 时 I → 1.168539（超出 +0.000000），
   而代价 D 把三个平面分开（642.58 / 656.08 / 914.69）⇒"受限传感的惩罚全在代价侧"。
 · 但那三个平面都是**坐标面**。"地板与 V 无关"如果只在坐标面成立，它是数值巧合不是定理。

⇒ 本轮改判据为**随机子空间**上的成/败，并且预注册（数字出来之前先写死）：

 [T1] 平局无关性：对每个可检测（ρ((I−LC)A)<1 在该尺度成立）且不稳定模态全可见的随机子空间 V，
      lim_{s→0} I(s) = R_exp。判据带**等号分支**（§66-E-2 的教训）：
        |I−R_exp| ≤ 1e-6 ⇒ 记"相等"；1e-6 < |·| ≤ 1e-3 ⇒ 记"近平等，须查"；> 1e-3 ⇒ 记"不等"。
 [T2] 代价发散：对同一批 V，D(s) 在 s→0 单调发散（D(1e-12)/D(1e-2) ≫ 1）。
 [T3] 前提外：不可检测 / 有不稳定模态不可见的 V 必须**排除在 T1 之外**，并单独报它的 I 极限
      （若 < R_exp，那正是 Th.1 前提被绕开的样子，不是反例）。
 [T4] 否证：存在可检测 V 使 lim I > R_exp+1e-6 ⇒ 无关性作废；
      存在可检测 V 使 lim I < R_exp−1e-6 ⇒ 禁止再把 R_exp 当**本车道设计**的下界引用，须自首。

 [U1] 指数普适性：地板附近 D ≈ J_c + κ(V)·(I−R_exp)^{−ν}。对每个 V 用
      ΔI ∈ [1e-3, 1e-1] 这段做 log-log 斜率 −ν 的最小二乘。
      若 detectable 样本间 ν 的极差 ≤ 5%（相对）⇒ 支持"指数普适、只有前因子 κ 依赖 V"。
 [U2] 否证阈值（沿用 §65-D-3 定下的 16%）：极差 > 16% ⇒ 不许写"普适指数"，只报实测分布。
"""
import sys
import numpy as np
from scipy.linalg import solve_discrete_are, schur

np.set_printoptions(precision=6, suppress=True, linewidth=170)
OUT = open('p0/e101_out.txt', 'w', encoding='utf-8')


def p(*args):
    s = ' '.join(str(a) for a in args)
    OUT.write(s + '\n')
    OUT.flush()


_src = open('p0/exp_c_audit.py', encoding='utf-8').read().split("print(r'== E53")[0]
_ns = {'__name__': 'p'}
exec(compile(_src, 'p0/exp_c_audit.py[preamble]', 'exec'), _ns)
A, W, TH, JC, n = _ns['A'], _ns['W'], _ns['TH'], _ns['JC'], _ns['n']
sym = _ns['sym']
ln2 = np.log(2.0)

ev = np.linalg.eigvals(A)
unst = np.sort(np.abs(ev[np.abs(ev) >= 1.0]))
R_exp = float(np.sum(np.log2(unst)))
w_, Vv_ = np.linalg.eig(A)
VU = np.array([np.real(Vv_[:, j]) if abs(np.imag(ev[j])) < 1e-12 else None
               for j in range(n)], dtype=object)

p('== (0) 锚点与地板 ==')
p(f' lambda(A) = {np.sort_complex(ev)}')
p(f' 不稳定模态 |λ|≥1: {unst}  → R_exp = {R_exp:.6f} bit   log2|det A| = '
  f'{np.log2(abs(np.linalg.det(A))):.6f}')
p(f' J_c = tr(W Pc) = {JC:.6f}   λ(TH) = {np.linalg.eigvalsh(sym(TH))}')
p(f' ⇒ TH ≻ 0：只要 Pp 在任何方向发散，D = J_c + tr(TH·Pp) 必发散（[T2] 的充分条件）')


def design(Z, s):
    """观测子空间 span(Z)（Z: n×r 正交基）+ 尺度 s：C = sqrt(s)·Zᵗ。返回 (rho, D, I, r)。"""
    r = Z.shape[1]
    C = np.sqrt(s) * Z.T
    try:
        Pm = sym(solve_discrete_are(A.T, C.T, W, np.eye(r)))
    except Exception:
        return np.nan, np.nan, np.nan, r
    Sm = C @ Pm @ C.T + np.eye(r)
    Lk = Pm @ C.T @ np.linalg.inv(Sm)
    Pp = sym(Pm - Lk @ Sm @ Lk.T)
    rho = float(np.abs(np.linalg.eigvals((np.eye(n) - Lk @ C) @ A)).max())
    D = JC + float(np.trace(TH @ Pp))
    I = 0.5 * (np.linalg.slogdet(Pm)[1] - np.linalg.slogdet(Pp)[1]) / ln2
    return rho, D, I, r


def visibility(Z):
    """每个不稳定模态在 span(Z) 里的可见度（PBH：C·v=0 ⇔ v⊥span(Z)）。"""
    out = []
    for j in range(n):
        if abs(ev[j]) < 1.0:
            continue
        v = Vv_[:, j]
        vr, vi = np.real(v), np.imag(v)
        out.append(float(max(np.linalg.norm(Z.T @ vr) / max(np.linalg.norm(vr), 1e-300),
                             np.linalg.norm(Z.T @ vi) / max(np.linalg.norm(vi), 1e-300))))
    return np.array(out)


GRID = np.array([1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-8, 1e-10, 1e-12])
FIT = 10.0 ** -np.arange(1.0, 9.1, 0.4)


def run(Z, tag):
    res = []
    for s in GRID:
        rho, D, I, r = design(Z, s)
        res.append((s, rho, D, I))
    det_all = all(np.isfinite(r) and r < 1.0 for _, r, _, _ in res)
    vis = visibility(Z)
    see_all = bool(np.all(vis > 1e-6))
    s_lo = GRID[-1]
    I_end = res[-1][3]
    D_end = res[-1][2]
    D_hi = res[0][2]
    dif = I_end - R_exp if np.isfinite(I_end) else np.nan
    verdict = ('相等' if abs(dif) <= 1e-6 else ('近平等' if abs(dif) <= 1e-3 else '不等')) if np.isfinite(dif) else '无值'
    # 指数拟合：ΔI ∈ [1e-3,1e-1]
    fit = [(s, design(Z, s)) for s in FIT]
    xs, ys = [], []
    for s, (rho, D, I, r) in fit:
        if np.isfinite(I) and np.isfinite(D) and 1e-3 <= I - R_exp <= 1e-1 and D > JC:
            xs.append(np.log(I - R_exp))
            ys.append(np.log(D - JC))
    nu = np.nan
    b_ = np.nan
    if len(xs) >= 4:
        m_, b_ = np.polyfit(np.array(xs), np.array(ys), 1)
        nu = -m_
    return dict(tag=tag, r=r, res=res, det_all=det_all, see_all=see_all, vis=vis, Z=Z.copy(),
                I_end=I_end, dif=dif, verdict=verdict,
                D_ratio=D_end / D_hi if np.isfinite(D_end) and D_hi > 0 else np.nan,
                nu=nu, kappa=float(np.exp(b_)) if np.isfinite(nu) else np.nan)


p('\n== (1) 对照组：坐标尾平面（§66-C 的四个） ==')
_, U = schur(A, output='real', sort=lambda a: abs(a) < 1.0)[:2]
coord = []
for k, nm in [(4, 'tail-4'), (3, 'tail-3'), (2, 'tail-2'), (1, 'tail-1')]:
    r_ = run(U[:, n - k:], nm)
    coord.append(r_)
    p(f" {nm:<8} r={r_['r']} 可检测={r_['det_all']} 全可见={r_['see_all']} "
      f"vis={np.round(r_['vis'], 4)} I(1e-12)={r_['I_end']:.6f} 超出{r_['dif']:+.6f} "
      f"[{r_['verdict']}] D比值={r_['D_ratio']:.3e} ν={r_['nu']:.4f} κ={r_['kappa']:.4g}")

p('\n== (2) 主实验：随机正交子空间，rank r=1,2,3 × 30 样本（seed 固定） ==')
rng = np.random.default_rng(20260930)
samples = []
for rk in (1, 2, 3):
    group = []
    for c in range(30):
        G = rng.standard_normal((n, rk))
        Z, _ = np.linalg.qr(G)
        group.append(run(Z, f'rand-r{rk}-{c:02d}'))
    samples.append((rk, group))
    d_ = [g for g in group if g['det_all'] and g['see_all']]
    und = [g for g in group if not (g['det_all'] and g['see_all'])]
    p(f'\n -- rank {rk}: 可检测且全可见 {len(d_)}/30，其余 {len(und)}/30')
    for g in d_:
        p(f"   {g['tag']:<15} vis={np.round(g['vis'], 4)} I_end={g['I_end']:.6f} "
          f"超出{g['dif']:+.6f} [{g['verdict']}] D比值={g['D_ratio']:.3e} ν={g['nu']:.4f}")
    for g in und:
        p(f"   {g['tag']:<15} [T3 排除] vis={np.round(g['vis'], 6)} 可检测={g['det_all']} "
          f"全可见={g['see_all']} I_end={g['I_end']} ρ(1e-12)={g['res'][-1][1]}")

p('\n== (3) 判据汇总 ==')
all_det = [g for rk, gr in samples for g in gr if g['det_all'] and g['see_all']]
alls_ = coord + all_det
gt = [g for g in all_det if g['dif'] > 1e-6]
lt = [g for g in all_det if g['dif'] < -1e-6]
near = [g for g in all_det if 1e-6 < abs(g['dif']) <= 1e-3]
p(f" 进入 [T1] 的可检测+全可见样本 = {len(all_det)}（随机）+ "
  f"{sum(1 for g in coord if g['det_all'] and g['see_all'])}（坐标对照）= {len(alls_)}")
p(f' [T1] 相等 {sum(1 for g in alls_ if g["verdict"] == "相等")}/{len(alls_)}，'
  f'近平等 {len(near)}，不等 {sum(1 for g in alls_ if g["verdict"] == "不等")}')
p(f' [T4 否证-上] lim I > R_exp+1e-6 的个数 = {len(gt)}  → '
  f'{"无关性作废" if gt else "无关性未被否证"}')
p(f' [T4 否证-下] lim I < R_exp−1e-6 的个数 = {len(lt)}  → '
  f'{"须自首：不能把 R_exp 当本车道设计下界" if lt else "R_exp 可继续作为下界引用"}')
if gt or lt:
    for g in (gt + lt):
        p(f'   例外：{g["tag"]} 超出 {g["dif"]:+.8f}')
p(f' [T2] D 发散比 = {np.nanmin([g["D_ratio"] for g in alls_]):.3e} … '
  f'{np.nanmax([g["D_ratio"] for g in alls_]):.3e}')
nu = np.array([g['nu'] for g in alls_], dtype=float)
nu = nu[np.isfinite(nu)]
p(f' [U1] ν 分布：n={len(nu)} 均值 {nu.mean():.4f} 极差 {nu.max()-nu.min():.4f} '
  f'相对极差 {(nu.max()-nu.min())/max(nu.min(),1e-12)*100:.1f}%')
p(f'      ν 的 12 个分位：{np.round(np.percentile(nu, [0, 10, 25, 50, 75, 90, 100]), 4)}')
p(f' [U2] 阈值 16%：相对极差 {(nu.max()-nu.min())/max(nu.min(),1e-12)*100:.1f}% → '
  f'{"支持普适指数" if (nu.max()-nu.min())/max(nu.min(),1e-12)*100 <= 16 else "不许写普适指数"}')
kap = np.array([g['kappa'] for g in alls_], dtype=float)
kap = kap[np.isfinite(kap)]
if len(kap):
    p(f'      前因子 κ 分布（应强烈依赖 V）：{np.nanmin(kap):.4g} … {np.nanmax(kap):.4g}'
      f'，跨 {np.nanmax(kap)/max(np.nanmin(kap),1e-300):.3g}×')

p('\n== (4) 一个可检测的 rank-1 随机子空间：它的整条 (I, D) 曲线 ==')
best = None
for g in all_det:
    if g['r'] == 1:
        best = g
        break
if best is None:
    p('  rank-1 组里没有可检测样本（预期不该发生：不可见条件是余维≥1 的零测集）')
else:
    p(f"  取 {best['tag']}，Zᵗ = {np.round(best['Z'].T, 4)}")
    for s, rho, D, I in best['res']:
        p(f'   s={s:8.0e}  ρ={rho:.5f}  I={I:.6f} (Δ={I-R_exp:+.6f})  D={D:.4f}')

OUT.close()
print('written p0/e101_out.txt')
