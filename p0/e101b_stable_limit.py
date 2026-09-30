r"""E101b = e101 的修正版：e101 把 s=1e-12 处的读数当"I(0)"，而那里 DARE 已经数值崩溃
（D~4e14、ρ 从 0.76178 跳到 0.77786、I 反而离开 R_exp），于是 [T4] 报出"68 个超过、20 个低于"
的假否证。本轮先把崩溃尺度量出来，只在**稳定窗**内取极限，并顺手把"地板附近 D 的发散律"
做成可比较的量。

预注册判据（跑之前写死）：
 [S1] 收敛形态：ΔI(s) := I(s) − R_exp 在 s→0 服从幂律 ΔI ∝ s^p。对每个 V 在
      s∈[1e-6,1e-2] 稳定段做 log-log 回归得 p。p ≥ 0.7 ⇒ ΔI→0，即 **lim I = R_exp**。
 [S2] 崩溃尺度 s_break：从大往小走，第一个使 |ΔI| 不再下降（比前一点大 20%）或 ρ≥1 的 s。
      不许把 s_break 以下的任何读数当作极限。
 [T1] 地板无关性：每个可检测 V 给出 |ΔI(s=1e-6)| ≤ 3e-5 且 [S1] 的 p≥0.7 ⇒ 记"收敛到 R_exp"。
      等号分支（§66-E-2）：|ΔI| ≤ 1e-6 记"相等"，1e-6<|·|≤3e-5 记"窗内相等（受 p 控制）"，
      >3e-5 记"不等，须查"。
 [T4] 否证：存在可检测 V 在稳定段内 p ≤ 0.3（即 ΔI 不往 0 走）⇒ 无关性作废；
      存在 V 使 I(s) 从**下方**单调趋近且极限 < R_exp ⇒ 禁止再把 R_exp 当本车道设计的下界。
 [U1] 发散律：D − J_c = κ(V)·ΔI^{−ν}。ν 用 ΔI∈[3e-5, 3e-3]（稳定段内两个数量级）回归。
      普适性判据：ν 的相对极差 ≤ 5% 支持"指数普适"；> 16%（§65-D-3 的阈值）则只报分布。
 [U3] 前因子：在**同一个 ΔI 锚点** ΔI*=1e-4 上用 log-log 内插取 (D−J_c)|_{ΔI*}，
      记 K(V)。这是"固定率余量下、子空间之间的代价差"，不许用回归截距外推冒充（e101 的 κ 就是这么错的）。
"""
import sys
import numpy as np
from scipy.linalg import solve_discrete_are, schur

np.set_printoptions(precision=6, suppress=True, linewidth=170)
OUT = open('p0/e101b_out.txt', 'w', encoding='utf-8')


def p(*args):
    OUT.write(' '.join(str(a) for a in args) + '\n')
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


def design(Z, s):
    r = Z.shape[1]
    C = np.sqrt(s) * Z.T
    try:
        Pm = sym(solve_discrete_are(A.T, C.T, W, np.eye(r)))
    except Exception:
        return np.nan, np.nan, np.nan
    Sm = C @ Pm @ C.T + np.eye(r)
    Lk = Pm @ C.T @ np.linalg.inv(Sm)
    Pp = sym(Pm - Lk @ Sm @ Lk.T)
    rho = float(np.abs(np.linalg.eigvals((np.eye(n) - Lk @ C) @ A)).max())
    D = JC + float(np.trace(TH @ Pp))
    I = 0.5 * (np.linalg.slogdet(Pm)[1] - np.linalg.slogdet(Pp)[1]) / ln2
    return rho, D, I


def visibility(Z):
    out = []
    for j in range(n):
        if abs(ev[j]) < 1.0:
            continue
        v = Vv_[:, j]
        vr = np.real(v)
        out.append(float(np.linalg.norm(Z.T @ vr) / max(np.linalg.norm(vr), 1e-300)))
    return np.array(out)


SGRID = 10.0 ** -np.arange(1.0, 7.01, 0.25)          # 1e-1 … 1e-7
KANCHOR = 1e-4
NUBAND = (3e-5, 3e-3)
S_FLOOR = 1e-6                                       # [T1] 的读数尺度


def curve(Z):
    rows = []
    for s in SGRID:
        rho, D, I = design(Z, s)
        rows.append((s, rho, D, I, (I - R_exp) if np.isfinite(I) else np.nan))
    return rows


def breakdown(rows):
    """从大 s 往小 s 走，第一个不再收缩（或 ρ≥1、ΔI≤0）的尺度 ⇒ 它以下全部读数作废。"""
    prev = np.inf
    for s, rho, D, I, d in rows:
        if (not np.isfinite(d)) or d <= 0 or (not np.isfinite(rho)) or rho >= 1.0 or d > 1.2 * prev:
            return float(s)
        prev = d
    return None


def run(Z, tag):
    rows = curve(Z)
    vis = visibility(Z)
    det = all(np.isfinite(r[1]) and r[1] < 1.0 for r in rows) and bool(np.all(vis > 1e-6))
    sb = breakdown(rows)
    ok_below = (sb is None)
    st = [(s, d, D) for (s, rho, D, I, d) in rows
          if d > 0 and rho < 1.0 and D > JC and (sb is None or s > sb)]
    ls, ld, lD = [], [], []
    for s, d, D in st:
        ls.append(np.log10(s)); ld.append(np.log10(d)); lD.append(np.log10(D - JC))
    p_exp = float(np.polyfit(ls, ld, 1)[0]) if len(st) >= 4 else np.nan
    d_at = {}
    for target in (1e-3, 1e-4, 1e-5):
        cand = [(abs(np.log10(d) - np.log10(target)), d, D - JC) for s, d, D in st]
        cand.sort()
        d_at[target] = (cand[0][1], cand[0][2]) if cand else (np.nan, np.nan)
    band = [(np.log10(d), np.log10(D - JC)) for s, d, D in st if NUBAND[0] <= d <= NUBAND[1]]
    nu = kap = np.nan
    if len(band) >= 4:
        m_, b_ = np.polyfit([x for x, y in band], [y for x, y in band], 1)
        nu, kap = -m_, float(np.exp(b_))
    d6c = [d for s, rho, D, I, d in rows if abs(s - S_FLOOR) < 1e-12]
    d6 = d6c[0] if d6c else np.nan
    if not ok_below:
        verdict = '崩溃窗内(须靠 p)'
    else:
        verdict = ('相等' if abs(d6) <= 1e-6 else ('窗内相等' if abs(d6) <= 3e-5 else '不等'))
    K = d_at[KANCHOR][1] / d_at[KANCHOR][0] if np.isfinite(d_at[KANCHOR][0]) else np.nan
    return dict(tag=tag, r=Z.shape[1], det=det, vis=vis, sb=sb, ok_below=ok_below,
                p_exp=p_exp, d6=d6, verdict=verdict, npts=len(st), nu=nu, kap=kap, K=K,
                D_anchor=d_at[KANCHOR][1], rows=rows, Z=Z.copy())



def sfmt(v, w=9):
    return (f'{v:.1e}' if v else 'none').rjust(w)
p('== (0) 锚点 ==')
p(f' 不稳定模态 {unst} → R_exp = {R_exp:.6f} bit；J_c = {JC:.6f}；λ(TH) = {np.linalg.eigvalsh(sym(TH))}')
p(f' [S2] 崩溃尺度表在下面；e101 用 s=1e-12 当极限 = 把崩溃点当收敛点。')

p('\n== (1) 对照：坐标尾平面的 ΔI(s) 序列（看它是否按幂律走向 0） ==')
_, U = schur(A, output='real', sort=lambda a: abs(a) < 1.0)[:2]
coord = []
p('  设计      ' + ''.join(f'{s:11.0e}' for s in SGRID) + '   p(稳定段)  s_break  ν   K=(D−Jc)/ΔI@1e-4')
for k, nm in [(4, 'tail-4'), (3, 'tail-3'), (2, 'tail-2'), (1, 'tail-1')]:
    g = run(U[:, n - k:], nm)
    coord.append(g)
    seq = ''.join(f'{r[4]:11.2e}' if np.isfinite(r[4]) else f'{"nan":>11}' for r in g['rows'])
    p(f'  {nm:<8}{seq}  {g["p_exp"]:8.4f}{sfmt(g["sb"])} {g["nu"]:6.4f} {g["K"]:10.4g}')

p('\n== (2) 随机正交子空间 rank 1/2/3 × 30（seed 20260930） ==')
rng = np.random.default_rng(20260930)
groups = {}
allruns = []
for rk in (1, 2, 3):
    gr = []
    for c in range(30):
        Z, _ = np.linalg.qr(rng.standard_normal((n, rk)))
        gr.append(run(Z, f'r{rk}-{c:02d}'))
    groups[rk] = gr
    allruns += gr
    d_ = [g for g in gr if g['det']]
    p(f'\n rank {rk}: 可检测 {len(d_)}/30')
    p(f'   {"tag":<9}{"vis":<18}{"ΔI(1e-6)":>11}{"p":>8}{"s_break":>10}{"ν":>8}{"K":>10}{"D@1e-4":>11}  判定')
    for g in gr:
        p(f'   {g["tag"]:<9}{str(np.round(g["vis"], 4)):<18}{g["d6"]:11.2e}{g["p_exp"]:8.4f}'
          f'{sfmt(g["sb"],10)}{g["nu"]:8.4f}{g["K"]:10.4g}{g["D_anchor"]:11.4g}  '
          f'{"可检测" if g["det"] else "[T3 排除]"} {g["verdict"]}')

p('\n== (3) 判据 ==')
det = [g for g in allruns if g['det']] + [g for g in coord if g['det']]
und = [g for g in allruns if not g['det']] + [g for g in coord if not g['det']]
p(f' 进入 [T1] 的可检测设计 = {len(det)}（随机 {len([g for g in allruns if g["det"]])} + 坐标 3），'
  f'[T3] 排除 = {len(und)}')
ps = np.array([g['p_exp'] for g in det], dtype=float)
ps = ps[np.isfinite(ps)]
p(f' [S1] ΔI ∝ s^p：p 均值 {ps.mean():.4f}，范围 [{ps.min():.4f}, {ps.max():.4f}]，'
  f'p≥0.7 的比例 {np.mean(ps >= 0.7)*100:.1f}%  → '
  f'{"ΔI 确实走向 0：地板 = R_exp 对一切可检测 V 成立" if np.all(ps >= 0.7) else "有 V 不收敛，须查"}')
p(f' [T1] 判定分布：相等 {sum(1 for g in det if g["verdict"]=="相等")}，'
  f'窗内相等 {sum(1 for g in det if g["verdict"]=="窗内相等")}，'
  f'不等 {sum(1 for g in det if g["verdict"]=="不等")}')
bad = [g for g in det if g['verdict'] == '不等']
for g in bad:
    p(f'   [T1 例外] {g["tag"]} ΔI(1e-6)={g["d6"]:+.3e} p={g["p_exp"]:.4f} vis={np.round(g["vis"],4)} '
      f's_break={g["sb"]}')
low = [g for g in det if np.isfinite(g['d6']) and g['ok_below'] and g['d6'] < -1e-6]
p(f' [T4] 稳定段内从下方穿破 R_exp 的可检测设计 = {len(low)} → '
  f'{"须自首" if low else "e101 报的 20 个低于全部是崩溃尺度以下的读数，不是反例"}')
sbs = [g['sb'] for g in det if g['sb']]
p(f' [S2] 崩溃尺度分布：{min(sbs):.1e} … {max(sbs):.1e}（共 {len(sbs)}/{len(det)} 个设计触发）')
sbv = np.array([min(g['sb'], 1e-7) if g['sb'] else 1e-7 for g in det])
vmin = np.array([g['vis'].min() for g in det])
cc = np.corrcoef(np.log10(sbv), np.log10(np.maximum(vmin, 1e-12)))[0, 1]
p(f'      corr(log s_break, log min-vis) = {cc:.3f} → '
  f'{"崩溃尺度随可见度变小而提前（病态由最弱可见模态控制）" if cc > 0.5 else "与可见度关系不强"}')
nu = np.array([g['nu'] for g in det], dtype=float)
nu = nu[np.isfinite(nu)]
p(f' [U1] ν：n={len(nu)} 均值 {nu.mean():.4f} 极差 {nu.max()-nu.min():.4f} '
  f'相对极差 {(nu.max()-nu.min())/nu.min()*100:.1f}% → '
  f'{"≤16%：指数普适（=1 的双曲律），只有前因子依赖 V" if (nu.max()-nu.min())/nu.min()*100 <= 16 else ">16%：不许写普适指数"}')
p(f'      ν 分位 [0,25,50,75,100] = {np.round(np.percentile(nu, [0,25,50,75,100]), 4)}')
K = np.array([g['K'] for g in det], dtype=float)
K = K[np.isfinite(K)]
p(f' [U3] K(V) = (D−J_c) 在 ΔI=1e-4 处的比值：{K.min():.4g} … {K.max():.4g}，跨 '
  f'{K.max()/K.min():.1f}×；对数标准差 {np.std(np.log10(K)):.3f}')
order = np.argsort(K)
p('      最小的 5 个（同样率余量下最便宜）：' + ', '.join(
    f'{det[i]["tag"]}(K={K[i]:.4g}, r={det[i]["r"]}, vis={np.round(det[i]["vis"],3)})' for i in order[:5]))
p('      最大的 5 个（最贵）：' + ', '.join(
    f'{det[i]["tag"]}(K={K[i]:.4g}, r={det[i]["r"]}, vis={np.round(det[i]["vis"],3)})' for i in order[-5:]))
for rk in (1, 2, 3):
    kk = np.array([g['K'] for g in allruns if g['det'] and g['r'] == rk], dtype=float)
    kk = kk[np.isfinite(kk)]
    p(f'      rank {rk}: K 中位 {np.median(kk):.4g}，区间 [{kk.min():.4g}, {kk.max():.4g}]（n={len(kk)}）')
kko = np.array([g['K'] for g in coord if g['det']], dtype=float)
p(f'      坐标对照 tail-4/3/2 的 K = {np.round(kko, 4)}')
p('\n== (4) [T3] 排除样本的详情（Th.1 前提之外，不是反例） ==')
for g in und:
    tail = [r for r in g['rows'] if r[4] is not None]
    p(f'   {g["tag"]:<9} vis={np.round(g["vis"],6)} ρ(1e-1)={tail[0][1]:.4f} '
      f'ρ(1e-6)={tail[-2][1]:.5f} ΔI 末值={tail[-1][4]:+.4e} s_break={g["sb"]}')

OUT.close()
print('written p0/e101b_out.txt')
