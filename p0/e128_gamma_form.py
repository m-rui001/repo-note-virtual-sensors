# -*- coding: utf-8 -*-
r"""R64/R65-D 实验 128：把 §78 从 2604.20369 eq(43) 反解出的形式当真拟基，检验
"§75-C 的限制句（K 只能排序、不能定价）是否只是基函数拟错造成的"。

闭式（已用 PDF 文本层逐字核对）：
  F(D) = [ log|a| + ½ log(1 + σ² m /(D − D_min))]_+ ，自然对数
反解（令 ΔI = I − log|a|/ln2）：
  D − D_min = K̃ / (e^{γ ΔI} − 1)，  K̃ = σ²m，  γ = 2 ln 2 = 1.3862944

对照我原来的律：F1 用 (D−J_c)·ΔI = K，即 γ=1 型双曲 + 地板写死 J_c。
两处偏离：① 指数（双曲 vs 指数减一）；② 地板（J_c vs 该设计自己的 D_min）。
本实验把这两处**分开**放开，回答"到底是哪一个在坏事"。

 [G1] 2 参数：放开 (K̃, γ)，地板写死 J_c。
 [G2] 3 参数：放开 (K̃, γ, D_min)。
 [G3] 1 参数：γ ≡ 2ln2 写死，地板写死 J_c（理论值直接用，不拟合）。
 [G4] 1 参数：γ ≡ 2ln2 写死，地板用该设计的 D_min（s→∞ 多解一次 DARE 得到）。
 [G5] 拟合窗**只用** ΔI ∈ [0.30, 1.00]（§78-D 写死的中间段），定价仍在 ΔI = 1.831461（I=3）。
      真值仍由 80 步二分给出，不许既拟合又当答案；窗内点数 <4 的设计逐个报出并剔除。
 [G6] 预注册判据（跑前写死，跑后不许改）：记 med|ẽ| 为四个口径里最好者的中位绝对误差。
      只有当 med|ẽ| ≤ ½ × 30.67% = 15.335% **且** rank 1/2/3 的有符号中位不再翻号，
      才允许把 §75-C / §76-E 的限制句换成新句；否则限制句原样保留，
      并把本实验的否证（哪个自由度不够）逐条记进 board。
"""
import sys
import numpy as np
from scipy.linalg import solve_discrete_are, eigh
from scipy.optimize import curve_fit
sys.stdout.reconfigure(encoding='utf-8')

_src = open('p0/exp_c_audit.py', encoding='utf-8').read().split("print(r'== E53")[0]
_ns = {'__name__': 'p'}
exec(compile(_src, 'p0/exp_c_audit.py[preamble]', 'exec'), _ns)
A, W, TH, JC, n = _ns['A'], _ns['W'], _ns['TH'], _ns['JC'], _ns['n']
sym = _ns['sym']
R_EXP = 1.168539
TARGET = 3.0
DIT = TARGET - R_EXP
GAMMA_TH = 2.0 * np.log(2.0)
WLO, WHI = 0.30, 1.00
SGRID = 10.0 ** np.linspace(-6.0, 2.5, 90)
SFLOOR = 10.0 ** 8.0
F1_MED = 0.3067  # e127 的 F1 med|e|，判据 [G6] 的基准


def curve(Z, s):
    r = Z.shape[1]
    Ir = np.eye(r)
    C = np.sqrt(s) * Z.T
    Pm = sym(solve_discrete_are(A.T, C.T, W, Ir))
    Sm = C @ Pm @ C.T + Ir
    Lk = Pm @ C.T @ np.linalg.inv(Sm)
    Pp = sym(Pm - Lk @ Sm @ Lk.T)
    I = 0.5 * np.log(np.linalg.det(Ir + C @ Pm @ C.T)) / np.log(2.0)
    return I, JC + np.trace(TH @ Pp)


def true_at(Z, target=TARGET):
    lo, hi = -7.0, 5.0
    if not (curve(Z, 10 ** lo)[0] < target < curve(Z, 10 ** hi)[0]):
        return None
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if curve(Z, 10 ** mid)[0] < target:
            lo = mid
        else:
            hi = mid
    return curve(Z, 10 ** (0.5 * (lo + hi)))[1]


def f_exp(di, K, gam):
    return K / (np.exp(gam * di) - 1.0)


def f_exp_free(di, K, gam, dmin):
    return dmin + K / (np.exp(gam * di) - 1.0)


def f_fixedgam(di, K):
    return K / (np.exp(GAMMA_TH * di) - 1.0)


def f_fixedgam_floor(di, K, dmin):
    return dmin + K / (np.exp(GAMMA_TH * di) - 1.0)


_, Vth = eigh(TH)
rng = np.random.default_rng(20260930)
pool = []
for rk in (1, 2, 3):
    for c in range(10):
        Z, _ = np.linalg.qr(rng.standard_normal((n, rk)))
        pool.append((f'r{rk}-{c:02d}', rk, Z))
    pool.append((f'THmin-{rk}', rk, Vth[:, :rk]))
    pool.append((f'THmax-{rk}', rk, Vth[:, n - rk:]))

print('== 基函数对照：拟合窗 ΔI∈[%.2f,%.2f]，定价 ΔI=%.6f（I=3） ==' % (WLO, WHI, DIT))
print('   γ 的理论值 2ln2=%.6f   F1(双曲+地板Jc) 的 med|e| 基准=%.2f%%' % (GAMMA_TH, 100 * F1_MED))
print('  design    rank nwin   Dfloor     D_true    e(G1)%   e(G2)%   e(G3)%   e(G4)%   gam_fit   K_fit')
rows = []
for name, rk, Z in pool:
    Dtrue = true_at(Z)
    if Dtrue is None:
        print(f'  {name:<9}{rk:>4}   (I=3 不在扫描范围)')
        continue
    Is, Ds = [], []
    for s in SGRID:
        try:
            I, D = curve(Z, s)
        except Exception:
            continue
        if np.isfinite(I) and np.isfinite(D):
            Is.append(I); Ds.append(D)
    Is = np.array(Is); Ds = np.array(Ds)
    di = Is - R_EXP
    ok = (di >= WLO) & (di <= WHI) & (Ds > JC)
    if ok.sum() < 4:
        print(f'  {name:<9}{rk:>4}   窗内点不足（{ok.sum()}），剔除')
        continue
    dif, Df = di[ok], Ds[ok] - JC
    dmin_floor = curve(Z, SFLOOR)[1] - JC
    out = {}
    try:
        p1, _ = curve_fit(f_exp, dif, Df, p0=[Df.mean() * dif.mean(), 1.0],
                          bounds=([1e-12, 1e-3], [np.inf, 20.0]), maxfev=40000)
        out['G1'] = (JC + f_exp(DIT, *p1)) / Dtrue - 1.0
        out['g1'] = p1
    except Exception as e:
        out['G1'] = None; out['g1'] = (np.nan, np.nan); print(f'   G1 fail {name}: {e}')
    try:
        p2, _ = curve_fit(f_exp_free, dif, Df, p0=[Df.mean() * dif.mean(), 1.0, 0.05 * Df.mean()],
                          bounds=([1e-12, 1e-3, 0.0], [np.inf, 20.0, np.inf]), maxfev=40000)
        out['G2'] = (JC + f_exp_free(DIT, *p2)) / Dtrue - 1.0
        out['g2'] = p2
    except Exception as e:
        out['G2'] = None; out['g2'] = (np.nan,) * 3; print(f'   G2 fail {name}: {e}')
    try:
        p3, _ = curve_fit(f_fixedgam, dif, Df, p0=[Df.mean() * dif.mean()],
                          bounds=([1e-12], [np.inf]), maxfev=40000)
        out['G3'] = (JC + f_fixedgam(DIT, *p3)) / Dtrue - 1.0
        out['g3'] = p3
    except Exception as e:
        out['G3'] = None; out['g3'] = (np.nan,)
    try:
        # 地板已知 ⇒ 只剩一个尺度参数，线性最小二乘即可（无需迭代）
        base = f_fixedgam(dif, 1.0)
        K4 = float(np.dot(base, Df - dmin_floor) / np.dot(base, base))
        out['G4'] = (JC + f_fixedgam_floor(DIT, K4, dmin_floor)) / Dtrue - 1.0
        out['g4'] = (K4,)
    except Exception as e:
        out['G4'] = None; out['g4'] = (np.nan,); print(f'   G4 fail {name}: {e}')

    def fmt(k):
        v = out[k]
        return '    n/a ' if v is None else f'{100 * v:>8.2f}'
    gam = out['g1'][1] if out['g1'] is not None else np.nan
    kk = out['g1'][0] if out['g1'] is not None else np.nan
    rows.append((name, rk, ok.sum(), dmin_floor, Dtrue, out))
    print(f'  {name:<9}{rk:>4}  {ok.sum():>4}  {JC + dmin_floor:>9.4f}  {Dtrue:>9.4f}  '
          f'{fmt("G1")} {fmt("G2")} {fmt("G3")} {fmt("G4")}  {gam:>7.3f}  {kk:>8.3f}')

print(f'\n== 汇总（{len(rows)} 个设计；每个窗内 4~{max(r[2] for r in rows)} 个点）==')
lab = {'G1': 'G1 γ自由·地板Jc (2参)', 'G2': 'G2 γ自由·地板也拟合 (3参)',
       'G3': 'G3 γ=2ln2写死·地板Jc (1参)', 'G4': 'G4 γ=2ln2写死·地板用DARE (1参+1解)'}
best = None
for k in ('G1', 'G2', 'G3', 'G4'):
    a = np.array([r[5][k] for r in rows if r[5][k] is not None])
    m = np.median(np.abs(a))
    print(f'  {lab[k]:<30} n={len(a):>3}  med|e|={100 * m:>7.2f}%  最差={100 * np.max(np.abs(a)):>8.2f}%  '
          f'med有符号={100 * np.median(a):>+8.2f}%  |e|<5%占比={np.mean(np.abs(a) < 0.05):.2f}')
    if best is None or m < best[1]:
        best = (k, m)
print(f'  最好口径 = {best[0]}，med|e|={100 * best[1]:.2f}%；判据要求 ≤ {100 * 0.5 * F1_MED:.3f}%（e127 F1 的一半）')

print('\n== rank 有符号中位（§75-C 的硬约束：不许翻号）==')
flip = {}
for k in ('G1', 'G2', 'G3', 'G4'):
    med_by_rank = []
    for rk in (1, 2, 3):
        a = np.array([r[5][k] for r in rows if r[1] == rk and r[5][k] is not None])
        if len(a) == 0:
            continue
        med_by_rank.append(100 * np.median(a))
        print(f'  {k} rank={rk}: n={len(a):>2}  med={100 * np.median(a):>+7.2f}%  '
              f'负/正={int(np.sum(a < 0))}/{int(np.sum(a > 0))}')
    signs = set(np.sign(v) for v in med_by_rank if v != 0)
    flip[k] = len(signs) > 1
    print(f'  {k}: rank 间是否翻号 = {"翻号（不合格）" if flip[k] else "同号"}')

gams = np.array([r[5]['g1'][1] for r in rows if r[5]['g1'] is not None])
print(f'\n== γ 的可识别性（[G6] 附报）==')
print(f'  γ_fit（G1）: med={np.median(gams):.4f}  IQR=[{np.percentile(gams, 25):.4f},{np.percentile(gams, 75):.4f}]'
      f'  全距=[{gams.min():.4f},{gams.max():.4f}]   理论 2ln2={GAMMA_TH:.4f}')
print(f'  γ_fit 落在 [1.0,2.0] 的占比={np.mean((gams > 1.0) & (gams < 2.0)):.2f}')

dmin_arr = np.array([r[3] for r in rows])
Dtrue_arr = np.array([r[4] for r in rows])
print(f'  该设计自身地板 D_floor−Jc：med={np.median(dmin_arr):.4f}  全距=[{dmin_arr.min():.4f},{dmin_arr.max():.4f}]')
print(f'  同池 D_true−Jc：med={np.median(Dtrue_arr - JC):.4f}  全距=[{(Dtrue_arr - JC).min():.4f},{(Dtrue_arr - JC).max():.4f}]'
      f'   ⇒ 地板占误差量的比例见 G2/G4 相对 G1/G3 的改善')

print('\n== 成本口径 ==')
print('  e127 F1：1 次 DARE 解（近地板取中位）。本实验 G1/G2/G3：窗内 %d~%d 次解 + 非线性最小二乘；'
      % (min(r[2] for r in rows), max(r[2] for r in rows)))
print('  G4：窗内解 + 1 次 s=1e8 的 DARE 解给地板（仍是"少数几次解"量级）。真值 80 次二分。')
print('\n判据 [G6]：')
print(f'  条件A med|e|(最好口径 {best[0]})={100 * best[1]:.2f}% ≤ 15.335% ？{"是" if best[1] <= 0.5 * F1_MED else "否"}')
print(f'  条件B {best[0]} 的 rank 有符号中位不翻号？{"是" if not flip[best[0]] else "否"}')
ok_A = best[1] <= 0.5 * F1_MED
ok_B = not flip[best[0]]
print(f'  ⇒ {"允许替换限制句" if (ok_A and ok_B) else "限制句原样保留，并记录否证"}')
