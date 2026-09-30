# -*- coding: utf-8 -*-
r"""实验 130：§79-E 留的新限制是"把 G2 救回来的那个自由地板 δ 不是 (43) 的 D_min"，
因为 G4（γ=2ln2 写死 + 真实地板）退回 15.42% 且 rank 翻号。但 G4 同时钉死了 γ，
所以它**不能**区分两种解释：
   (H1) 地板的**身份**错了（拟合出的 δ 只是个截距，和 s→∞ 的真地板不是一回事）；
   (H2) 地板身份没错，错在 γ 被钉死（真地板 + 自由 γ 就够）。
本实验只放开 γ、把地板钉在**实测真值**上，直接判 H1/H2，并把 δ_fit 与真地板逐设计对齐比较。

 [Z1] 口径：拟合窗 $\Delta I\in[0.30,1.00]$（同 e128），定价 $\Delta I=1.831461$，真值 80 步二分。
      TF2 = 地板取实测 $D_{\rm floor}$（$s=10^{8}$ 一次 DARE）、$(K,\gamma)$ 自由拟合（2 参）。
      EF3 = 地板也自由（3 参，= e128 的 G2）。
      HF2 = 地板取 $D_{\rm floor}$、幂律 $K/\Delta I^\nu$（2 参）—— 检验"是否只有指数式挑食"。
      另加 FIX：$\gamma\equiv2\ln2$ + 真地板（= e128 的 G4，复现用，应当 $=15.42\%$）。
 [Z2] 收敛性自证：$D_{\rm floor}$ 用 $s=10^{6},10^{8},10^{12}$ 三档算，逐设计报最大相对漂移；
      漂移 $>1\%$ 的真地板不可用于判决（说明极限没取到，H1/H2 都判不了）。
 [Z3] $\delta_{\rm fit}$ vs $D_{\rm floor}$：逐设计比值，报 rank 分层的中位比值与散布。
 [Z4] 判据（跑前写死）：
      若 TF2 的 $\mathrm{med}\|e\|\le 8\%$ **且** rank 有符号中位同号 ⇒ 判 **H2**：
         "闭式的 $D_{\min}$ 就是缺失那一块"可写进正文，§79-E 的不可解释限制撤销；
      否则若 $\delta_{\rm fit}/D_{\rm floor}$ 的中位偏离 1 超过 30% ⇒ 判 **H1**：
         $\delta$ 只是截距，正文必须写 "with $\delta$ left free"（§79-F 现措辞保持）；
      两者都成立则报"身份成立但需同时放开 $\gamma$"；都不成立 ⇒ 明写"两问都否证，另找机制"。
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
R_EXP, TARGET = 1.168539, 3.0
DIT = TARGET - R_EXP
GT = 2.0 * np.log(2.0)
WLO, WHI = 0.30, 1.00
SGRID = 10.0 ** np.linspace(-6.0, 2.5, 90)
FLOORS = [10.0 ** 6, 10.0 ** 8, 10.0 ** 12]


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


_, Vth = eigh(TH)
rng = np.random.default_rng(20260930)
pool = []
for rk in (1, 2, 3):
    for c in range(10):
        Z, _ = np.linalg.qr(rng.standard_normal((n, rk)))
        pool.append((f'r{rk}-{c:02d}', rk, Z))
    pool.append((f'THmin-{rk}', rk, Vth[:, :rk]))
    pool.append((f'THmax-{rk}', rk, Vth[:, n - rk:]))


def f_exp(di, K, g):
    return K / (np.exp(g * di) - 1.0)


def f_exp_free(di, K, g, d):
    return d + K / (np.exp(g * di) - 1.0)


def f_pow_floor(di, K, nu, d):
    return d + K / di ** nu


print('== Z2 真地板的收敛性（三档 s）==')
rows = []
for name, rk, Z in pool:
    fl = []
    for s in FLOORS:
        try:
            fl.append(curve(Z, s)[1] - JC)
        except Exception:
            fl.append(np.nan)
    fl = np.array(fl, dtype=float)
    drift = (np.nanmax(fl) - np.nanmin(fl)) / max(np.nanmin(fl), 1e-12)
    okc = drift <= 0.01
    print(f'  {name:<9} rk={rk}  Dfloor(1e6/1e8/1e12)-Jc='
          f'{fl[0]:>10.4f} {fl[1]:>10.4f} {fl[2]:>10.4f}  相对漂移={drift:>7.4f}'
          f'  {"可用" if okc else "**不可用**（>1%）"}')
    rows.append([name, rk, fl[1], okc, None, None, None, None, None])

print('\n== Z1/Z3 拟合（窗 ΔI∈[%.2f,%.2f]，定价 %.6f）==' % (WLO, WHI, DIT))
print('  design    rk  nwin   D_true    Dfloor    δ_fit    δ/Dfl   gam_fit   e(TF2)%  e(EF3)%  e(HF2)%  e(FIX)%')
for row in rows:
    name, rk, dfloor, okc = row[0], row[1], row[2], row[3]
    Z = [Z for nm, _, Z in pool if nm == name][0]
    Dtrue = true_at(Z)
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
    ok = (di >= WLO) & (di <= WHI)
    dif, Df = di[ok], Ds[ok] - JC
    tf2 = ef3 = hf2 = fix = np.nan
    delta_fit = g_fit = g_tf2 = row_ratio = np.nan
    if ok.sum() >= 4 and okc:
        try:
            p, _ = curve_fit(f_exp, dif, Df - dfloor, p0=[Df.mean() * dif.mean(), GT],
                             bounds=([1e-12, 1e-3], [np.inf, 20.0]), maxfev=60000)
            tf2 = (JC + dfloor + f_exp(DIT, *p)) / Dtrue - 1.0
            g_tf2 = p[1]
        except Exception:
            pass
        try:
            p, _ = curve_fit(f_exp_free, dif, Df, p0=[Df.mean() * dif.mean(), GT, 0.05 * Df.mean()],
                             bounds=([1e-12, 1e-3, 0.0], [np.inf, 20.0, np.inf]), maxfev=60000)
            ef3 = (JC + f_exp_free(DIT, *p)) / Dtrue - 1.0
            delta_fit, g_fit = p[2], p[1]
        except Exception:
            pass
        try:
            fpow = lambda di, K, nu: dfloor + K / di ** nu
            p, _ = curve_fit(fpow, dif, Df, p0=[Df.mean() * dif.mean(), 1.0],
                             bounds=([1e-12, 0.05], [np.inf, 5.0]), maxfev=60000)
            hf2 = (JC + fpow(DIT, *p)) / Dtrue - 1.0
        except Exception:
            pass
        try:
            # 复现 e128 的 G4：γ≡2ln2 写死 + 真实地板（地板必须在定价时加回来）
            bfun = lambda di, K: K / (np.exp(GT * di) - 1.0)
            base = bfun(dif, 1.0)
            Kfix = float(np.dot(base, Df - dfloor) / np.dot(base, base))
            fix = (JC + dfloor + bfun(DIT, Kfix)) / Dtrue - 1.0
        except Exception:
            pass
    row[4:] = [delta_fit, g_fit, tf2, ef3, hf2]
    row.append(g_tf2)
    row.append(fix)
    row_ratio = delta_fit / dfloor if (delta_fit == delta_fit and dfloor > 0) else np.nan
    row.append(row_ratio)
    print(f'  {name:<9} {rk:>3}  {ok.sum():>4}  {Dtrue:>9.3f} {JC + dfloor:>9.3f} {delta_fit:>8.3f} '
          f'{row_ratio:>7.3f} {g_fit:>8.3f} {100*tf2:>8.2f} {100*ef3:>9.2f} {100*hf2:>8.2f} {100*fix:>7.2f}')

use = [r for r in rows if r[3]]
print(f'\n== 汇总（可用地板的设计 {len(use)}/{len(rows)}）==')


def col(i):
    return np.array([r[i] for r in use if r[i] == r[i]])


names = {4: 'δ_fit', 5: 'γ_fit(EF3)', 6: 'e(TF2) 真地板+自由γ', 7: 'e(EF3) 地板也自由',
         8: 'e(HF2) 真地板+幂律', 10: 'e(FIX)', 11: 'δ/Dfloor'}
for i in (6, 7, 8):
    a = col(i)
    if len(a):
        print(f'  {names[i]:<26} n={len(a):>3}  med|e|={100*np.median(np.abs(a)):>7.2f}%  '
              f'最差={100*np.max(np.abs(a)):>8.2f}%  有符号med={100*np.median(a):>+7.2f}%')
a = col(10)
if len(a):
    print(f'  {"e(FIX) γ=2ln2+真地板 (=e128 G4 复现)":<40} n={len(a):>3}  '
          f'med|e|={100*np.median(np.abs(a)):>7.2f}%（e128 G4 报 15.42%）')
rat = col(11)
if len(rat):
    print(f'  δ_fit/D_floor: med={np.median(rat):.3f}  IQR=[{np.percentile(rat,25):.3f},{np.percentile(rat,75):.3f}]'
          f'  全距=[{rat.min():.3f},{rat.max():.3f}]  |ln 比|>ln1.3 占比={np.mean(np.abs(np.log(rat))>np.log(1.3)):.2f}')
gam = col(5)
if len(gam):
    print(f'  EF3 γ_fit: med={np.median(gam):.4f}  IQR=[{np.percentile(gam,25):.4f},{np.percentile(gam,75):.4f}]  2ln2={GT:.4f}')
gt2 = col(9)
if len(gt2):
    print(f'  TF2 γ_fit: med={np.median(gt2):.4f}  IQR=[{np.percentile(gt2,25):.4f},{np.percentile(gt2,75):.4f}]'
          f'  全距=[{gt2.min():.4f},{gt2.max():.4f}]  贴下界(<0.01)占比={np.mean(gt2 < 0.01):.2f}'
          f'  落在[1.0,2.0]占比={np.mean((gt2 > 1.0) & (gt2 < 2.0)):.2f}')
    for rk in (1, 2, 3):
        a = np.array([r[9] for r in use if r[1] == rk and r[9] == r[9]])
        if len(a):
            print(f'    TF2 γ rank={rk}: med={np.median(a):.4f}  全距=[{a.min():.4f},{a.max():.4f}]')

print('\n== rank 有符号中位（H1/H2 判决的同号要求）==')
for i, lab in ((6, 'TF2'), (7, 'EF3'), (10, 'FIX')):
    meds = []
    for rk in (1, 2, 3):
        a = np.array([r[i] for r in use if r[1] == rk and r[i] == r[i]])
        if len(a) == 0:
            continue
        meds.append(100 * np.median(a))
        print(f'  {lab} rank={rk}: n={len(a):>2} med={100*np.median(a):>+7.2f}%  负/正={int(np.sum(a<0))}/{int(np.sum(a>0))}')
    signs = set(np.sign(v) for v in meds if v != 0)
    print(f'  {lab}: 同号={"是" if len(signs)<=1 else "否（翻号）"}')

print('\n== [Z4] 判决 ==')
tf2m = np.median(np.abs(col(6))) if len(col(6)) else np.nan
ef3m = np.median(np.abs(col(7))) if len(col(7)) else np.nan
med_by = [np.median(np.array([r[6] for r in use if r[1] == rk and r[6] == r[6]]))
          for rk in (1, 2, 3)]
med_by = [m for m in med_by if m == m]
same_sign = len(set(np.sign(m) for m in med_by if m != 0)) <= 1
H2 = (tf2m == tf2m) and tf2m <= 0.08 and same_sign
rr = np.median(np.abs(col(11))) if len(col(11)) else np.nan
H1 = (rr == rr) and abs(rr - 1.0) > 0.30
print(f'  TF2 med|e|={100*tf2m:.2f}%（判据 ≤8%？{"是" if H2 else "否"}）  rank 同号？{"是" if same_sign else "否"}')
print(f'  δ_fit/D_floor 中位={rr:.3f}（偏离 1 超 30%？{"是" if H1 else "否"}）')
print(f'  ⇒ H2（真地板可用、错在钉死 $\gamma$）：{"成立" if H2 else "不成立"}；'
      ' H1：delta_fit 是 D_floor 的好估计子？' + ('否' if H1 else '是') + f'（中位比 {rr:.3f}）')
if H2 and not H1:
    print('  ⇒ 只判 H2。')
elif H2 and H1:
    print('  ⇒ **合成读法**：地板的**身份成立**（用 $s\\to\\infty$ 的真地板比拟合地板更好：2.07% vs 4.35%），'
          '但 $\delta_{\rm fit}$ **不是**真地板的好估计子（中位比 1.389，64% 偏差 $>30\\%$）'
          '⇒ 正文必须写"代入可计算的真地板 + 放开 $\gamma$"，而不是"拟合一个自由截距"。')
elif not H2 and H1:
    print('  ⇒ 判 H1：$\delta$ 只是截距，$\gamma$ 也救不回来。')
else:
    print('  ⇒ 两问都否证：另找机制（下一嫌疑：秩亏把 $D_{\min}$ 变成矩阵值，标量地板假设本身失效）')
print(f'  附加（$\gamma$ 随 rank）：rank-1 med={np.median(np.array([r[9] for r in use if r[1]==1])):.4f}  '
      f'rank-2 med={np.median(np.array([r[9] for r in use if r[1]==2])):.4f}  '
      f'rank-3 med={np.median(np.array([r[9] for r in use if r[1]==3])):.4f}  （$2\ln2=1.3863$）')
tf2s = np.array([r[6] for r in use if r[6] == r[6]])
print(f'  反向偏置警示：TF2 有符号 med={100*np.median(tf2s):+.2f}%、最差 {100*tf2s.min():+.2f}% ⇒ '
      f'预测代价**低于**真值的比例即"非保守"方向，若要当下界证书用须另加裕度')
