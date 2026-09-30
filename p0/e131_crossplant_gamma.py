# -*- coding: utf-8 -*-
r"""实验 131（§80-G 预注册）：80-C 的 "$\gamma$ 随 rank 缩水" 只在**一棵锚点植物**上测过。
本实验把它交给跨植物复现，并且**允许它被否证**。

口径与 e130 逐字相同（同设计族、同窗 $\Delta I\in[0.30,1.00]$、同定价 $\Delta I_{\rm tgt}=3-R_{\exp}$、
同真值 80 步二分、同 $D_{\rm floor}= $ s 取 $10^6/10^8/10^{12}$ 三档且要求相对漂移 $\le1\%$），
只换植物。$\Theta,J_c$ 一律由 e128/e130 用的同一个 $\mathrm{ctrl}(Q{=}I,R{=}I)$ 构造重算（不抄别人数）。

 [Y1] 植物集合：锚点（`1510.04214` §V）$+$ 4 株随机 $n{=}4$ 植物（seed 1–4）$+$ 1 株 $n{=}6$、
      含强不稳定实极点与一对不稳定复极点（专门测"$R_{\exp}$ 大"时定律还成不成立）。
      随机植物的 unstable 极点个数**不做筛选**，但逐株印出 $\lambda(A)$ 与 $R_{\exp}$，无不稳定极点的株要明写。
 [Y2] 每株 rank 1/2/3 各 6 个随机正交 $Z$ $+$ $\Theta$min/$\Theta$max 两个极端方向 = 18 设计。
 [Y3] 三列口径都跑：TF2（真地板 $+$ 自由 $\gamma$）、FIX（真地板 $+$ 钉死 $2\ln2$）、HF2（真地板 $+$ 幂律 $\nu$）。
 [Y4] 预注册判据（跑前写死）：记一株植物"支持 80-C" ⟺
      rank-1 的 $\gamma_{\rm TF2}$ 中位 $\in[1.1,1.8]$ **且** rank-3 的中位**严格低于** rank-1 的中位。
      在 4 株新植物里支持数 $\ge3$ ⇒ 80-C 升为"跨植物实测律"；$\le2$ ⇒ 退回单锚点陈述并记为否证。
      另外无论 Y4 结果如何，都要报 TF2 的 $\mathrm{med}\|e\|$ 是否 $\le8\%$（定价器本身是否跨植物成立）。
 [Y5] 不许用"某株植物点数不足"悄悄减样本：逐株报有效设计数，有效数 $<12$ 的株不参与 Y4 计数并明写。
"""
import sys
import numpy as np
from scipy.linalg import solve_discrete_are, eigh
from scipy.optimize import curve_fit
sys.stdout.reconfigure(encoding='utf-8')

_src = open('p0/exp_c_audit.py', encoding='utf-8').read().split("print(r'== E53")[0]
_ns = {'__name__': 'p'}
exec(compile(_src, 'p0/exp_c_audit.py[preamble]', 'exec'), _ns)
A0, B0, W0, ctrl, sym = _ns['A'], _ns['B'], _ns['W'], _ns['ctrl'], _ns['sym']
R = np.eye(4)
R_EXP0, TARGET, GT = 1.168539, 3.0, 2.0 * np.log(2.0)
WLO, WHI = 0.30, 1.00
SGRID = 10.0 ** np.linspace(-6.0, 2.5, 90)
FLOORS = [10.0 ** 6, 10.0 ** 8, 10.0 ** 12]


def make_plant(kind, seed):
    rng = np.random.default_rng(seed)
    if kind == 'anchor':
        A, B, W, m = A0, B0, W0, 4
    elif kind == 'rand4':
        m = 4
        A = rng.normal(0.0, 1.05, (m, m))
        B = rng.normal(0.0, 1.15, (m, m))
        M = rng.normal(0.0, 1.0, (m, m))
        W = sym(M @ M.T + 0.6 * np.eye(m))
    elif kind == 'big6':
        m = 6
        # 块对角：不稳定实极点 2.1、不稳定复极点 1.35e^{±0.9i}、两个稳定块
        bl = [np.array([[2.1]]),
              1.35 * np.array([[np.cos(0.9), -np.sin(0.9)], [np.sin(0.9), np.cos(0.9)]]),
              np.array([[-0.55, 0.3], [-0.2, 0.22]]), np.array([[0.34]])]
        A = np.zeros((m, m)); pos = 0
        for blk in bl:
            k = blk.shape[0]; A[pos:pos + k, pos:pos + k] = blk; pos += k
        A = A + 0.16 * rng.normal(0, 1, (m, m))
        A = sym(A) if False else A
        B = rng.normal(0.0, 1.2, (m, m))
        M = rng.normal(0.0, 1.0, (m, m))
        W = sym(M @ M.T + 0.6 * np.eye(m))
    RR = np.eye(m)
    Pc, K, TH, JC = ctrl_full(A, B, W, np.eye(m), RR)
    return A, B, W, TH, JC, m


def ctrl_full(Ax, Bx, Wx, Qt, Rx):
    Pc = sym(solve_discrete_are(Ax, Bx, sym(Qt), Rx))
    K = np.linalg.solve(Rx + Bx.T @ Pc @ Bx, Bx.T @ Pc @ Ax)
    return Pc, K, sym(K.T @ (Rx + Bx.T @ Pc @ Bx) @ K), float(np.trace(Wx @ Pc))


def r_exp(A):
    mod = np.abs(np.linalg.eigvals(A))
    return float(np.log(mod[mod > 1.0]).sum() / np.log(2.0))


PLANTS = [('anchor', '锚点 1510.04214 §V', ('anchor', 0)),
          ('rand-1', '随机 n=4 seed1', ('rand4', 1)),
          ('rand-2', '随机 n=4 seed2', ('rand4', 2)),
          ('rand-3', '随机 n=4 seed3', ('rand4', 3)),
          ('rand-4', '随机 n=4 seed4', ('rand4', 4)),
          ('big-6', 'n=6 强不稳定(2.1, 1.35e±0.9i)', ('big6', 11))]


def curve(Z, s, A, C_W, TH, JC):
    r = Z.shape[1]
    Ir = np.eye(r)
    C = np.sqrt(s) * Z.T
    Pm = sym(solve_discrete_are(A.T, C.T, C_W, Ir))
    Sm = C @ Pm @ C.T + Ir
    Lk = Pm @ C.T @ np.linalg.inv(Sm)
    Pp = sym(Pm - Lk @ Sm @ Lk.T)
    I = 0.5 * np.log(np.linalg.det(Ir + C @ Pm @ C.T)) / np.log(2.0)
    return I, JC + np.trace(TH @ Pp)


def true_at(Z, A, W, TH, JC, rexp, target):
    lo, hi = -7.0, 5.0
    if not (curve(Z, 10 ** lo, A, W, TH, JC)[0] < target < curve(Z, 10 ** hi, A, W, TH, JC)[0]):
        return None
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if curve(Z, 10 ** mid, A, W, TH, JC)[0] < target:
            lo = mid
        else:
            hi = mid
    return curve(Z, 10 ** (0.5 * (lo + hi)), A, W, TH, JC)[1]


f_exp = lambda di, K, g: K / (np.exp(g * di) - 1.0)

print('== 植物一览（逐株印 $\lambda(A)$ 与 $R_{\exp}$，不做筛选）==')
for lab, desc, spec in PLANTS:
    A, B, W, TH, JC, m = make_plant(*spec)
    mod = np.sort(np.abs(np.linalg.eigvals(A)))[::-1]
    print(f'  {lab:<8} {desc:<30} n={m}  |lam| 排序={np.round(mod,4)}  Rexp={r_exp(A):.6f}  '
          f'JC={JC:.4f}  cond(W)={np.linalg.cond(W):.2f}')

print('\n== 每株：TF2（真地板+自由γ） / FIX（钉2ln2） / HF2（幂律） ==')
verdict = {}
for pi, (lab, desc, spec) in enumerate(PLANTS):
    A, B, W, TH, JC, m = make_plant(*spec)
    REXP = r_exp(A)
    if REXP <= 0:
        REXP = 0.0
    DIT = TARGET - REXP
    _, Vth = eigh(TH)
    rng = np.random.default_rng(20260930 + 7*pi)
    pool = []
    for rk in (1, 2, min(3, m - 1)):
        for c in range(6):
            Z, _ = np.linalg.qr(rng.standard_normal((m, rk)))
            pool.append((rk, Z))
        pool.append((rk, Vth[:, :rk]))
        pool.append((rk, Vth[:, m - rk:]))
    rec = {1: [], 2: [], 3: []}
    err = {'TF2': [], 'FIX': [], 'HF2': []}
    gam = {1: [], 2: [], 3: []}
    nun, bad_floor = 0, 0
    for rk, Z in pool:
        fl = [curve(Z, s, A, W, TH, JC)[1] - JC for s in FLOORS]
        if (max(fl) - min(fl)) / max(min(fl), 1e-12) > 0.01:
            bad_floor += 1
            continue
        dfloor = fl[1]
        Dt = true_at(Z, A, W, TH, JC, REXP, TARGET)
        if Dt is None:
            nun += 1
            continue
        Is, Ds = [], []
        for s in SGRID:
            try:
                I, D = curve(Z, s, A, W, TH, JC)
            except Exception:
                continue
            if np.isfinite(I) and np.isfinite(D):
                Is.append(I); Ds.append(D)
        Is = np.array(Is); Ds = np.array(Ds)
        di = Is - REXP
        ok = (di >= WLO) & (di <= WHI)
        if ok.sum() < 4:
            nun += 1
            continue
        dif, Df = di[ok], Ds[ok] - JC
        try:
            p, _ = curve_fit(f_exp, dif, Df - dfloor, p0=[Df.mean() * dif.mean(), GT],
                             bounds=([1e-12, 1e-3], [np.inf, 20.0]), maxfev=60000)
            err['TF2'].append((JC + dfloor + f_exp(DIT, *p)) / Dt - 1.0)
            gam[rk].append(p[1])
        except Exception:
            continue
        bfun = lambda di, K: K / (np.exp(GT * di) - 1.0)
        base = bfun(dif, 1.0)
        Kfix = float(np.dot(base, Df - dfloor) / np.dot(base, base))
        err['FIX'].append((JC + dfloor + bfun(DIT, Kfix)) / Dt - 1.0)
        try:
            fpow = lambda di, K, nu: dfloor + K / di ** nu
            p, _ = curve_fit(fpow, dif, Df, p0=[Df.mean() * dif.mean(), 1.0],
                             bounds=([1e-12, 0.05], [np.inf, 5.0]), maxfev=60000)
            err['HF2'].append((JC + fpow(DIT, *p)) / Dt - 1.0)
        except Exception:
            pass
    nvalid = len(err['TF2'])
    print(f'\n  [{lab}] 有效设计 {nvalid}/18（I=3 不可达或点不足 {nun}，地板不收敛 {bad_floor}）  ΔI_tgt={DIT:.4f}')
    for k in ('TF2', 'FIX', 'HF2'):
        a = np.array(err[k])
        if len(a):
            print(f'    {k:<5} med|e|={100*np.median(np.abs(a)):>7.2f}%  最差={100*np.max(np.abs(a)):>7.2f}%  '
                  f'有符号med={100*np.median(a):>+7.2f}%')
    gmed = {}
    for rk in (1, 2, 3):
        if gam[rk]:
            gmed[rk] = float(np.median(gam[rk]))
            print(f'    γ_TF2 rank={rk}: n={len(gam[rk])} med={gmed[rk]:.4f} 全距=[{min(gam[rk]):.4f},{max(gam[rk]):.4f}]')
    if nvalid >= 12 and 1 in gmed and 3 in gmed:
        sup = (1.1 <= gmed[1] <= 1.8) and (gmed[3] < gmed[1])
        verdict[lab] = (sup, nvalid, gmed, np.median(np.abs(err['TF2'])))
        print(f'    [Y4] 有效数≥12 ⇒ 支持 80-C？{"是" if sup else "否"}'
              f'（rank1 {gmed[1]:.4f} 落在[1.1,1.8]？{1.1 <= gmed[1] <= 1.8}；rank3 {gmed[3]:.4f} < rank1？{gmed[3] < gmed[1]}）')
    else:
        verdict[lab] = None
        print(f'    [Y5] 有效数 {nvalid}<12 或 rank 不全 ⇒ 不参与 Y4 计数')

print('\n== [Y4] 跨植物判决 ==')
newp = [k for k in verdict if k != 'anchor' and verdict[k] is not None]
sup = [k for k in newp if verdict[k][0]]
print(f'  新植物（排除锚点）参与计数 {len(newp)} 株：{newp}')
print(f'  支持 80-C 的株数 = {len(sup)} / {len(newp)} ⇒ '
      f'{"≥3 ⇒ 升为跨植物实测律" if len(sup) >= 3 else "≤2 ⇒ 退回单锚点陈述，并记为否证"}')
allm = [verdict[k][3] for k in newp if verdict[k] is not None]
print(f'  TF2 各株 med|e|：{[round(100*v,2) for v in allm]}  ⇒ 全部 ≤8%？'
      f'{"是 ⇒ 定价器跨植物成立" if allm and max(allm) <= 0.08 else "否 ⇒ 定价器本身也要打植物限定"}')
