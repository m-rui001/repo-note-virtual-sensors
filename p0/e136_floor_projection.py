# -*- coding: utf-8 -*-
r"""实验 136（§85-F-2 / §86-G-2 预注册）：我的地板是 $s=10^8$ 的 DARE 值，C 的 c73 用的是**投影极限**
（$P\leftarrow P-PQ(Q^\top PQ)^{-1}Q^\top P$，$P\leftarrow APA^\top+W$ 的不动点，$\mathrm{tr}(\Theta P)$）。
两口径的差会不会改我 §81-B/§83-B 的读数？本轮只回答这一问。
投影极限是**我自己重写**的实现（只读过 C 的定义，未运行 C 的脚本）。

 [X1] 三口径同池并列：`dfl_8`（$s=10^8$，我现口径）、`dfl_12`（$s=10^{12}$）、`dfl_proj`（投影极限）。
      逐设计报三值、相对散布，以及 $(\text{proj}-8)/\text{proj}$。池与 e131 逐字相同（6 株 $\times$ 18 设计）。
 [X2] 判据（跑前写死）：记 $\rho=$ 全体设计 $(\text{proj}-8)/\text{proj}$ 的中位绝对值。
      $\rho\le1\%$ ⇒ 我的三档 $s$ 门够用，**不需要**重跑 §81-B/§83-B，只需在正文写明两口径差；
      $\rho>1\%$ ⇒ **必须**用投影地板重跑 TF2/FIX/HF2 全表并逐列报漂移（本轮直接一并跑出来）。
 [X3] 双向自证：投影极限必须与 $s=10^{12}$ 一致（相对差中位 $\le0.5\%$）。不一致 ⇒ 报**谁偏了**、
      偏多少，并明写"要么我的 $s$ 档没取到极限，要么定点迭代没收敛"——不许默认自己是对的。
 [X4] 顺带对照（不是新判据）：换地板后 e131 的 $\gamma_r/\gamma_1>1/r$（§85-D 六路不等式）还成不成立。
      若 6/6 仍成立 ⇒ 不等式对地板口径稳健；若有株翻转 ⇒ §85-D 的英文句必须加口径限定。
 [X5] 不许悄悄改判据、不许用"某设计点数不足"悄悄减样本：逐株报有效数与剔除原因。
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
TARGET, GT = 3.0, 2.0 * np.log(2.0)
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
        bl = [np.array([[2.1]]),
              1.35 * np.array([[np.cos(0.9), -np.sin(0.9)], [np.sin(0.9), np.cos(0.9)]]),
              np.array([[-0.55, 0.3], [-0.2, 0.22]]), np.array([[0.34]])]
        A = np.zeros((m, m)); pos = 0
        for blk in bl:
            k = blk.shape[0]; A[pos:pos + k, pos:pos + k] = blk; pos += k
        A = A + 0.16 * rng.normal(0, 1, (m, m))
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


def dare_floor(Z, s, A, W, TH, JC):
    return curve(Z, s, A, W, TH, JC)[1] - JC


def proj_floor(Z, A, W, TH, itmax=60000, tol=1e-16):
    """s→∞ 极限的直接不动点：永远以无穷精度测 span(Z)。
    返回 tr(Θ P_post)（**投影后**，与 curve() 里 Pp 同口径）；发散返回 inf。"""
    Q, _ = np.linalg.qr(Z)
    Pt = sym(W.copy())
    P = Pt
    for _ in range(itmax):
        P = sym(Pt - Pt @ Q @ np.linalg.solve(Q.T @ Pt @ Q, Q.T @ Pt))
        Pn = sym(A @ P @ A.T + W)
        done = np.max(np.abs(Pn - Pt)) < tol * max(1.0, np.max(np.abs(Pn)))
        Pt = Pn
        if done:
            break
        if not np.all(np.isfinite(Pt)) or np.max(np.abs(Pt)) > 1e14:
            return np.inf
    return float(np.trace(TH @ P))


def true_at(Z, A, W, TH, JC, target):
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
bfun = lambda di, K: K / (np.exp(GT * di) - 1.0)
fpow = lambda di, K, nu: K / di ** nu

allrel = []           # [X1] (proj-8)/proj 全体设计
rows = []             # 逐设计记录（供汇总）
gam_tbl = {}          # [X4] {lab: {rk: {'old': [γ], 'new': [γ]}}}
drift = {'TF2': [], 'FIX': [], 'HF2': []}   # [X2] 逐设计读数漂移（绝对百分点）

for pi, (lab, desc, spec) in enumerate(PLANTS):
    A, B, W, TH, JC, m = make_plant(*spec)
    REXP = max(r_exp(A), 0.0)
    DIT = TARGET - REXP
    _, Vth = eigh(TH)
    rng = np.random.default_rng(20260930 + 7 * pi)
    pool = []
    for rk in (1, 2, min(3, m - 1)):
        for c in range(6):
            Z, _ = np.linalg.qr(rng.standard_normal((m, rk)))
            pool.append((rk, Z))
        pool.append((rk, Vth[:, :rk]))
        pool.append((rk, Vth[:, m - rk:]))
    gam_tbl[lab] = {1: {'old': [], 'new': []}, 2: {'old': [], 'new': []}, 3: {'old': [], 'new': []}}
    nun = 0
    ndiv = 0
    badfl = 0
    print('\n== [%s] ΔI_tgt=%.6f  地板三口径（相对差 = (proj-8)/proj）==' % (lab, DIT))
    print('  rk  dfl_8        dfl_12       dfl_proj     proj-8/proj   12-8/12    收敛')
    for rk, Z in pool:
        f6, f8, f12 = [dare_floor(Z, s, A, W, TH, JC) for s in FLOORS]
        fp = proj_floor(Z, A, W, TH)
        conv = np.isfinite(fp)
        if not conv:
            ndiv += 1
        rel_p8 = (fp - f8) / fp if conv else np.nan
        rel_128 = (f12 - f8) / f12
        if conv:
            allrel.append(abs(rel_p8))
        if (max(f6, f8, f12) - min(f6, f8, f12)) / max(min(f6, f8, f12), 1e-12) > 0.01:
            badfl += 1
        rows.append((lab, rk, f8, f12, fp, rel_p8))
        print('  %d   %-12.8f %-12.8f %-12.8f %+11.3f%%  %+9.3f%%  %s'
              % (rk, f8, f12, fp if conv else float('nan'), 100 * rel_p8, 100 * rel_128,
                 'yes' if conv else 'DIVERGE'))
        if not conv:
            continue
        Dt = true_at(Z, A, W, TH, JC, TARGET)
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
        two = []
        for tag, dfloor in (('old', f8), ('new', fp)):
            try:
                p, _ = curve_fit(f_exp, dif, Df - dfloor, p0=[Df.mean() * dif.mean(), GT],
                                 bounds=([1e-12, 1e-3], [np.inf, 20.0]), maxfev=60000)
                e_tf2 = (JC + dfloor + f_exp(DIT, *p)) / Dt - 1.0
            except Exception:
                two.append(None); continue
            Kfix = float(np.dot(bfun(dif, 1.0), Df - dfloor) / np.dot(bfun(dif, 1.0), bfun(dif, 1.0)))
            e_fix = (JC + dfloor + bfun(DIT, Kfix)) / Dt - 1.0
            try:
                pp, _ = curve_fit(lambda x, K, nu: dfloor + K / x ** nu, dif, Df,
                                  p0=[Df.mean() * dif.mean(), 1.0],
                                  bounds=([1e-12, 0.05], [np.inf, 5.0]), maxfev=60000)
                e_hf2 = (JC + fpow(DIT, *pp)) / Dt - 1.0
            except Exception:
                e_hf2 = None
            two.append((p[1], e_tf2, e_fix, e_hf2))
            gam_tbl[lab][rk][tag].append(p[1])
        if two[0] and two[1]:
            for k, idx in (('TF2', 1), ('FIX', 2), ('HF2', 3)):
                if two[0][idx] is not None and two[1][idx] is not None:
                    drift[k].append(100 * abs(two[1][idx] - two[0][idx]))
    print('  [X5] %s：投影发散 %d / 真值或点数不足 %d / e131 地板门会剔 %d / 本轮计价设计 %d'
          % (lab, ndiv, nun, badfl,
             len(gam_tbl[lab][1]['old']) + len(gam_tbl[lab][2]['old']) + len(gam_tbl[lab][3]['old'])))

ar = np.array(allrel)
print('\n== [X1][X2] 两口径差的总账（全 %d 设计）==' % len(ar))
print('  |(proj-8)/proj|  中位=%.4f%%  p90=%.4f%%  最差=%.4f%%  >1%% 的株数=%d'
      % (100 * np.median(ar), 100 * np.percentile(ar, 90), 100 * np.max(ar), int((ar > 0.01).sum())))
r12 = np.array([abs(x) for x in [100 * (rows[i][4] - rows[i][3]) / rows[i][3] if np.isfinite(rows[i][3]) and np.isfinite(rows[i][4]) and rows[i][3] != 0 else np.nan for i in range(len(rows))] if np.isfinite(x)])
print('  [X3] 反向核对 (proj-12)/12：中位=%.4f%% 最差=%.4f%%  n=%d'
      % (np.median(r12), np.max(r12), len(r12)))
print('  [X3] 判读：投影 vs $10^{12}$ 中位差 >0.5%% ⇒ %s'
      % ('**两口径不一致，要报谁偏了**' if np.median(r12) > 0.5 else '两口径一致 ⇒ 我的 s 档取到了极限'))
gate = np.median(ar) <= 0.01
print('  [X2] 判决：ρ=%.4f%% ⇒ %s' % (100 * np.median(ar),
      'ρ≤1% ⇒ 三档 s 门够用，§81-B/§83-B 读数不重跑（下表仅作对照）' if gate
      else 'ρ>1% ⇒ 地板口径会改读数，必须按投影地板重报（下表即为重报值）'))
for k in ('TF2', 'FIX', 'HF2'):
    a = np.array(drift[k])
    if len(a):
        print('    %-4s 逐设计读数漂移 |新-旧|：中位=%.3f 百分点  最差=%.3f 百分点' % (k, np.median(a), np.max(a)))

print('\n== [X4] 换地板后 γ_r/γ_1 > 1/r 还成立吗（§85-D 的稳健性对照，非新判据）==')
print('  株        γ1旧    γ1新    γ2旧/γ1旧 γ2新/γ1新   >1/2?   γ3旧/γ1旧 γ3新/γ1新   >1/3?')
hold_old, hold_new = [], []
for lab, _, _ in PLANTS:
    g = gam_tbl[lab]
    if not all(len(g[rk]['old']) for rk in (1, 2, 3)):
        print('  %-8s 秩不全（n=%d/%d/%d）⇒ 跳过并明写' % (lab, len(g[1]['old']), len(g[2]['old']), len(g[3]['old'])))
        continue
    o = {rk: float(np.median(g[rk]['old'])) for rk in (1, 2, 3)}
    nn = {rk: float(np.median(g[rk]['new'])) for rk in (1, 2, 3)}
    ro2, rn2 = o[2] / o[1], nn[2] / nn[1]
    ro3, rn3 = o[3] / o[1], nn[3] / nn[1]
    hold_old.append(ro2 > 0.5 and ro3 > 1 / 3)
    hold_new.append(rn2 > 0.5 and rn3 > 1 / 3)
    print('  %-8s %7.4f %7.4f   %7.4f %7.4f  %-6s   %7.4f %7.4f  %-6s'
          % (lab, o[1], nn[1], ro2, rn2, 'yes' if ro2 > .5 and rn2 > .5 else 'NO',
             ro3, rn3, 'yes' if ro3 > 1 / 3 and rn3 > 1 / 3 else 'NO'))
print('  旧口径成立 %d/%d；新口径成立 %d/%d ⇒ %s'
      % (sum(hold_old), len(hold_old), sum(hold_new), len(hold_new),
         '两口径都 6/6 ⇒ §85-D 不等式对地板口径稳健'
         if sum(hold_new) == len(hold_new) == len(hold_old) == sum(hold_old)
         else '**有翻转 ⇒ §85-D 的英文句必须加口径限定**'))
print('EXIT=0')
