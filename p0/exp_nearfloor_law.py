"""E23：近地板发散到底是幂律还是对数——把三条互斥读数一次讲清。

背景：#9 我按 E13 的最左两点拟合出 p≈0.53 并判它"网格伪值"；C 的 c16 稠密网格在同一 x 窗口给
−0.516；两边各自的全局拟合又都朝 −1.0（C: −1.052/−1.118；我: −1.003）。三个数不能同时是极限指数。

E22 的机制（本脚本要独立复验并推广）：$\min_{\mathcal K_F}J=$ 地板只在 $\Gamma\to\infty$ 的极限里达到
（E20：锥上只看 J 的下确界精确等于 $D_{\min}^{unc}+\Phi$），所以 $D\downarrow$ 地板必然 $\gamma_j\to\infty$、
$\Delta\to\infty$。若后验误差按 $P(\gamma)-P_\infty=\sum_j c_j/\gamma_j+o(1/\gamma_{\min})$ 收敛（$c_j>0$），
则
    $\min\ \tfrac12\sum_j\log_2\gamma_j\ \ \mathrm{s.t.}\ \ \sum_j c_j/\gamma_j\le x$
的 KKT 解是 $\gamma_j=c_j\,r/x$（所有通道一起发散，比例固定），代回得

    **Prop E（对数律）** $\Delta(D)=\dfrac{r}{2}\log_2\dfrac1x+O(1)$，$x=D-D_{\min}(\mathcal K_F)$。

系数 $r/2$ 只依赖 $\mathrm{rank}F=r$，与 $A,Q,F$ 的取向、$\Phi$ 都无关；这些都进 $O(1)$ 截距。
可检验的两件事：
  (A) 各向同性分支（$\Gamma=\gamma I_r$，无需优化）：F2 应给 $3.3219$ bit/decade，F3 应给 $4.9829$。
  (B) F2 真最优（$\Gamma$ 自由）：斜率应与 (A) 同（同系数），只差一个常数截距；且最优 $\Gamma$ 的两个
      特征值之比应随 $x\to0$ 趋于常数（幂次相同、比例固定）。
若 (B) 的斜率显著小于 $3.32$，说明有通道 $c_j=0$（该方向免费），有效系数按活跃通道数取。

数值边界：$\gamma\gtrsim10^{7}$ 后 `ss_filter` 的双精度地板让 $x$ 掉头回升（E22 已见），拟合窗口一律
限制在 $x$ 单调段的 $[10^{-5},0.3]$。
"""
import sys
import numpy as np
from scipy.linalg import schur
from scipy.optimize import minimize
sys.stdout.reconfigure(encoding='utf-8')
from p0.p0_replicate_letter import (A, B, W, n, sym, ln2, ctrl, rate_cost, unconstrained_sdp)
from p0.exp_authoritative import noiseless

Pc, K, Th = ctrl(np.eye(n))
UNC_FLOOR = float(np.trace(W @ Pc))
Tsch, U, sdim = schur(A, output='real', sort=lambda a: np.abs(a) < 1.0)
nu = n - sdim
iT = list(range(n - nu, n))
Fs = {'F2 尾部2(封闭, r=2)': U[:, iT].T,
      'F3 尾部3(不封闭, r=3)': U[:, [1, 2, 3]].T}
Z = np.ascontiguousarray(U[:, :sdim])


def floor_of(F):
    P, res, it = noiseless(F, A, W, iters=80000, tol=1e-15)
    return UNC_FLOOR + float(np.trace(Th @ P)), float(np.trace(Th @ P)), res


# ---------- 参数化：Gamma = L L^T，L 下三角、对角写成 log，覆盖全部 SPD 且不限制量级 ----------
# p 的前 r 个 = log diag(L)，其后 = L 的严格下三角（按行）。r=2 时即 [log L00, log L11, L10]，与旧版逐字一致。
def _rank_of(p):
    return int(round((np.sqrt(1 + 8 * len(p)) - 1) / 2))


def unpack(p):
    r = _rank_of(p)
    M = np.zeros((r, r))
    M[np.diag_indices(r)] = np.exp(p[:r])
    M[np.tril(np.ones((r, r), dtype=bool), -1)] = p[r:]
    return sym(M @ M.T)


def pack(G):
    r = G.shape[0]
    L = np.linalg.cholesky(G + 1e-12 * np.eye(r))
    return np.concatenate([np.log(np.diag(L)), L[np.tril(np.ones((r, r), dtype=bool), -1)]])


def ij(F, p):
    try:
        I, J, _, _ = rate_cost(F, unpack(p), Th, Pc)
    except Exception:
        return 1e6, 1e6
    if not np.isfinite(I) or not np.isfinite(J):
        return 1e6, 1e6
    return I, J


def solve_cone(F, D, jmin, starts, budget=4000):
    """min I s.t. J<=D，两段式（罚函数 Nelder-Mead → SLSQP 硬约束）。返回 (I,J,Gamma)。
    jmin = 地板：低于它的一定是数值垃圾（E23 首轮在 x=3e-4 抓到 J=-7.6e9）。"""
    best = None
    for p0 in starts:
        pen = lambda p: (lambda I, J: I + 2000. * max(0., (J - D) / D) ** 2
                         + 200. * max(0., (J - D) / D))(*ij(F, p))
        r1 = minimize(pen, p0, method='Nelder-Mead',
                      options=dict(maxiter=budget, xatol=1e-10, fatol=1e-12))
        r2 = minimize(lambda p: ij(F, p)[0], r1.x, method='SLSQP',
                      constraints=[dict(type='ineq', fun=lambda p: D - ij(F, p)[1])],
                      options=dict(maxiter=500, ftol=1e-12))
        for cand in (r1.x, r2.x):
            I, J = ij(F, cand)
            if J <= D * (1 + 1e-9) and J > jmin - 1e-6 and (best is None or I < best[0]):
                best = (I, J, unpack(cand))
    return best


def gamma_iso(F, D, lo=-2., hi=13.):
    """二分求各向同性 $\gamma$ 使 $J(\gamma I_r)=D$；给出一个保证可行的起点。"""
    f = lambda lg: rate_cost(F, np.exp(lg) * np.eye(F.shape[0]), Th, Pc)[1] - D
    flo, fhi = f(lo), f(hi)
    if flo * fhi > 0:
        return None
    for _ in range(200):
        mid = .5 * (lo + hi)
        if f(mid) * flo <= 0:
            hi = mid
        else:
            lo = mid
    return .5 * (lo + hi)


def logfit(xs, ds, tag, pred):
    xs, ds = np.array(xs), np.array(ds)
    lx = np.log10(1. / xs)
    a, b = np.polyfit(lx, ds, 1)                       # Δ = a·log₁₀(1/x) + b
    ld = np.log2(ds)
    a2, b2 = np.polyfit(lx, ld, 1)                     # log₂Δ = a2·log₁₀(1/x)+b2 → 幂律指数
    resid = np.max(np.abs(ds - (a * lx + b)))
    print('  %-22s 对数律拟合 Δ = a·log₁₀(1/x)+b：a=%+.4f (预言 %+.4f)  b=%+.4f  最大残差=%.3f  点=%d'
          % (tag, a, pred, b, resid, len(xs)))
    print('  %-22s 同窗口幂律拟合 log₂Δ 对 log₁₀(1/x) 斜率=%+.4f（≈0 即对数律成立，≈p 即幂律）'
          % ('', a2))
    return a, b


if __name__ == '__main__':
    print('D_min^unc=%.4f  |λ|=%s  sdim=%d' % (UNC_FLOOR, np.abs(np.linalg.eigvals(A)).round(4), sdim))
    floors = {}
    for tag, F in Fs.items():
        fl, phi, res = floor_of(F)
        floors[tag] = fl
        print('[%s] 地板=%.5f  Φ=%.5f  定点残差=%.1e' % (tag, fl, phi, res))

    # ---------------- (A) 各向同性分支：无需优化，直接扫 γ ----------------
    print('\n===== (A) 各向同性分支 Γ=γI_r：对数律系数是否只依赖 r =====')
    for tag, F in Fs.items():
        r = F.shape[0]
        rows = []
        for lg in np.arange(-0.5, 12.6, 0.25):
            G = np.exp(lg) * np.eye(r)
            I, J, _, _ = rate_cost(F, G, Th, Pc)
            x = J - floors[tag]
            if x > 1e-7:
                rows.append((lg, x, J, I))
        xs = [b[1] for b in rows]
        turn = next((i for i in range(len(xs) - 1) if xs[i + 1] >= xs[i]), None)
        rows = rows[:turn] if turn else rows          # 掉头之后是数值垃圾，整段丢弃
        win = [b for b in rows if 1e-5 <= b[1] <= 0.3]
        print('\n[%s] x 单调段到 γ=%.1e（其后掉头=双精度底），拟合窗口 %d 点 x∈[%.1e,%.1e] = %.2f 个十倍程'
              % (tag, np.exp(rows[turn - 1][0]) if turn else np.inf, len(win),
                 win[-1][1], win[0][1], np.log10(win[0][1] / win[-1][1])))
        pairs = [(b[1], b[3] - unconstrained_sdp(b[2], np.eye(n), Th, Pc)['I_true']) for b in win]
        logfit([p for p, _ in pairs], [q for _, q in pairs],
               '各向同性 r=%d' % r, (r / 2.) * np.log2(10.))

    # ---------------- (B) F2 真最优：Γ 自由，热启动链 ----------------
    print('\n===== (B) F2 真最优（Γ 自由）：最优分支与各向同性分支同斜率、只差常数截距？ =====')
    print('  质检：各向同性点在锥上必然可行 ⟹ Δ_opt ≤ Δ_iso 是硬上界。违反即未收敛，剔除后再拟合。')
    F = Fs['F2 尾部2(封闭, r=2)']
    fl = floors['F2 尾部2(封闭, r=2)']
    prev_x, prev_G = None, None
    good = []
    for x in (2.0, 1.4, 0.8, 0.377, 0.2, 0.1, 0.05, 0.03, 0.02, 0.01, 0.005, 0.003, 0.002, 0.001):
        D = fl + x
        Iu = unconstrained_sdp(D, np.eye(n), Th, Pc)['I_true']
        lg = gamma_iso(F, D)
        starts, I_iso = [], np.nan
        if lg is not None:
            g = np.exp(lg)
            I_iso = rate_cost(F, g * np.eye(2), Th, Pc)[0] - Iu
            starts += [pack(np.diag([g, g]))]
        if prev_G is not None:                      # KKT 结构热启动：γ_j = c_j·r/x ⟹ Γ ∝ 1/x
            s = prev_x / x
            G0 = prev_G * s
            starts += [pack(G0)]
            for d in (.3, 3.):
                Rd = np.diag([np.sqrt(d), np.sqrt(1. / d)])   # 合同变换只改两通道比，不动 det 量级
                starts.append(pack(sym(Rd @ G0 @ Rd)))
        res = solve_cone(F, D, fl, starts)
        if res is None:
            print('x=%9.4g  求解失败（无可信可行点）' % x)
            continue
        It, Jt, Gt = res
        ev = np.sort(np.linalg.eigvalsh(Gt))[::-1]
        dt = It - Iu
        ok = np.isnan(I_iso) or dt <= I_iso + 1e-6
        print('x=%9.4g  Δ_opt=%9.4f  Δ_iso=%9s  差=%7s  I_opt=%9.4f  J=%.5f  eig(Γ)=[%.3e %.3e] 比=%8.2f  %s'
              % (x, dt, '%.4f' % I_iso if not np.isnan(I_iso) else '—',
                 '%+.4f' % (I_iso - dt) if not np.isnan(I_iso) else '—',
                 It, Jt, ev[0], ev[1], ev[0] / ev[1], '' if ok else '×未收敛'))
        if ok:
            good.append((x, dt, ev))
            prev_x, prev_G = x, Gt
    if len(good) >= 4:
        print('\n  通过质检的 %d 点（x 从 %.3g 到 %.3g）：' % (len(good), good[0][0], good[-1][0]))
        # 第二段质检：KKT 预言 eig(Γ) 之比 → 常数 c₁/c₂。比值塌向 1 = 解退化到各向同性 = 没收敛。
        rat = np.array([g[2][0] / g[2][1] for g in good])
        keep = [g for g, rt in zip(good, rat) if rt > 0.5 * np.max(rat)]
        print('  eig(Γ) 之比：%s' % ' '.join('%.1f' % rt for rt in rat))
        print('  比值可信段（>%.1f）取 %d/%d 点参与拟合' % (0.5 * np.max(rat), len(keep), len(good)))
        logfit([a for a, _, _ in keep], [b for _, b, _ in keep],
               '真最优 F2 r=2', 3.3219)
        print('  逐段局部指数 + 对数律不变量 Δ̄·|指数| →> (r/2)log₂e = %.4f（幂律则该值恒 = p）：'
              % (1.0 * np.log2(np.e)))
        print('%12s %12s %10s %14s' % ('x 段', 'Δ 段', '局部指数', 'Δ̄·|指数|'))
        for i in range(len(keep) - 1):
            x1, d1 = keep[i][0], keep[i][1]
            x2, d2 = keep[i + 1][0], keep[i + 1][1]
            sl = np.log(d2 / d1) / np.log(x2 / x1)
            print('%5.4g→%6.4g %6.3f→%6.3f %+10.3f %14.3f'
                  % (x1, x2, d1, d2, sl, .5 * (d1 + d2) * abs(sl)))
