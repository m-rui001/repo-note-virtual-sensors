r"""E72 = 把 E71 的下界缺口**归因**：是真对偶间隙，还是我用常数替掉先验项的损失？

E71 的判定：蓝平面 I0=3 上 sup L = 33.014，远低于验收线 43.80，也低于 Powell 上界 45.4537。
合法（没有任何一行越过已知可达值）但太弱。弱在哪有两个候选：
 (a) **真对偶间隙**：前沿 g(I) 在 I0 处没有支撑线 => sup_mu[p(mu)-mu I0] < g(I0)，
     那么"认证下界"这条路本身就够不到验收线，跟实现无关。
 (b) **我的凸代理的松弛**：E71 把 (mu/2)logdet G(X) 用常数 logdet(Z'WZ) 替掉，
     损失 = (mu/2ln2)(logdet G - logdet Z'WZ) **正比于 mu**，而对偶恰恰需要大 mu
     （无约束前沿在 I=3 的左斜率 ~ -16，右斜率 ~ -5.5，所以 mu* ~ 10 量级）。
     若是 (b)，这条路能修好；若是 (a)，Remark 2 只能继续降档。

p(mu) := inf over **物理设计** M>=0 的 [D(M) + mu I(M)]，用精确 Riccati（sdare）逐点算，
多起点 Powell 全局搜。那么
    sup_mu [p(mu) - mu I0]  <=  g(I0)   （对每个 mu 都成立，见 E71 §1）
等号 <=> 前沿在 I0 处凸/可支撑。**这一条不需要凸松弛，是"真"拉格朗日对偶的数值实现。**
它同时给出 g(I0) 的第二个独立估计（从下方），直接对着 43.80 验收线判。

判据（写在前）：
 [0] 自由平面：E63f 的前沿点 g(2)=58.1931, g(3)=42.0427, g(4)=36.5386。
     sup[p-mu*3] 必须 <= 42.0427（越界 = 搜索口径或率算错）。若接近 42.04 => 无对偶间隙 => (b)。
 [1] 蓝平面：sup[p-mu*3] 与 43.80 比。> 43.80 => Remark 2 有升级的希望（等认证）；
     <= 43.711 => 蓝线在平面内可能可行，Remark 2 彻底降档。
 [2] 一致性：sup[p-mu*3] 也必须 <= Powell 上界 45.4537。
 [3] 方向声明（不要再犯第二次）：p_hat(mu) 是多重起点局部搜出来的**上界**，
     所以 L_hat = p_hat - mu*I0 是上界 on the dual，**不是证书**。它的用处只在
     (i) 定位 mu* 和前沿的下沿估计，(ii) 告诉我认证需要补哪条腿（见 E73 的割线外逼近）。
 首跑（DI 里没有稳定性护栏）在 mu=6(自由)/mu=2(蓝) 出过 D=-9.8e6 / -4.6e14 的垃圾行：
 sdare 在不可检测的 (A,F) 上返回非稳定化解。已补 E69 JIfast 同一条闭环谱半径 guard（debt #6）。
"""
import sys, time
sys.stdout.reconfigure(encoding='utf-8')
import numpy as np
from scipy.linalg import schur, solve_discrete_are as sdare
from scipy.optimize import minimize

np.set_printoptions(precision=4, suppress=True, linewidth=170)
_src = open('p0/exp_c_audit.py', encoding='utf-8').read().split("print(r'== E53")[0]
_ns = {'__name__': 'p'}
exec(compile(_src, 'p0/exp_c_audit.py[preamble]', 'exec'), _ns)
A, W, TH, JC, n, sym = _ns['A'], _ns['W'], _ns['TH'], _ns['JC'], _ns['n'], _ns['sym']
ln2 = np.log(2.0)
EPS = 1e-12
print('== E72：真拉格朗日对偶（在物理设计集上直接罚）+ 缺口归因 ==')


def DI(S):
    """精确：S -> (D, I)，信息型定点经 sdare；S 奇异/退化 -> inf。"""
    ev, EV = np.linalg.eigh(sym(S))
    keep = ev > max(1e-11 * max(ev.max(), 1e-30), 1e-300)
    if not keep.any():
        return np.inf, np.inf
    F = (EV[:, keep] * np.sqrt(ev[keep])).T
    k = F.shape[0]
    try:
        Pt = sym(sdare(A.T, F.T, sym(W), np.eye(k)))
    except Exception:
        return np.inf, np.inf
    if not np.all(np.isfinite(Pt)):
        return np.inf, np.inf
    FPt = F @ Pt
    try:
        R = np.linalg.solve(FPt @ F.T + np.eye(k), FPt)
    except Exception:
        return np.inf, np.inf
    P = sym(Pt - FPt.T @ R)
    # 稳定性护栏：sdare 在 (A,F) 不可检测时会给非稳定化解，代价可为负数（E72 首跑就是这么脏的）。
    # 与 E69 的 JIfast 逐字同一条 guard（debt #6）。
    cl = A - A @ Pt @ F.T @ np.linalg.solve(FPt @ F.T + np.eye(k), F)
    if np.max(np.abs(np.linalg.eigvals(cl))) >= 1.0 - 1e-9:
        return np.inf, np.inf
    if np.min(np.linalg.eigvalsh(P)) <= 0 or np.min(np.linalg.eigvalsh(Pt)) <= 0:
        return np.inf, np.inf
    ld = np.linalg.slogdet(Pt)[1] - np.linalg.slogdet(P)[1]
    if not np.isfinite(ld) or ld <= 0:
        return np.inf, np.inf
    return JC + float(np.trace(TH @ P)), 0.5 * ld / ln2


Ts, U = schur(A, output='real', sort=lambda a: abs(a) < 1.0)[:2]
Z3 = U[:, n - 3:]
PLANES = {'free_R4': U[:, :], 'blue_tail3': Z3}
I0 = 3.0
MUS = [1.0, 2.0, 4.0, 6.0, 8.0, 11.0, 15.0, 20.0, 28.0]
G3 = {'free_R4': 42.0427, 'blue_tail3': 45.4537}


def penal(Z, mu, nz):
    """inf_M [D + mu I]，M = L L^t（下三角指数参数化，覆盖秩亏边界）。"""
    r = Z.shape[1]

    def unpack(x):
        Lm = np.zeros((r, r)); k = 0
        for i in range(r):
            for j in range(i + 1):
                Lm[i, j] = np.exp(x[k]) if i == j else x[k]
                k += 1
        return sym(Lm @ Lm.T)

    def f(x):
        M = unpack(x)
        d, I = DI(sym(Z @ M @ Z.T))
        if not np.isfinite(d):
            return 1e7 + 1e3 * float(np.sum(np.abs(x)))
        return d + mu * I

    rng = np.random.default_rng(7)
    wT = np.clip(np.linalg.eigvalsh(TH), 1e-9, None)[::-1]
    starts = [np.full(nz, -1.0), np.full(nz, 0.5), np.full(nz, 1.5), np.full(nz, 3.0),
              np.concatenate([np.log(np.sqrt(wT))[:r], np.zeros(nz - r)])]
    starts += [rng.normal(0.5, 1.5, nz) for _ in range(6)]
    bb = (np.inf, None)
    for x0 in starts:
        rq = minimize(f, x0, method='Powell', options=dict(maxfev=1200, xtol=1e-5, ftol=1e-10))
        if rq.fun < bb[0]:
            bb = (float(rq.fun), rq.x)
    M = unpack(bb[1])
    d, I = DI(sym(Z @ M @ Z.T))
    return bb[0], d, I


t0 = time.time()
print('\n[0]/[1] p(mu) 与 L_true(mu) = p(mu) - mu*I0，I0=%.1f' % I0)
for name, Z in PLANES.items():
    r = Z.shape[1]
    nz = r * (r + 1) // 2
    rows, bestL = [], (-np.inf, None)
    for mu in MUS:
        p, d, I = penal(Z, mu, nz)
        L = p - mu * I0
        rows.append((mu, p, L, d, I))
        if L > bestL[0]:
            bestL = (L, mu)
    print('  %s (r=%d)  已知 g(3)>=%.4f' % (name, r, G3[name]))
    for mu, p, L, d, I in rows:
        print('     mu=%6.1f  p=%10.4f  L=%10.4f  | 胜者设计 D=%9.4f I=%7.4f' % (mu, p, L, d, I))
    over = [L for _, _, L, _, _ in rows if L > G3[name] + 1e-3]
    print('     sup L = %.5f 在 mu=%.1f | 越过前沿上界的行数 %d %s'
          % (bestL[0], bestL[1], len(over), '**口径错**' if over else '合法'))
    if name == 'free_R4':
        print('     => 与可达上界 42.0427 差 %.4f。**不要读成"真对偶间隙"**：p_hat 是 p 的上界，'
              '所以 L_hat 也是对偶值的上界，这个差混了"罚问题没搜满"与"对偶间隙"两个来源，分不开。'
              '可说的是反面：没有任何一行越过可达值 => 与"零间隙 + 前沿在 I0 凸"一致（不是证明）。'
              '认证版下界见 E74。' % (42.0427 - bestL[0]))
    else:
        print('     => 对验收线 43.80：%s | 对 43.711：%s（上界方向，只用来定位 mu* 与量级；'
              '认证见 E74）' % ('超过' if bestL[0] > 43.80 else '未超过',
                                '超过' if bestL[0] > 43.711 else '未超过'))

print('\n[2] 归因：E71 蓝平面 sup L=33.014 vs 本脚本 sup L=?; 两者之差 = 凸代理（常数先验）的损失')
print('  耗时 %.0f s' % (time.time() - t0))
