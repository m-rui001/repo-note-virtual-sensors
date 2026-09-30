r"""E96c = 精确率侧前沿的自由解 vs **平面限制**解：Remark 2 那一格的最终数字。

背景（留言板 §64-F-3 / §63-E）：Tanaka 的式 (18) 是被证明**精确**的率侧机器。E96b 已复现
原文 §V 三个发表锚点（补上目标里被我漏掉的 $\\frac12\\log\\det W$ 常数后：
$D{=}33\\to6.1727$ bit、$40\\to3.2719$、$80\\to1.6031$，秩 $3/2/1$ 全对）。
本脚本要的是两个数：
 (1) $D^\\*_{free}(3\\ \\text{bit})$ —— 用 SDP 反演，用来审我车道的 42.0427；
 (2) $D_V(3\\ \\text{bit})$ —— tail-3 平面上的最小可达代价。

平面限制为什么不能交给 SDP：E96 的 [K4] 已用数字否证了凸性（两个精确平面可行点的
中点残差 $1.7\\times10^{-4}$，端点残差 $2\\times10^{-16}$）。但平面族本身有一个**干净
的参数化**：任何 range $\\subseteq V=\\operatorname{span}(Z)$ 的信息矩阵都唯一写成
$\\Sigma=Z\\,T\\,Z^t$（$T\\succ0$，$k\\times k$，$k=\\dim V$ 个自由度 $=6$），
因为 $C^tV_{obs}^{-1}C$ 里 $C=MZ^t$ 只通过 $T=(M^tV_{obs}^{-1}M)^{-1}$ 起作用。
于是平面上的搜索是 6 维全局搜索，不需要凸性。

== 判据（跑前写死，跑后不许改）==
 [V1] 同一条搜索通道在 $k=4$（自由）必须复现 (1) 的 SDP 值，误差 $\\le0.05$；
      复现不了说明搜索不可信，(2) 的数字一律不引用。
 [V2] (1) 与 42.0427 比较：若 SDP 值 $>42.0427+0.05$，则我车道"无约束锥 $I=3$ 已知
      最优 42.0427"与已发表精确前沿**矛盾**，必须自首（它是可达点，不是最优）。
 [V3] (2) 与 $43.711\\pm0.007$（数字化蓝线）比较：若 (2) $<43.711$ ⇒ 蓝线在该平面**可行**，
      Remark 2 的"蓝线可疑"当场塌掉，我的 45.4537 上界同时作废；
      若 (2) $\\ge43.80$ ⇒ 只说明"搜不到更便宜的设计"，**不构成任何认证**（搜索给上界）。
 [V4] 附带审 45.4537：若 (2) $<45.4537-10^{-4}$ ⇒ 我车道印在注记里的"平面锥最优"
      不是最优，注记 Remark 2 的数值要改。
"""
import sys
import time as _t
sys.stdout.reconfigure(encoding='utf-8')
import numpy as np
import cvxpy as cp
from scipy.linalg import schur, solve_discrete_are
from scipy.optimize import minimize, brentq

np.set_printoptions(precision=5, suppress=True, linewidth=170)
_src = open('p0/exp_c_audit.py', encoding='utf-8').read().split("print(r'== E53")[0]
_ns = {'__name__': 'p'}
exec(compile(_src, 'p0/exp_c_audit.py[preamble]', 'exec'), _ns)
A, W, TH, JC, n, sym = _ns['A'], _ns['W'], _ns['TH'], _ns['JC'], _ns['n'], _ns['sym']
ln2 = np.log(2.0)
LOGDET_W = float(np.linalg.slogdet(W)[1])
Ts, U = schur(A, output='real', sort=lambda a: abs(a) < 1.0)[:2]


def Ptil(P):
    return sym(A @ P @ A.T + W)


def free_sdp(D):
    """式 (18)，目标 = -½logdet Π + ½logdet W  （E96 漏了后者，这里补上）。"""
    P = cp.Variable((n, n), symmetric=True)
    Pi = cp.Variable((n, n), symmetric=True)
    cons = [P >> 0, Pi >> 0, P << Ptil(P), cp.trace(TH @ P) + JC <= D,
            cp.bmat([[P - Pi, P @ A.T], [A @ P, Ptil(P)]]) >> 0]
    prob = cp.Problem(cp.Minimize(-0.5 * cp.log_det(Pi) + 0.5 * LOGDET_W), cons)
    prob.solve(solver=cp.CLARABEL)
    if P.value is None:
        return np.nan, None, prob.status
    return (prob.value / ln2), sym(P.value), prob.status


def design(T, Z, check=False, rtol=1e-8):
    """平面信息矩阵 $\\Sigma=ZTZ^t\\succeq0$。取等效实现 $C=\\Sigma^{1/2}$（即 $C^tC=\\Sigma$、$V=I$，
    零方向直接丢掉），先验由滤波 DARE 解，后验 $P=P^--P^tC^t(CP^tC^t+I)^{-1}CP^-$。
    额外物理门：$\\rho((I-LC)A)<1$（可检测），否则设计不可实现。"""
    Sg = Z @ sym(T) @ Z.T
    wv, Vv = np.linalg.eigh(sym(Sg))
    mx = max(wv.max(), 1.0)
    if wv.min() < -1e-6 * mx:
        return np.inf, np.inf, 'T_not_PSD'
    idx = np.where(wv > rtol * max(wv.max(), 1e-300))[0] if wv.max() > 0 else np.array([], int)
    if len(idx) == 0:
        return np.inf, np.inf, 'no_info'
    C = (np.sqrt(wv[idx])[:, None] * Vv[:, idx].T)      # r x n
    r = len(idx)
    try:
        Pm = sym(solve_discrete_are(A.T, C.T, W, np.eye(r)))
    except Exception as e:
        return np.inf, np.inf, f'dare:{type(e).__name__}'
    Ir = np.eye(r)
    Sm = C @ Pm @ C.T + Ir
    Lk = Pm @ C.T @ np.linalg.inv(Sm)                  # n x r
    Pp = sym(Pm - Lk @ Sm @ Lk.T)
    rho = float(np.abs(np.linalg.eigvals((np.eye(n) - Lk @ C) @ A)).max())
    if rho >= 1.0:
        return np.inf, np.inf, f'undetectable:rho={rho:.4f}'
    D = JC + float(np.trace(TH @ Pp))
    I = 0.5 * (np.linalg.slogdet(Pm)[1] - np.linalg.slogdet(Pp)[1]) / ln2
    if not np.isfinite(D) or not np.isfinite(I):
        return np.inf, np.inf, 'nan'
    if check:
        Pq = W.copy()
        stt = 'maxit'
        for _ in range(40000):
            Pn = sym(np.linalg.inv(np.linalg.inv(Ptil(Pq)) + Sg))
            if np.abs(Pn - Pq).max() > 1e14:
                stt = 'blow'
                break
            if np.abs(Pn - Pq).max() < 1e-13 * max(1.0, np.abs(Pn).max()):
                stt = 'conv'
                Pq = Pn
                break
            Pq = Pn
        Itr = 0.5 * (np.linalg.slogdet(Ptil(Pq))[1] - np.linalg.slogdet(Pq)[1]) / ln2
        Dtr = JC + float(np.trace(TH @ Pq))
        print(f'   [自检 DARE vs 定点] |D|={abs(D-Dtr):.2e} |I|={abs(I-Itr):.2e}  定点={stt}'
              f'  rho((I-LC)A)={rho:.5f}  r={r}')
    return D, I, 'ok'


def scale_to(T, Z, I0):
    """沿射线 $sT$ 找出 $I=I_0$ 的尺度。$I(s)$ 随 $s$ 单调增，所以先在小范围内找变号区间，
    再用 brentq。候选点本身已由罚函数带到 $I_0$ 附近，所以窗口不需要大；
    遇到不可实现的尺度（DARE 失败/不可检测）就停止扩张，不做外推。"""
    def at(ls):
        D, I, st = design(np.exp(ls) * T, Z)
        return (I, D) if st == 'ok' else None

    base = at(0.0)
    if base is None:
        return np.inf, np.inf, None, 'base_fail'
    lo = hi = None
    for step in (0.02, 0.1, 0.4, 1.5, 5.0, 12.0):
        if lo is None:
            r = at(-step)
            if r is not None and r[0] < I0:
                lo = -step
            elif r is None:
                lo = -999.0        # 更小的尺度不可实现，视为已到边界
        if hi is None:
            r = at(step)
            if r is not None and r[0] > I0:
                hi = step
            elif r is None:
                hi = 999.0
        if lo is not None and hi is not None:
            break
    if lo is None or hi is None or not (-100 < lo < 100 and -100 < hi < 100):
        return np.inf, base[0], T, 'nofix'
    f = lambda ls: (at(ls)[0] - I0) if at(ls) is not None else np.nan
    ls = brentq(f, lo, hi, xtol=1e-11, rtol=1e-13, maxiter=200)
    Tt = np.exp(ls) * T
    D, I, st = design(Tt, Z)
    return D, I, Tt, st


def search(Z, I0, nstart=60, seed=0, lam=300.0, budget=3000, mode='full'):
    """min D s.t. I = I0（率约束必然取等：D 随 I 单调减）。
    坐标：mode='full' 取 $T=LL^t$ 的上三角元素（$k(k+1)/2$ 维，无冗余维）；
    mode='diag' 只取 Schur 基下对角 $T=\\mathrm{diag}(e^{x})$（$k$ 维，即"逐本征方向配功率"子族，
    用来审 45.4537 是否只是受限族的 optimum）。
    两阶段：(i) 双侧罚 $D+\\lambda|I-I_0|$ 的 Nelder-Mead 多起点；(ii) 对每个候选沿射线
    精确二分尺度，保证报告的每个点都是**严格可行**的设计（因此是上界，不是认证）。"""
    k = Z.shape[1]
    iu = np.triu_indices(k)
    m = len(iu[0]) if mode == 'full' else k
    rng = np.random.default_rng(seed)

    def toT(x):
        if mode == 'full':
            Lm = np.zeros((k, k))
            Lm[iu] = x
            return Lm @ Lm.T
        return np.diag(np.exp(x))

    def pen_obj(x):
        D, I, st = design(toT(x), Z)
        if st != 'ok':
            return 1e7 + float(np.dot(x, x))
        return D + lam * abs(I - I0)

    cands = []
    t0 = _t.time()
    for t in range(nstart):
        if mode == 'full':
            x0 = (0.5 * np.eye(k))[iu] if t == 0 else rng.normal(0.0, 0.7, m)
        else:
            x0 = np.full(k, -0.5) if t == 0 else rng.normal(-0.5, 1.5, m)
        for _ in range(2):
            x0 = minimize(pen_obj, x0, method='Nelder-Mead',
                          options=dict(maxiter=budget, maxfev=budget, xatol=1e-11, fatol=1e-13)).x
        cands.append(x0)
        if (t + 1) % 10 == 0 or t + 1 == nstart:
            Dp, Ip, stp = design(toT(x0), Z)
            print(f'   [{k}维/{mode}] start {t+1}/{nstart} 用时 {_t.time()-t0:.0f}s '
                  f'罚值目标 {pen_obj(x0):.4f}  (D,I,状态)=({Dp:.4f},{Ip:.4f},{stp})', flush=True)
    best = (np.inf, None, np.inf, None)
    feas = []
    for x in cands:
        D, I, Tt, st = scale_to(toT(x), Z, I0)
        if st == 'ok':
            feas.append(D)
            if D < best[0]:
                best = (D, Tt, I, design(Tt, Z))
    feas = np.sort(np.array(feas)) if feas else np.array([np.inf])
    return best, feas


print('== (1) 自由前沿：SDP 网格 + 3 bit 反演 ==')
for D in (33.0, 40.0, 42.0427, 46.0, 80.0):
    v, Pv, st = free_sdp(D)
    print(f' D={D:8.4f}  DI={v:.4f} bit  status={st}')
Dstar = brentq(lambda D: free_sdp(D)[0] - 3.0, 34.0, 90.0, xtol=1e-5)
v, Pstar, st = free_sdp(Dstar)
print(f' [V2] D*_free(3 bit) = {Dstar:.4f}   （核验 DI={v:.5f}）  与 42.0427 差 {Dstar-42.0427:+.4f}')
SNR = sym(np.linalg.inv(Pstar) - np.linalg.inv(Ptil(Pstar)))
sv = np.sort(np.linalg.svd(SNR, compute_uv=False))[::-1]
print(f'      rank(0.1%)={int((sv>1e-3*sv.max()).sum())}  sv={sv}')
Zb = U[:, n - 3:]
Qm = np.eye(n) - Zb @ Zb.T
w_, Vv = np.linalg.eigh(sym(SNR))
print('      自由最优 SNR 逐本征方向：权重 / 落在 tail-3 平面内的比例 ||Z^t v||')
for i in np.argsort(w_)[::-1]:
    if w_[i] > 1e-6 * max(w_.max(), 1.0):
        vv = Vv[:, i]
        print(f'        w={w_[i]:10.5f}   ||Z^t v||={np.linalg.norm(Zb.T @ vv):.6f}'
              f'   在 U 基下的坐标 {np.abs(U.T @ vv)}')
print(f'      该 P* 的平面残差 ||Q Sigma||/||Sigma|| = {np.linalg.norm(Qm @ SNR)/np.linalg.norm(SNR):.3e}')

print('\n== (2) 同一条搜索通道在 k=4 上复核 [V1] ==')
Z4 = U[:, :4]
# 2a 先把 SDP 解塞回 design()：T = U^t Sigma U 必须复现 (D*, 3 bit)
Tchk = sym(Z4.T @ SNR @ Z4)
tw, tv = np.linalg.eigh(Tchk)
Tchk = tv @ np.diag(np.maximum(tw, 0.0)) @ tv.T
Dchk, Ichk, stchk = design(Tchk, Z4, check=True)
print(f' design(SDP 解的 T) = D {Dchk:.4f} (SDP {Dstar:.4f}), I {Ichk:.5f} (SDP {v:.5f})  状态 {stchk}')
b4, f4 = search(Z4, 3.0, nstart=25, seed=1)
print(f' 搜索 best D={b4[0]:.4f}  I={b4[2]:.5f}   SDP D*={Dstar:.4f}  差 {b4[0]-Dstar:+.4f}'
      f'   [V1] {"通过" if abs(b4[0]-Dstar)<0.05 else "失败——搜索不可信"}')
print(f'   可行起点 {int(np.isfinite(f4).sum())}/25，D 分位 min/med/max = {f4[0]:.4f}/{np.median(f4):.4f}/{f4[-1]:.4f}')

print('\n== (3) tail-3 平面上的 D_V(3 bit) ==')
b3, f3 = search(Zb, 3.0, nstart=80, seed=2)
print(f' 搜索 best D={b3[0]:.4f}  I={b3[2]:.5f}   可行起点 {int(np.isfinite(f3).sum())}/80'
      f'   前 6 个 D: {f3[:6]}')
print(f' [V3/V4] 对照：蓝线数字化 43.711±0.007 / 验收线 43.80 / 我车道 Powell 上界 45.4537'
      f' / 自由精确下界 {Dstar:.4f}')
print(f'   D_plane - 43.711 = {b3[0]-43.711:+.4f}   D_plane - 45.4537 = {b3[0]-45.4537:+.4f}')
if b3[1] is not None:
    Tm = b3[1]
    print(f'   T = {Tm}')
    print(f'   eig(T) = {np.linalg.eigvalsh(Tm)}')
    Sn = Zb @ Tm @ Zb.T
    Dv, Iv, stv = design(Tm, Zb, check=True)
    print(f'   复核 design(T) -> D={Dv:.5f} I={Iv:.6f} 状态 {stv}'
          f'   平面残差(构造级) {np.linalg.norm(Qm @ Sn):.3e}')
    print(f'   rank(0.1%)={int((np.linalg.svd(Sn, compute_uv=False) > 1e-3*np.linalg.svd(Sn, compute_uv=False).max()).sum())}')

print('\n== (4) 对角子族（Schur 基下 T 对角）—— 审 45.4537 的口径 ==')
Ddiag, Tdiag, Idiag, _ = search(Zb, 3.0, nstart=20, seed=7, mode='diag')[0]
print(f' 对角 T 的平面最优 D = {Ddiag:.4f}（I={Idiag:.5f}）  与 45.4537 差 {Ddiag-45.4537:+.4f}')
print(f' eig(T 对角族) = {np.linalg.eigvalsh(Tdiag)}')
Z2 = U[:, n - 2:]
b2, f2 = search(Z2, 3.0, nstart=40, seed=3)
d2 = search(Z2, 3.0, nstart=20, seed=8, mode='diag')[0]
print(f' tail-2 平面：全 3 维搜索 D={b2[0]:.4f}（起点可行 {int(np.isfinite(f2).sum())}/40），'
      f'对角子族 D={d2[0]:.4f}  对照我车道 50.6023（差 {b2[0]-50.6023:+.4f} / {d2[0]-50.6023:+.4f}）')
