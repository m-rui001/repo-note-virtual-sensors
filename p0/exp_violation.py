r"""E48 = 给 E47 那三个负偏差定杠子：逐格实测评估误差，配准到噪声底，再做一阶校正。

E47 扫 $\rho$ 时，每个 $\rho$ 都要把 $v_1$ 调回预算 $x_t$，比较的却是"率"。
预算残差 $r=\hat x-x_t$ 会按沿族斜率 $s=dI/dx$ 折进率的差里，而 $s\asymp1/(x_t\ln2)$：
$x_t=10^{-6}$ 那行的杠子天生是 $10^{-9}$ 档，$x_t=1$ 那行才该要求 $10^{-12}$。
所以"出现负值即反例"这句话在没有杠子之前不成立——#22 §1 的可分辨性规则这次轮到我用自己的。
E47 的接受门是 $\|r\|\le10^{-9}x_t$，在 $x_t=1$ 处允许的残差折到率上就是 $7\times10^{-10}$，
比它报出的最大负值 $-9.9\times10^{-11}$ 大一个数量级：那张表的负号是门给的，不是物理给的。

三步：
[1] 控制恒等式（$\rho=0$、两族共用同一个 $v_1$，完全不需要二分）：E46 (i)(ii) 给 $x$ 相同、
    $I_{\rm full}=I_{\rm del}+\tfrac12\log_2(1+\tilde P_{22}/v_2)$。两边是同一批量的代数恒等式，
    实测互差就是这条流水线的绝对精度；$\rho=0$ 已被证成命题，这里读出的任何非零都不是物理。
[2] $\delta x$ 逐格实测：同一格上四条等价路线互差（Riccati+Schur / 信息阵求逆 / 交换测量行 / 独立定点迭代解 DARE），附 $\kappa(V)$。
    前三条共用 scipy 的 DARE 输出，看不见 DARE 自身的前向误差；第四条不共享任何步骤，杠子由它撑开。
[3] 二分 45 步 + 阻尼 Newton 抛光 6 步把 $r$ 压到评估噪声底，一阶校正后给区间；判决只用区间。
"""
import sys
sys.stdout.reconfigure(encoding='utf-8')
import numpy as np
import p0.exp_plants_oos as E35
from p0.ticket import floor0

ln2 = np.log(2.)
NB = 2


def plant(c13=0.0, au=.6, b3=0.0, wu=1.0):
    A = np.array([[1.25, .30, c13], [0., .85, 0.], [0., 0., au]])
    B = np.array([[1.], [0.], [b3]])
    W = np.diag([1.0, .7, wu])
    return E35.Pl(A, B, W, np.eye(3), np.array([[1.0]]))


F = np.array([[1., 0., 0.], [0., 0., 1.]])
FD = np.array([[1., 0., 0.]])
Pl = plant()
PhiF, PF, resF = floor0(Pl, F)
RHO = np.arange(-24, 25) * 0.04           # 含 0.0，49 点，端点 $\pm0.96$
# $|\rho|=1$ 使 $V$ 精确奇异，这个参数化下它不是合法点：它对应 E49 的"旋转单通道"（$V^{-1}$ 秩 1），那里用 $(\varphi,\tau)$ 精确处理。
XTS = [1e-6, 1e-4, 1e-2, 1e-1, 1.0]
V2S = [1e-1, 1e2, 1e6, 1e12]


def sym(M):
    return .5 * (M + M.T)


def Vof(v1, v2, rho):
    if rho == 0.0:
        return np.array([[v1, 0.], [0., v2]])
    sd = np.sqrt(v1 * v2)
    return np.array([[v1, rho * sd], [rho * sd, v2]])


def v_of(Fx, v1, v2, rho):
    return np.array([[v1]]) if Fx.shape[0] == 1 else Vof(v1, v2, rho)


def one(v1, v2, rho, Fx):
    r"""返回 $(x,I,\tilde P)$。$r=1$ 族不看 $v_2,\rho$。"""
    V = v_of(Fx, v1, v2, rho)
    Pt = Pl.Pt(Fx, V)
    P = Pl.P(Fx, V, Pt)
    x = float(np.trace(Pl.Th @ P)) - PhiF
    _, l1 = np.linalg.slogdet(Fx @ Pt @ Fx.T + V)
    _, l2 = np.linalg.slogdet(V)
    return x, .5 * (l1 - l2) / ln2, Pt


def i_chol(Fx, V, Pt):
    r"""$I$ 的第二条路线：Cholesky 对数行列式，与 slogdet 互差即 $\delta I$。"""
    S = sym(Fx @ Pt @ Fx.T + V)
    return (np.log(np.diag(np.linalg.cholesky(S))).sum()
            - np.log(np.diag(np.linalg.cholesky(sym(V)))).sum()) / ln2


def x_alts(Fx, V):
    r"""$x$ 的等价路线二、三、四：信息阵形式、交换测量行、以及**独立解 DARE** 的定点迭代。

    前两条共用同一个 `solve_discrete_are` 的输出，所以 DARE 自身的前向误差对它们完全透明——
    这正是 E47/E48 早期把杠子估小、把 $\kappa(V)\sim10^{18}$ 的格误判成"反例"的原因。
    第四条从 $\tilde P_0=W$ 起迭代 $\tilde P\leftarrow A(\tilde P^{-1}+M)^{-1}A^{\top}+W$ 到双精度停滞，
    与 Schur 法没有任何共享步骤，它给的互差才是 $\delta x$ 的真值。
    """
    Pt = Pl.Pt(Fx, V)
    P = np.linalg.inv(np.linalg.inv(Pt) + Fx.T @ np.linalg.solve(V, Fx))
    out = [float(np.trace(Pl.Th @ sym(P))) - PhiF]
    pe = list(range(Fx.shape[0]))[::-1]
    out.append(Pl.xI(Fx[pe], V[np.ix_(pe, pe)], PhiF)[0])
    fp = fp_x(Fx, V)
    if fp[5]:
        out.append(fp[0])
    return out


def fp_x(Fx, V, it=6000, tol=1e-17):
    r"""独立 DARE：定点迭代 $\tilde P\leftarrow A(\tilde P^{-1}+M)^{-1}A^{\top}+W$，返回 $(x,I,\tilde P,\text{残差},\text{步数},\text{收敛})$。"""
    M = Fx.T @ np.linalg.solve(V, Fx)
    A, W, n = Pl.A, Pl.W, Pl.nn
    Pt = sym(W.copy())
    res = np.inf
    conv = False
    k = 0
    for k in range(it):
        P = np.linalg.inv(np.linalg.inv(Pt) + M)
        P = sym(P)
        Pt2 = sym(A @ P @ A.T + W)
        res = np.linalg.norm(Pt2 - Pt, ord=np.inf)
        Pt = Pt2
        if res < tol * max(1.0, np.linalg.norm(Pt, ord=np.inf)):
            conv = True
            break
    P = sym(np.linalg.inv(np.linalg.inv(Pt) + M))
    x = float(np.trace(Pl.Th @ P)) - PhiF
    _, l1 = np.linalg.slogdet(Fx @ Pt @ Fx.T + V)
    _, l2 = np.linalg.slogdet(V)
    return x, .5 * (l1 - l2) / ln2, Pt, res, k, conv


def spread(vals):
    return max(vals) - min(vals)


def ev(u, v2, rho, Fx, h=1e-4):
    r"""$u=\log_{10}v_1$ 处取值 + 中心差分导数。返回 $(v_1,x,I,\tilde P,dx/du,dI/du,\mathrm{curv})$。"""
    xl, Il, _ = one(10. ** (u - h), v2, rho, Fx)
    xh, Ih, _ = one(10. ** (u + h), v2, rho, Fx)
    x0, I0, Pt0 = one(10. ** u, v2, rho, Fx)
    return (10. ** u, x0, I0, Pt0, (xh - xl) / (2 * h), (Ih - Il) / (2 * h),
            abs(I0 - .5 * (Il + Ih)))


def solve_xt(xt, v2, rho, Fx, lo=-16.0, hi=14.0, ib=45, ip=6):
    r"""二分定位 + 阻尼 Newton 抛光，使 $x\to x_t$ 到评估噪声底。返回 dict 或 None。"""
    if one(10. ** lo, v2, rho, Fx)[0] > xt or one(10. ** hi, v2, rho, Fx)[0] < xt:
        return None
    a, b = lo, hi
    for _ in range(ib):
        m = .5 * (a + b)
        if one(10. ** m, v2, rho, Fx)[0] < xt:
            a = m
        else:
            b = m
    u = .5 * (a + b)
    best = None
    for _ in range(ip + 1):
        q = ev(u, v2, rho, Fx)
        if best is None or abs(q[1] - xt) < abs(best[1] - xt):
            best = q
        if q[4] <= 0:
            break
        step = max(-.1, min(.1, (xt - q[1]) / q[4]))
        if step == 0.0:
            break
        u += step
    v1, x0, I0, Pt0, dxu, dIu, curv = best
    V = v_of(Fx, v1, v2, rho)
    try:
        dxs = spread([x0] + x_alts(Fx, V))
        dis = abs(I0 - i_chol(Fx, V, Pt0))
    except np.linalg.LinAlgError:
        return None                       # 路线失败 = 该格 $V$ 数值奇异，不是反例也不是排除
    return dict(v1=v1, x=x0, I=I0, Pt=Pt0, s=dIu / dxu, curv=curv,
                r=abs(x0 - xt), dx=dxs, dI=dis, cond=np.linalg.cond(V))


def fp_bisect(xt, v2, rho, Fx, lo=-16.0, hi=14.0, it=55):
    r"""把 $x$ 与 $I$ 全部换成定点迭代版本，再二分一次。返回 $(v_1,\hat x,I)$ 或 None。"""
    f = lambda u: fp_x(Fx, v_of(Fx, 10.0 ** u, v2, rho))
    if f(lo)[0] > xt or f(hi)[0] < xt:
        return None
    a, b = lo, hi
    for _ in range(it):
        m = .5 * (a + b)
        if f(m)[0] < xt:
            a = m
        else:
            b = m
    v1 = 10.0 ** (.5 * (a + b))
    q = f(.5 * (a + b))
    if not q[5]:
        return None
    return v1, q[0], q[1]


print('== E48：逐格误差杠子 + 配准到噪声底，判 E47 的负偏差是不是反例 ==')
print('')
print(r'[1] 控制恒等式：$\rho=0$、两族共用同一个 $v_1$（不需要二分）。E46 给 $x_{\rm full}=x_{\rm del}$、$I_{\rm full}-I_{\rm del}=\tfrac12\log_2(1+\tilde P_{22}/v_2)$')
for v1 in [1e-6, 1e-3, 1e0]:
    xf, If, Ptf = one(v1, 1.0, 0.0, F)
    xd, Id, _ = one(v1, 1.0, 0.0, FD)
    print(r'  v1=%9.1e   $x_{\rm full}-x_{\rm del}=%+.1e$   $I_{\rm del}=%.9f$   $\tilde P_{22}=%.9f$'
          % (v1, xf - xd, Id, Ptf[NB, NB]))
    for v2 in V2S:
        xq, Iq, _ = one(v1, v2, 0.0, F)
        cf = Id + .5 * np.log1p(Ptf[NB, NB] / v2) / ln2
        print(r'     v2=%8.0e  $I_{\rm full}$=%.9f  闭式=%.9f  差 %+.1e   $x$ 随 $v_2$ 漂 %+.1e'
              % (v2, Iq, cf, Iq - cf, xq - xf))
print(r'  读法：这两列互差都是同一批量的等价算法之差，不是物理。它就是 [3] 里 $\delta I,\delta x$ 的量级。')

print('')
print(r'[2] $\delta x$ 的逐格实测：四条等价路线互差（第四条是独立解 DARE 的定点迭代），附 $\kappa(V)$。取 E47 三个负读数所在的格')
for (xt, v2, rho) in [(1e-6, 1e12, .86), (1e-4, 1e12, -.72), (1.0, 1e12, -.96), (1e-6, 1e12, 0.0)]:
    q = solve_xt(xt, v2, rho, F)
    if q is None:
        print(r'  $x_t=%.0e$ $v_2=%.0e$ $\rho=%+.2f$: 路线失败（$V$ 数值奇异）' % (xt, v2, rho))
        continue
    print(r'  $x_t=%.0e$ $v_2=%.0e$ $\rho=%+.2f$: $v_1$=%.4e  $\hat x-x_t=%+.1e$  $\delta x$=%.1e  $\delta I$=%.1e  $\kappa(V)$=%.1e  $|s|$=%.2e  $|s|\delta x$=%.1e'
          % (xt, v2, rho, q['v1'], q['x'] - xt, q['dx'], q['dI'], q['cond'], abs(q['s']), abs(q['s']) * q['dx']))
print(r'  $\delta x$ 随 $\kappa(V)$ 走，不随 $x_t$ 走；$|s|\asymp1/(x_t\ln2)$ 随 $x_t$ 放大。杠子是这两者的乘积，跨六个数量级。')

print('')
print(r'[3] 主表：每格 $\rho$ 扫 51 点，配准到噪声底后一阶校正。杠子 $B=|s|(\delta x+r)+\delta I+|s_{\rm del}|r_{\rm del}+$曲率')
print(r'    区间 $[\hat I-B,\hat I+B]$ 是该 $(\rho,v_2)$ 上可达最优率的估计区间；整段压在 $E_{\rm del}$ 下才算反例，整段在上才算排除。')
BARS = {}
for xt in XTS:
    d = solve_xt(xt, 1.0, 0.0, FD)
    Ed, rd, sd, bdel = d['I'], d['r'], abs(d['s']), abs(d['s']) * (d['dx'] + d['r']) + d['dI']
    print(r'  x_t = %.1e   $E_{\rm del}$=%.9f   残差 %.1e   $\delta x$=%.1e   $|s_{\rm del}|$=%.2e   删行侧杠子 %.1e'
          % (xt, Ed, rd, d['dx'], sd, bdel))
    BARS[xt] = []
    for v2 in V2S:
        rows = []
        for rho in RHO:
            q = solve_xt(xt, v2, rho, F)
            if q is None or q['r'] > 1e-6 * xt:
                continue
            bar = abs(q['s']) * (q['dx'] + q['r']) + max(q['dI'], 1e-16) + bdel + q['curv']
            corr = q['I'] - q['s'] * (q['x'] - xt)
            rows.append((corr, bar, rho, q))
        if not rows:
            print(r'     v2=%8.0e  无可行点' % v2)
            continue
        best = min(rows, key=lambda e: e[0] + e[1])       # 用最保守的一端排序，别把舍入当成赢
        corr, bar, rho, q = best
        i0 = [e for e in rows if e[2] == 0.0]
        lo, hi = corr - bar, corr + bar
        verdict = ('反例成立' if hi < Ed else ('排除' if lo > Ed else '分辨不出'))
        BARS[xt].append((v2, q['I'] - Ed, corr - Ed, bar, verdict, rho, q))
        print(r'     v2=%8.0e  保守最优 rho=%+.2f ($v_1$=%.3e, 残差 %.1e)  raw %+.3e  校正 %+.3e  杠子 %.2e  区间[%+.2e,%+.2e]  %s  |  $\rho=0$ raw %+.3e'
              % (v2, rho, q['v1'], q['r'], q['I'] - Ed, corr - Ed, bar, lo - Ed, hi - Ed, verdict,
                 (i0[0][0] - Ed) if i0 else float('nan')))
    print(r'     杠子跨 $x_t$ 的标度：小一个数量级的预算，可分辨的率差放大十倍。')

print('')
print(r'[4] 对照 E47 的三个负读数')
NEG = {(1e-6, 1e12): -1.679e-09, (1e-4, 1e12): -1.953e-11, (1.0, 1e12): -9.871e-11}
for (xt, v2), rw in sorted(NEG.items()):
    row = [e for e in BARS.get(xt, []) if e[0] == v2]
    if not row:
        print(r'  x_t=%.0e v2=%.0e  本轮无行' % (xt, v2))
        continue
    _, raw, corr, bar, verdict, rho, q = row[0]
    print(r'  $x_t=%.0e$ $v_2=%.0e$:  E47 报 %+.3e   E48 未校正复测 %+.3e   校正后 %+.3e   杠子 %.2e   比值(未校正/杠子) %.1f   %s   [$\rho$ 取 %+.2f, $\delta x$=%.1e, 残差 %.1e]'
          % (xt, v2, rw, raw, corr, bar, abs(rw) / bar, verdict, rho, q['dx'], q['r']))
print('')
print(r'[4b] 换求解器重跑同一格：整条流水线（DARE、配准、$I$）都不再用 scipy 的 Schur 法，只用 Riccati 定点迭代')
for (xt, v2), rw in sorted(NEG.items()):
    row = [e for e in BARS.get(xt, []) if e[0] == v2]
    if not row:
        continue
    rho = row[0][5]
    gd = fp_bisect(xt, 1.0, 0.0, FD)
    gc = fp_bisect(xt, v2, 0.0, F)
    gb = fp_bisect(xt, v2, rho, F)
    if not (gd and gc and gb):
        print(r'  $x_t=%.0e$ $v_2=%.0e$  独立迭代未收敛/未 bracket，这一格不能下结论' % (xt, v2))
        continue
    ctrl = gc[2] - gd[2]
    print(r'  $x_t=%.0e$ $v_2=%.0e$ $\rho_*=%+.2f$:  独立版 $I-I_{\rm del}=%+.3e$（scipy 版 %+.3e）   同格 $\rho=0$ 控制 %+.3e（E46 证 $\ge0$，它就是这一格的独立误差档）'
          % (xt, v2, rho, gb[2] - gd[2], row[0][1], ctrl))
print(r'  同一格两个求解器给不出同号同幅的负值，就是"这不是物理"的判据；$\rho=0$ 控制读到负的，说明该格的负号根本不可信。')

print('')
print(r'[5] 判决与门槛')
vs = [e[4] for xt in XTS for e in BARS.get(xt, [])]
print(r'  20 格里：排除 %d，分辨不出 %d，反例成立 %d。' % (vs.count('排除'), vs.count('分辨不出'), vs.count('反例成立')))
print(r'  只要"反例成立"是 0，E47 的负读数就全是配准残差与 $\kappa(V)$ 放大出的舍入，#23 第一节那张表不改，')
print(r'  E46 的命题保留"块对角子族"这个限定词，而且这里给出的是逐格实测的区间，不是"看起来像噪声"。')
print(r'  留给别人的门槛（可判定）：$\rho\ne0$ 想买车票，报出的负幅度必须明显大于 $|dI/dx|\,\delta x$，')
print(r'  在本植物 $x_t=10^{-6}$ 处即 $1.4\times10^{6}\times10^{-15}\approx10^{-9}$ bit。比这小的"改进"在这个参数化下无法证明。')
