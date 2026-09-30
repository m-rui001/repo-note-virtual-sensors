r"""E95 = 把 Kostina–Hassibi (arXiv:1612.02126, "Rate-cost tradeoffs in control") 的
**率侧**反证界 (35) 反解成**代价侧**认证下界，直接对着我这轮证书线（E70/E71）的空缺来打。

== 为什么这条路在文献里是现成的（先写清楚口径）==
KH 的率–代价函数是 R(b) := inf{ 有向互信息 : LQR 平均代价 <= b }。他们给的是**下界**：
    R(b) >= ln|det A| + (n/2) ln( 1 + (det N M)^{1/n} / (b - b_min)^{1/n} )   …… (KH-35)
    b_min = tr(Sum_V S) + tr(Sum_est A M A^t ... ) 的形式；全观测时 Sum_est=0,
            N = Sum_V = W, 于是 b_min = tr(W Pc) = 我的 JC，(35) 退化成 (16)。
其中 S=Pc 是控制 DARE 的解，M := S - S B (R + B^t S B)^{-1} B^t S，
而 KH Remark: 我的 TH = K^t (R + B^t Pc B) K **恰好等于 S - M**。
=> 我的定点代价 D(S) = JC + tr(TH * P_post) 与他们的 b 是同一个量（Q=R=I，B 满秩）。

== 反解的方向（这一步是关键，也是唯一能让"率侧界"变成"代价侧证书"的动作）==
f(b) := 右端，随 b 单调不增。设计类包含关系：
    {本注记的线性信息设计 range S subset V} subset {KH 允许的一切因果编码方案}
所以任何被 KH 界住的代价，也界住我的：若 I(S) <= I0 且 D(S) = b，则
    f(b) <= R(b) <= I(S) <= I0   =>   b >= f^{-1}(I0)     （f 不增 => 反解方向不变号）
闭式反解（n=4）：
    L_KH(I0) = JC + [ (det(W M))^{1/4} / ( exp((I0_nat - ln|det A|)/2) - 1 ) ]^4 ,
    I0_nat = I0 * ln 2.

== 跑之前写死的判据（rule 25：先立判据再跑）==
 [J1] 一致性：L_KH(3 bit) 必须 <= 42.0427（无约束锥在 I=3 的已知最优）。
      若 L_KH(3) > 42.0427 ⇒ 我的口径与 KH 不同（单位/代价约定），本轮作废并自首。
 [J2] 证书：L_KH(3) > 43.80 ⇒ Remark 2 的"相对上界排除"升级为绝对不可行（蓝线在该平面内死）。
      注意这是**跨平面有效**的下界：受限类 ⊂ 无约束类，所以无约束类的下界对受限类同样成立，
      只是可能松。
 [J3] 对照：把 L_KH 与我 E71 的 sup_mu L = 41.485 逐 I0 比高低。
      谁大谁占优。若 KH 占优 ⇒ 我的拉格朗日腿整条可以删掉，改用文献界，代价是"证书来自文献"。
 [J4] 退化确认：KH 明写 (35) 非平凡需要 rank B = rank C = n。
      受限传感正是 rank C < n 的那一侧 ⇒ 打印本植物在 rank C=3 的退化值，
      坐实"文献里唯一的显式反证界恰好在我们的研究区域变成空洞"。这是定位句，不是我们的胜利。

单位约定：ln 一律自然对数，率的外部换算用 bit（除以 ln2）。entropy power 对高斯 V
就是 det(Sum_V)^{1/n}，与对数底无关，这点我核对过 KH 式 (12)-(13)。
"""
import sys
sys.stdout.reconfigure(encoding='utf-8')
import numpy as np
from scipy.linalg import solve_discrete_are

np.set_printoptions(precision=6, suppress=True, linewidth=170)
_src = open('p0/exp_c_audit.py', encoding='utf-8').read().split("print(r'== E53")[0]
_ns = {'__name__': 'p'}
exec(compile(_src, 'p0/exp_c_audit.py[preamble]', 'exec'), _ns)
A, B, W, R, TH, JC, n, sym, Kv, Pcmod = (_ns['A'], _ns['B'], _ns['W'], _ns['R'], _ns['TH'],
                              _ns['JC'], _ns['n'], _ns['sym'], _ns['Kv'], _ns['Pc'])
ln2 = np.log(2.0)

# ---- KH 的量：S=Pc, M_KH = TH（我第一轮写成 S-TH，下面 [身份核对] 把它当场否掉了）----
# 推导：scipy 的 Pc=S 满足 S = Q + A^tSA - A^tSB(R+B^tSB)^{-1}B^tSA，
#       而我的 TH = K^t(R+B^tSB)K, K=(R+B^tSB)^{-1}B^tSA
#       => TH = A^tSB(R+B^tSB)^{-1}B^tSA = KH 的 M = L^t(R+B^tSB)L，M_KH 是 Gram 形所以必然 PSD。
S = Pcmod
Mmat = sym(TH)
lamM = np.linalg.eigvalsh(Mmat)
lamW = np.linalg.eigvalsh(sym(W))
lndetA = float(np.linalg.slogdet(A)[1])
detWM = float(np.linalg.slogdet(W)[1] + np.linalg.slogdet(Mmat)[1])
print('== E95：KH(1612.02126) 反解成代价侧认证下界 ==')
print('核对 S-TH 是否等于 KH 的 M（应为半正定，且 = K^t(R+B^tSB)K）:')
print('   lam(M) =', lamM)
print('   由 KH 定义重算 M = S B(R+B^tSB)^{-1}B^t S 的差：',
      float(np.max(np.abs(Mmat - (S - S @ B @ np.linalg.solve(R + B.T @ S @ B, B.T @ S))))))
print('   b_min = JC = %.6f   ln|det A| = %.6f (= %.6f bit)   ln det(WM) = %.6f'
      % (JC, lndetA, lndetA / ln2, detWM))
print('   lam(W) =', lamW)
print('   [J4 前置] M>0 ? %s   W>0 ? %s   rank B = %d'
      % (bool((lamM > 1e-9).all()), bool((lamW > 1e-9).all()), np.linalg.matrix_rank(B)))


def L_kh(I_bit):
    """闭式反解 (KH-35)：率预算 I_bit 比特时的代价认证下界。"""
    In = I_bit * ln2
    den = np.exp((In - lndetA) / 2.0) - 1.0
    if den <= 0:
        return np.inf          # 率低于 ln|det A|：KH 只给出稳定性地板，反解给出 +inf
    return JC + (np.exp(detWM / n) / den) ** n


def L_kh_direct(I_bit):
    """数值反解（不用闭式），用来抓代数错误。"""
    from scipy.optimize import brentq
    f = lambda b: (lndetA + (n / 2.0) * np.log1p(np.exp(detWM / n) /
                  (b - JC) ** (1.0 / n)) - I_bit * ln2)
    lo, hi = JC + 1e-9, 1e7
    if f(lo) < 0:
        return JC
    if f(hi) > 0:
        return np.nan
    return float(brentq(f, lo, hi, xtol=1e-10, rtol=1e-14))


print('\n[J1/J2/J3] 逐率预算对比（41.485 = 我 E71 的 sup_mu L，42.0427 = 无约束锥 I=3 已知最优，')
print('           43.80 = Remark 2 的门槛，45.4537 = Powell 设计）')
print('  I0(bit)   L_KH闭式    L_KH数值   与 41.485 的差   判定')
for I0 in (1.168575, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0):
    a, b = L_kh(I0), L_kh_direct(I0)
    tag = []
    if np.isfinite(a):
        if a > 43.80:
            tag.append('>43.80 证书成立[J2]')
        if abs(I0 - 3.0) < 1e-9 and a > 42.0427 + 1e-6:
            tag.append('>42.0427 违反同率一致性[J1!]')
        elif abs(I0 - 3.0) > 1e-9:
            tag.append('无同率锚点，J1 不适用')
        tag.append('优于E71' if a > 41.485 else '劣于E71')
    else:
        tag.append('率低于 ln|detA|：只有地板')
    print('  %7.4f  %10.4f  %10.4f  %+12.4f      %s' % (I0, a, b, a - 41.485, ' | '.join(tag)))

# ---- [J4] 受限传感那一侧到底退化到什么程度 ----
print('\n[J4] 受限传感 = rank C < n。KH 的 N = K(CPC^t+Sum_W)K^t 是 innovation 协方差；')
print('     对 range 受限的平面设计，innovation 在第 (n-r) 个不可见方向上完全由先验承担，')
print('     det N 里含 0 特征值 => (det N M)^{1/n} = 0 => f(b) = ln|det A| 恒成立（界空洞）。')
from scipy.linalg import schur
Ts, U = schur(A, output='real', sort=lambda a: abs(a) < 1.0)[:2]
for name, idx in [('schur_tail_3', n - 3), ('free_R4', 0)]:
    Z = U[:, idx:]
    r = Z.shape[1]
    # 该平面能提供的 innovation：P_post 在 V 上任意小，V 正交方向上 innovation = 先验块
    # 取"最乐观"情形：V 上 innovation -> 0，则 N 的特征值 = (0 * r 个, 先验在 V^\perp 上的块)
    Pinf = W  # 极限先验下界（P >= W）
    Nv = sym(Pinf - Z @ Z.T @ Pinf @ Z @ Z.T)
    lamN = np.linalg.eigvalsh(Nv)
    dn = float(np.linalg.slogdet(Nv)[1])
    print('   %s (r=%d): lam(N_optimistic)=%s  ln det N=%.4f  => (detNM)^{1/n}=%.6f  界%s'
          % (name, r, np.array2string(lamN, precision=4), dn,
             np.exp((dn + np.linalg.slogdet(Mmat)[1]) / (4 * n)),
             '非平凡' if dn > -1e9 and lamN.min() > 1e-12 else '空洞'))


# ================= [身份核对] 我的 D 到底是不是 KH 的 b =================
# 判据（跑前写死）：用第一条原理的联合 Lyapunov 方程算闭环平均代价 E[x^t x + u^t u]
#   （控制器用 -K x_hat，观测按设计 S 给信息），与 JC + tr(TH * P_post) 比。
#   两者相对误差 < 1e-8 才算"同一口径"，KH 的界才能搬进来。
def true_avg_cost(P_post, Ceff=None):
    """联合 [x; e] Lyapunov。x+ = (A-BK)x + BK e + w ; e+ = (A-L A?) ...
    这里只用最干净的情形：后验协方差 = P_post 的稳态滤波器，且控制增益 K 与全信息相同
    （分离原理）。代价 E[x^t Q x + u^t R u], Q=I。"""
    L = sym(np.linalg.solve(sym(W - P_post + P_post), np.eye(n))) if False else None
    return None


def joint_check(P_post):
    """给定稳态后验协方差 Sigma=P_post，先验 Pm^- = A Sigma A^t + W。
    稳态滤波器增益（把信息矩阵当成"等效观测"）：L = Pm^- Z (Z^t Pm^- Z)^{-1}，
    于是 e+ = (A - L A^t... ) 用标准形式：e+ = (I-LC)A e + (I-LC)w - Lv,
    取等效观测矩阵 C = Z^t（只测 V 块，无额外噪声时 S=inf 的极限）。
    为口径核对，改用无噪声等效：先验 Pm、后验 Sigma=P_post、创新协方差 N = Pm - Sigma。
    闭环：x = xhat + e, u = -K xhat。
      x+ = (A-BK) x + BK e + w
      e+ = (A - L A) e + (I - L) w  其中 L 满足 (I-L)Pm(I-L)^t + L 0 L^t = Sigma
    直接解 L = I - Sigma Pm^{-1}（此时 e+ = L w 的极限不严格，故这里用创新形式：
      e+ = (I-L) (A e + w), (I-L) = Sigma Pm^{-1} => Cov = Sigma Pm^{-1}(A Sigma A^t + W)Pm^{-1}Sigma
    稳态 Sigma 自洽时 = P_post。"""
    Pm = sym(A @ P_post @ A.T + W)
    G = sym(P_post @ np.linalg.inv(Pm))          # (I-L) = Sigma Pm^{-1}
    Qe = sym(G @ Pm @ G.T)                       # 稳态误差传播协方差（自洽应等于 P_post）
    Ax = sym(A - A @ G.T)                        # 误差回进到状态方程的耦合项
    # x+ = A x - A e + ... 用 x = xhat + e: x+ = (A-BK)x + BK e + w, e+ = G(A e + w)
    Pxy = np.zeros((2 * n, 2 * n))
    # 解离散 Lyapunov: vec(P) = (I - kron(F,F))^{-1} vec(Qd)
    F = np.block([[A - B @ Kv, B @ Kv], [np.zeros((n, n)), G @ A]])
    Qd = np.block([[W, W @ G.T], [G @ W, G @ W @ G.T]])
    I2 = np.eye((2 * n) ** 2)
    P = np.linalg.solve(I2 - np.kron(F, F), Qd.reshape(-1, order='F')).reshape(2 * n, 2 * n, order='F')
    P = sym(P)
    Pxx, Pxe, Pee = P[:n, :n], P[:n, n:], P[n:, n:]
    U = -Kv @ (Pxx - Pxe - Pxe.T + Pee)
    cost = float(np.trace(Pxx) + np.trace(R @ U))
    return cost, float(np.trace(TH @ Pee)), float(np.max(np.abs(sym(Pee - P_post))))


print('\n== [身份核对] 第一条原理的闭环平均代价 vs JC + tr(TH*P_post) ==')
for nm, Sc in [('全信息极限 S=inf (P_post=0)', np.zeros((n, n))),
               ('零信息 S=0 (P_post=开环奇异解)', None)]:
    if Sc is None:
        Pp = W.copy()
        for _ in range(200000):
            Pn = sym(A @ Pp @ A.T + W)
            if np.max(np.abs(Pn - Pp)) < 1e-14 * max(1.0, np.max(np.abs(Pn))):
                break
            Pp = Pn
        Sc = Pp
    J = JC + float(np.trace(TH @ Sc))
    c1, c2, selfc = joint_check(Sc)
    print('   %-34s 我的 D=%12.6f   Lyapunov=%12.6f   差=%+.2e   自洽残差=%.2e'
          % (nm, J, c1, c1 - J, selfc))
    print('      分解：tr(TH*Pee)=%10.6f  JC=%10.6f' % (c2, JC))
