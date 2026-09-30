r"""E75 = 前沿的**闭式下界**：AM-GM 一步推出 C(V)·2^{-2I/r}，直接顶验收线。

背景：E71/E73/E74 三条凸松弛/割线腿全部不够用（E74 第一版还踩了"1 维特征值分格不是覆盖"
的真 bug —— 硬门当场抓住：L 冲到 58.57 > 可达值 42.04）。E72 显示真对偶间隙几乎为零，
所以瓶颈是"怎么把 +mu/2 logdet G 这个凹项处理干净"。换一条不需要凸化的路。

== 命题（冻结先验 X 上的精确分解 + AM-GM）==
Z 是 V 的正交基，G := Z^t X Z，设计的信息 M >= 0（S = Z M Z^t），Y := (M^{-1}+G)^{-1}，
T := G^{1/2} Y G^{1/2}，则 0 <= T <= I，且
    D = j_c + tr(Theta X) - tr(Gamma T),   Gamma := G^{-1/2}(Z^t X Theta X Z)G^{-1/2} >= 0,
    I = -1/2 log_2 det(I - T).
证明要点：P = (X^{-1}+S)^{-1} = X - X Z Y Z^t X（求逆引理），所以
    tr(Theta P) = tr(Theta X) - tr(Z^t X Theta X Z Y) = tr(Theta X) - tr(Gamma T)；
    det P = det X * det(I - G Y) = det X * det(I - T)（行列式引理），所以率 = -1/2 log2 det(I-T)。
令 R := I - T，则 tr(Gamma R) >= r * det(Gamma R)^{1/r}（AM-GM 用在 Gamma^{1/2} R Gamma^{1/2} 的
特征值上），det R >= 2^{-2I}（率约束），det(Gamma R) = det Gamma * det R，于是
    D >= j_c + tr(Theta X) - tr Gamma + r * det(Gamma)^{1/r} * 2^{-2I/r}.
再把前两项认出来：X Z I^{-1} Z^t X 就是"T=I（V 上完全传感）"的后验，所以
    j_c + tr(Theta X) - tr Gamma = j_c + tr(Theta P_wall(X)) =: W_X，
而 X -> P_wall(X) = (X^{-1} + Z Z^t)^{-1} 单调，可达先验都满足 X >= X_w（墙的先验），故
    W_X >= j_c + tr(Theta P_wall(X_w)) =: D_min(V)      【墙】
给出**前沿的闭式下界**
    D(I) >= D_min(V) + C(V) * 2^{-2I/r},   C(V) := r * det(Gamma(X_w))^{1/r}.        (*)
最后一步用到 det Gamma(X) >= det Gamma(X_w)（X >= X_w 时）—— 这条单调性我**没有证明**，
本脚本的职责就是判它：在所有实测前沿设计上验 (*) 是否成立、以及 det Gamma 是否真的随 X 单调。
若 (*) 在所有点成立且 (X, detGamma) 同向，则 (*) 是可发表的定理雏形；若出现反例，则 (*) 作废，
只能保留"冻结 X 版本"（对每个设计自身成立，无需单调性）。

判据（写在前）：
 [0] 墙复现：D_min(red)=46.1231、D_min(blue)=36.6032（E68）。对不上=我这套 X_w 迭代有 bug。
 [1] 逐个前沿点验 (*) 的红/蓝两侧：bound(I) <= 实测前沿值（越界=命题或单调性错）。
 [2] 验收线：蓝平面 I=3 的 (*) 界 > 43.80 => Remark 2 从"相对上界排除"升级为**闭式认证不可行**；
     落在 43.711~43.80 => 在读图带宽内，符号仍不可读；<= 43.711 => 不成立。
 [3] det Gamma 单调性抽检：沿可达先验（各前沿设计的 X）与墙先验比 det Gamma。
"""
import sys, time
sys.stdout.reconfigure(encoding='utf-8')
import numpy as np
from scipy.linalg import schur

np.set_printoptions(precision=4, suppress=True, linewidth=170)
_src = open('p0/exp_c_audit.py', encoding='utf-8').read().split("print(r'== E53")[0]
_ns = {'__name__': 'p'}
exec(compile(_src, 'p0/exp_c_audit.py[preamble]', 'exec'), _ns)
A, W, TH, JC, n, sym = _ns['A'], _ns['W'], _ns['TH'], _ns['JC'], _ns['n'], _ns['sym']
ln2 = np.log(2.0)
print('== E75：AM-GM 闭式前沿下界 ==')


def stat(S, maxit=60000):
    """设计的定点 (X_prior, P_post)；不收敛返回 None。"""
    Pt = W.copy()
    for _ in range(maxit):
        Pm = sym(np.linalg.inv(np.linalg.inv(Pt) + S))
        Pn = sym(A @ Pm @ A.T + W)
        if not np.all(np.isfinite(Pn)) or np.max(np.abs(Pn)) > 1e14:
            return None
        if np.max(np.abs(Pn - Pt)) < 1e-13 * max(1.0, np.max(np.abs(Pn))):
            return Pt, Pm
        Pt = Pn
    return None


def wall_pair(Z, maxit=60000):
    """V 上完全传感的极限：X_{k+1} = A (X_k - X_k Z G^-1 Z^t X_k) A^t + W。"""
    X = W.copy()
    for _ in range(maxit):
        G = sym(Z.T @ X @ Z)
        P = sym(X - X @ Z @ np.linalg.solve(G, Z.T @ X))
        Xn = sym(A @ P @ A.T + W)
        if not np.isfinite(Xn).all():
            return None, None
        if np.max(np.abs(Xn - X)) < 1e-13 * max(1.0, np.max(np.abs(Xn))):
            X = Xn
            break
        X = Xn
    G = sym(Z.T @ X @ Z)
    P = sym(X - X @ Z @ np.linalg.solve(G, Z.T @ X))
    return X, P


def gamma_of(Z, X):
    G = sym(Z.T @ X @ Z)
    Lam = sym(Z.T @ X @ TH @ X @ Z)
    w, V = np.linalg.eigh(G)
    Gi2 = (V / np.sqrt(np.clip(w, 1e-300, None))) @ V.T
    return sym(Gi2 @ Lam @ Gi2), sym(Lam), G


Ts, U = schur(A, output='real', sort=lambda a: abs(a) < 1.0)[:2]
Z2 = U[:, n - 2:]
Z3 = U[:, n - 3:]
RHO_TH = {2.0: (0.0, 138.0), 2.5: (0.0, 132.0), 3.0: (0.0, 130.0), 4.0: (0.01, 128.0),
          4.90: (0.1, 128.0)}
RED = {2.0: 62.8975, 2.5: 54.2924, 3.0: 50.6023, 4.0: 47.9396, 4.90: 47.1006}
BLUE = {2.0: 59.5611, 3.0: 45.4537, 4.0: 40.8072}
WALLS = {'red': 46.1231, 'blue': 36.6032}


def shape(Z, rho, th):
    Rv = np.array([[np.cos(th), -np.sin(th)], [np.sin(th), np.cos(th)]]) if Z.shape[1] == 2 \
        else np.eye(3)
    M = sym(Rv @ np.diag([1.0, rho]) @ Rv.T) if Z.shape[1] == 2 else np.eye(3)
    return sym(Z @ M @ Z.T)


print('\n[0] 墙复现（对着 E68）')
Xw = {}
for tag, Z in (('red', Z2), ('blue', Z3)):
    Xw[tag], Pw = wall_pair(Z)
    Dm = JC + float(np.trace(TH @ Pw))
    print('  %-5s D_min(本次)=%10.5f  E68=%9.4f  差 %+.5f | lam(X_w)=%s'
          % (tag, Dm, WALLS[tag], Dm - WALLS[tag],
             np.array2string(np.linalg.eigvalsh(Xw[tag])[:4], precision=3)))

print('\n[3] det(Gamma) 单调性抽检（红平面各前沿设计的 X vs 墙先验）')
Gw_red, _, _ = gamma_of(Z2, Xw['red'])
detGw = float(np.linalg.det(Gw_red))
print('  det Gamma(X_wall_red) = %.6e' % detGw)
for It in RED:
    rho, th_deg = RHO_TH[It]
    S0 = shape(Z2, rho, np.radians(th_deg))
    lo, hi = -3.0, 13.0
    for _ in range(46):
        m = 0.5 * (lo + hi)
        sp = stat((10.0 ** m) * S0)
        if sp is None:
            hi = m; continue
        _, Pm = sp
        I = 0.5 * (np.linalg.slogdet(sp[0])[1] - np.linalg.slogdet(Pm)[1]) / ln2
        if I >= It:
            hi = m
        else:
            lo = m
    Sp = (10.0 ** hi) * S0
    sp = stat(Sp)
    Xd = sp[0]
    Gd, _, _ = gamma_of(Z2, Xd)
    detd = float(np.linalg.det(Gd))
    D = JC + float(np.trace(TH @ sp[1]))
    I = 0.5 * (np.linalg.slogdet(Xd)[1] - np.linalg.slogdet(sp[1])[1]) / ln2
    r = 2
    C = r * detd ** (1.0 / r)
    Cw = r * detGw ** (1.0 / r)
    print('  I=%4.2f: 实测 D=%9.4f (I实=%.4f) | detGamma 设计/墙 = %.3e / %.3e %s'
          % (It, D, I, detd, detGw, '设计>=墙 OK' if detd >= detGw * (1 - 1e-6) else '**反例**'))
    print('        冻结X版界 D >= %9.4f | 墙版(*) C=%7.4f -> %9.4f (+%.4f)'
          % (JC + float(np.trace(TH @ (Xd - Xd @ Z2 @ np.linalg.solve(Z2.T @ Xd @ Z2, Z2.T @ Xd))))
             + C * 2.0 ** (-2.0 * I / r), Cw, WALLS['red'] + Cw * 2.0 ** (-2.0 * It / r),
             RED[It] - (WALLS['red'] + Cw * 2.0 ** (-2.0 * It / r))))

print('\n[1]/[2] 蓝平面 I=3 对验收线')
Gw_blue, _, _ = gamma_of(Z3, Xw['blue'])
Cw = 3.0 * float(np.linalg.det(Gw_blue)) ** (1.0 / 3.0)
b3 = WALLS['blue'] + Cw * 2.0 ** (-2.0 * 3.0 / 3.0)
print('  C(blue) = %.5f -> bound(I=3) = %.5f | 实测蓝上界 45.4537 | vs 43.711 %s | vs 43.80 %s'
      % (Cw, b3, '超过' if b3 > 43.711 else '未超过',
         '超过 => 闭式认证不可行' if b3 > 43.80 else '未超过'))
for It in BLUE:
    print('  I=%4.2f: 界 %9.4f | 可达上界 %9.4f | %s' % (It, WALLS['blue'] + Cw * 2.0 ** (-2 * It / 3),
                                                          BLUE[It],
                                                          'OK' if WALLS['blue'] + Cw * 2.0 ** (-2 * It / 3)
                                                          <= BLUE[It] + 1e-6 else '**越界：命题错**'))
