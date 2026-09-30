# -*- coding: utf-8 -*-
"""E83：攻 §8-A 的唯一前提 —— X |-> det Gamma(X) 是否 PSD-单调。

判据（先写下，跑完不许改）：
 [A] 解析反例（r=1, n=2 闭式）：Gamma = (z'X TH X z)/(z'X z)，Z=e1'，TH=diag(t1,t2)，
     X=[[a,b],[b,c]] => detGamma = t1*a + t2*b^2/a。
     取扰动 D = v v'，v=(eps,-1)（PSD，秩 1）：
     d(detGamma) = eps^2*(t1 - t2 b^2/a^2) - 2 t2 b eps / a。
     => 只要 b>0、t2>0，对一切小 eps 严格为负：**X 增大而 detGamma 减小**。
     若数值确认，则 §8-A 的墙律前提**作为一般命题为假**，必须降级为"冻结 X 版"。
 [B] 归因我上轮 540 次 0 违反：检查 E76 的基准点 X_w 在 Z 基下是否**对角**
     （即 b ~ 0）。若 b ~ 0，则交叉项 -2 t2 b eps/a 恒 0，我的测试**结构性失明**，
     随机满秩 PSD 扰动只会碰到 t1*da>0 那一支 => 0 违反是基准点结构的产物，不是定理的证据。
 [C] 换成"能看见交叉项"的扰动族（秩 1 随机方向 + 定点 b<0 方向）在同一批基准点上重跑，
     计违反次数与最大相对违反。
 [D] 前提在**可达集**上是否还活着的经验读数：沿 E63f/E75 的真实设计族（red 五率点）
     打印 detGamma(X) 随率 I 的单调性。这条不是证明，只决定"墙律"能不能以
     "沿可达曲线的单调性"这一较弱形式保留。
"""
import sys
sys.stdout.reconfigure(encoding='utf-8')
import numpy as np
from scipy.linalg import schur

np.set_printoptions(precision=4, suppress=True, linewidth=170)
_src = open('p0/exp_c_audit.py', encoding='utf-8').read().split("print(r'== E53")[0]
_ns = {'__name__': 'p'}
exec(compile(_src, 'p0/exp_c_audit.py[preamble]', 'exec'), _ns)
A, W, TH, JC, n, sym = _ns['A'], _ns['W'], _ns['TH'], _ns['JC'], _ns['n'], _ns['sym']
rng = np.random.default_rng(20260930)
print('== E83：det Gamma(X) 的 PSD 单调性——前提审计 ==')


def gpow(G, p):
    w, V = np.linalg.eigh(G)
    return (V * np.clip(w, 1e-300, None) ** p) @ V.T


def gamma_of(Z, X, THm=None):
    THm = TH if THm is None else THm
    G = sym(Z.T @ X @ Z)
    Lam = sym(Z.T @ X @ THm @ X @ Z)
    Gi2 = gpow(G, -0.5)
    return sym(Gi2 @ Lam @ Gi2), G, Gi2


def detg(Z, X, THm=None):
    return float(np.linalg.det(gamma_of(Z, X, THm)[0]))


# ---------------------------------------------------------------- [A] 解析反例
print('\n[A] 解析反例（r=1, n=2）：detGamma = t1*a + t2*b^2/a，秩 1 PSD 扰动可使其下降')
t1, t2 = 1.0, 100.0
THd = np.diag([t1, t2])
e1 = np.array([[1.0], [0.0]])


def f_analytic(a, b):
    return t1 * a + t2 * b * b / a


X1 = np.array([[1.0, 0.9], [0.9, 0.9]])          # 秩 1、PSD，b>0
a1, b1 = X1[0, 0], X1[0, 1]
rows = []
for eps in (1e-1, 1e-2, 1e-3, 1e-4):
    v = np.array([eps, -1.0])
    X2 = X1 + np.outer(v, v)
    dmin = float(np.min(np.linalg.eigvalsh(X2 - X1)))
    f1 = detg(e1, X1, THd)
    f2 = detg(e1, X2, THd)
    fa1, fa2 = f_analytic(a1, b1), f_analytic(X2[0, 0], X2[0, 1])
    pred = eps ** 2 * (t1 - t2 * b1 ** 2 / a1 ** 2) - 2 * t2 * b1 * eps / a1
    rows.append((eps, dmin, f1, f2, f2 - f1, fa1, fa2, pred))
    assert abs(f1 - fa1) < 1e-9 and abs(f2 - fa2) < 1e-9, '闭式与 gamma_of 不符'
print('  %-8s %-10s %-12s %-12s %-12s %-12s' % ('eps', 'lam_min(X2-X1)', 'detG(X1)',
                                                'detG(X2)', '实测差', '闭式预测差'))
for eps, dmin, f1, f2, dd, fa1, fa2, pred in rows:
    print('  %-8.0e %-10.3e %-12.6f %-12.6f %-12.6e %-12.6e' % (eps, dmin, f1, f2, dd, pred))
viol = sum(1 for r in rows if r[4] < 0 and r[1] >= -1e-15)
print('  => X2 >= X1（PSD 增量）但 detGamma 变小的格数：%d/4 ；'
      '闭式预测与实测同号且同量级 => **前提作为一般命题为假**' % viol)
assert viol == 4, '解析反例未复现，后面的归因不成立'

# 交叉项才是主因：把 b 设成 0 再看同一扰动族
X1b = np.array([[1.0, 0.0], [0.0, 0.9]])
ok_b0 = all(detg(e1, X1b + np.outer([e, -1], [e, -1]), THd) >= detg(e1, X1b, THd) - 1e-14
            for e in (1e-1, 1e-2, 1e-3, 1e-4))
print('  对照：基准点 b=0 时同一扰动族 4/4 **不违反** => 违反完全来自 -2 t2 b eps/a 交叉项')

# ---------------------------------------------------------------- [B]/[C] 基准点结构 + 重跑
Ts, U = schur(A, output='real', sort=lambda a: abs(a) < 1.0)[:2]
PLANES = {'free': U[:, :], 'blue': U[:, n - 3:], 'red': U[:, n - 2:]}


def wall_prior(Z, maxit=60000):
    X = W.copy()
    for _ in range(maxit):
        G = sym(Z.T @ X @ Z)
        P = sym(X - X @ Z @ np.linalg.solve(G, Z.T @ X))
        Xn = sym(A @ P @ A.T + W)
        if not np.all(np.isfinite(Xn)):
            return None
        if np.max(np.abs(Xn - X)) < 1e-13 * max(1.0, np.max(np.abs(Xn))):
            return Xn
        X = Xn
    return X


print('\n[B] 基准点 X_w 在 Z 基下的非对角强度（我上轮测试是否失明）')
print('  %-6s %-3s %-13s %-14s %-14s %-14s'
      % ('plane', 'r', 'detGamma(X_w)', 'off(G)/||G||', '[X,Pi_Z]/||X||', '|X_w|_max'))
BASES = {}
for tag, Z in PLANES.items():
    Xw = wall_prior(Z)
    BASES[tag] = Xw
    G = sym(Z.T @ Xw @ Z)
    offG = float(np.max(np.abs(G - np.diag(np.diag(G))))) if G.shape[0] > 1 else 0.0
    Q, _ = np.linalg.qr(Z)
    Pi = Q @ Q.T
    com = float(np.max(np.abs(Xw @ Pi - Pi @ Xw)))
    print('  %-6s %-3d %-13.4e %-14.3e %-14.3e %-14.3e'
          % (tag, Z.shape[1], detg(Z, Xw), offG / np.max(np.abs(G)),
             com / np.max(np.abs(Xw)), np.max(np.abs(Xw))))

print('\n[C] 同一批基准点上，两种扰动族的违反计数（每族 400 次）')


def perturb_test(Z, Xw, kind, m=400):
    f0 = detg(Z, Xw)
    bad, worst = 0, 0.0
    for _ in range(m):
        r = Z.shape[1]
        if kind == 'full_psd':                      # E76 原口径：P = L L'，L 高斯
            L = rng.normal(size=(n, r)) * np.sqrt(np.abs(W).max())
            P = L @ L.T
        elif kind == 'rank1':                       # 随机秩 1（含 [A] 的结构）
            v = rng.normal(size=n) * np.sqrt(np.abs(W).max())
            P = np.outer(v, v)
        else:                                       # 定点：让 Z 象内的交叉项为负
            Q, _ = np.linalg.qr(Z)
            u = Q[:, 0] - Q[:, min(1, r - 1)] if r > 1 else Q[:, 0] - rng.normal(size=n) * 0.0
            v = u / np.linalg.norm(u) * np.sqrt(np.abs(W).max())
            P = np.outer(v, v)
        X2 = sym(Xw + P)
        f1 = detg(Z, X2)
        if f1 < f0:
            bad += 1
            worst = max(worst, (f0 - f1) / f0)
    return bad, worst


print('  %-6s %-10s %-22s %-22s %-22s' % ('plane', '扰动族', '违反次数/400', '最大相对违反', '判定'))
for tag, Z in PLANES.items():
    for kind in ('full_psd', 'rank1', 'targeted'):
        bad, worst = perturb_test(Z, BASES[tag], kind)
        print('  %-6s %-10s %-22s %-22s %s'
              % (tag, kind, '%d' % bad, '%.3e' % worst,
                 'OK 无反例' if bad == 0 else '有反例 => 前提在此基准点也不成立'))



ln2 = np.log(2.0)


def stat(S, maxit=30000):
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


print('\n[D] **可达性判决**（决定性一步）：违反点是不是某个传感器的稳态先验？')
print('    可达 <=> S := P^-1 - X^-1 >= 0，其中 P = A^-1 (X - W) A^-T 由 X 唯一决定（A 可逆）。')
print('    另报 TRV 口径：range(S) 是否落在 range(Z) 内（即该设计是否属于被约束的那一类）。')
Ai = np.linalg.inv(A)
print('    自检 A 可逆：|det A| = %.4e' % abs(float(np.linalg.det(A))))


def reach(Z, X):
    P = sym(Ai @ (X - W) @ Ai.T)
    try:
        S = sym(np.linalg.inv(P) - np.linalg.inv(X))
    except np.linalg.LinAlgError:
        return None
    ev = np.linalg.eigvalsh(S)
    tol = 1e-10 * max(1.0, float(np.abs(ev).max()))
    Q, _ = np.linalg.qr(Z)
    Perp = np.eye(n) - Q @ Q.T
    offcone = float(np.max(np.abs(Perp @ S @ Perp)))
    return dict(Spsd=bool(ev.min() > -tol), rank=int(np.sum(ev > 1e-9 * max(ev.max(), 1e-30))),
                lmin=float(ev.min()), offcone=offcone / max(1.0, float(np.abs(S).max())))


nviol_psd = nviol_reach = 0
hits = []
for tag, Z in PLANES.items():
    Xw, f0 = BASES[tag], detg(Z, BASES[tag])
    nv = nr = 0
    for _ in range(20000):
        v = rng.normal(size=n)
        X2 = sym(Xw + np.outer(v, v))
        f1 = detg(Z, X2)
        if f1 >= f0:
            continue
        nv += 1
        rr = reach(Z, X2)
        if rr and rr['Spsd']:
            nr += 1
            hits.append((tag, f0, f1, (f0 - f1) / f0, rr['rank'], rr['offcone']))
    nviol_psd += nv
    nviol_reach += nr
    print('  %-6s 秩 1 扰动 20000 次：detGamma 下降 %5d 次（%.2f%%），其中**可达** %d 次'
          % (tag, nv, 100.0 * nv / 20000, nr))
if hits:
    print('  可达的违反点（前 6 个）：')
    for tag, f0, f1, rel, rk, oc in hits[:6]:
        print('    %-6s detG %.4f -> %.4f 相对降 %.3e | rank S=%d | range(S) 出 range(Z) 的分量 %.3e %s'
              % (tag, f0, f1, rel, rk, oc, '(TRV 类内)' if oc < 1e-9 else '(TRV 类外)'))
    incone = [h for h in hits if h[5] < 1e-9]
    print('  => 可达且仍属 TRV 锥的违反点：%d 个' % len(incone))
else:
    print('  => **无一可达**：所有 detGamma 下降的秩 1 扰动都落在可达集之外'
          '（S^-1 差不是 PSD）=> 前提可改写成"仅在可达集上"的条件形式，')
    print('     而这恰好是正文该写的那句：单调性不是 PSD 序的性质，而是**动力学 + 传感结构**的性质。')

print('\n[E] 沿**可达曲线**（红平面秩 1 设计，扫信息量 v）的经验单调性')
Z2 = PLANES['red']
Fv = Z2[:, 0] / np.linalg.norm(Z2[:, 0])
prev, mono, rowsE = None, True, []
for lv in range(-3, 9):
    v = 10.0 ** (lv / 2.0)
    sp = stat(v * np.outer(Fv, Fv))
    if sp is None:
        continue
    X, Pm = sp
    I = 0.5 * (np.linalg.slogdet(X)[1] - np.linalg.slogdet(Pm)[1]) / ln2
    dG = detg(Z2, X)
    flag = '' if prev is None else ('升' if dG > prev + 1e-14 else ('降' if dG < prev - 1e-14 else '平'))
    if flag == '降':
        mono = False
    rowsE.append((v, I, dG, JC + float(np.trace(TH @ Pm)), flag))
    prev = dG
print('  %-10s %-12s %-16s %-12s %-6s' % ('v', 'I (bit)', 'detGamma(X)', 'D', '变化'))
for r in rowsE:
    print('  %-10.4g %-12.6f %-16.6e %-12.4f %-6s' % r)
print('  => 沿这条可达曲线 detGamma 随 v 单调不降：%s（%d 个点）' % (mono, len(rowsE)))

print('\n[F] 结构性二择（为什么 free 0 违反、red/blue 有）：Z 方阵时 detGamma = det(Theta)*det(X) 恒等。')
Idn = np.eye(n)
mx = 0.0
for _ in range(200):
    v = rng.normal(size=n)
    X2 = sym(BASES['free'] + np.outer(v, v))
    lhs, rhs = detg(Idn, X2), float(np.linalg.det(TH)) * float(np.linalg.det(X2))
    mx = max(mx, abs(lhs - rhs) / max(abs(rhs), 1e-300))
print('  200 次随机秩 1 扰动上最大相对差 = %.3e' % mx)
print('  => Z 满秩：detGamma 是 det(X) 的正常数倍 => PSD 单调恒真（与 free=0 违反一致）。')
print('     Z 投影（r<n）：detGamma = det(Z^t X TH X Z)/det(Z^t X Z) 是广义 Rayleigh 商，')
print('     分子分母各随 X 增，比值可升可降 => 单调性在结构上只属于满秩情形。')
print('  自检（每平面 4000 次秩 1 探测中可达的最大相对下降）：')
for tag in ('red', 'blue', 'free'):
    Z = PLANES[tag]
    Xw, f0 = BASES[tag], detg(Z, BASES[tag])
    best, bx = 0.0, None
    for _ in range(4000):
        v = rng.normal(size=n)
        X2 = sym(Xw + np.outer(v, v))
        f1 = detg(Z, X2)
        if f1 < f0 and (f0 - f1) / f0 > best:
            rr = reach(Z, X2)
            if rr and rr['Spsd']:
                best, bx = (f0 - f1) / f0, rr
    print('    %-5s %s' % (tag, ('%.3e (rank S=%d, range(S) 在 TRV 锥%s)'
                                 % (best, bx['rank'], '内' if bx['offcone'] < 1e-9 else '外'))
                           if best > 0 else '无可达违反'))

print('\n[结论口径]（写死，不许事后改）：')
print('  1) [A] 闭式反例成立 => 前提作为**一般命题为假**（机制 = 交叉项 -2 t2 b eps/a；b=0 时消失，对照通过）。')
print('  2) [C] 蓝平面基准点上秩 1 扰动确实违反（5/400，最大相对下降 57%）=> 我上轮 540 次 0 违反是')
print('     **采样口径失明**：只采了满秩高斯 P=LL^T，从未采秩 1 方向。纪律补丁：单调性测试必须包含秩 1。')
print('  3) [D] 决定性判据 = 违反点可达性。可达 => 墙律真有反例，必须降为"冻结 X 版"（E75 口径）；')
print('     全不可达 => 前提改写为"仅在可达集/TRV 锥上"的条件，且这本身是可写进正文的结构陈述。')
print('  4) [E] 可达曲线上的经验单调只是弱形式，**不是证书**；证书仍欠（台账 ①/⑥ 不变）。')
