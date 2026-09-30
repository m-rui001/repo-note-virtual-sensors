r"""E95b = 第一条原理核对"D(S) = JC + tr(TH * P_post)"到底是不是真实的闭环 LQR 平均代价。

为什么必须做：E95 想把 Kostina-Hassibi 的率侧反证界反解成代价侧下界，那一步成立的**唯一**
前提是"我车道量的 b 和 KH 式 (8)/(30) 里的 b 是同一个数"。E95 里我随手写的联合 Lyapunov
把 E[u^t R u] 写成了 tr(R*E[xhat xhat^t])（少了 K），全信息极限当场对不上
（22.23 vs 31.48）。这个脚本用干净的实现重做一遍，并且顺便把注记正文的主 measurable
从第一性原理坐实。

== 实现（不借任何道岔里的公式）==
给定设计 S >= 0、range S ⊆ span(Z)：写 S = Z T Z^t（T = Z^t S Z，r x r 正定），
等效传感实现取 C = Z^t（只测 V 块），观测噪声 V_obs = T^{-1}，则 C^t V_obs^{-1} C = S 成立。
稳态 Kalman 增益 L = Pm C^t (C Pm C^t + V_obs)^{-1}，后验 Σ 与先验 Pm = A Σ A^t + W 自洽。
误差动态（后验误差）： e+ = (I - L C)(A e + w) - L v
状态动态（u = -K xhat，x = xhat + e）： x+ = (A - B K) x + B K e + w
=> 联合 F = [[A-BK, BK],[0, (I-LC)A]]，Qd = [[W, W(I-LC)^t],[... , (I-LC)W(I-LC)^t + L V_obs L^t]]
解 DLYAP（kron 形式），代价 = tr(Pxx) + tr(K^t R K * P_xhat_xhat)，其中
P_xhat_xhat = Pxx - Pxe - Pxe^t + Pee。

== 跑前写死的判据 ==
 [C1] 三种设计（全信息极限 / V=R^4 的等方向射线 @3bit / Schur-tail-3 平面 @3bit）都必须
      |D_lyap - (JC + tr(TH*P_post))| / D < 1e-8。
 [C2] 若 [C1] 失败，则 E95 的整张对照表作废（口径不同不能比），并在留言板上自首。
 [C3] 附带核对：率的两条式子 ½log det(Pm/Σ) 与 ½log det(I + V_obs^{-1} C Pm C^t) 必须相等。
"""
import sys
sys.stdout.reconfigure(encoding='utf-8')
import numpy as np
from scipy.linalg import schur

np.set_printoptions(precision=6, suppress=True, linewidth=170)
_src = open('p0/exp_c_audit.py', encoding='utf-8').read().split("print(r'== E53")[0]
_ns = {'__name__': 'p'}
exec(compile(_src, 'p0/exp_c_audit.py[preamble]', 'exec'), _ns)
A, B, W, R, TH, JC, n, sym, Kv = (_ns['A'], _ns['B'], _ns['W'], _ns['R'], _ns['TH'],
                                  _ns['JC'], _ns['n'], _ns['sym'], _ns['Kv'])
ln2 = np.log(2.0)
Ts, U = schur(A, output='real', sort=lambda a: abs(a) < 1.0)[:2]


def are_post(C, Vobs, tol=1e-14, maxit=400000):
    """自洽后验 Σ：Σ = (I-LC)(AΣA^t+W)(I-LC)^t + L Vobs L^t, L = PmC^t(CPmC^t+Vobs)^-1"""
    Sig = W.copy()
    for _ in range(maxit):
        Pm = sym(A @ Sig @ A.T + W)
        L = Pm @ C.T @ np.linalg.inv(sym(C @ Pm @ C.T + Vobs))
        M = np.eye(n) - L @ C
        SigN = sym(M @ Pm @ M.T + L @ Vobs @ L.T)
        if np.max(np.abs(SigN - Sig)) < tol * max(1.0, np.max(np.abs(SigN))):
            return SigN, Pm, L
        Sig = SigN
    return None, None, None


def dlyap(F, Qd):
    m = F.shape[0]
    P = np.linalg.solve(np.eye(m * m) - np.kron(F, F), Qd.reshape(-1, order='F'))
    return sym(P.reshape(m, m, order='F'))


def cost_lyap(Sig, L, C):
    M = np.eye(n) - L @ C
    F = np.block([[A - B @ Kv, B @ Kv], [np.zeros((n, n)), M @ A]])
    Qd = sym(np.block([[W, W @ M.T], [M @ W, M @ W @ M.T + L @ Vobs_glob @ L.T]]))
    rad = max(abs(np.linalg.eigvals(F)))
    P = dlyap(F, Qd)
    Pxx, Pxe, Pee = P[:n, :n], P[:n, n:], P[n:, n:]
    Uxx = sym(Pxx - Pxe - Pxe.T + Pee)
    return float(np.trace(Pxx) + np.trace(Kv.T @ R @ Kv @ Uxx)), float(rad), float(np.max(np.abs(sym(Pee - Sig))))


print('== E95b：闭环代价的第一条原理核对 ==')
cases = []
cases.append(('全信息极限 C=I, Vobs=1e-8 I', np.eye(n), 1e-8 * np.eye(n)))
for I0 in (3.0, 2.5):
    # 等方向射线：S = s Z Z^t，用二分找 s 使率 = I0（率随 s 单调增）
    for nm, Z in [('free_R4', U[:, :]), ('schur_tail_3', U[:, n - 3:])]:
        r = Z.shape[1]
        lo, hi = 1e-8, 1e8
        for _ in range(200):
            s = np.sqrt(lo * hi)
            T = s * np.eye(r)
            Sig, Pm, L = are_post(Z.T, sym(np.linalg.inv(T)))
            if Sig is None:
                hi = s; continue
            I_form = 0.5 * (np.linalg.slogdet(Pm)[1] - np.linalg.slogdet(Sig)[1]) / ln2
            if I_form < I0:
                lo = s
            else:
                hi = s
        T = s * np.eye(r)
        cases.append(('%s @I0=%.1f (s=%.6f)' % (nm, I0, s), Z.T, sym(np.linalg.inv(T))))

print('  设计                                JC+tr(TH*Sig)     Lyapunov        相对差      收敛半径  Sig自洽残差')
ok = True
for nm, C, Vobs in cases:
    Vobs_glob = Vobs
    Sig, Pm, L = are_post(C, Vobs)
    if Sig is None:
        print('  %-34s 发散（无稳态后验）' % nm); ok = False; continue
    D_form = JC + float(np.trace(TH @ Sig))
    D_ly, rad, res = cost_lyap(Sig, L, C)
    rel = abs(D_ly - D_form) / max(1.0, abs(D_form))
    I_a = 0.5 * (np.linalg.slogdet(Pm)[1] - np.linalg.slogdet(Sig)[1]) / ln2
    I_b = 0.5 * float(np.linalg.slogdet(np.eye(C.shape[0])
                  + np.linalg.solve(Vobs, C @ Pm @ C.T))[1]) / ln2
    print('  %-34s %13.6f %14.6f  %9.2e  %8.4f  %9.2e   率: %.6f vs %.6f'
          % (nm, D_form, D_ly, rel, rad, res, I_a, I_b))
    if rel > 1e-8:
        ok = False
print('\n[C1] 全部一致 ? %s   => E95 的口径前提%s' % (ok, '成立，KH 的界可直接同轴对照' if ok else '不成立，E95 表作废'))
