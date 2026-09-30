# -*- coding: utf-8 -*-
r"""e151（回应 §67-B 的"e150 是玩具 / 真实对象就是 $\Theta$ 谱上的水填"）：
两条独立审计，都在我车道自算，不读 C 的解。

 [W1] **真空检验**：主角度对"满秩"这一档是否携带信息。
      取锚点 $\Theta$，造两个满秩 SPD：一个与 $\Theta$ 可对易（在 $\Theta$ 基里对角），
      一个随机不可对易。若"$\Theta$ 前 $k{=}4$ 特征子空间的主角度"给两者同一个数（$0.0°$），
      那个读数就**不是**判据。同时打印两个**能**分开它们的量：
      对易子相对范数 $|[S,\Theta]|_F/(|S|_F|\Theta|_F)$、$\Theta$ 基里的非对角质量 $|\mathrm{offdiag}(U^\top SU)|_F/|S|_F$。
 [W2] **文献条件清点**（1912.07640 Proposition 3 + Theorem 5 的字面前提）：
      该线把"结构简化到逆水填"写成**条件定理**，三选一：
      (i) $A=\alpha I_p$ 且 $\bar\Sigma\succeq0$；(ii) $A$ 实对称且 $\bar\Sigma=\sigma^2I_p$；(iii) $A=\bar\Sigma\succ0$
      $\\Rightarrow$ "只有满足其一，才有 $(A,\Sigma_\xi,\bar\Sigma)$ 两两对易，进而 $(\Sigma_\xi,\Pi_\xi)$ 对易"（式 (60) 之后，其 Thm 5 明写 "Suppose that one of the conditions of Proposition 3 hold"）。
      对我六株逐一**打印**这三条成不成立（$\bar\Sigma$ 取本车道代价侧的 $W$，$\Theta$ 取 $K^\top(R+B^\top P_cB)K$）。
 [W3] 打印每株 $\Theta$ 的谱与相邻间隙，供"$6\text{–}8°$ 偏离是否落在近似退化子空间内"判断。

预注册判据：
 [Z1] 若 W1 里两个满秩矩阵的主角度相同而对易子不同 $\\Rightarrow$ §66-C/§67-B 的"$k{=}4$ 主角度 $0.0°$ $\\Rightarrow$ 严格对角"
      这一步**不成立**（真空判据），必须换成对易子或非对角质量。
 [Z2] 若六株里满足 Prop 3 任一条件的株数 $=0$ $\\Rightarrow$ "水填在 $\Theta$ 谱上"**不能**从这条文献借来，
      C 侧要么自己证（新结构定理），要么降级为"$S^*$ 的支撑与 $\Theta$ 前 $k$ 子空间夹角小"。
 [Z3] 若某株满足任一条件 $\\Rightarrow$ 点名是哪株、哪条，并只在该株上谈"文献已知"。
"""
import sys
import numpy as np
from scipy.linalg import solve_discrete_are
sys.stdout.reconfigure(encoding='utf-8')

_src = open('p0/exp_c_audit.py', encoding='utf-8').read().split("print(r'== E53")[0]
_ns = {'__name__': 'p'}
exec(compile(_src, 'p0/exp_c_audit.py[preamble]', 'exec'), _ns)
A0, B0, W0 = _ns['A'], _ns['B'], _ns['W']


def sym(X):
    return 0.5 * (X + X.T)


def ctrl_full(Ax, Bx, Wx, Qt, Rx):
    Pc = sym(solve_discrete_are(Ax, Bx, sym(Qt), Rx))
    K = np.linalg.solve(Rx + Bx.T @ Pc @ Bx, Bx.T @ Pc @ Ax)
    return Pc, K, sym(K.T @ (Rx + Bx.T @ Pc @ Bx) @ K), float(np.trace(Wx @ Pc))


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
    Pc, K, TH, JC = ctrl_full(A, B, W, np.eye(m), np.eye(m))
    return A, B, W, TH, JC, m


def angles_to(S, U, k):
    """range(S) 与 U 前 k 个特征向量张成子空间的主角度（度）。"""
    Qs, _ = np.linalg.qr(S)                      # range(S) 的正交基（S 满秩时是 R^n 的任意基）
    r = np.linalg.matrix_rank(S, tol=1e-9 * np.linalg.norm(S))
    Qs = Qs[:, :r]
    Sv = U[:, :k]
    c = np.linalg.svd(Qs.T @ Sv, compute_uv=False)
    c = np.clip(c, -1.0, 1.0)[:min(r, k)]
    return np.degrees(np.arccos(c))


def comm_rel(S, TH):
    n = np.linalg.norm(S @ TH - TH @ S, 'fro')
    return n / (np.linalg.norm(S, 'fro') * np.linalg.norm(TH, 'fro'))


def offdiag_rel(S, U):
    M = U.T @ S @ U
    off = M - np.diag(np.diag(M))
    return np.linalg.norm(off, 'fro') / np.linalg.norm(S, 'fro')


print('== [W1] 真空检验（锚点 Theta，满秩两例） ==')
A, B, W, TH, JC, n = make_plant('anchor', 0)
w, U = np.linalg.eigh(TH)
w = w[::-1]; U = U[:, ::-1]                      # 降序
rng = np.random.default_rng(7)
G = rng.normal(0, 1, (n, n))
S_rand = sym(G @ G.T + 0.2 * np.eye(n))          # 一般位置：不与 Theta 对易
d = np.array([118.5, 17.58, 4.563, 8.73e-2])     # 用 C 在 D=32.00 报的谱，构造"对角版"
S_diag = U @ np.diag(d) @ U.T
for name, S in [('S_diag (Theta 基里对角)', S_diag), ('S_rand (一般位置满秩)', S_rand)]:
    a4 = angles_to(sym(S + 1e-12 * np.eye(n)), U, n)
    print('%-28s rank=%d  主角度(前%d)=%s  max=%.4f deg | 对易子=%.4e | 非对角质量=%.4f'
          % (name, np.linalg.matrix_rank(S), n, np.array2string(a4, precision=4), a4.max(),
             comm_rel(S, TH), offdiag_rel(S, U)))
a_d = angles_to(S_diag, U, n); a_r = angles_to(S_rand, U, n)
print('两例主角度向量的最大差 = %.3e deg（两者 max 均为 %.3e deg）；对易子之比 = %.3e 倍'
      % (np.abs(a_d - a_r).max(), max(a_d.max(), a_r.max()),
         comm_rel(S_rand, TH) / max(comm_rel(S_diag, TH), 1e-300)))
print('[Z1] 判据：两例 max 主角度均 <1e-3 deg 而对易子相差 >1e10 倍 $\\Rightarrow$ 满秩档的主角度不是判据。'
      ' 结论 = %s'
      % ('成立（真空判据）' if max(a_d.max(), a_r.max()) < 1e-3 and
         comm_rel(S_rand, TH) / max(comm_rel(S_diag, TH), 1e-300) > 1e10 else '不成立'))

print('\n== [W2] 1912.07640 Prop 3 三条件在我六株上的清点 ==')
PLANTS = [('anchor', ('anchor', 0)), ('rand-1', ('rand4', 1)), ('rand-2', ('rand4', 2)),
          ('rand-3', ('rand4', 3)), ('rand-4', ('rand4', 4)), ('big-6', ('big6', 11))]
hit = []
for nm, spec in PLANTS:
    Ai, Bi, Wi, THi, JCi, mi = make_plant(*spec)
    # (i) A = alpha*I
    a_i = np.linalg.norm(Ai - (np.trace(Ai) / mi) * np.eye(mi)) / np.linalg.norm(Ai)
    c1 = a_i < 1e-9
    # (ii) A 实对称 且 Sigma_bar = sigma^2 I（本车道 Sigma_bar 取过程噪声 W）
    symA = np.linalg.norm(Ai - Ai.T) / np.linalg.norm(Ai)
    sc = np.linalg.norm(Wi - (np.trace(Wi) / mi) * np.eye(mi)) / np.linalg.norm(Wi)
    c2 = (symA < 1e-9) and (sc < 1e-9)
    # (iii) A = Sigma_bar
    c3 = np.linalg.norm(Ai - Wi) / max(np.linalg.norm(Ai), 1e-30) < 1e-9
    # 附加：Theta 是否标量阵（对角水填在 Theta 基之外的另一类退化）
    tc = np.linalg.norm(THi - (np.trace(THi) / mi) * np.eye(mi)) / np.linalg.norm(THi)
    print('%-8s n=%d |A-aI|/|A|=%.3e |A-A^T|/|A|=%.3e |W-sI|/|W|=%.3e |Theta-tI|/|Theta|=%.3e '
          '=> (i)%s (ii)%s (iii)%s' % (nm, mi, a_i, symA, sc, tc, c1, c2, c3))
    if c1 or c2 or c3:
        hit.append(nm)
print('满足 Prop 3 任一条件的株：%s（共 %d/6）' % (hit if hit else '无', len(hit)))

print('\n== [W3] 各株 Theta 谱与相邻间隙（判 6-8 度偏离是否落在近退化子空间） ==')
for nm, spec in PLANTS:
    Ai, Bi, Wi, THi, JCi, mi = make_plant(*spec)
    ww = np.linalg.eigvalsh(THi)[::-1]
    gaps = (ww[:-1] - ww[1:]) / ww[:-1]
    print('%-8s spec(Theta)=%s 相对间隙=%s' % (nm, np.array2string(ww, precision=4),
                                               np.array2string(gaps, precision=4)))
print('\n== 附：8 度倾斜的能量漏出 ==')
for degv in [3.0, 3.4, 5.0, 6.0, 6.8, 8.0]:
    print('theta=%.1f deg -> sin^2=%.5f（即该方向上有 %.3f%% 的能量落在被砍子空间）'
          % (degv, np.sin(np.radians(degv)) ** 2, 100 * np.sin(np.radians(degv)) ** 2))
