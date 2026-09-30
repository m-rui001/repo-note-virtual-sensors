r"""E98 = 把 2606.31396（新普查捞出的旗舰风险条）的判据搬到锚点上，逐条对表。

普查（`p0/e97_out.txt` §3）在旧"12 条全集"之外新捞出两条 2026 论文，标题与我们的核心 framing
直接撞车：`2606.31396`（How Much Sensing Information Is Needed to Control an Unstable Linear
System?）与 `2601.12782`（Sensing-Limited Control ... Under Nonlinear Observations）。
2606.31396 的 Theorem 1 给的是 **对任意观测律 p(y_t|x_t)（可非线性、可非高斯、可秩亏）** 的必要条件
    I(z→y) ≥ R_exp := Σ_{|λ_i(A)|≥1} log2|λ_i|
Proposition 1 把它延到加性过程噪声；Corollary 1 在**线性高斯**观测下把可达率算成
    I(x→y) = ½ log2 det(I + CΣCᵗ R⁻¹)，Σ 为稳态 Kalman 预测协方差（DARE 解）
—— 看着就是我车道 `design()` 里那条率公式。

⇒ 用数字回答四件事：
 [P1] 锚点 `R_exp` 多少 bit？与 Tanaka 无约束率地板 `1.1685 bit`、KH 秩亏退化右端
      `log2|det A| = 0.354774 bit` 谁强。若 `R_exp < 1.1685` ⇒ 他们的必要条件严格弱于 [4]
      的精确前沿，**撞不掉** §65-A 的证书腿。
 [P2] 他们的定理能否给出我们的"墙"（固定 F 下 $(j_c,j_c+\Phi_0(F))$ 整段不可行）？
      逐条列作用域，不许用"类似"混过去。
 [P3] 他们的前提 `(A,C)` detectable 是**假设**不是结论 ⇒ 对秩亏传感，只有过滤收敛的平面
      才落在 Th.1 内。用 ρ((I−LC)A) 门把"哪些平面在前提之外"量化。
 [P4] Cor.1 的率公式与我车道 `I` 是否同一件事（数值对表到机器精度）。
"""
import sys
sys.stdout.reconfigure(encoding='utf-8')
import numpy as np
from scipy.linalg import solve_discrete_are, schur
np.set_printoptions(precision=6, suppress=True, linewidth=170)

_src = open('p0/exp_c_audit.py', encoding='utf-8').read().split("print(r'== E53")[0]
_ns = {'__name__': 'p'}
exec(compile(_src, 'p0/exp_c_audit.py[preamble]', 'exec'), _ns)
A, W, TH, JC, n, sym = _ns['A'], _ns['W'], _ns['TH'], _ns['JC'], _ns['n'], _ns['sym']
Ts, U = schur(A, output='real', sort=lambda a: abs(a) < 1.0)[:2]
ln2 = np.log(2.0)
TANAKA_FLOOR = 1.1685

ev = np.linalg.eigvals(A)
print('== (0) 锚点 A 的谱 ==')
print(' lambda(A) =', np.sort_complex(ev))
unst = np.sort(np.abs(ev[np.abs(ev) >= 1.0]))
R_exp = float(np.sum(np.log2(unst)))
detA = float(np.log2(abs(np.linalg.det(A))))
print(f' n_u = {len(unst)}  |lambda_u| = {unst}')
print(f' [P1] R_exp = Σ log2|λ_u| = {R_exp:.6f} bit/sample     ← 2606.31396 Th.1 的必要条件')
print(f'      Tanaka [4] 无约束水平渐近线（文献印出） = {TANAKA_FLOOR} bit')
print(f'      log2|det A|（含稳定模态，KH Th.5 退化右端） = {detA:.6f} bit')
print(f'      比值 {TANAKA_FLOOR/max(R_exp,1e-12):.3f}×，差 {TANAKA_FLOOR-R_exp:+.4f} bit')
print(f'      ⇒ R_exp {"<" if R_exp < TANAKA_FLOOR else ">="} 1.1685：'
      f'{"[4] 的精确前沿更强，Th.1 撞不掉 §65-A 的证书腿" if R_exp < TANAKA_FLOOR else "Th.1 反超，须自首"}')

print('\n== (1) 逐不稳定模态：谁负责多少 bit ==')
for lam in unst:
    print(f'   |λ|={lam:.6f}  → log2|λ| = {np.log2(lam):.6f} bit')
print(f'   合计 {R_exp:.6f}   （对照：整植物 log2|det A| = {detA:.6f}，稳定模态把它拉低）')


def eval_design(T, Z):
    """与 e96c `design()` 同一条特征分解通路，但**不**用 ρ 门筛掉设计——本轮要问的正是
    "哪些设计落在 2606.31396 的可检测性前提之外"。"""
    Sg = Z @ sym(T) @ Z.T
    wv, Vv = np.linalg.eigh(sym(Sg))
    idx = np.where(wv > 1e-8 * max(wv.max(), 1e-300))[0]
    if len(idx) == 0:
        return np.nan, np.inf, np.inf, 0, np.inf
    C = np.sqrt(wv[idx])[:, None] * Vv[:, idx].T
    r = len(idx)
    try:
        Pm = sym(solve_discrete_are(A.T, C.T, W, np.eye(r)))
    except Exception:
        return np.nan, np.nan, np.nan, r, np.nan
    Sm = C @ Pm @ C.T + np.eye(r)
    Lk = Pm @ C.T @ np.linalg.inv(Sm)
    Pp = sym(Pm - Lk @ Sm @ Lk.T)
    rho = float(np.abs(np.linalg.eigvals((np.eye(n) - Lk @ C) @ A)).max())
    D = JC + float(np.trace(TH @ Pp))
    I_lane = 0.5 * (np.linalg.slogdet(Pm)[1] - np.linalg.slogdet(Pp)[1]) / ln2
    I_cor = float(np.linalg.slogdet(np.eye(r) + C @ Pm @ C.T)[1]) / (2 * ln2)
    return rho, D, I_lane, r, I_cor


print('\n== (2) [P3/P4] 候选平面 × 信息强度：ρ、可达代价、两条率公式对表 ==')
print(f" {'候选':<28}{'ρ((I−LC)A)':>11}  可检测  Th.1在前提内     D         I_车道      I_Cor.1   差")
for k, nm in [(4, 'tail-4 全信息'), (3, 'tail-3 平面'), (2, 'tail-2 平面'), (1, 'tail-1 平面')]:
    Z = U[:, n - k:]
    for mult, lab in [(20.0, '×20 强'), (1.0, '×1'), (0.01, '×0.01 弱')]:
        rho, D, I, r, Icor = eval_design(mult * np.eye(k), Z)
        ok = rho < 1.0
        Ds = f'{D:9.4f}' if ok else '     无效'
        Is = f'{I:10.5f}' if ok else '      无效'
        Cs = f'{Icor:10.5f}' if ok else '      无效'
        print(f' {(nm+" "+lab):<28}{rho:11.5f}   {"是" if ok else "否":<5}  '
              f'{"是" if ok else "否":<12} {Ds} {Is} {Cs} '
              f'{(abs(I-Icor) if ok else np.nan):9.2e}  r={r}')
print('\n（ρ≥1 的行：DARE 给的是反稳定解，D/I 一律标"无效"——不是数值失败，是设计本身不可实现。）')

print('\n== (2b) 每个平面的"率地板"：尺度 → 0 时的 I，与 R_exp 对表 ==')
for k, nm in [(4, 'tail-4 全信息'), (3, 'tail-3 平面'), (2, 'tail-2 平面'), (1, 'tail-1 平面')]:
    Z = U[:, n - k:]
    vals = []
    for mult in (1e-4, 1e-6, 1e-8):
        rho, D, I, r, Icor = eval_design(mult * np.eye(k), Z)
        vals.append((mult, rho, I))
    rho_l, I_l = vals[-1][1], vals[-1][2]
    if np.isfinite(rho_l) and rho_l < 1.0:
        print(f' {nm:<14} ρ={rho_l:.5f}  I(1e-4)={vals[0][2]:.6f}  I(1e-6)={vals[1][2]:.6f}  '
              f'I(1e-8)={I_l:.6f}   vs R_exp={R_exp:.6f}  超出 {I_l-R_exp:+.6f} bit')
    else:
        print(f' {nm:<14} ρ={rho_l}  不可检测 ⇒ 代价对任何率不可达'
              f'（这面墙由可检测性给出＝L-CSS Thm 2.1，不是本轮新内容）')

print('\n== (2c) 哪个平面漏看哪根不稳定模态 ==')
w_, Vv_ = np.linalg.eig(A)
for k, nm in [(4, 'tail-4'), (3, 'tail-3'), (2, 'tail-2'), (1, 'tail-1')]:
    Z = U[:, n - k:]
    parts = []
    for j, lam in enumerate(w_):
        v = Vv_[:, j]
        seen = np.linalg.norm(Z.T @ np.real(v)) + np.linalg.norm(Z.T @ np.imag(v))
        if abs(lam) >= 1.0:
            parts.append(f'λ={lam.real:+.4f} |λ|={abs(lam):.4f} 可见度 {seen:.4f}'
                         f'{"  ← 看不见" if seen < 1e-6 else ""}')
    print(f' {nm:<8} ' + ' | '.join(parts))

print('\n ⇒ ρ≥1 的行：过滤不收敛 ⇒ 2606.31396 的 Th.1 前提不满足，他们的定理对那一行**不作任何陈述**；')
print('   而我们那面墙要说的正是这些行（代价对任何率都不可达）。这就是作用域差别，不是"类似"。')
print('\n== (3) 已知构造点：各自**本身**的率 vs R_exp（注意口径，别把反推率当设计率）==')
for lab, Dv, I_own in [('无约束构造 42.0427（E71/E74）', 42.0427, 3.0000),
                       ('tail-3 平面 Powell 45.4537', 45.4537, 3.0000),
                       ('tail-2 平面 Powell 50.6023', 50.6023, 3.0000)]:
    print(f'   {lab:<34} 该设计自身的率 {I_own:.4f} bit  vs R_exp={R_exp:.4f} bit → '
          f'{"满足 Th.1（必要条件下界对这些设计不构成额外约束）" if I_own >= R_exp else "违反 Th.1！须查"}')
print('   注：§65-A 的 2.6678/2.3261 bit 是**无约束前沿**在这些代价上反推的最小率，')
print('       不是这些平面设计自身花的率（它们都被构造成 I=3.0000）。两个口径不许混。')
