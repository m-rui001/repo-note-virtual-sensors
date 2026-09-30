r"""E104 = e103 暴露出来的真正结构：地板不是"逐模态"的，而是**标度 SNR 行列式的乘积不变量**。

e103 实测：rank-1 子空间的标度预测方差是单个 θ*=4.05（不是 {0.7232, 1.9322} 里任何一个），
而 ½log2(1+4.0528)=1.16855=R_exp；rank-2/3 的 θ 集合各不相同（r2-00=[0.008,4.015]、
r2-01=[0.243,3.067]、tail-2=[0.710,1.954]），但 Σ½log2(1+θ_i) 全部落到 1.1685。
⇒ 猜测（本轮要判的）：
 [X1] 对任何**可检测**的 V：det(I_r + s·ZᵗPm Z) → |det A_u|²（s→0），与 V 无关；
      各通道怎么分这个乘积由 V 决定 ⇒ 地板是**行列式级**不变量，不是逐通道不变量。
 [X2] 对**不可检测**的 V（漏看某些不稳定模态）：极限 = 只含"可见"不稳定模态的乘积，
      严格小于 |det A_u|² ⇒ 2606.31396 Th.1 的地板在那里**不成立**（前提失效），
      而 1510.04214 图 5 的水平渐近线 1.169 只数了可见的那部分。
 [X3] 判据带等号分支（§66-E-2）：相对偏离 |Π/|detA_u|² − 1| ≤ 1e-4 记"相等"，
      ≤1e-3 记"近平等（受 O(s) 修正支配，用外推裁决）"，>1e-3 记"不等，猜测作废"。
外推：S(s) 的一阶修正是 O(s)（e101b 的 p≈1 就是这件事），所以用最小三个 s 做 S(0)=S(s)−s·S'(s)
线性外推；同时**报原始值**，不许只报外推。
"""
import numpy as np
from scipy.linalg import solve_discrete_are, schur

np.set_printoptions(precision=6, suppress=True, linewidth=170)
OUT = open('p0/e104_out.txt', 'w', encoding='utf-8')


def p(*a):
    OUT.write(' '.join(str(x) for x in a) + '\n')
    OUT.flush()


_src = open('p0/exp_c_audit.py', encoding='utf-8').read().split("print(r'== E53")[0]
_ns = {'__name__': 'p'}
exec(compile(_src, 'p0/exp_c_audit.py[preamble]', 'exec'), _ns)
A, W, TH, JC, n = _ns['A'], _ns['W'], _ns['TH'], _ns['JC'], _ns['n']
sym = _ns['sym']
ln2 = np.log(2.0)
ev = np.linalg.eigvals(A)
u_idx = np.where(np.abs(ev) >= 1.0)[0]
DET_U2 = float(np.prod(np.abs(ev[u_idx]) ** 2))
R_exp = float(np.sum(np.log2(np.abs(ev[np.abs(ev) >= 1.0]))))
w_, Vv_ = np.linalg.eig(A)
p('== (0) 目标 ==')
p(f' 不稳定模态 {np.round(np.abs(ev[u_idx]),5)} → |det A_u|² = {DET_U2:.6f}，'
  f'  ½log2(该值) = {0.5*np.log2(DET_U2):.6f} = R_exp = {R_exp:.6f}')


def S_of(Z, s):
    r = Z.shape[1]
    C = np.sqrt(s) * Z.T
    Pm = sym(solve_discrete_are(A.T, C.T, W, np.eye(r)))
    rho_ok = True
    Sm = C @ Pm @ C.T + np.eye(r)
    Lk = Pm @ C.T @ np.linalg.inv(Sm)
    rho = float(np.abs(np.linalg.eigvals((np.eye(n) - Lk @ C) @ A)).max())
    Sv = float(np.linalg.det(np.eye(r) + s * (Z.T @ Pm @ Z)))
    th = np.sort(np.linalg.eigvalsh(sym(Z.T @ (s * Pm) @ Z)))
    return rho, Sv, th, r


def visible(Z):
    out = []
    for j in u_idx:
        vr = np.real(Vv_[:, j])
        out.append(float(np.linalg.norm(Z.T @ vr) / max(np.linalg.norm(vr), 1e-300)))
    return np.array(out)


GR = 10.0 ** -np.arange(2.0, 7.01, 1.0)   # 1e-2 … 1e-7
_, U = schur(A, output='real', sort=lambda a: abs(a) < 1.0)[:2]
rng = np.random.default_rng(20260932)
designs = [(U[:, n - k:], f'tail-{k}') for k in (4, 3, 2, 1)]
for rk in (1, 2, 3):
    for c in range(3):
        Z, _ = np.linalg.qr(rng.standard_normal((n, rk)))
        designs.append((Z, f'r{rk}-{c}'))
# 一个"只漏一根不稳定模态"的构造 rank-2：与 v1 正交、含 v2
v1, v2 = np.real(Vv_[u_idx[0]]), np.real(Vv_[u_idx[1]])
v1 = v1 / np.linalg.norm(v1)
P1 = np.eye(n) - np.outer(v1, v1)
Q1, _ = np.linalg.qr(P1)
Zx = Q1[:, :2]
designs.append((Zx, '漏λ1(rank2)'))

p('\n== (1) [X1/X2/X3] 逐设计：S(s)=det(I+s ZᵗPmZ) 的原始序列与外推 ==')
p(f" {'设计':<13}{'r':>2} {'ρ(1e-7)':>9} 可见度 " + ''.join(f'{s:>11.0e}' for s in GR)
  + f"  {'S(0)外推':>12}{'相对偏离':>11}  判定")
for Z, nm in designs:
    rows = []
    for s in GR:
        try:
            rows.append((s,) + S_of(Z, s)[:2])
        except Exception:
            rows.append((s, np.nan, np.nan))
    sv = [r[2] for r in rows]
    rho_end = rows[-1][1]
    last3 = [(r[0], r[2]) for r in rows[-3:] if np.isfinite(r[2])]
    ext = np.nan
    if len(last3) == 3:
        x = np.array([a for a, b in last3]); y = np.array([b for a, b in last3])
        ext = float(np.polyfit(x, y, 1)[1])
    rel = (ext / DET_U2 - 1) if np.isfinite(ext) else np.nan
    visv = np.round(visible(Z), 4)
    det_ok = np.all(visible(Z) > 1e-6)
    lab = ('相等' if abs(rel) <= 1e-4 else ('近平等' if abs(rel) <= 1e-3 else '不等')) if np.isfinite(rel) else '无值'
    if not det_ok:
        lab += '（[X2]：不可检测，本来就该小）'
    p(f' {nm:<13}{Z.shape[1]:>2} {rho_end:>9.5f} {str(visv):<20}'
      + ''.join(f'{v:>11.5f}' if np.isfinite(v) else f'{"nan":>11}' for v in sv)
      + f'  {ext:>12.5f}{rel*100:>10.3f}%  {lab}')

p('\n== (2) 各设计的 θ 集合（同一个乘积、不同的分配） ==')
for Z, nm in designs:
    if not np.all(visible(Z) > 1e-6):
        continue
    rho, Sv, th, r = S_of(Z, 1e-6)
    prod = float(np.prod(1 + np.maximum(th, 0)))
    p(f' {nm:<12} θ(1e-6)={np.round(th,4)}  Π(1+θ)={prod:.5f}  vs |det A_u|²={DET_U2:.5f} '
      f'（相对 {(prod/DET_U2-1)*100:+.3f}%）  各通道率 {np.round(0.5*np.log2(1+np.maximum(th,0)),4)}')

p('\n== (3) [X2] 专门那一行：漏看 λ1 的 rank-2 子空间，它的极限应当是 |λ2|² 而非 |det A_u|² ==')
rho, Sv, th, r = S_of(Zx, 1e-7)
lam2 = float(abs(ev[u_idx[1]]) ** 2)
p(f' 可见度 = {np.round(visible(Zx), 6)}  ρ={rho:.5f}（应 ≥1：不可检测）')
p(f' S(1e-7)={Sv:.6f}   |λ2|²={lam2:.6f}   |det A_u|²={DET_U2:.6f}')
p(f' ⇒ 若 S 贴住 |λ2|²，则"地板=可见不稳定模态的乘积"，不可见的那根**不进地板**；')
p(f'   这正是 2606.31396 Th.1 的前提（均方可观测）失效的那一侧，他们的 R_exp 在那里不约束任何设计。')
OUT.close()
print('written p0/e104_out.txt')
