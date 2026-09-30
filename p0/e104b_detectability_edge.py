r"""E104b = 补 e104 的 [X2] 那一行。e104 的"漏 λ1"构造**失败且原因是我的 bug**：
我用 ev=eigvals(A) 的下标去取 w_,Vv_=eig(A) 的特征向量，两个 LAPACK 驱动的特征序不保证一致
⇒ 投影掉的可能根本不是 λ1 的方向（实测可见度 0.953，设计照样可检测、S 照样等于 |det A_u|²）。
本轮改用 eig() 自己返回的 (w_, Vv_) 成对索引，并把"确实漏看"作为**前置断言**（可见度 < 1e-8 才算构造成功）。

 [X2a] 断言先过：构造出的 Z 对某根不稳定模态可见度 <1e-8，且 ρ((I−LC)A) ≥ 1（不可检测）。
 [X2b] 判据（带等号分支）：该设计的 S(0) = 只含可见不稳定模态的 |λ|² 乘积，相对偏离 ≤1e-4。
       若 S(0) 反而等于全乘积 |det A_u|² ⇒ [X2] 作废，地板与可见性无关（须当场自首）。
"""
import numpy as np
from scipy.linalg import solve_discrete_are

OUT = open('p0/e104b_out.txt', 'w', encoding='utf-8')


def p(*a):
    OUT.write(' '.join(str(x) for x in a) + '\n')
    OUT.flush()


_src = open('p0/exp_c_audit.py', encoding='utf-8').read().split("print(r'== E53")[0]
_ns = {'__name__': 'p'}
exec(compile(_src, 'p0/exp_c_audit.py[preamble]', 'exec'), _ns)
A, W, TH, JC, n = _ns['A'], _ns['W'], _ns['TH'], _ns['JC'], _ns['n']
sym = _ns['sym']
ln2 = np.log(2.0)
w_, Vv = np.linalg.eig(A)          # 成对，不再与 eigvals 混用
u_idx = np.where(np.abs(w_) >= 1.0)[0]
p('== (0) 成对谱 ==')
for j in range(n):
    p(f'   idx{j}  λ={w_[j]:+.6f}  |λ|={abs(w_[j]):.6f}  {"不稳定" if abs(w_[j]) >= 1 else "稳定"}')
DET_U2 = float(np.prod(np.abs(w_[u_idx]) ** 2))
p(f' |det A_u|² = {DET_U2:.6f}   R_exp = {np.sum(np.log2(np.abs(w_[u_idx]))):.6f} bit')


def S_of(Z, s):
    r = Z.shape[1]
    C = np.sqrt(s) * Z.T
    Pm = sym(solve_discrete_are(A.T, C.T, W, np.eye(r)))
    Sm = C @ Pm @ C.T + np.eye(r)
    Lk = Pm @ C.T @ np.linalg.inv(Sm)
    rho = float(np.abs(np.linalg.eigvals((np.eye(n) - Lk @ C) @ A)).max())
    Sv = float(np.linalg.det(np.eye(r) + s * (Z.T @ Pm @ Z)))
    return rho, Sv


def vis(Z):
    out = []
    for j in u_idx:
        vr, vi = np.real(Vv[:, j]), np.imag(Vv[:, j])
        out.append(max(np.linalg.norm(Z.T @ vr) / max(np.linalg.norm(vr), 1e-300),
                       np.linalg.norm(Z.T @ vi) / max(np.linalg.norm(vi), 1e-300)))
    return np.array(out)


p('\n== (1) 逐根不稳定模态构造"漏看它"的子空间 ==')
for drop in u_idx:
    v = np.real(Vv[:, drop])
    if np.linalg.norm(v) < 1e-12:
        v = np.imag(Vv[:, drop])
    v = v / np.linalg.norm(v)
    P = np.eye(n) - np.outer(v, v)
    Q, Rr = np.linalg.qr(P)
    for rk in (1, 2, 3):
        Z = Q[:, :rk]
        vv = vis(Z)
        seen = bool(np.all(vv > 1e-8))
        line = f' 丢掉 λ={w_[drop]:+.5f}(|λ|={abs(w_[drop]):.5f}) 的 rank-{rk}: 可见度={np.round(vv,6)}'
        if not seen:
            try:
                rho, Sv = S_of(Z, 1e-7)
            except Exception as e:
                p(line + f'  DARE 失败 {type(e).__name__}')
                continue
            keep = [j for j in u_idx if j != drop]
            tgt = float(np.prod(np.abs(w_[keep]) ** 2)) if keep else 1.0
            rel = Sv / tgt - 1
            p(line + f' ⇒ [X2a] 通过（确实漏看），ρ={rho:.5f} '
              f'{"不可检测" if rho >= 1 else "居然可检测(须查)"}；S(1e-7)={Sv:.6f} vs 可见乘积 {tgt:.6f} '
              f'相对 {rel*100:+.3f}% → {"[X2b] 成立：地板只数可见的不稳定模态" if abs(rel) <= 1e-2 else "[X2b] 不成立"}')
        else:
            p(line + ' ⇒ [X2a] 断言失败：该 Z 仍然看见所有不稳定模态（正交补的前 rk 列没躲开它），不作 [X2b] 判定')

p('\n== (2) 对照：一个可检测 rank-1 的 S 序列（确认 [X1] 与 bug 无关） ==')
rng = np.random.default_rng(7)
for c in range(3):
    Z, _ = np.linalg.qr(rng.standard_normal((n, 1)))
    seq = []
    for s in 10.0 ** -np.arange(3.0, 7.01, 1.0):
        try:
            seq.append(S_of(Z, s)[1])
        except Exception:
            seq.append(np.nan)
    p(f'  r1-{c} vis={np.round(vis(Z),4)} S=' + ' '.join(f'{v:.6f}' for v in seq)
      + f'  vs |det A_u|²={DET_U2:.6f}（末点相对 {seq[-1]/DET_U2*100-100:+.4f}%）')
OUT.close()
print('written p0/e104b_out.txt')
