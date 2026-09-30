r"""E59 = 对 C 的 §10（R28-C6，Grassmannian 下降）做独立复核。

他交付了 `.work3/c24_bestF_r{1,2,3}.npy` 并点名让我"接自己的 SCA/定点管线交叉核对"（§10.5）。
我这边已有的路线是 E53：多起点 + 流形上 Powell 精修，标量定点迭代 + 奇异 DARE 双路算 $\Phi$。
这一跑四件事：
  [1] 载入他的三个 $F$，用我的两条独立管线算 $\Phi$，与他的 10.169224 / 0.735143 / 0.158701 对表；
  [2] 复算他 §10.3 整张候选规则表（$\Theta$ 前 $r$、标准基、$W$ 前 $r$、$A$ 的不稳模、不稳 Schur 尾）；
  [3] 复算他 §10.4 的主角（最优子空间 vs $\Theta$ 前 $r$ 特征子空间）；
  [4] **不复用他的代码**，我自己写 $r=1$ 的球面投影梯度下降，24 个 Haar 起点，检验"24/24 命中同一个值"。
"""
import sys
sys.stdout.reconfigure(encoding='utf-8')
import numpy as np

_src = open('p0/exp_c_audit.py', encoding='utf-8').read().split("print(r'== E53")[0]
_ns = {'__name__': 'e53_preamble'}
exec(compile(_src, 'p0/exp_c_audit.py[preamble]', 'exec'), _ns)
A, W, TH, JC, n, sym = _ns['A'], _ns['W'], _ns['TH'], _ns['JC'], _ns['n'], _ns['sym']
phi_iter, phi_dare = _ns['phi_iter'], _ns['phi_dare']
print(r'== E59：C 的 §10 交付物独立复核 ==')
print(r'  $j_c=\operatorname{tr}(WP_c)=%.4f$   百分数口径：$\Phi/j_c$（+46.5%% $\Leftrightarrow$ $\Phi=14.6398$）' % JC)


def phi_both(F):
    """两条独立管线：定点迭代 vs 奇异 DARE（都作用在 $F$ 的**行空间**上，$F$ 为 $(r,n)$）。
    返回 (定点值, DARE值, 相对分歧)。"""
    F = np.atleast_2d(F)
    v1 = phi_iter(F)[0]
    v2 = phi_dare(F)[0]
    return v1, v2, abs(v1 - v2) / max(1e-300, abs(v1))


def pangle(Fa, Fb):
    """子空间主角（度）。Fa,(n,ra) Fa,(n,rb) 正交标架。"""
    sv = np.linalg.svd(Fa.T @ Fb, compute_uv=False)
    return np.degrees(np.arccos(np.clip(sv, -1, 1)))


# ---------------- [1] 他的三个 F ----------------
print('\n[1] 载入 `.work3/c24_bestF_r{1,2,3}.npy`，用我的两条管线算 $\Phi$')
C_VALS = {1: 10.169224, 2: 0.735143, 3: 0.158701}
E53_VALS = {1: 10.169224, 2: 0.735143, 3: 0.158701}
print('  %-3s %-12s %-12s %-12s %-12s %-12s %-10s' %
      ('r', '形状', '正交性', '我的定点', '我的DARE', 'C 公布', '相对差'))
BEST = {}
for r in (1, 2, 3):
    F = np.load('.work3/c24_bestF_r%d.npy' % r)
    if F.shape[0] != n:                      # 他可能存成 (r,n)
        F = F.T
    BEST[r] = np.ascontiguousarray(F)        # (n,r) 标架；管线要 (r,n)
    ortho = np.abs(F.T @ F - np.eye(r)).max()
    v1, v2, dv = phi_both(F.T)
    C_VAL = C_VALS[r]
    print('  %-3d %-12s %-12.2e %-12.6f %-12.6f %-12.6f %-10.2e  两路分歧 %.1e' %
          (r, str(F.shape), ortho, v1, v2, C_VAL, abs(v1 - C_VAL) / C_VAL, dv))

# ---------------- [2] 候选规则表 ----------------
print('\n[2] 复算 C 的 §10.3 候选规则表（每格 = 我的定点 $\Phi$ / 占 $j_c$ / 与他值之差）')
wT, vT = np.linalg.eigh(TH)
order = np.argsort(-wT)
wT, vT = wT[order], vT[:, order]
wA, vA = np.linalg.eig(A)
UNST = [k for k in range(n) if abs(wA[k]) >= 1.0]
VU = np.array([np.real(vA[:, k]) for k in UNST]).T
VU = VU / np.linalg.norm(VU, axis=0, keepdims=True)
wW, vW = np.linalg.eigh(W)
ordW = np.argsort(-wW)
wW, vW = wW[ordW], vW[:, ordW]
STAB = [k for k in range(n) if abs(wA[k]) < 1.0]
VS = np.array([np.real(vA[:, k]) for k in STAB]).T          # 稳定不变子空间 (n, 2)
VS = np.linalg.qr(VS)[0]


def frame(B, r):
    """任意列组 $\to$ 前 $r$ 个正交方向；不足则用与已有列正交的随机方向补齐。"""
    B = np.asarray(B, float)
    Q, R = np.linalg.qr(B)
    if Q.shape[1] < r:
        rngc = np.random.default_rng(1000 + r)
        while Q.shape[1] < r:
            v = rngc.standard_normal(n)
            v = v - Q @ (Q.T @ v)
            nv = np.linalg.norm(v)
            if nv > 1e-6:
                Q = np.hstack([Q, (v / nv)[:, None]])
    return Q[:, :r]


def phi_of_basis(Qn_r):
    return phi_both(Qn_r.T)


C_T = {(1, 'theta'): 12.2224, (2, 'theta'): 0.7918, (3, 'theta'): 0.1903,
       (1, 'std'): 30.9810, (2, 'std'): 18.6282, (3, 'std'): 7.8631,
       (1, 'W'): 520.3206, (2, 'W'): 8.3182, (3, 'W'): 1.0567,
       (1, 'A'): 179.5693, (2, 'A'): 78.6010,
       (1, 'desc'): 10.169224, (2, 'desc'): 0.735143, (3, 'desc'): 0.158701,
       (2, 'schur'): 14.6398, (3, 'schur'): 14.6398}
LBL = {'theta': r'$\Theta$ 前 r', 'std': '标准基', 'W': '$W$ 前 r',
       'A': r'$A$ 不稳模(2维)', 'stab': r'$A$ 稳定模(2维)',
       'desc': r'下降最优', 'schur': r'不稳 Schur 尾(原文)'}
for r in (1, 2, 3):
    cands = {
        'desc': BEST[r],
        'theta': vT[:, :r],
        'std': np.eye(n)[:, :r],
        'W': vW[:, :r],
        'A': VU[:, :min(r, 2)],
        'stab': VS[:, :min(r, 2)],
    }
    if r >= 2:
        cands['schur'] = VU[:, :2]                 # 实 Schur：不稳块=span(两条实本征向量)
    for key, B in cands.items():
        Q = frame(B, r)
        v1, v2, dv = phi_of_basis(Q)
        ref = C_T.get((r, key))
        tag = '' if ref is None else ('   C 表 %+8.4f → 差 %+.2e' % (ref, (v1 - ref) / abs(ref)))
        print('  r=%d  %-22s  $\\Phi$=%12.4f  占 $j_c$=%+8.2f%%  两路分歧 %.1e%s'
              % (r, LBL[key], v1, 100 * v1 / JC, dv, tag))

# ---------------- [3] 主角 ----------------
print('\n[3] 他 §10.4 的主角（最优子空间 vs $\Theta$ 前 $r$ 特征子空间，单位：度）')
for r in (1, 2, 3):
    pa = pangle(BEST[r], vT[:, :r])
    print('  r=%d  %s   （C 写：r=2 含首向量 0.06 度，r=3 含前两 0.00/0.00，r=1 距首 10.59 度）'
          % (r, np.array2string(pa, precision=2)))
    print('      与 $A$ 不稳模子空间的夹角 %s 度；到最近坏超平面距离 d=%.4f'
          % (np.array2string(pangle(BEST[r], VU), precision=2),
             np.abs(BEST[r].T @ VU).min()))

# ---------------- [4] 我自己的 r=1 下降 ----------------
print('\n[4] 我自己的 $r=1$ 下降（球面投影梯度 + 回溯线搜索，24 个 Haar 起点；不用他的代码）')
rng = np.random.default_rng(771)


def phi_f(f):
    f = f / np.linalg.norm(f)
    F = f.reshape(1, n)
    p = phi_iter(F)
    return p[0] if isinstance(p, tuple) else p


def grad_f(f, h=1e-5):
    f = f / np.linalg.norm(f)
    Q = np.eye(n) - np.outer(f, f)
    U, s, _ = np.linalg.svd(Q)
    basis = U[:, s > 0.5][:, :n - 1]            # 切空间正交基 (n, n-1)
    g = np.zeros(n - 1)
    f0 = phi_f(f)
    for i in range(n - 1):
        d = basis[:, i] * h
        g[i] = (phi_f(f + d) - phi_f(f - d)) / (2 * h)
    return f0, basis @ g, np.linalg.norm(g)


NSTART = 24
vals, its = [], []
for k in range(NSTART):
    f = rng.standard_normal(n)
    f /= np.linalg.norm(f)
    L = 1.0
    it = 0
    for it in range(400):
        val, g, gn = grad_f(f)
        if gn < 1e-7:
            break
        step = L / (1.0 + gn)
        improved = False
        for _ in range(30):                      # 回溯
            fn = f - step * g
            fn /= np.linalg.norm(fn)
            vn = phi_f(fn)
            if vn < val - 1e-12:
                f, L = fn, min(step * 3.0, 1e4)
                improved = True
                break
            step *= 0.5
        if not improved:
            L *= 0.25
            if L < 1e-14:
                break
    vals.append(phi_f(f))
    its.append(it)
vals = np.array(vals)
print('  24 起点收敛值：最小 %.6f  最大 %.6f  唯一值 %d 个' % (vals.min(), vals.max(), len(np.unique(np.round(vals, 6)))))
print('  迭代数：min %d  max %d  median %.0f' % (min(its), max(its), np.median(its)))
print('  与 C 的 10.169224 的相对差：min %.2e  max %.2e' %
      (np.abs(vals - 10.169224).min() / 10.169224, np.abs(vals - 10.169224).max() / 10.169224))
print('  命中 %.2f%%/起点（若全部等于最优则他的 24/24 在我这条独立路线上复现）' % (100.0 * (np.abs(vals - 10.169224) < 1e-4).mean()))
print('\n  判读：[1][2][3] 是核对，[4] 是复现。若 [4] 出现非最优收敛值，则"任意随机起点 100% 命中"需要加限定词。')
