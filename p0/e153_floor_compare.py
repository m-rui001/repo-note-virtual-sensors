# -*- coding: utf-8 -*-
r"""e153（结 §98-B 挂的账："核完 $D^{min}$ 的构成式之前不许写'同一地板'"）。

1912.07640 附录 C（抽取文本第 2880–2900 行）给的是
 $$D^{min}_{[0,\\infty]}=\\mathrm{trace}(\\bar\\Sigma),\\qquad \\bar\\Sigma=\\lim_{t\\to\\infty}\\Sigma_t,$$
即**无权**迹、且是**给定数据的常数**（$\\Pi^\\xi_{t|t-1}=A\\Pi^\\xi_{t-1|t-1}A^{\\mathsf T}+\\bar\\Sigma_t$ 里的噪声序列极限）。
本站的 $D_{\\rm floor}$ 是 $\\lim_{s\\to\\infty}\\big(J_C+\\mathrm{tr}(\\Theta P_p)\\big)$ $-$ $-$ **带任务权 $\\Theta$、且随设计 $Z$ 变**。
本实验把两者在同一批植物上**并排印出来**，判"是不是同一个对象"。

 [F1] 每株印：$\\mathrm{tr}(W)$（他们形状的无权地板类比）、$J_C$、以及我侧 $s=10^7$ 端的
      $\\Theta$-加权地板 $\\min_Z\\mathrm{tr}(\\Theta P_p)$（$r=1..n$，$Z$ 取 $\\Theta$ 轴与随机正交阵各若干）。
 [F2] 同一档再印**无权**版本 $\\min_Z\\mathrm{tr}(W^{1/2\\cdot}P_p)$ 一类的对照 $\\mathrm{tr}(P_p)$，
      用来量"加权"这一项本身贡献多少 $-$ $-$ 若加权与无权同序，则差别在常数；若不同序，差别在泛函。
 [F3] 预注册（跑前写死，两问分开）：
      (a) 按 $r$ 的地板序列在两口径下是否同序 $-$ $-$ 这条**弱**（$r$ 单调几乎必然同序），只作体检；
      (b) 真正的判据：**同一 $r$ 下**加权 argmin 子空间与无权 argmin 子空间是否同一个。
          主角度 $\\ge30°$（只取 $r<n$ 的档 $-$ $-$ $r=n$ 时两侧都是全空间、角度恒 $0$，属 §99-B 我刚刚点名过的真空判据）
          $\\Rightarrow$ 两个"地板"是**不同泛函**，"同一地板"判否；若全部 $<5°$ 才可写"同族不同常数"。

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
    return 0.0 if X is None else 0.5 * (X + X.T)


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


def post_at(Z, A, W, s):
    r = Z.shape[1]
    Ir = np.eye(r)
    C = np.sqrt(s) * Z.T
    Pm = sym(solve_discrete_are(A.T, C.T, W, Ir))
    Sm = C @ Pm @ C.T + Ir
    Lk = Pm @ C.T @ np.linalg.inv(Sm)
    return sym(Pm - Lk @ Sm @ Lk.T)


def floors(A, W, TH, JC, m, s=1e7, nrand=400):
    """返回按 r 的 (加权地板, 无权地板)：min over Z（Theta 轴补集 + 随机正交阵）。"""
    wth, Uth = np.linalg.eigh(TH)
    order = np.argsort(-wth)
    Ww, Wu = [], []
    Zw, Zu = [], []
    for r in range(1, m + 1):
        cands = [Uth[:, order[:r]]]
        rng = np.random.default_rng(1000 + r)
        for _ in range(nrand):
            Q, Rr = np.linalg.qr(rng.normal(0, 1, (m, r)))
            cands.append(Q * np.sign(np.diag(Rr)))
        bw = bu = np.inf
        zb = zt = cands[0]
        for Z in cands:
            try:
                Pp = post_at(Z, A, W, s)
            except Exception:
                continue
            v = float(np.trace(TH @ Pp))
            if v < bw:
                bw, zb = v, Z
            v2 = float(np.trace(Pp))
            if v2 < bu:
                bu, zt = v2, Z
        Ww.append(bw); Wu.append(bu); Zw.append(zb); Zu.append(zt)
    return np.array(Ww), np.array(Wu), np.array(wth[::-1]), Zw, Zu


def subspace_angle(Za, Zb):
    sv = np.clip(np.linalg.svd(Za.T @ Zb, compute_uv=False), -1.0, 1.0)
    return np.degrees(np.arccos(sv))



print('== [F1]/[F2] 加权 vs 无权地板，以及他们形状的常数 trace(W) ==')
ORDER_DIFF = []
for nm, spec in [('anchor', ('anchor', 0)), ('rand-1', ('rand4', 1)), ('rand-2', ('rand4', 2)),
                 ('rand-3', ('rand4', 3)), ('rand-4', ('rand4', 4)), ('big-6', ('big6', 11))]:
    A, B, W, TH, JC, m = make_plant(*spec)
    fw, fu, ev, Zw, Zu = floors(A, W, TH, JC, m)
    so_w = tuple(np.argsort(-fw)); so_u = tuple(np.argsort(-fu))
    same = so_w == so_u
    ORDER_DIFF.append(same)
    ang = [subspace_angle(Zw[i], Zu[i]) for i in range(m)]
    print('%-8s n=%d tr(W)=%9.4f JC=%9.4f | 加权地板(r=1..n)=%s' % (nm, m, float(np.trace(W)), JC,
                                                                    np.array2string(fw, precision=3)))
    print('%-8s          无权地板(r=1..n)=%s | 加权/无权逐档比=%s | 两口径按 r 的排序相同？%s'
          % ('', np.array2string(fu, precision=3), np.array2string(fw / np.maximum(fu, 1e-30), precision=2), same))
    print('%-8s          同一 r 下"加权 argmin 子空间 vs 无权 argmin 子空间"主角度(max，只取 r<n)：%s deg'
          % ('', np.array2string(np.array([a.max() for a in ang[:m - 1]]), precision=1)))
    print('%-8s          tr(W)/加权地板(r=n)=%.3f 倍；tr(W)/加权地板(r=1)=%.3f 倍'
          % ('', float(np.trace(W)) / fw[-1], float(np.trace(W)) / fw[0]))

print('\n== [F3] 结算 ==')
print('六株两口径按 r 的排序全同？%s；不同序的株数 %d' % (all(ORDER_DIFF), sum(1 for x in ORDER_DIFF if not x)))
print('判读：真正的"同一泛函"检验是**同一 r 下两侧 argmin 子空间是否同一个**（上一节每株的 max 主角度）。'
      '角度大 $\\Rightarrow$ 加权与无权选的是不同设计 $\\Rightarrow$ 两个"地板"不是同一泛函在不同坐标下的读数。')

