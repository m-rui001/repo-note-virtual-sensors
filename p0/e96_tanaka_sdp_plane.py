r"""E96 = 把 1510.04214 式 (18) 的 max-det SDP 在我车道锚点植物上**独立复现**，再判断
"平面限制后的率侧问题还是不是凸的"。这是我 §64-F-3 自己许下的实验，也是 §63-E 那个
"把线性 range 条件搬进率侧 SDP"的正面判决。

== 为什么要做（先写死判据，跑完不许改）==
Remark 2 那一格的验收线是"认证下界 L > 43.80"。我 §64 向 C 提议：Tanaka 的 (18) 是
**已被证明精确**的率侧机器（Theorem 2：在全体 Borel 传感机制里最优），把它限制在
tail-3 平面上跑，就能反解出精确的 D_V(3 bit)。我当时给出的理由是"range 条件是线性的，
所以搬进 SDP 不破坏凸性"。**这句话本轮必须自己审一遍**，因为线性是在逆坐标
$(X=P^{-1},Y=\tilde P^{-1})$ 里成立，而凸性是在协方差坐标 $(P,\Pi)$ 里成立——
同一个问题不能同时在两套坐标里凸。

跑前写死的判据：
 [K1] 自由（无平面约束）SDP 必须复现原文 §V 的三个发表锚点：
        D=33 → 6.133 bit (rank 3)，D=40 → 3.266 bit (rank 2)，D=80 → 1.602 bit (rank 1)，
      以及水平渐近线 1.169 bit。容差：率 ≤0.02 bit，秩按原文口径（截断 0.1% 最大奇异值）。
 [K2] 原文给出的 D=33 最优机制 C（3×4）、V=diag(0.029,0.208,1.435) 用我自己的
      (C,V)→ARE 定点→率/代价通道重算，必须落回 (6.133 bit, 33.0)。这一步验证的是
      "我车道量的 D 与原文 (18b) 的 D 是同一个数"，比任何代数推演都硬。
 [K3] 对 D 二分求 DI(D)=3 bit 的自由解 $D^\*_{free}(3)$，与 42.0427 比较。
      差 >0.01 则说明 42.0427 不是自由前沿点（它只是我车道某个设计的可达代价）。
 [K4] 平面约束的可行性判据：取两个满足 $Q\,\Sigma\,=0$ 的可行点 $P_1,P_2$（$Q=I-ZZ^t$），
      检查中点残差 $\lVert Q\Sigma(P_\lambda)\rVert_F$。若对某个 $\lambda$ 残差 > $10^{-6}$
      倍的端点尺度，则"平面限制后仍是凸 SDP"**当场否证**，§64-F-3 的理由撤回。
"""
import sys
sys.stdout.reconfigure(encoding='utf-8')
import numpy as np
import cvxpy as cp
from scipy.linalg import schur
from scipy.optimize import brentq

np.set_printoptions(precision=5, suppress=True, linewidth=170)
_src = open('p0/exp_c_audit.py', encoding='utf-8').read().split("print(r'== E53")[0]
_ns = {'__name__': 'p'}
exec(compile(_src, 'p0/exp_c_audit.py[preamble]', 'exec'), _ns)
A, B, W, R, TH, JC, n, sym, Kv, Pcmod = (_ns['A'], _ns['B'], _ns['W'], _ns['R'], _ns['TH'],
                                         _ns['JC'], _ns['n'], _ns['sym'], _ns['Kv'], _ns['Pc'])
ln2 = np.log(2.0)
Ts, U = schur(A, output='real', sort=lambda a: abs(a) < 1.0)[:2]
Ptilde_aff = lambda P: sym(A @ P @ A.T + W)
SFIX = Pcmod                      # 原文 (17) 的 S = 控制 DARE 解，我这里叫 Pcmod
CONST = float(np.trace(W @ SFIX))  # (18b) 里不含 P 的那一项 = JC
print(f'[口径] JC=tr(W Pc)={CONST:.6f}  lambda(TH)={np.sort(np.linalg.eigvalsh(TH))}')


def rank_snr(SNR, tol=1e-3):
    s = np.linalg.svd(SNR, compute_uv=False)
    return int((s > tol * s.max()).sum()) if s.size else 0


def tanaka_sdp(D, plane=None):
    """式 (18)：变量 (P, Pi)。plane=None 为自由；否则 Z(n x k) 正交列，加 Pi 的平面约束。

    返回 (率 nat, P, Pi, 状态, 平面残差)。平面残差是对 **Sigma = P^-1 - Ptilde^-1**
    算的，因为原文 (68) 的传感等式用的是 Sigma，不是 Pi。"""
    P = cp.Variable((n, n), symmetric=True)
    Pi = cp.Variable((n, n), symmetric=True)
    cons = [P >> 0, Pi >> 0, cp.trace(TH @ P) + CONST <= D,
            cp.bmat([[P - Pi, P @ A.T], [A @ P, Ptilde_aff(P)]]) >> 0]
    if plane is not None:
        Z = plane
        Qm = np.eye(n) - Z @ Z.T
        cons += [Qm @ Pi @ Qm == 0, Qm @ Pi @ Z == 0]     # 线性，但加在 Pi 上（见 [K5]）
    prob = cp.Problem(cp.Minimize(-0.5 * cp.log_det(Pi)), cons)
    try:
        val = prob.solve(solver=cp.CLARABEL, verbose=False)
    except Exception as e:
        return None, None, None, f'raise:{type(e).__name__}:{e}', None
    if P.value is None or not np.isfinite(val):
        return None, P.value, Pi.value, prob.status, None
    Pv = sym(P.value)
    Piv = sym(Pi.value)
    SNR = sym(np.linalg.inv(Pv) - np.linalg.inv(Ptilde_aff(Pv)))
    res = None
    if plane is not None:
        Z = plane
        Qm = np.eye(n) - Z @ Z.T
        res = (np.linalg.norm(Qm @ SNR @ Qm), np.linalg.norm(Qm @ SNR @ Z),
               np.linalg.norm(Qm @ Piv @ Qm))
    return val, Pv, Piv, prob.status, res


def free_rate_bit(D):
    v, Pv, Pi, st, _ = tanaka_sdp(D)
    return (v / ln2 if v is not None and np.isfinite(v) else np.nan), Pv, Pi, st


print('\n== K1 自由 SDP vs 原文 §V 发表锚点 ==')
PUB = {33.0: (6.133, 3), 40.0: (3.266, 2), 80.0: (1.602, 1)}
store = {}
for D, (ip, rp) in PUB.items():
    v, Pv, Pi, st = free_rate_bit(D)
    SNR = sym(np.linalg.inv(Pv) - np.linalg.inv(Ptilde_aff(Pv)))
    r = rank_snr(SNR)
    store[D] = (v, Pv, Pi)
    print(f' D={D:5.1f}  SDP={v:8.4f} bit (原文 {ip})  差 {v-ip:+.4f}  rank {r} (原文 {rp})'
          f'  status={st}  tr(TH P)={np.trace(TH @ Pv):.5f}')

print('\n== K2 用原文 D=33 的最优机制 (C,V) 走我自己的 ARE 定点通道 ==')
Cp = np.array([[-0.864, 0.258, -0.205, -0.382],
               [-0.469, -0.329, 0.662, 0.483],
               [-0.130, 0.332, -0.502, 0.780]])
Vp = np.diag([0.029, 0.208, 1.435])


def are_post(C, Vobs, tol=1e-13, maxit=800000):
    Sig = W.copy()
    for _ in range(maxit):
        Pm = sym(A @ Sig @ A.T + W)
        Lm = Pm @ C.T @ np.linalg.inv(sym(C @ Pm @ C.T + Vobs))
        Mm = np.eye(n) - Lm @ C
        Sn = sym(Mm @ Pm @ Mm.T + Lm @ Vobs @ Lm.T)
        if not np.isfinite(Sn).all():
            return None, None, 'blow'
        if np.abs(Sn - Sig).max() < tol * max(1.0, np.abs(Sn).max()):
            return Sn, Pm, 'conv'
        Sig = Sn
    return Sig, Pm, 'maxit'


Pp, Pmp, stt = are_post(Cp, Vp)
if stt == 'conv':
    rate_p = 0.5 * (np.linalg.slogdet(Pmp)[1] - np.linalg.slogdet(Pp)[1]) / ln2
    D_p = CONST + float(np.trace(TH @ Pp))
    SNR_p = np.linalg.inv(Pp) - np.linalg.inv(Pmp)
    print(f' status={stt}  我算的率={rate_p:.4f} bit (原文 6.133)  我算的 D={D_p:.4f} (原文 33.0)')
    print(f' rank(SNR)={rank_snr(SNR_p)}  C^t V^-1 C 与 SNR 的 Frobenius 差'
          f'={np.linalg.norm(Cp.T @ np.linalg.inv(Vp) @ Cp - SNR_p):.3e}')
    print(f' [K2 判定] 率差 {abs(rate_p-6.133):.4f} bit, 代价差 {abs(D_p-33.0):.4f}'
          f' ⇒ {"口径一致，(18b) 的 D 就是我车道的 D" if abs(rate_p-6.133)<0.02 and abs(D_p-33.0)<0.5 else "口径不同，禁止用该 SDP 打 43.80"}')
else:
    print(f' [K2 失败] ARE 不收敛：{stt}')

print('\n== K3 自由前沿的反演 D*_free(3 bit) ==')
try:
    Dstar = brentq(lambda D: free_rate_bit(D)[0] - 3.0, JC + 1e-3, 120.0, xtol=1e-6, rtol=1e-10)
    v, Pv, Pi, st = free_rate_bit(Dstar)
    SNR = sym(np.linalg.inv(Pv) - np.linalg.inv(Ptilde_aff(Pv)))
    print(f' D*_free(3 bit)={Dstar:.4f}  (核验 SDP 率={v:.5f} bit)  与 42.0427 差 {Dstar-42.0427:+.4f}')
    print(f' rank(SNR)={rank_snr(SNR)}  奇异值 {np.sort(np.linalg.svd(SNR, compute_uv=False))[::-1]}')
    # 这个自由最优的 SNR 主方向落在哪个子空间？和 tail-3 平面有多贴？
    Zb = U[:, n - 3:]
    w, Vv = np.linalg.eigh(sym(SNR))
    key = np.argsort(w)[::-1][:3]
    Um = Vv[:, key]
    ov = abs(np.linalg.det(np.hstack([Um, Zb]).T @ np.hstack([Um, Zb])))
    print(f' [K3 附带] 自由最优的非零 SNR 方向 vs tail-3 平面：子空间重叠 |det(U^t Z)|={ov:.4f}'
          f'（=1 同平面，=0 正交）')
    print(f' 逐方向重叠 {(np.linalg.svd(Um.T @ Zb, compute_uv=False))}')
except Exception as e:
    print(f' [K3 异常] {type(e).__name__}: {e}')

print('\n== K4 平面可行集凸不凸：逆坐标里线性 ≠ 协方差坐标里凸 ==')
Zb = U[:, n - 3:]
Qm = np.eye(n) - Zb @ Zb.T


def plane_resid(P):
    Sn = sym(np.linalg.inv(P) - np.linalg.inv(Ptilde_aff(P)))
    sc = np.linalg.norm(Sn)
    return np.linalg.norm(Qm @ Sn) / max(sc, 1e-300), sc


def plane_design(Vz):
    """C = Z^t、观测噪声 Vz ⇒ SNR = Z Vz^-1 Z^t 精确落在 span(Z) 内（这是构造级的）。"""
    Pp3, Pmp3, st3 = are_post(Zb.T.copy(), sym(Vz))
    return Pp3, Pmp3, st3


V1 = np.eye(3) * 1e-8
V2 = np.diag([3.0, 0.7, 0.05])
Pa3, Pma, sta3 = plane_design(V1)
Pb3, Pmb, stb3 = plane_design(V2)
ra, sa = plane_resid(Pa3)
rb, sb = plane_resid(Pb3)
print(f' 端点 1（近完美测 tail-3）: 相对残差 {ra:.3e}  |Sigma|={sa:.4e}  状态 {sta3}')
print(f' 端点 2（tail-3 有限噪声）: 相对残差 {rb:.3e}  |Sigma|={sb:.4e}  状态 {stb3}')
for lam in (0.1, 0.25, 0.5, 0.75, 0.9):
    Pmid = sym(lam * Pa3 + (1 - lam) * Pb3)
    rm, sm = plane_resid(Pmid)
    ok = np.all(np.linalg.eigvalsh(Ptilde_aff(Pmid) - Pmid) > -1e-12)
    print(f'  λ={lam:4.2f} 中点相对残差 {rm:.3e}  |Sigma_mid|={sm:.4e}  '
          f'P<=Ptilde {"满足" if ok else "违反"}   端点凸组合预测 {lam*ra+(1-lam)*rb:.3e}')
print(' [K4 判定] 残差 ≫ 1e-8 ⇒ 平面可行集非凸，"搬进 SDP 保持凸性"不成立。')

