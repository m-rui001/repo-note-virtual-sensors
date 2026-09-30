r"""E96b = 诊断：我的 (18) 复现在 D=33 与原文差 0.03 bit（可接受），D=40/80 差 0.48/0.58 bit。
要么我少写了约束，要么原文那两个数不是同一个问题。这个脚本把 SDP 解拿回到**独立通道**
（K2 那条：SNR→(C,V)→ARE 定点→率/代价）复核，并打印 Π 的松紧、SNR 谱、以及
"如果按原文口径把 Σ 的秩截到 r 会得到什么"。
"""
import sys
sys.stdout.reconfigure(encoding='utf-8')
import numpy as np
import cvxpy as cp
from scipy.linalg import schur

np.set_printoptions(precision=5, suppress=True, linewidth=170)
_src = open('p0/exp_c_audit.py', encoding='utf-8').read().split("print(r'== E53")[0]
_ns = {'__name__': 'p'}
exec(compile(_src, 'p0/exp_c_audit.py[preamble]', 'exec'), _ns)
A, B, W, R, TH, JC, n, sym, Kv, Pcmod = (_ns['A'], _ns['B'], _ns['W'], _ns['R'], _ns['TH'],
                                         _ns['JC'], _ns['n'], _ns['sym'], _ns['Kv'], _ns['Pc'])
ln2 = np.log(2.0)
LOGDET_W = float(np.linalg.slogdet(W)[1])
Ptilde_aff = lambda P: sym(A @ P @ A.T + W)


def solve(D, rankmax=None, solver=cp.CLARABEL):
    P = cp.Variable((n, n), symmetric=True)
    Pi = cp.Variable((n, n), symmetric=True)
    cons = [P >> 1e-12 * np.eye(n), Pi >> 0, P << Ptilde_aff(P),
            cp.trace(TH @ P) + JC <= D,
            cp.bmat([[P - Pi, P @ A.T], [A @ P, Ptilde_aff(P)]]) >> 0]
    prob = cp.Problem(cp.Minimize(-0.5 * cp.log_det(Pi)), cons)
    prob.solve(solver=solver)
    return prob.value, sym(P.value), sym(Pi.value), prob.status


def indep(P):
    """把 P 当作设计，走 SNR→(C,V)→ARE 定点，返回 (率 bit, 代价, 后验, 定点状态)。"""
    Pt = Ptilde_aff(P)
    SNR = sym(np.linalg.inv(P) - np.linalg.inv(Pt))
    w, Vv = np.linalg.eigh(SNR)
    idx = np.where(w > 1e-9 * max(w.max(), 1.0))[0]
    k = len(idx)
    C = np.sqrt(np.maximum(w[idx], 0))[:, None] * Vv[:, idx].T     # C^t V^-1 C = SNR, V = I
    Sig = W.copy()
    st = 'maxit'
    for _ in range(400000):
        Pm = sym(A @ Sig @ A.T + W)
        L = Pm @ C.T @ np.linalg.inv(C @ Pm @ C.T + np.eye(k))
        M = np.eye(n) - L @ C
        Sn = sym(M @ Pm @ M.T + L @ L.T)
        if np.abs(Sn - Sig).max() < 1e-14 * max(1.0, np.abs(Sn).max()):
            Sig = Sn
            st = 'conv'
            break
        Sig = Sn
    Pe = np.linalg.eigvals((np.eye(n) - L @ C) @ A)
    rate = 0.5 * (np.linalg.slogdet(Pm)[1] - np.linalg.slogdet(Sig)[1]) / ln2
    return rate, JC + float(np.trace(TH @ Sig)), Sig, st, np.abs(Sig - P).max(), \
        float(np.abs(Pe).max()), k, SNR


for D in (33.0, 40.0, 80.0):
    val, Pv, Piv, st = solve(D)
    sdp_rate = (val + 0.5 * LOGDET_W) / ln2
    Pt = Ptilde_aff(Pv)
    Sch = sym(Pv - Pv @ A.T @ np.linalg.inv(Pt) @ A @ Pv)
    lamPi = np.linalg.eigvalsh(Piv)
    lamSch = np.linalg.eigvalsh(Sch)
    r2, c2, Sig2, stt, gap, rho, k2, SNR = indep(Pv)
    sv = np.sort(np.linalg.svd(SNR, compute_uv=False))[::-1]
    print(f'\n=== D={D} ===')
    print(f' SDP: status={st}  目标(补上 ½logdet W 后)={sdp_rate:.4f} bit')
    print(f' Π 特征值 {lamPi}   Schur 上界特征值 {lamSch}')
    print(f' Π 与上界的最大差 {np.abs(lamPi - lamSch).max():.3e}  ⇒ LMI {"紧" if np.abs(lamPi-lamSch).max()<1e-6 else "不紧"}')
    print(f' 独立通道: 率={r2:.4f} bit  代价={c2:.4f}  |P_sdp - P_ARE|={gap:.3e}  定点={stt}  rank={k2}')
    print(f' SNR 奇异值 {sv}   (0.1% 截断秩 = {int((sv > 1e-3*sv.max()).sum())})')
    print(f' 误差动态谱半径 rho((I-LC)A) = {rho:.6f}')

print('\n== 对照：原文三条发表机制各自走独立通道（K2 口径） ==')
PUBMECH = {
    33.0: (np.array([[-0.864, 0.258, -0.205, -0.382],
                     [-0.469, -0.329, 0.662, 0.483],
                     [-0.130, 0.332, -0.502, 0.780]]), np.diag([0.029, 0.208, 1.435]), 6.133),
    40.0: (np.array([[-0.886, 0.241, -0.170, -0.359],
                     [-0.431, -0.350, 0.647, 0.523]]), np.diag([0.208, 2.413]), 3.266),
    80.0: (np.array([[-0.876, 0.271, -0.169, -0.362]]), np.array([[1.775]]), 1.602),
}


def are_post(Dc, Vv_, Cx):
    Sig = W.copy()
    Lm = None
    Pm = None
    for _ in range(Dc):
        Pm = sym(A @ Sig @ A.T + W)
        Lm = Pm @ Cx.T @ np.linalg.inv(sym(Cx @ Pm @ Cx.T + Vv_))
        Mm = np.eye(n) - Lm @ Cx
        Sn = sym(Mm @ Pm @ Mm.T + Lm @ Vv_ @ Lm.T)
        if not np.isfinite(Sn).all():
            return None, None, None, 'blow'
        cv = np.abs(Sn - Sig).max() < 1e-13 * max(1.0, np.abs(Sn).max())
        Sig = Sn
        if cv:
            return Sig, Pm, Lm, 'conv'
    return Sig, Pm, Lm, 'maxit'


for D, (Cc, Vobs_, ip) in PUBMECH.items():
    Ps, Pms, Ls, stt = are_post(800000, Vobs_, Cc)
    if stt != 'conv':
        print(f' D={D}: 定点 {stt}')
        continue
    rate = 0.5 * (np.linalg.slogdet(Pms)[1] - np.linalg.slogdet(Ps)[1]) / ln2
    cost = JC + float(np.trace(TH @ Ps))
    SNRp = np.linalg.inv(Ps) - np.linalg.inv(Pms)
    rho = float(np.abs(np.linalg.eigvals((np.eye(n) - Ls @ Cc) @ A)).max())
    sv = np.sort(np.linalg.svd(SNRp, compute_uv=False))[::-1]
    print(f' D={D}: 我算的率={rate:.4f} bit (原文 {ip})  我算的代价={cost:.4f}'
          f'  rho((I-LC)A)={rho:.4f}  sv={sv}'
          f'  |C^tV^-1C-SNR|={np.linalg.norm(Cc.T @ np.linalg.inv(Vobs_) @ Cc - SNRp):.2e}'
          f'  定点={stt}')

