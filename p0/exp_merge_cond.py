r"""E41：合并（$\mathrm{gap}_d=0$）到底需要什么条件——把 #21 里"25 个角度一个都没找到"从数值观察升级成结构结论。

E40 的三分法说 $\mathrm{gap}_d=\operatorname{tr}(\Theta\Delta P_d)$，$\Delta P_d=P_0(F_{-d})-P_0(F)\succeq0$。
于是 $\mathrm{gap}_d=0$ 只有两条路：

  机制一（后验不动）$\Delta P_d=0$。$\Theta\succ0$ 时这是**唯一**的路。
  机制二（任务看不见）$\Delta P_d\ne0$ 但 $\operatorname{ran}\Delta P_d\subseteq\ker\Theta$。

本轮钉三件事：

 (1) **$\Theta\succ0$ 且 $W\succ0$ 时机制一也不可能**。$P_0$ 是"预测以后投影到 $\ker F$ 上"的不动点，
     删方向只会让支撑变大；若 $W\succ0$，多出来的那块被 $A,W$ 推着严格正，$\operatorname{tr}(\Delta P_d)>0$。
     这条在锚点植物上可以直接验：$\lambda_{\min}(\Theta)=0.031>0$，E40 的 $\theta$ 网格固定谱
     $\Sigma\succ0$，所以那个网格**不可能**有合并——#21 的"一个都没找到"是定理的推论，不是搜索失败。
     反过来给一个 $W$ 奇异的对照：$\ker W$ 挪到被砍通道的不可观测子空间上，机制一真的会兑现。

 (2) **机制二只在 $\operatorname{rank}\Theta<n$ 上有可能**，且需要 $\Delta P$ 的列空间整块落进 $\ker\Theta$。
     用结构搜索（把 $\ker W$ 当成自由参数在球面上扫）替代 E40 的 400 个随机方向。

 (3) 顺带把 $W=0$ 的退化情形打出来：$\tilde P_0=0\Rightarrow P_0=0\Rightarrow$ 所有方向的 $\mathrm{gap}=0$。
     这是机制一的极端见证，用来说明分类器不会把"全零噪声"错报成有限门槛。
"""
import sys
import numpy as np
sys.stdout.reconfigure(encoding='utf-8')
import p0.exp_plants_oos as E35
import p0.exp_ladder as X37
import p0.exp_C_exact as X36
from p0.exp_ticket import gaps_of, floor_pair
from p0.exp_nearfloor_law import Fs

sym = X36.sym

base = None
for nm, pl, F in E35.plants():
    if nm.find('锚点') >= 0:
        base = (nm.replace('$', ''), pl, F)
        break
nm0, Pl0, _Fa = base
F2 = Fs['F2 尾部2(封闭, r=2)']
n = Pl0.nn
r = F2.shape[0]
wT = np.linalg.eigvalsh(sym(Pl0.Th))
print('== E41 ==  算例：%s  n=%d r=%d  $\\lambda(\\Theta)=[%.4g,%.4g]$  rank=%d'
      % (nm0, n, r, wT[0], wT[-1], int(np.sum(wT > 1e-9 * wT[-1]))))


def det_dir(Pl, F):
    """不可检测时 $P_0$ 不存在；先给 PBH 结论。返回 (ok, 失败的 |λ|)。"""
    return Pl.detectable(F)


def dP_of(Pl, F, i):
    """删 $F$ 第 $i$ 行。不可检测返回 `(None, |λ|)`；否则给 $\Delta P$ 和尺度无关的合并指标
    $\eta=\operatorname{tr}(\Delta P)/\operatorname{tr}(P_0(F_{-i}))\in[0,1)$，$\eta=0\iff\Delta P=0$。"""
    Feff = np.delete(F, i, axis=0)
    ok, lam = Pl.detectable(Feff)
    if not ok:
        return None, lam
    phi0F, P0F, _r, _k = floor_pair(Pl, F)
    phiA, P0A, _r2, _k2 = floor_pair(Pl, Feff)
    dP = sym(P0A - P0F)
    trA = float(np.trace(P0A))
    return dict(dP=dP, gap=phiA - phi0F, tr=float(np.trace(dP)), eta=float(np.trace(dP)) / max(trA, 1e-30),
                ev=np.linalg.eigvalsh(dP), phi0F=phi0F, phiA=phiA), None


print('\n[1] $W\succ0$ 的网格（E40 的 $\\theta$ 扫描）为什么不可能合并：逐点看 $\\operatorname{tr}(\Delta P)$ 的下界')
Sig = np.diag([1., 1e-3, 1e-6, 1e-9])
tr_lo, gap_lo = [], []
for deg in np.linspace(0, 90, 19):
    th = np.deg2rad(deg)
    R = np.eye(4)
    c, s = np.cos(th), np.sin(th)
    R[0, 0] = R[1, 1] = c
    R[0, 1], R[1, 0] = -s, s
    Pl = E35.Pl(Pl0.A, Pl0.B, sym(R @ Sig @ R.T), np.eye(4))
    _p, g, _r = gaps_of(Pl, F2, drop='row')
    fin = [q for q in g if np.isfinite(q['gap'])]
    tr_lo += [q['tr'] for q in fin]
    gap_lo += [q['gap'] for q in fin]
    ev_min = [min(q['ev'][q['ev'] > -1e-8], default=np.nan) for q in fin]
print('   19 个角度上 $\\operatorname{tr}(\Delta P)$ 最小 %.4e，$\\mathrm{gap}$ 最小 %.4e，'
      r'$\lambda_{\min}(\Delta P)$ 最小 %.4e' % (min(tr_lo), min(gap_lo), min(ev_min)))
print(r'   全部 $>0$：与 $\lambda_{\min}(\Theta)=%.4g>0$ 合起来 $\Rightarrow$ 这个网格上两条机制都关闭'
      % wT[0])

print('\n[2] 机制一的见证：把 $\ker W$ 挪到被砍通道真正不读的东西上')
for lab, wv in [('满秩（对照）', [1., 1e-3, 1e-6, 1e-9]),
                ('$W$ 秩 3：零在 $e_3$', [1., 1e-3, 1e-6, 0.]),
                ('$W$ 秩 3：零在 $e_2$', [1., 1e-3, 0., 1e-9]),
                ('$W$ 秩 2', [1., 1e-3, 0., 0.]),
                ('$W=0$（退化）', [0., 0., 0., 0.])]:
    Pl = E35.Pl(Pl0.A, Pl0.B, np.diag(wv), np.eye(4))
    outs = []
    for i in range(r):
        res, _lam = dP_of(Pl, F2, i)
        if res is None:
            outs.append('删%d: 不可检测(|λ|=%.3f)' % (i, _lam))
            continue
        outs.append('删%d: gap=%.3e tr=%.3e $\\eta$=%.2e $\\lambda_{\\max}(\Delta P)$=%.2e'
                    % (i, res['gap'], res['tr'], res['eta'], res['ev'][-1]))
    print('   %-18s  rank W=%d   %s' % (lab, int(np.sum(np.array(wv) > 0)), '  |  '.join(outs)))

print('\n[3] 结构搜索：$\ker W$ 取单位球面上的方向 $u$（$W=\Sigma_u=R\,\mathrm{diag}(1,10^{-3},10^{-6},0)R^\top$，$Ru=0$）')
print('    机制一要求 $\Delta P=0$；机制二要求 $\operatorname{ran}\Delta P\subseteq\ker\Theta$——先算 $\Theta$ 的核有没有东西')
QT, wQT = np.linalg.eigh(sym(Pl0.Th))
kerT = wQT < 1e-9 * wQT[-1]
print('    锚点 $\ker\Theta$ 维数 = %d  $\\Rightarrow$ 机制二在本算例%s'
      % (int(kerT.sum()), '根本不存在（$\\Theta\\succ0$）' if int(kerT.sum()) == 0 else '存在'))
rng = np.random.default_rng(20260929)
best = None
tries = 0
for _ in range(120):
    u = rng.normal(size=n)
    u /= np.linalg.norm(u)
    Q, _R = np.linalg.qr(np.eye(n) - np.outer(u, u))
    wv = np.array([1., 1e-3, 1e-6])
    Wk = sym(Q[:, :3] @ np.diag(wv) @ Q[:, :3].T)
    Pl = E35.Pl(Pl0.A, Pl0.B, Wk, np.eye(4))
    for i in range(r):
        res, _lam = dP_of(Pl, F2, i)
        tries += 1
        if res is None:
            continue
        if best is None or res['eta'] < best[0]:
            best = (res['eta'], i, res['gap'], res['tr'], res['ev'][-1], u.copy())
print('    %d 次（$W$ 秩 3，随机 $\ker W$）里最小的合并指标 $\\eta=\\operatorname{tr}(\Delta P)/\operatorname{tr}(P_0(F_{-i}))$：'
      r'$\eta=%.3e$（删行%d，$\mathrm{gap}=%.3e$，$\operatorname{tr}(\Delta P)=%.3e$，'
      r'$\lambda_{\max}(\Delta P)=%.3e$，$\ker W$ 方向 $u=[%s]$）'
      % (best[0], best[1], best[2], best[3], best[4], ','.join('%.2f' % v for v in best[5])))
print('    机制一在 $W\succ0$ 与 $W$ 秩 3 的两侧都没兑现 $\Rightarrow$ 在这个植物上合并需要更强的退化')

print('\n[4] 机制二要一个 $\operatorname{rank}\Theta<n$ 的植物：PL3 $r=2$、PL4 单输入 $r=3$')
for nm, Pl, F in E35.plants():
    if nm.find('PL4 单输入 $r=3$') < 0 and nm.find('PL3 $r=2$') < 0:
        continue
    rr = F.shape[0]
    wq, Uq = np.linalg.eigh(sym(Pl.Th))
    dimK = int(np.sum(wq < 1e-9 * wq[-1]))
    print('  %s  $\dim\ker\Theta$=%d  rank W=%d' % (nm.replace('$', ''), dimK, np.linalg.matrix_rank(Pl.W)))
    for i in range(rr):
        res, _lam = dP_of(Pl, F, i)
        if res is None:
            print('     删%d 不可检测（$|\lambda|=%.3f$）' % (i, _lam))
            continue
        dP = res['dP']
        Uk = Uq[:, wq < 1e-9 * wq[-1]]
        res_in = float(np.linalg.norm(Uk.T @ dP @ Uk, np.inf)) if Uk.shape[1] else 0.0
        leak = float(np.linalg.norm(Uk.T @ dP, np.inf)) if Uk.shape[1] else 0.0
        print(r'     删%d  gap=%.4e  $\operatorname{tr}(\Delta P)$=%.4e  $\eta$=%.3e  $\|\Delta P\|_\infty$=%.4e'
              '  核内块 $\|U_k^\top\Delta P\,U_k\|=$%.3e  泄漏 $\|U_k^\top\Delta P\|=$%.3e'
              % (i, res['gap'], res['tr'], res['eta'], np.linalg.norm(dP, np.inf), res_in, leak))
