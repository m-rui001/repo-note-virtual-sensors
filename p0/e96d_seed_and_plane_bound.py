r"""E96d = 把 E96c 的 [V1] 失败**定量归因**，并给出 tail-3 平面的**合法上界**（不依赖全局搜索）。

E96c 的判据 [V1] 已经失败：同一条 Nelder-Mead 多起点通道在**自由** 10 维族里只做到
$D=43.6682$，而式 (18) 的 SDP 精确值是 $D^*_{free}(3\,\text{bit})=42.0402$ ⇒ 搜索偏差
$+1.6281$。按我自己写死的规则，平面数字一律不得当作"最优"引用。本脚本要区分两件事：

 (A) **评价通路**（$\Sigma=ZTZ^t\to$ 滤波 DARE $\to(D,I)$）是否正确。
     检验法：把 SDP 的解直接当起点塞进同一条通道。E96c 的 (2a) 已给 $D\;42.0402$，
     这里要求它在**带起点的搜索**里仍成立（seeded 搜索 best 必须 $=D^*$，偏差 $\le10^{-3}$）。
 (B) **搜索**是否只是掉进局部极小。若是，(A) 通过 + (B) 失败 ⇒ 平面值只能作为**上界**报告。

平面族的合法上界根本不需要全局搜索：把自由最优的信息矩阵**投影**进平面
（$T=\mathcal Z^t\Sigma^*\mathcal Z$，再沿射线精确二分到 $I=3$）就得到一个**严格可行**的
平面设计，于是
$$D_{\text{plane}}(3\,\text{bit})\le D_{\text{seed}}.$$
这条上界与"45.4537 是我车道 Powell 上界"是同一级别的陈述；它的用途是把 Remark 2 的
数字挤进一个**有上界、有下界**的区间：下界 $D^*_{free}=42.0402$（SDP 全局最优，任何受限族都不可能在更低价）。

== 判据（跑前写死）==
 [W1] seeded 自由搜索的 best 与 $D^*_{free}$ 差 $\le 10^{-3}$。否则评价通路有 bug，
      本轮全部平面数字作废（包括 $D_{\text{seed}}$）。
 [W2] $D_{\text{seed}}$ 与随机搜索 best（E96c 第 (3) 节）取小者为当前平面**上界**；
      若该上界 $<43.711-0.007$ ⇒ 蓝线在该平面可行，Remark 2 的"蓝线可疑"塌；
      若 $\ge 43.80$ ⇒ 只报"上界尚未降到蓝线以下"，**不构成任何认证**。
 [W3] 报出"搜索偏差度量"：自由族里 (随机搜索 best $-$ SDP 精确值) $=+1.6281$，
      平面族里用 (对角子族 best $-$ 全平面随机 best) 作为同一通道的族内偏差参照。
"""
import sys
import time as _t
sys.stdout.reconfigure(encoding='utf-8')
import numpy as np
from scipy.optimize import minimize, brentq

np.set_printoptions(precision=5, suppress=True, linewidth=170)
_src = open('p0/e96c_plane_frontier.py', encoding='utf-8').read().split("print('== (1)")[0]
_ns = {'__name__': 'e96d'}
exec(compile(_src, 'p0/e96c_plane_frontier.py[helpers]', 'exec'), _ns)
globals().update({k: _ns[k] for k in
                  ('A', 'W', 'TH', 'JC', 'n', 'sym', 'Ts', 'U', 'ln2', 'free_sdp',
                   'design', 'scale_to', 'Ptil')})


def coords_full(T, k):
    r"""把 PSD 的 $T$ 写成上三角因子 $L$（$T=LL^t$）的上三角坐标。
    反序共轭恒等式：$L=P\,\mathrm{chol}(P T P)\,P$，$P$ 为反序置换（numpy 的 [::-1,::-1]）。"""
    rev = lambda X: X[::-1, ::-1]
    for jit in (1e-14, 1e-11, 1e-8, 1e-6):
        try:
            S = np.linalg.cholesky(rev(T + jit * np.eye(k)))
            break
        except np.linalg.LinAlgError:
            continue
    L = rev(S)
    iu = np.triu_indices(k)
    Chk = np.zeros((k, k))
    Chk[iu] = L[iu]
    return Chk, Chk[iu]


def search2(Z, I0, T_seeds=(), nstart=25, seed=0, lam=300.0, budget=3000):
    """与 E96c 的 search 同一条通道，只多一件事：$T_{seeds}$ 的坐标直接进候选池。"""
    k = Z.shape[1]
    iu = np.triu_indices(k)
    m = len(iu[0])
    rng = np.random.default_rng(seed)

    def toT(x):
        Lm = np.zeros((k, k))
        Lm[iu] = x
        return Lm @ Lm.T

    def pen_obj(x):
        D, I, st = design(toT(x), Z)
        if st != 'ok':
            return 1e7 + float(np.dot(x, x))
        return D + lam * abs(I - I0)

    cands = []
    for Ts in T_seeds:
        Chk, x = coords_full(sym(Ts), k)
        rec = np.linalg.norm(Chk @ Chk.T - Ts)
        print(f'   种子坐标重建误差 ||LL^t-T||={rec:.3e}', flush=True)
        for _ in range(2):
            x = minimize(pen_obj, x, method='Nelder-Mead',
                         options=dict(maxiter=budget, maxfev=budget, xatol=1e-11, fatol=1e-13)).x
        cands.append(x)
    t0 = _t.time()
    for t in range(nstart):
        x0 = rng.normal(0.0, 0.7, m)
        for _ in range(2):
            x0 = minimize(pen_obj, x0, method='Nelder-Mead',
                          options=dict(maxiter=budget, maxfev=budget, xatol=1e-11, fatol=1e-13)).x
        cands.append(x0)
        if (t + 1) % 10 == 0:
            print(f'   [{k}维] 随机起点 {t+1}/{nstart} 用时 {_t.time()-t0:.0f}s', flush=True)
    best = (np.inf, None, np.inf, None)
    feas = []
    seen = []
    for i, x in enumerate(cands):
        D, I, Tt, st = scale_to(toT(x), Z, I0)
        tag = 'SEED' if i < len(T_seeds) else f'rnd{i-len(T_seeds)}'
        if st == 'ok':
            feas.append(D)
            seen.append((D, tag))
            if D < best[0]:
                best = (D, Tt, I, tag)
    seen.sort()
    grp = (' '.join(f'{d:.4f}:{t}' for d, t in seen[:6]))
    print(f'   前 6 个可行点（D:来源）{grp}', flush=True)
    rnd = [d for d, t in seen if t != 'SEED']
    return best, np.sort(np.array(feas)) if feas else np.array([np.inf]), (min(rnd) if rnd else np.inf)


print('== (0) 重算自由精确前沿上的 3 bit 点 ==')
Dstar = brentq(lambda D: free_sdp(D)[0] - 3.0, 34.0, 90.0, xtol=1e-5)
v, Pstar, st = free_sdp(Dstar)
SNR = sym(np.linalg.inv(Pstar) - np.linalg.inv(Ptil(Pstar)))
print(f' D*_free(3 bit) = {Dstar:.4f}（DI={v:.5f}, status={st}）  与 42.0427 差 {Dstar-42.0427:+.4f}')

Z4 = U[:, :4]
Zb = U[:, n - 3:]
T_seed_free = sym(Z4.T @ SNR @ Z4)
T_seed_pl = sym(Zb.T @ SNR @ Zb)          # 把自由最优**投影**进 tail-3 平面
tw, tv = np.linalg.eigh(T_seed_pl)
T_seed_pl = tv @ np.diag(np.maximum(tw, 0.0)) @ tv.T

print('\n== (1) 投影种子的平面可行点（不需要搜索，直接是上界）==')
Ds_raw, Is_raw, st_raw = design(T_seed_pl, Zb, check=True)
print(f' 投影未重标定：D={Ds_raw:.4f} I={Is_raw:.5f} 状态 {st_raw}')
Ds, Is, Tt_s, st_s = scale_to(T_seed_pl, Zb, 3.0)
print(f' [W2] 沿射线二分到 I=3：D_seed = {Ds:.4f}  I={Is:.6f}  状态 {st_s}'
      f'   ⇒ 平面 3 bit 的上界 {Ds:.4f}')
print(f'      与蓝线 43.711 差 {Ds-43.711:+.4f}   与验收线 43.80 差 {Ds-43.80:+.4f}'
      f'   与我车道 Powell 上界 45.4537 差 {Ds-45.4537:+.4f}')
print(f'      投影丢掉的率：||Q Sigma||/||Sigma||={np.linalg.norm((np.eye(n)-Zb@Zb.T)@SNR)/np.linalg.norm(SNR):.3e}')

print('\n== (2) seeded 自由搜索：评价通路必须在种子处复现精确值 [W1] ==')
b4, f4, r4 = search2(Z4, 3.0, T_seeds=(T_seed_free,), nstart=15, seed=11)
print(f' seeded 自由 best D={b4[0]:.5f}（来源 {b4[3]}）  SDP D*={Dstar:.5f}  差 {b4[0]-Dstar:+.5f}'
      f'   [W1] {"通过" if abs(b4[0]-Dstar)<1e-3 else "失败——评价通路有 bug，本轮平面数字全废"}')
print(f'   自由族：随机起点 best {r4:.4f}  随机 best - SDP = {r4-Dstar:+.4f}（E96c 25 起点时 +1.6281）'
      f'   全部可行 {int(np.isfinite(f4).sum())}/16')

print('\n== (3) seeded tail-3 搜索：上界能否压到蓝线以下 [W2] ==')
b3, f3, r3 = search2(Zb, 3.0, T_seeds=(T_seed_pl,), nstart=25, seed=12)
print(f' tail-3 best D={b3[0]:.4f}（来源 {b3[3]}） I={b3[2]:.6f}   可行 {int(np.isfinite(f3).sum())}/{len(f3)}')
print(f' [W3] 族内偏差参照：平面随机 best {r3:.4f} - 上界 {min(Ds,b3[0]):.4f} = {r3-min(Ds,b3[0]):+.4f}')
print(f' 投影种子 {Ds:.4f}  搜索 best {b3[0]:.4f}   平面 3 bit 上界 = {min(Ds, b3[0]):.4f}')
print(f' 合法区间：[{Dstar:.4f}, {min(Ds, b3[0]):.4f}]   宽度 {min(Ds, b3[0])-Dstar:.4f}')
print(f' [W3] 族内偏差参照：对角子族 vs 全平面 —— 见 E96c 第 (4) 节')
np.save('p0/e96d_T_seed_plane.npy', T_seed_pl if Tt_s is None else Tt_s)
np.save('p0/e96d_SNR_free.npy', SNR)
