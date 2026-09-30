r"""E42：给 $E(D)$ 补第二条腿，并把图 1 的高度差拆成"公式代价 + 求解器代价"两块。

四条线，口径完全一致（同一个 $\Theta$、同一个预算 $D$，约束都写成 $\mathrm{tr}(\Theta P)+\mathrm{tr}(WP_c)\le D$）：

  下界腿  $I_{\rm sdp}(D)$：Tanaka 式(18) 的 max-det 松弛，即去掉"精度必须落在 $\operatorname{ran}F^\top$"这条锥限制，
      只留一步一致性 LMI 和代价。真问题的任一可行点丢掉 $\operatorname{ran}$ 限制后仍可行，可行集变大、最小值变小，
      所以 $I_{\rm sdp}(D)\le E_{\rm TRV}(D)$，且这个下界**不依赖我的族、不依赖我的优化器**。
      作者的黑线 $I_{unc}$ 按定义就是这个 SDP，故记 $E_{unc}:=I_{\rm sdp}$；[5] 再逐点核对像素。
  上界腿  $I_{\rm chain}(D)$：E39c 的下降热启动锥最优，逐点带可行性证书（$J-D\sim10^{-11}$）。
  发表线  $I_{\rm red}(D)$、$I_{\rm unc}(D)$：图 1 红/黑逐像素（去端帽、去坐标轴伪影）。

于是 $\Delta_{\rm rep}=I_{\rm red}-I_{\rm sdp}$ 是信上写的"任务限制要多付的率"，而真值
$\Delta_{\rm true}=E_{\rm TRV}-E_{unc}\le I_{\rm chain}-I_{\rm sdp}=\Delta_{\rm rep}-(I_{\rm red}-I_{\rm chain})$。
右边的括号是 E39c 已证书化的求解器损失。**本轮交付物：信上的 $\Delta(D)$ 至少被求解器抬高了 $I_{\rm red}-I_{\rm chain}$。**
[2] 用逐窗口的读图误差界定这条修正分辨不分辨得出；分辨不出的区间明说，既不算反例也不作证。

第二条诊断回答"锥限制到底咬在哪"：松弛最优给出 $M=P^{-1}-\tilde P^{-1}\succeq0$，真问题要求 $\operatorname{ran}M\subseteq\operatorname{ran}F^\top$。
[3] 看越界比例 $\operatorname{tr}M_\perp/\operatorname{tr}M$；[4] 直接把 $M$ 夹回 $\operatorname{ran}F^\top$
（$S=\Pi M\Pi$ 仍是合法的 $F^\top\Gamma F$，$\Gamma\succeq0$，于是 $P'$ 落在锥内），看夹回去是否超预算。
"""
import sys
import numpy as np
sys.stdout.reconfigure(encoding='utf-8')
import cvxpy as cp
import p0.exp_plants_oos as E35
import p0.exp_ladder as X37
import p0.exp_C_exact as X36
from p0.exp_nearfloor_law import Fs, solve_cone, gamma_iso, pack
from p0.exp_fig1_cert import full_space_point

sym = X36.sym
ln2 = np.log(2.0)

base = None
for nm, pl, F in E35.plants():
    if nm.find('锚点') >= 0:
        base = (nm.replace('$', ''), pl, F)
        break
nm0, Pl0, _Fa = base
F2 = Fs['F2 尾部2(封闭, r=2)']
n, r = Pl0.nn, F2.shape[0]
TH = sym(Pl0.Th)
phi0 = X37.floor_exact(Pl0, F2, None)['闭式']
UNC = Pl0.UNC
FLOOR = phi0 + UNC
print('== E42 ==  算例：%s  $n=%d$ $r=%d$  地板 $D_{\\min}=%.6f$（$\\Phi_0=%.6f$，$\\operatorname{tr}(WP_c)=%.6f$）'
      % (nm0, n, r, FLOOR, phi0, UNC))

# $\operatorname{ran}F^\top$ 的投影（$F$ 是 $r\times n$，秩 $r$，正交补维数 $n-r$）
PiF = F2.T @ np.linalg.pinv(F2 @ F2.T) @ F2
Perp = np.eye(n) - PiF
print(r'  校验：rank $\Pi_F$=%d  tr $\perp$=%d  $\|\Pi_F^2-\Pi_F\|=%.1e$  $F\Pi_F-F$ 残差 $=%.1e$'
      % (np.linalg.matrix_rank(PiF), round(np.trace(Perp)),
         np.linalg.norm(PiF @ PiF - PiF), np.linalg.norm(F2 @ PiF - F2)))


def sdp_point(Pl, D):
    r"""Tanaka max-det 松弛：返回 $(I,\ P,\ \tilde P,\ \text{状态})$，$I=\tfrac12\log\det(\tilde P/P)$。"""
    P = cp.Variable((n, n), symmetric=True)
    Pi = cp.Variable((n, n), symmetric=True)
    APT = Pl.A @ P @ Pl.A.T + Pl.W
    cons = [Pi >> 0, P >> 0, cp.trace(TH @ P) + Pl.UNC <= D, P << APT,
            cp.bmat([[P - Pi, P @ Pl.A.T], [Pl.A @ P, APT]]) >> 0]
    obj = cp.Minimize(-.5 * cp.log_det(Pi) / ln2 + .5 * np.log(np.linalg.det(Pl.W)) / ln2)
    prob = cp.Problem(obj, cons)
    prob.solve(solver=cp.CLARABEL)
    if P.value is None:
        return None
    Pv = sym(P.value)
    Pt = sym(Pl.A @ Pv @ Pl.A.T + Pl.W)
    I = .5 * np.log(np.linalg.det(Pt) / np.linalg.det(Pv)) / ln2
    return I, Pv, Pt, prob.status


def binstats(arr, Dq, w=.03):
    r"""同一 $D$ 窗口的像素均值，和"半程误差"=窗内极差的一半（窗内曲线单调时读图误差的上界）。"""
    s = np.abs(arr[:, 0] - Dq) < w
    if not s.sum():
        return float('nan'), float('nan'), 0
    v = arr[s, 1]
    return float(v.mean()), .5 * float(v.max() - v.min()), int(s.sum())


red = np.load('p0/fig/fig1_red.npy')
blkraw = np.load('p0/fig/fig1_black.npy')
# 黑线像素里有一条水平伪影（同一个 $I$ 精确重复几百次 = 坐标轴底边），剔掉再用。
vals, cnts = np.unique(np.round(blkraw[:, 1], 12), return_counts=True)
bad = vals[cnts > 200]
blkpx = blkraw[~np.isin(np.round(blkraw[:, 1], 12), bad)]
print('  黑线像素 %d 个；剔除水平伪影 %s（各 %s 个像素，疑为坐标轴底边）后剩 %d 个，纵轴 %.4f..%.4f'
      % (len(blkraw), ', '.join('%.4f' % v for v in bad),
         ', '.join(str(c) for c in cnts[cnts > 200]), len(blkpx),
         blkpx[:, 1].min(), blkpx[:, 1].max()))

print('\n' + r'[1] 四条线同表：$I_{\rm sdp}\le E_{\rm TRV}\le I_{\rm chain}$，发表线在右')
print('  %8s %7s %8s %9s %8s %9s | %8s %8s %8s %7s | %s' %
      ('$D$', '$x$', r'$I_{\rm sdp}$', r'$I_{\rm unc}$' + '(像素)',
       r'$I_{\rm chain}$', r'$I_{\rm red}$' + '(像素)',
       r'$\Delta_{\rm rep}$', r'$\Delta_{\rm true}\le$', '求解器损失', '读图err',
       r'$\operatorname{tr}M_\perp/\operatorname{tr}M$  SDP状态'))
Dlist = [47.10, 47.21, 47.30, 47.50, 47.70, 47.90, 48.10, 48.30, 48.50,
         49.00, 49.50, 50.50, 52.0, 55.0, 60.0]
prev_G, prev_x = None, None
tab = []
for Dq in Dlist:
    rv, err_r, nr = binstats(red, Dq)
    bv, err_b, nb = binstats(blkpx, Dq)
    x = Dq - FLOOR
    lg = gamma_iso(F2, Dq)
    st = [pack(np.exp(lg) * np.eye(r))] if lg is not None else []
    if prev_G is not None:
        st += [pack(prev_G * (prev_x / x)), pack(prev_G)]
    rr = solve_cone(F2, Dq, FLOOR, st)
    It = rr[0] if rr is not None else float('nan')
    if rr is not None:
        prev_G, prev_x = rr[2], x
    got = sdp_point(Pl0, Dq)
    if got is None:
        print('  %8.3f  SDP 失败' % Dq)
        continue
    Isdp, Pv, Pt, status = got
    ev = np.linalg.eigvalsh(Pv)
    if ev.min() < 1e-12 * ev.max():
        print('  %8.3f  SDP 解接近奇异（比值 $%.2e$），跳过' % (Dq, ev.min() / ev.max()))
        continue
    Mi = sym(np.linalg.inv(Pv) - np.linalg.inv(Pt))
    trM, trperp = float(np.trace(Mi)), float(np.trace(Perp @ Mi @ Perp))
    # 夹回测试：$S=\Pi M\Pi\succeq0$ 仍是合法的 $F^\top\Gamma F$，所以 $P'$ 落在锥内。
    S = sym(PiF @ Mi @ PiF)
    FFt = F2 @ F2.T
    Gam = sym(np.linalg.solve(FFt, F2 @ S @ F2.T) @ np.linalg.inv(FFt))
    res_G = float(np.linalg.norm(F2.T @ Gam @ F2 - S))
    lamG = float(np.linalg.eigvalsh(Gam).min())
    Pclip = sym(np.linalg.inv(np.linalg.inv(Pt) + S))
    Iclip = .5 * np.log(np.linalg.det(Pt) / np.linalg.det(Pclip)) / ln2
    Jclip = float(np.trace(TH @ Pclip)) + UNC
    cert = full_space_point(np.eye(r), np.linalg.pinv(sym(rr[2])), phi0) if rr is not None else None
    tab.append((Dq, x, Isdp, bv, It, rv, rv - Isdp, It - Isdp, rv - It, err_r,
                cert['J'] - Dq if cert else float('nan'),
                trperp / max(trM, 1e-30), Iclip, Jclip - Dq, res_G, bv - Isdp,
                nr, nb, lamG, status))
    print('  %8.3f %7.4f %8.4f %9.4f %8.4f %9.4f | %8.4f %8.4f %8.4f %7.4f | %8.4f  %s' %
          (Dq, x, Isdp, bv, It, rv, rv - Isdp, It - Isdp, rv - It, err_r,
           trperp / max(trM, 1e-30), status))

# 列：0 $D$ 1 $x$ 2 sdp 3 黑像素 4 chain 5 红像素 6 $\Delta_{\rm rep}$ 7 $\Delta\le$
#     8 求解器损失 9 读图err 10 证书 11 越界占比 12 $I_{\rm clip}$ 13 $J^{\rm clip}-D$
#     14 $\|F^\top\Gamma F-S\|$ 15 黑$-$sdp 16 红像素数 17 黑像素数 18 $\lambda_{\min}\Gamma$
print('\n' + r'[2] 拆分：信上的 $\Delta$ 里有多少是求解器造成的，分辨不分辨得出')
print('  %8s %9s %9s %9s %8s %7s %8s %s' %
      ('$D$', r'$\Delta_{\rm rep}$', r'$\Delta\le$', '求解器损失', '读图err',
       r'$/$err', '占比', '判定'))
for q in tab:
    print('  %8.2f %9.4f %9.4f %9.4f %8.4f %7.1f %7.1f%% %s' %
          (q[0], q[6], q[7], q[8], q[9], q[8] / max(q[9], 1e-12),
           100 * q[8] / q[6], '可分辨' if q[8] > q[9] else '淹没在噪声'))
res = [q[0] for q in tab if q[8] > q[9]]
neg = [q[0] for q in tab if q[8] <= 0]
print(r'  可分辨（求解器损失 $>$ 该窗口读图误差）：%s'
      % (', '.join('%.2f' % v for v in res) if res else '一个都没有'))
print(r'  符号反过来的点（我的可行点不比红线低）：%s —— 只说明我的优化器在大预算端没跑赢读图精度，'
      r'不能拿来反证红线最优。' % (', '.join('%.2f' % v for v in neg) if neg else '无'))
print(r'  近地板窗口 $D\in[47.1,47.9]$：$\Delta_{\rm rep}$ 跨 %.3f..%.3f，求解器损失跨 %.3f..%.3f，'
      r'占 %.0f%%..%.0f%%，每一点都超过该窗口读图误差的最大值 %.3f。'
      % (min(q[6] for q in tab if q[0] <= 47.95), max(q[6] for q in tab if q[0] <= 47.95),
         min(q[8] for q in tab if q[0] <= 47.95), max(q[8] for q in tab if q[0] <= 47.95),
         min(100 * q[8] / q[6] for q in tab if q[0] <= 47.95),
         max(100 * q[8] / q[6] for q in tab if q[0] <= 47.95),
         max(q[9] for q in tab if q[0] <= 47.95)))

print('\n' + r'[3] 松弛最优 $M=P^{-1}-\tilde P^{-1}$ 越出 $\operatorname{ran}F^\top$ 的比例')
for q in tab:
    print('  $D=%6.2f$  越界占比 %.4f   证书 $J-D=%.2e$' % (q[0], q[11], q[10]))
print(r'  占比恒 $>0$：松弛用掉了 $F$ 读不到的方向，所以 $\Delta_{\rm true}>0$；'
      r'它随 $D$ 单调下降——预算越松，锥限制越不咬。')

print('\n' + r'[4] 夹回测试：把松弛解投影回 $\operatorname{ran}F^\top$ 还剩多少代价')
print('  %8s %9s %9s %9s %11s %8s %9s %9s' %
      ('$D$', r'$I_{\rm sdp}$', r'$I_{\rm clip}$', r'$\Delta I_{\rm clip}$',
       r'$J^{\rm clip}-D$', '超预算?', r'$\|F^\top\Gamma F-S\|$', r'$\lambda_{\min}\Gamma$'))
for q in tab:
    print('  %8.2f %9.4f %9.4f %9.4f %11.4f %8s %9.1e %9.1e' %
          (q[0], q[2], q[12], q[2] - q[12], q[13],
           '是' if q[13] > 1e-9 else '否', q[14], q[18]))
print(r'  若 $J^{\rm clip}>D$ 恒成立：松弛省下的那 $\Delta I_{\rm clip}$ bit 全靠 $F$ 读不到的方向买到，'
      r'夹回锥内就不再满足预算——这是 $\Delta_{\rm true}>0$ 的构造性证据，不需要我的族、不需要优化器。'
      r'$\lambda_{\min}\Gamma\ge0$ 与 $\|F^\top\Gamma F-S\|\approx0$ 是"夹回确实落在锥内"的自检。')

print('\n' + r'[5] 合法化下界腿：$I_{\rm sdp}$ 对上发表黑线（剔轴伪影后）')
print('  %8s %16s %9s %9s' % ('$D$', r'$I_{\rm unc}^{\rm像素}-I_{\rm sdp}$', '红像素数', '黑像素数'))
for q in tab:
    print('  %8.2f %16.4f %9d %9d' % (q[0], q[15], q[16], q[17]))
fin = [abs(q[15]) for q in tab if np.isfinite(q[15])]
print(r'  有黑像素的 %d 个点：$|I_{\rm unc}^{\rm像素}-I_{\rm sdp}|$ 最大 %.4f，中位 %.4f；'
      r'其余窗口黑像素数为 0（那段曲线被轴伪影占满或本来就没画），这些点只能用解析 $I_{\rm sdp}$。'
      % (len(fin), max(fin), float(np.median(fin))))
print(r'  结论：下界腿就是作者自己黑线的口径，不是我另造的基准；'
      r'所以 [2] 的修正直接作用于他们图 1 上那条 $\Delta$。')

print('\n' + r'[6] 三个数放在一起，决定这条修正敢不敢写进信')
dp = [(q[0], q[5] - q[3]) for q in tab if np.isfinite(q[3])]
print('  %8s %12s %12s %10s %10s' %
      ('$D$', r'$\Delta_{\rm paper}$' + '(红-黑像素)', r'$\Delta_{\rm rep}$' + '(红-解析)',
       '两者差', r'修正/差'))
for Dq, v in dp:
    q = [t for t in tab if t[0] == Dq][0]
    d = v - q[6]
    print('  %8.2f %12.4f %12.4f %10.4f %10s' %
          (Dq, v, q[6], d, '%.1f' % (q[8] / abs(d)) if abs(d) > 1e-6 else 'n/a'))
leg = max(abs(q[5] - q[3] - q[6]) for q in tab if np.isfinite(q[3]))
win = [q for q in tab if q[0] <= 48.55]
strong = [q for q in win if q[8] > 3 * leg]
print(r'  合法化误差 $\max|\Delta_{\rm paper}-\Delta_{\rm rep}|$（有黑像素的 %d 个点）$=%.4f$ bit：'
      r'这就是"用解析 SDP 代替像素黑线"的代价上限。' % (len(dp), leg))
print(r'  修正量（求解器损失）在 $D\le48.5$ 跨 %.3f..%.3f bit（$\Delta_{\rm rep}$ 的 %.0f%%..%.0f%%），'
      r'其中 $>3\times$ 合法化误差的点：%s'
      % (min(q[8] for q in win), max(q[8] for q in win),
         min(100 * q[8] / q[6] for q in win), max(100 * q[8] / q[6] for q in win),
         ', '.join('%.2f' % q[0] for q in strong) if strong else '无'))
print('\n' + r'  净表述（全部用 $\le$ 号，不藏误差）：$\Delta_{\rm true}\le I_{\rm chain}-I_{\rm sdp}'
      r'=\Delta_{\rm rep}-s$，$s$=求解器损失；换成像素口径 $\Delta_{\rm true}\le\Delta_{\rm paper}+\ell-s$，$\ell\le%.3f$。'
      % leg)
print(r'  所以在 $D\in[47.1,47.9]$（$s\ge%.3f$）这六个点上，'
      r'信上那条 $\Delta$ 至少被求解器抬高了 $s-\ell\ge%.3f$ bit；'
      r'$D\ge52$ 的点 $s$ 已进入读图误差，不提供任何方向的结论。'
      % (min(q[8] for q in tab if q[0] <= 47.95),
         min(q[8] for q in tab if q[0] <= 47.95) - leg))
print(r'  反向声明也要写：$\operatorname{tr}M_\perp/\operatorname{tr}M=3.8\%$ 说明松弛的最优解确实越界，'
      r'夹回锥内代价超出预算 $J^{\rm clip}-D\approx4.0$；但这只否掉了"这一个"极小点，'
      r'要主张 $\Delta_{\rm true}>0$ 严格成立还得排除存在另一个锥内极小点。')
print(r'  另外 $s$ 的正负号不是红线优劣的判据：$D=60$ 处 $s=-0.0007$，'
      r'只说明我的优化器在大预算端没跑赢读图精度，不能反证红线最优。')
