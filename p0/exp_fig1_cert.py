r"""E39c：把"原文图 1 的红曲线不是最优"变成**可审计的证书**，并定位弱通道何时被关掉。

E39b 的读数：在 $D\in[47.2,48.5]$ 上，逐像素红线比我的解析可行点高 $0.08\sim0.13$ bit，
而 $D\ge51$ 之后两者差到 $<0.005$ bit。这句话要成立，必须回答两个问题：
 (1) **可行性**。$E(D)=\min_A I_A(x)$ 用的是降阶传感器 $F_A=U_A^\top F$，也就是"只传输 $C$ 的一个方向"，
     精度矩阵 $\Gamma=U_A V_A^{-1}U_A^\top$ 是**奇异**的。它是不是原锥问题（传感器 $F$、精度 $\Gamma\succeq0$）
     的合法设计点？回答是：$P=(\tilde P^{-1}+F^\top\Gamma F)^{-1}$ 与 $\Gamma$ 的秩无关，而
     $I_{dir}=\tfrac12\log\det(I+\Gamma^{1/2}F\tilde PF^\top\Gamma^{1/2})$ 在零精度方向上贡献 0，
     正好等于降阶公式。这一轮把这个等式在**全空间**上重算一遍（先验 DARE 残差 + 代价 + 率），
     不看降阶版本。
 (2) **轴标定**。同一张图的黑/品红曲线已经对到 $\pm0.02$ bit（E39a 的 [1]），所以 $0.1$ bit 的差
     不是像素误差；但红线最左 3 列是折线的端点帽（E39b 的 [4]），必须排除。
顺带：数值最优的 $V^\ast$ 在 $C$ 本征基下的两个方差，谁在发散——这决定"$D\gtrsim48$ 之后 F2 其实只用一个通道"
这句话能不能说成是被证明的。
"""
import sys
import numpy as np
sys.stdout.reconfigure(encoding='utf-8')
from scipy.linalg import solve_discrete_are
import p0.p0_replicate_letter as base
import p0.exp_plants_oos as E35
import p0.exp_ladder as X37
from p0.exp_nearfloor_law import Fs, solve_cone, gamma_iso, pack, unpack
from p0.p0_replicate_letter import n

ln2 = np.log(2.0)
sym = X36sym = X37.sym
F2 = Fs['F2 尾部2(封闭, r=2)']
Pl = E35.Pl(base.A, base.B, base.W, np.eye(4))
Th, Pc, UNC = Pl.Th, Pl.Pc, Pl.UNC


def full_space_point(U_A, V_A, phi0):
    r"""在全空间上重算一个秩 deficient 设计点：$\Gamma=U_A V_A^{-1}U_A^\top$，传感器仍是 $F$。"""
    r = F2.shape[0]
    Gam = sym(U_A @ np.linalg.pinv(V_A) @ U_A.T)
    Feff = U_A.T @ F2
    Pt = sym(solve_discrete_are(Pl.A.T, Feff.T, Pl.W, V_A))
    for k in range(4000):
        Pk = sym(Pt - Pt @ Feff.T @ np.linalg.solve(Feff @ Pt @ Feff.T + V_A, Feff @ Pt))
        new = sym(Pl.A @ Pk @ Pl.A.T + Pl.W)
        d = np.linalg.norm(new - Pt)
        Pt = new
        if k > 5 and d <= 3e-17 * max(np.linalg.norm(Pt), 1e-30):
            break
    # 全空间核验：后验用奇异 $\Gamma$ 的信息形式，先验残差用 $F$（不是 $F_A$）
    S = Pt
    Pfull = sym(np.linalg.inv(np.linalg.inv(S) + F2.T @ Gam @ F2))
    res = float(np.linalg.norm(S - (Pl.A @ Pfull @ Pl.A.T + Pl.W)) / np.linalg.norm(S))
    J = float(np.trace(Th @ Pfull)) + UNC
    w, Q = np.linalg.eigh(Gam)
    sup = w > 1e-12 * w.max()
    Qs, ws = Q[:, sup], w[sup]
    M = np.eye(int(sup.sum())) + np.sqrt(ws) * (Qs.T @ (F2 @ S @ F2.T) @ Qs) * np.sqrt(ws)
    I_full = .5 * np.log(np.linalg.det(M)) / ln2
    return dict(J=J, I=I_full, res=res, P=Pfull, S=S, rank=int(sup.sum()),
                cond_Gam=float(w.max() / w[sup].min()))


if __name__ == '__main__':
    fx = X37.floor_exact(Pl, F2, None)
    phi0 = fx['闭式']
    C = X37.C_of(Pl, F2)
    fams, e, U, rho = X37.fams_of(Pl, F2, C, phi0)
    FLOOR = phi0 + UNC
    ffull = [f for f in fams if f['a'] == rho][0]
    f0 = [f for f in fams if f['name'] == '0'][0]
    print(r'  地板 $D_{\min}=%.6f$，$\lambda(C)=%s$，$\mathrm{gap}_{\{0\}}=%.4f$'
          % (FLOOR, np.array2string(e[:rho], precision=4), f0['gap']))

    # ---------- (1) 族切换点精扫 ----------
    print(r'\n===== [1] 包络的主动集切换：$\{0,1\}$ 射线 vs $\{0\}$ 射线（$E$ 口径） =====')

    def ray(fm, x):
        xA = x - fm['gap']
        if xA <= 30 * fm['res']:
            return None
        lo, hi = -15.0, 3.0
        for _ in range(46):
            mid = .5 * (lo + hi)
            try:
                xF = X37.pv(Pl, fm['F'], 10 ** mid * np.diag(1.0 / fm['w']), phi0)[0]
            except Exception:
                return None
            if xF - fm['gap'] < xA:
                lo = mid
            else:
                hi = mid
        t = 10 ** (.5 * (lo + hi))
        p = X37.pv(Pl, fm['F'], t * np.diag(1.0 / fm['w']), phi0)
        return p[1], t * np.diag(1.0 / fm['w'])

    last, sw = None, []
    for x in np.arange(1.60, 2.61, .02):
        a_, b_ = ray(ffull, x), ray(f0, x)
        if a_ is None or b_ is None:
            continue
        who = '0+1' if a_[0] <= b_[0] else '0'
        if last is not None and who != last:
            sw.append(x)
        if abs(x - round(x, 1)) < 1e-9 or who != last:
            print('    $x=%5.2f$  $D=%7.3f$  全支撑 %.4f  单通道 %.4f  差 %+.4f  胜 %s'
                  % (x, x + FLOOR, a_[0], b_[0], a_[0] - b_[0], who))
        last = who
    print(r'  切换点 $x^*=%s$  →  $D^*=%s$（而 $\mathrm{gap}_{\{0\}}=%.4f$ 只是 $\{0\}$ 族的**可行门槛**）'
          % (r'$\approx%.2f$' % sw[0] if sw else '未测到',
             r'$\approx%.2f$' % (sw[0] + FLOOR) if sw else '—', f0['gap']))

    # ---------- (2) 全空间可行性证书 ----------
    print('\n===== [2] 秩 deficient 设计点的全空间证书（不看降阶版本） =====')
    print('  %8s %10s %10s %10s %10s %8s %10s' %
          ('D', 'x', '$I_{全空间}$', '$J$', '$J-D$', 'rank$\\Gamma$', 'DARE 残差'))
    for x in (1.0, f0['gap'] + .02, 1.5, 2.0, 2.3769, 3.3769):
        for fm, nm in ((ffull, '0+1'), (f0, '0')):
            rr = ray(fm, x)
            if rr is None:
                continue
            Iv, Va = rr
            UA = U[:, fm['A']]
            pt = full_space_point(UA, Va, phi0)
            print('  %8.4f %10.4f %10.4f %10.4f %+10.2e %8d %10.1e   （族 %s，降阶率 %.4f，代价差 %+.2e）'
                  % (x + FLOOR, x, pt['I'], pt['J'], pt['J'] - (x + FLOOR), pt['rank'], pt['res'],
                     nm, Iv, pt['J'] - FLOOR - x))

    # ---------- (3) 红曲线 vs 证书：谁高谁低 ----------
    print('\n===== [3] 逐像素红线（去端帽）对可行点上界（锥最优用下降热启动链） =====')
    red = np.load('p0/fig/fig1_red.npy')
    capD = red[:, 0].min() + 3 * 0.0277
    Dlist = sorted([47.10, 47.21, 47.30, 47.50, 47.70, 47.90, 48.10, 48.30, 48.50,
                    49.00, 49.50, 50.50, 52.0, 55.0, 60.0], reverse=True)
    rows = []
    prev_G, prev_x = None, None
    for Dq in Dlist:
        if Dq <= capD:
            continue
        sel = np.abs(red[:, 0] - Dq) < .03
        if not sel.sum():
            continue
        rv = red[sel, 1].mean()
        x = Dq - FLOOR
        cand = []
        for fm in fams:
            q = ray(fm, x)
            if q is not None:
                cand.append((q[0], fm['name']))
        if not cand:
            continue
        Ev, nm = min(cand)
        lg = gamma_iso(F2, Dq)
        st = [pack(np.exp(lg) * np.eye(2))] if lg is not None else []
        if prev_G is not None:
            st += [pack(prev_G * (prev_x / x)), pack(prev_G)]
        r_ = solve_cone(F2, Dq, FLOOR, st)
        It = float('nan')
        cert = None
        if r_ is not None:
            It, prev_G, prev_x = r_[0], r_[2], x
            Vq = np.linalg.pinv(sym(r_[2]))
            cert = full_space_point(np.eye(F2.shape[0]), Vq, phi0)
        rows.append((Dq, x, rv, Ev, It, nm, cert))
    print('  %8s %8s %10s %10s %10s %10s %10s %10s %8s %s' %
          ('D', 'x', '红像素', 'E(D)解析', '红-E', '锥链最优', '红-锥', '$J-D$', 'DARE残差', '$I_{审计}$'))
    for Dq, x, rv, Ev, It, nm, cert in sorted(rows):
        print('  %8.3f %8.4f %10.4f %10.4f %+10.4f %10.4f %+10.4f %10.2e %8.1e %10.4f  最优族 %s'
              % (Dq, x, rv, Ev, rv - Ev, It, rv - It,
                 cert['J'] - Dq if cert else float('nan'),
                 cert['res'] if cert else float('nan'),
                 cert['I'] if cert else float('nan'), nm))

    # ---------- (4) $V^\ast$ 沿 $C$ 本征基的发散方向 ----------
    print('\n===== [4] 数值最优的后验噪声方差在 $C$ 本征基下：谁被关掉？ =====')
    print(r'  %8s %14s %14s %10s %14s %10s' %
          ('D', r'$v_{强}(\lambda_C=%.3f)$' % e[0], r'$v_{弱}(\lambda_C=%.3f)$' % e[1],
           r'弱/强', r'$\cos$(发散向-强)', '非对角/几何均值'))
    prev_G, prev_x = None, None
    for Dq in [47.0, 47.5, 48.0, 48.5, 49.5, 51.0, 55.0, 60.0]:
        x = Dq - FLOOR
        lg = gamma_iso(F2, Dq)
        st = [pack(np.exp(lg) * np.eye(2))] if lg is not None else []
        if prev_G is not None:
            st += [pack(prev_G * (prev_x / x))]
        r_ = solve_cone(F2, Dq, FLOOR, st)
        if r_ is None:
            continue
        prev_G, prev_x = r_[2], x
        V = np.linalg.pinv(sym(r_[2]))
        Vc = U[:, :rho].T @ V @ U[:, :rho]
        wq, Aq = np.linalg.eigh(sym(Vc))
        ang = abs(Aq[:, -1] @ np.array([1., 0.]))        # 最大方差方向与强 $C$ 方向的夹角余弦
        print('  %8.2f %14.4g %14.4g %10.3g %14.4f %10.2g'
              % (Dq, Vc[0, 0], Vc[1, 1], Vc[1, 1] / max(Vc[0, 0], 1e-30), ang,
                 Vc[0, 1] / np.sqrt(Vc[0, 0] * Vc[1, 1])))
