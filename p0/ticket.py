r"""ticket.py — 门票定理的可调用实现：输入 $(A,B,W,Q,R,F)$，输出每个"删一个传输方向"的分支、门槛与判据残差。

为什么要有这个文件：#21/#22 的全部结论都是三分法
    $\mathrm{gap}_d=+\infty \iff (A,F_{-d})$ 失去可检测性（PBH 秩检验）；
    $\mathrm{gap}_d=0 \iff \operatorname{ran}\Delta P_d\subseteq\ker\Theta$；
    否则 $\lambda_{\min}(\Theta)\operatorname{tr}(\Delta P_d)\le\mathrm{gap}_d\le\lambda_{\max}(\Theta)\operatorname{tr}(\Delta P_d)$，
而我在 #21 里亲手把一个 $\lambda_{\max}$ 的下标写错、制造了 22 个假违例。所以这里把三件事绑在一起做：
算门槛、验两界、报残差，任何一项不过就打印 `!!`，不让"只会亮绿的检查"再进第二次。

全程闭式：奇异 DARE（`solve_discrete_are(Aᵀ,Fᵀ,W,0)`）+ 定点抛光，没有 SDP、没有优化器。
两种删除口径都给：`C` 的本征基方向（精确，先解一次 $C$）和直接删 $F$ 的行（便宜，#21 §2 已证不可用，留着只为对照）。
"""
import sys
import numpy as np
sys.stdout.reconfigure(encoding='utf-8')
from scipy.linalg import solve_discrete_are
import p0.exp_plants_oos as E35
import p0.exp_ladder as X37
import p0.exp_C_exact as X36

sym = X36.sym
ORT = 1e-9


def floor0(Pl, F):
    r"""$(\Phi_0,P_0,\text{DARE 残差})$，$P_0=\lim_{V\to0}P(V)$。"""
    r = F.shape[0]
    Pt0 = sym(solve_discrete_are(Pl.A.T, F.T, Pl.W, np.zeros((r, r))))
    Pt0, Pk, res, _ = X37.polish_Pt0(Pl, F, Pt0)
    return float(np.trace(Pl.Th @ Pk)), Pk, res


def ker_theta(Pl):
    Th = sym(Pl.Th)
    wT, QT = np.linalg.eigh(Th)
    sup = wT > ORT * max(wT.max(), 1e-30)
    return wT, QT[:, ~sup], QT[:, sup], wT[sup]


def floor0_safe(Pl, F):
    """奇异 DARE 在无（可检测性意义的）有限解时抛 LinAlgError；普查里这不是崩溃而是数据。"""
    try:
        return floor0(Pl, F)
    except np.linalg.LinAlgError:
        return None


def branch(Pl, F, Feff, label, base=None):
    r"""一个删除动作的完整判定。返回 dict，键见 `_print`。

    两界这里给**两条**，因为 #21 只写了粗的那条、并且被我实现错了：
      粗界   $\lambda_{\min}(\Theta)\operatorname{tr}\Delta P\le\mathrm{gap}\le\lambda_{\max}(\Theta)\operatorname{tr}\Delta P$，
             $\lambda_{\min}$ 取 $\Theta$ 的**全体**特征值（奇异时就是 $0$，界变平凡但不假）。
      细界   在支撑投影 $\Pi$ 上：$\lambda_{\min}^+(\Theta)\operatorname{tr}(\Pi\Delta P\Pi)\le\mathrm{gap}\le\lambda_{\max}(\Theta)\operatorname{tr}(\Pi\Delta P\Pi)$，
             即"把 $\ker\Theta$ 藏起来的那部分迹"根本不该出现在门槛里 —— mech2 的 $0$ 门票正是 $\operatorname{tr}(\Pi\Delta P\Pi)=0$。
    """
    ok, lam = Pl.detectable(Feff)
    out = dict(label=label, detectable=ok, pbh=lam)
    if not ok:
        out.update(gap=float('inf'), branch='infeasible')
        return out
    if base is None:
        base = floor0_safe(Pl, F)
    if base is None:
        out.update(gap=float('inf'), branch='infeasible')
        return out
    phiF, PF, resF = base
    got = floor0_safe(Pl, Feff)
    if got is None:
        out.update(gap=float('inf'), branch='infeasible', phiF=phiF)
        return out
    phiD, PD, resD = got
    dP = sym(PD - PF)
    wT, ker, sup, wsup = ker_theta(Pl)
    Z = np.hstack([sup, ker])
    nb = sup.shape[1]
    Db = sym(Z.T @ dP @ Z)
    Dss = Db[:nb, :nb]
    Dsk = Db[:nb, nb:]
    Dkk = Db[nb:, nb:]
    # 定理的前提要验，不能默认：ΔP 必须真的 PSD（否则"零对角块 ⇒ 整行为零"不成立）
    psd_min = float(np.linalg.eigvalsh(dP).min())
    vis = float(np.trace(Dss))                     # 支撑上看得见的迹
    hid = float(np.trace(Dkk))                     # 藏进 ker Θ 的迹
    gap = phiD - phiF
    trdP = vis + hid
    exact = abs(gap - float(np.trace(np.diag(wsup) @ Dss)))
    crude_lo = max(wT.min(), 0.0) * trdP if len(wT) else 0.0
    crude_hi = (wT.max() if len(wT) else 0.0) * trdP
    ref_lo = (wsup.min() if len(wsup) else 0.0) * vis
    ref_hi = (wsup.max() if len(wsup) else 0.0) * vis
    tol = 1e-9 * max(abs(crude_hi), 1.0)
    out.update(phiF=phiF, phiD=phiD, resF=resF, resD=resD, gap=gap, trdP=trdP,
               psd_min=psd_min, vis=vis, hid=hid, exact=exact,
               cross=float(np.linalg.norm(Dsk)),
               rankdP=int(np.linalg.matrix_rank(dP, tol=1e-9)),
               lo=crude_lo, hi=crude_hi, rlo=ref_lo, rhi=ref_hi,
               inbnd=(crude_lo - tol <= gap <= crude_hi + tol
                      and ref_lo - tol <= gap <= ref_hi + tol),
               trcheck=abs(gap - float(np.trace(sym(Pl.Th) @ dP))),
               branch=('mech1' if abs(gap) <= 1e-12 and abs(trdP) <= 1e-12 else
                       ('mech2' if abs(gap) <= 1e-12 else 'active')))
    return out


def ticket(Pl, F, C=None, verbose=True):
    r"""全部删除方向（$C$ 本征基）+ 全部删行口径。门票 = 有限正 gap 的最小值。"""
    r = F.shape[0]
    wT, ker, sup, wsup = ker_theta(Pl)
    rows = []
    base = floor0_safe(Pl, F) if Pl.detectable(F)[0] else None
    if base is None:
        if verbose:
            print(r'  整株在任何口径下都没有有限基线（$F$ 自身不可检测，或奇异 DARE 无有限解）—— 门票无定义。')
        return dict(rows=[], entry=float('nan'), rankC=0, elam=np.array([]),
                    lamTh=wT, dimKer=ker.shape[1], n1=0, n2=0)
    if C is None:
        C = X37.C_of(Pl, F)
    rho, elam, U = X36.rank_C(C)
    for i in range(rho):
        d = U[:, i]
        wv, Vv = np.linalg.eigh(np.eye(r) - np.outer(d, d))
        Feff = Vv[:, wv > .5].T @ F
        rows.append(branch(Pl, F, Feff, r'$-u_%d$（$\lambda_C$=%.4g）' % (i, elam[i]), base=base))
    for i in range(r):
        rows.append(branch(Pl, F, np.delete(F, i, axis=0), '删行 %d' % i, base=base))
    fin = [q['gap'] for q in rows if q['branch'] in ('active', 'mech1', 'mech2')]
    entry = min([g for g in fin if g > 1e-12], default=(0.0 if fin else float('inf')))
    n1 = sum(1 for q in rows if q['branch'] == 'mech1')
    n2 = sum(1 for q in rows if q['branch'] == 'mech2')
    if verbose:
        hdr = '  rank $\\Theta$=%d  $\\dim\\ker\\Theta$=%d  $\\lambda(\\Theta)$=%s' % \
              (len(wsup), ker.shape[1], ', '.join('%.4g' % v for v in wT))
        print(hdr)
        print('  %-26s %9s %10s %10s %10s %10s %10s %10s %s' %
              ('方向', '分支', '$\\Phi_0(F)$', '$\\Phi_0(F_{-d})$', 'gap',
               r'$\operatorname{tr}\Delta P$', r'可见迹 $\operatorname{tr}\Pi\Delta P\Pi$',
               r'隐藏迹', '两界 / 残差'))
        for q in rows:
            if q['branch'] == 'infeasible':
                print('  %-26s %9s  PBH 在 $|\\lambda|=%.4f$ 失败 —— 全族在任何预算下不可行'
                      % (q['label'], q['branch'], abs(q['pbh'])))
                continue
            flag = 'ok' if q['inbnd'] else '!! 两界违例'
            print('  %-26s %9s %10.6f %10.6f %10.3e %10.3e %10.3e %10.3e  粗 %.2e..%.2e 细 %.2e..%.2e %s'
                  '  判据残差 $\\Delta_{\\rm 界}$=%.1e $\\lambda_{\\min}(\\Delta P)$=%.1e $\\Delta_{\\rm 交叉}$=%.1e'
                  % (q['label'], q['branch'], q['phiF'], q['phiD'], q['gap'], q['trdP'],
                     q['vis'], q['hid'], q['lo'], q['hi'], q['rlo'], q['rhi'], flag,
                     q['exact'], q['psd_min'], q['cross']))
        print('  门票（有限正 gap 的最小值）$=%.6e$   mech1（$\\Delta P=0$，真冗余）%d 个   '
              'mech2（$\\Delta P\\ne0$ 但全藏进 $\\ker\\Theta$，非冗余的免费删）%d 个'
              % (entry, n1, n2))
    return dict(rows=rows, entry=entry, rankC=rho, elam=elam, lamTh=wT,
                dimKer=ker.shape[1], n1=n1, n2=n2)


if __name__ == '__main__':
    print('== ticket.py 自检 ==')
    print('\n' + r'[1] 信上的锚点（$n=4$，$\operatorname{rank}\Theta=4$ 满秩）：#21 的两个门槛应逐位复现')
    base = None
    for nm, pl, F in E35.plants():
        if nm.find('锚点') >= 0:
            base = (nm.replace('$', ''), pl, F)
            break
    nm0, Pl0, F0 = base
    t0 = ticket(Pl0, F0)
    got = {q['label']: q['gap'] for q in t0['rows'] if q['branch'] != 'infeasible'}
    strong = [v for k, v in got.items() if k.find('u_0') >= 0]
    weak = [v for k, v in got.items() if k.find('u_1') >= 0]
    print('  复现检查：$-$u_0 门槛 %.6f（#21 记的是 74.282449）  $-$u_1 门槛 %.6f（#21 记的是 1.084392）'
          % (strong[0] if strong else float('nan'), weak[0] if weak else float('nan')))
    print('  满秩 $\\Theta$ 下不应出现 mech2 支：%s'
          % ('确认（没有 mech2）' if all(q['branch'] != 'mech2' for q in t0['rows']) else '!! 出现了 mech2'))

    print('\n' + r'[2] E43 的机制二植物（$n=3$，$\operatorname{rank}\Theta=1$）：应出现 mech2 支')

    def mech_plant(c13=0.0, au=.6, b3=0.0, wu=1.0):
        A = np.array([[1.25, .30, c13], [0., .85, 0.], [0., 0., au]])
        B = np.array([[1.], [0.], [b3]])
        W = np.diag([1.0, .7, wu])
        return E35.Pl(A, B, W, np.eye(3), np.array([[1.0]]))

    Fm2 = np.array([[1., 0., 0.], [0., 0., 1.]])
    Plm = mech_plant()
    t1 = ticket(Plm, Fm2)
    fr = [q for q in t1['rows'] if q['branch'] == 'mech2']
    print('  mech2 支个数=%d  门票（只数有限正 gap）$=%.6e$' % (len(fr), t1['entry']))
    print('  两种删除口径枚举的是**不同的动作集合**：$C$ 本征基给 %d 个方向，删行给 %d 个方向，'
          '本例没有同一个动作的两种写法，所以"两口径结论冲突"读不出来 —— 收回 #21 那句口号式断言，'
          '换成本文件 `[4]` 的可判定版本。' % (t1['rankC'], Fm2.shape[0]))

    print('\n' + r'[3] 断言（不过就抛错，别只打印）')
    for nm_, t_ in (('锚点', t0), ('机制二', t1)):
        for q in t_['rows']:
            if q['branch'] == 'infeasible':
                continue
            assert q['inbnd'], '%s：%s 两界违例（粗 %.3e..%.3e 细 %.3e..%.3e gap=%.3e）' \
                               % (nm_, q['label'], q['lo'], q['hi'], q['rlo'], q['rhi'], q['gap'])
            assert q['exact'] < 1e-9, '%s：%s gap 与支撑块闭式不一致' % (nm_, q['label'])
            assert q['trcheck'] < 1e-9, '%s：%s gap 与 tr(ΘΔP) 不一致' % (nm_, q['label'])
            assert q['psd_min'] > -1e-8, '%s：%s ΔP 不半正，定理前提失效' % (nm_, q['label'])
            if q['branch'] == 'mech2':
                assert abs(q['gap']) <= 1e-12 and q['trdP'] > 1e-6, '%s：mech2 判据自相矛盾' % nm_
                assert q['vis'] < 1e-9 and q['cross'] < 1e-9, '%s：mech2 却仍有可见分量' % nm_
                assert abs(q['hid'] - q['trdP']) < 1e-9, '%s：mech2 的迹没有全部藏进 ker Θ' % nm_
    print('  全部通过：两界（粗+细）、两条独立复核、ΔP 半正、mech2 的三条子判据都成立。')

    print('\n' + r'[4] 这轮真正的收获：粗界实现错了，细界才是能用的那条')
    print(r'  #21 写的区间里 $\lambda_{\min}$ 是**全体**特征值，$\Theta$ 奇异时它 $=0$，')
    print(r'  下界退化成 $0\le\mathrm{gap}$：对 mech2 既不错也不狠。我第一版把下界写成')
    print(r'  "最小正特征值乘 $\operatorname{tr}\Delta P$"，于是 mech2 那行必然报违例——')
    print(r'  违的是我的实现，不是定理。换成正交分块把 $\Delta P$ 拆成可见/隐藏两块之后：')
    for nm_, t_ in (('锚点', t0), ('机制二', t1)):
        for q in t_['rows']:
            if q['branch'] == 'infeasible':
                continue
            print(r'    %-8s %-26s tr dP=%.6e  可见 %.6e  隐藏 %.6e  粗 [%.3e,%.3e] 细 [%.3e,%.3e]  gap=%.6e'
                  % (nm_, q['label'], q['trdP'], q['vis'], q['hid'],
                     q['lo'], q['hi'], q['rlo'], q['rhi'], q['gap']))
    full = [q for q in t0['rows'] if q['branch'] != 'infeasible']
    worst = max(abs(q['rlo'] - q['lo']) / max(abs(q['hi']), 1e-30)
                + abs(q['rhi'] - q['hi']) / max(abs(q['hi']), 1e-30) for q in full)
    print(r'  满秩时两界重合（锚点最大相对差 %.2e），奇异时只有细界还能定位 gap：'
          % worst)
    print(r'  mech2 的 0 门票被细界读成"可见迹为 0"，粗界只会读成"落在 $[0,\lambda_{\max}\operatorname{tr}\Delta P]$ 里"。')
    print(r'  （机制二植物的 $\lambda_{\min}(\Theta)$ 数值上是 $-5.6e^{-17}$，是 $\Theta\succeq0$ 的舍入，粗界下界按 $0$ 截断。）')
    print(r'  可直接用的推论（本轮新增）：门槛只对可见块 $\Pi\Delta P\Pi$ 收费。')
    print(r'  于是任何拿 $\operatorname{tr}\Delta P$ 或 $\|\Delta P\|$ 给删除动作排序的启发式，'
          r'在 $\dim\ker\Theta>0$ 的植物上会把免费动作排到贵动作之后，')
    print(r'  而且排序误差可以任意大——E43 的 $a_u$ 扫描里 $\operatorname{tr}\Delta P$ 从 $1.04$ 走到 $50.3$，门票恒为 $0$。')
    print('\n' + r'[5] 全植物普查：E35.plants() 的每一株、两种口径的每一个删除动作')
    print(r'  这是 #22 §4 那条猜测的可判定形式："免费删是不是 generic" —— 现在有数了。')
    print(r'  区分两支：mech1 是 $\Delta P=0$（这条传输本来就被别的通道覆盖，删了等于没删）；')
    print(r'  mech2 才是非平凡的那支——后验确实变了，但变的方向全在 $\ker\Theta$ 里，控制器不为此收费。')
    print('  %-30s %3s %3s %5s %5s %6s %6s %6s %8s %10s %10s' %
          ('植物', 'n', 'r', r'$\dim\ker$', '动作数', '不可行', 'mech1', 'mech2', '门票',
           '正gap最小', '正gap最大'))
    zoo = list(E35.plants())
    nrows_all = ninf = n1a = n2a = 0
    pl1 = pl2 = pl0 = 0
    for nm, pl, Fq in zoo:
        tq = ticket(pl, Fq, verbose=False)
        cnt = dict((k, sum(1 for q in tq['rows'] if q['branch'] == k))
                   for k in ('infeasible', 'mech1', 'mech2', 'active'))
        pos = [q['gap'] for q in tq['rows'] if q['branch'] == 'active']
        nrows_all += len(tq['rows'])
        ninf += cnt['infeasible']
        n1a += cnt['mech1']
        n2a += cnt['mech2']
        pl1 += 1 if cnt['mech1'] else 0
        pl2 += 1 if cnt['mech2'] else 0
        pl0 += 1 if tq['entry'] == 0.0 else 0
        print('  %-30s %3d %3d %5d %5d %6d %6d %6d %8.4f %s %s'
              % (nm.replace('$', '')[:30], pl.A.shape[0], Fq.shape[0], tq['dimKer'],
                 len(tq['rows']), cnt['infeasible'], cnt['mech1'], cnt['mech2'], tq['entry'],
                 '%.2e' % min(pos) if pos else '   —  ', '%.2e' % max(pos) if pos else '   —  '))
        for q in tq['rows']:
            if q['branch'] == 'infeasible':
                continue
            assert q['inbnd'], '%s：%s 两界违例' % (nm[:12], q['label'])
            assert q['exact'] < 1e-8 and q['trcheck'] < 1e-8, '%s：%s 复核不过' % (nm[:12], q['label'])
            assert q['psd_min'] > -1e-7, '%s：%s ΔP 不半正' % (nm[:12], q['label'])
            if q['branch'] == 'mech2':
                assert q['vis'] < 1e-9 and q['trdP'] > 1e-6, '%s：%s mech2 子判据不成立' % (nm[:12], q['label'])
    print(r'  普查规模 %d 个删除动作 / %d 株植物：不可行 %d、mech1 %d、mech2 %d。'
          % (nrows_all, len(zoo), ninf, n1a, n2a))
    print(r'  至少出现一个 mech1 的植物 %d 株，出现 mech2 的 %d 株，门票为 $0$ 的 %d 株。'
          % (pl1, pl2, pl0))
    print(r'  读法：#22 说"机制二是构造出来的、随机扫不到"。' + r'在信上与我的旧扫描里造的那批植物上，'
          r'mech2 是否为 $0$ 就是这条猜测的判决。')
    print(r'  如果 mech2 恒为 $0$ 而 mech1 常见，那"免费删"在日常植物上其实只有"本来就没传"这一种平凡形态，')
    print(r'  我的 E43 就该被读成**边界情形的手术刀**而不是普遍现象——这一点必须由合作者推翻或确认。')
    print(r'  两界断言在全普查上没有一处违例：细界没有反例，粗界在 $\dim\ker\Theta>0$ 处只是松弛。')
