# -*- coding: utf-8 -*-
r"""实验 150：**反向水填的纯算术**能否解释我们正文的 $\gamma_r=2\ln2/r$？（自我降级用）

动机（全文实读 1912.07640 的 Thm 3 式 (50)：$R^{G}_{[0,n],in}(D-D^{min}_{[0,n]})=\frac12\sum_t\log(\sigma^2_{\xi_{t|t-1}}/\sigma^2_{\xi_{t|t}})$
$=$ "率写成**地板上方失真**的函数"，与我的 $\Delta D=D-D_{\rm floor}$ 对 $\Delta I$ 同形；以及 1810.00298 的 Thm 4 式 (33)）。
于是必须问：$\gamma_r=2\ln2/r$ 究竟是我们测出来的，还是水填算术一给就有的？

跑前写死的口径：
 [W0] 只用对角源 $\\lambda=(5,3,1.2,0.4)$ 与教科书反向水填 $\\delta_i=\\min(\\nu,\\lambda_i)$、
      $D(\\nu)=\\sum_i\\delta_i$、$R(\\nu)=\\frac12\\sum_i\\log_2(\\lambda_i/\\delta_i)$。**不碰 DARE、不碰我车道任何代码。**
 [W1] 对每个"活跃维数" $r$（$\\nu$ 小于最小的 $r$ 个 $\\lambda$ 才叫活跃），在 $\\nu$ 网格上回归 $\\ln D$ 对 $R$，
      斜率应等于 $-2\\ln2/r$；判据 [W1-a] 相对差 $\\le1\\%$ ⇒ "常数由算术给出"成立。
 [W2] 同时打印活跃维数随 $\\nu$ 下降的**阶次**（零率维被剔除），核对 "floor 上方" 的说法在水填侧也成立。
 [W3] 把我们**实测**的锚株 $\\Delta I=1.00$ 三元组 $1.860/1.420/1.281$（`e145_out.txt` 已落）除以 $2\\ln2/r$，
      印出超出的百分比 ⇒ 这才是没被算术解释掉的那部分。
"""
import sys
import numpy as np
sys.stdout.reconfigure(encoding='utf-8')

lam = np.array([5.0, 3.0, 1.2, 0.4])
ln2 = np.log(2.0)


def wf(nu):
    d = np.minimum(nu, lam)
    D = d.sum()
    R = 0.5 * np.sum(np.log2(np.where(d > 0, lam / d, 1.0)))
    return D, R, d


print('== [W0] 教科书反向水填：$\\delta_i=\\min(\\nu,\\lambda_i)$，$\\lambda=$', lam)
print('\n== [W1] 活跃维数 $r$ 上的 $\\ln D$ 对 $R$ 斜率 vs $-2\\ln2/r$ ==')
for r in (1, 2, 3, 4):
    # 让最后 r 个特征值活跃：nu 必须 < lam 的第 (p-r) 大值，且 > 0
    desc = np.sort(lam)[::-1]
    hi = desc[r - 1] * 0.999
    lo = desc[r] * 1.001 if r < len(lam) else 1e-5
    nus = np.geomspace(hi, lo, 400)
    Rs, Ds, Fs = [], [], []
    for nu in nus:
        D, R, d = wf(nu)
        act = d < lam - 1e-12
        if np.sum(act) != r:
            continue
        Rs.append(R); Ds.append(D); Fs.append(lam[~act].sum())
    Rs = np.array(Rs); Ds = np.array(Ds); Fs = np.array(Fs)
    A_ = np.vstack([Rs, np.ones_like(Rs)]).T
    # 关键：水填里的 D 含"零率维自带的方差"这一地板，所以只有 $D-D_{\rm floor}$ 才该拿来回归
    slope, inter = np.linalg.lstsq(A_, np.log(Ds - Fs), rcond=None)[0]
    slope_tot, _ = np.linalg.lstsq(A_, np.log(Ds), rcond=None)[0]
    pred = -2.0 * ln2 / r
    rel = abs(slope / pred - 1.0)
    print('  r=%d  点数 %3d  R 范围 [%.3f, %.3f]  地板 $D_{\\rm floor}$=%.3f'
          % (r, len(Rs), Rs.min(), Rs.max(), np.median(Fs)))
    print('        $\\ln(D-D_{\\rm floor})$ 对 $R$ 斜率 %+.6f  vs 预测 %+.6f  相对差 %.2e  $\\Rightarrow$ [W1-a] %s'
          % (slope, pred, rel, '算术即给出' if rel <= 0.01 else '不成立'))
    print('        对照：$\\ln D$（不扣地板）斜率 %+.6f $=$ 预测的 %.2f 倍 ⇒ 不扣地板必读错'
          % (slope_tot, slope_tot / pred))

print('\n== [W2] 活跃维数随 $\\nu$ 的阶次（零率维被剔除 ⇒ "地板上方"在水填侧同样成立） ==')
for nu in (6.0, 1.3, 0.9, 0.5, 0.2, 0.05):
    D, R, d = wf(nu)
    act = int(np.sum(d < lam - 1e-12))
    print('  $\\nu=%5.2f$  D=%.4f  R=%.4f bits  活跃维数 r=%d  剔除 %d 个零率维'
          % (nu, D, R, act, len(lam) - act))

print('\n== [W3] 实测局域斜率除以算术渐近值：没被解释掉的部分（锚株 $\\Delta I{=}1.00$，见 e145_out.txt） ==')
meas = {1: 1.860, 2: 1.420, 3: 1.281}
for r in (1, 2, 3):
    pred = 2.0 * ln2 / r
    print('  r=%d  实测 $\\gamma_{\\rm loc}=%.3f$  渐近 $2\\ln2/r=%.4f$  $\\Rightarrow$ 高出 %+.1f%%'
          % (r, meas[r], pred, 100.0 * (meas[r] / pred - 1.0)))
print('  结论：$\\gamma_\\infty=2\\ln2/r$ 这一**常数本身**由反向水填算术给出（[W1]），'
      '不属于我们的发现；属于我们的是有限率端高出 %s 这一段与"哪个 $Z$、开关阈值 $I^*(k)$、墙先理处的 $C(V)$"。'
      % ('/'.join('%+.0f%%' % (100.0 * (meas[r] / (2.0 * ln2 / r) - 1.0)) for r in (1, 2, 3))))
print('EXIT=0')
