r"""E99 = 把 1912.03799 的发表证书 α̂（Theorem 5 式 (35)）搬到我们的锚点上算成一个数。

为什么值得搬：§64-G 留的问题"Ye 的『常数因子近似都不存在』与我们的连续锥之间要不要插一层松弛"，
在 1912.03799 里**文献自己已经插过了**——
 Remark 3: "it is impossible to obtain a non-trivial, universal (system-independent) bound on α̂ unless P = NP"
            （其 [16] = **Ye–Roy–Sundaram, ACC 2018, pp. 5049–5054**，即 2003.11951 的会议前身），
 并且 "α̂ → 0 as σ²_v → 0，即无噪声输出时基于近似超模性的保证**变得空泛**"。
 ⇒ 系统相关的近优证书（他们的 Th.5）与系统无关的不可近似（Ye）并不冲突：硬实例住在 σ_v→0 那一侧。
 ⇒ 这条共存是**已发表的**，我们不许把它当自己的观察卖。

本脚本要一个数（不是复现他们的定理）：在锚点植物 + 标准基候选传感器 + $V=I$ 的口径下，
$\hat\alpha=\lambda_{\min}(\Pi)/\lambda_{\max}(\Pi+\sum_u V_u)$ 到底是多大，
以及它隐含的证书比 $(1-e^{-\hat\alpha r/s})^{-1}$ 与实测贪心差距（C 的 1.08–1.20×）差多少个量级。

== 判据（跑前写死）==
 [A1] 噪声比读数：$\sigma^2_v/\sigma^2_w$（$R=I$ 对 $W$ 的本征值）。他们的 Fig.2 非平凡区是 $>10$。
      锚点若落在该比值之下 ⇒ 他们的证书在我们工作点**本来就不该期待有用**，这是定位句不是失败。
 [A2] 若 $\hat\alpha$ 隐含的证书比 $>5$，正文/注记**不许**写"greedy 的近优性由 [新引] 保证"，
      只能写"发表证书在锚点上为 $X\times$，实测 $1.08$–$1.20\times$ 比它紧"。
 [A3] $\Pi$ 的口径必须写死并用两版对照（全信息基线 / 去掉一根传感器），防止"挑一个好看的 $\Pi$"。
"""
import sys
sys.stdout.reconfigure(encoding='utf-8')
import numpy as np
from scipy.linalg import solve_discrete_are, schur
sys.stdout.reconfigure(encoding='utf-8')
np.set_printoptions(precision=5, suppress=True, linewidth=170)

_src = open('p0/exp_c_audit.py', encoding='utf-8').read().split("print(r'== E53")[0]
_ns = {'__name__': 'p'}
exec(compile(_src, 'p0/exp_c_audit.py[preamble]', 'exec'), _ns)
A, W, TH, JC, n, sym = _ns['A'], _ns['W'], _ns['TH'], _ns['JC'], _ns['n'], _ns['sym']
ln2 = np.log(2.0)

wv = np.linalg.eigvalsh(sym(W))
print('== (0) [A1] 噪声比：他们的非平凡区是 σ²_v > 10 σ²_w ==')
print(f' λ(W) = {wv}   σ_w = sqrt λ ∈ {np.sqrt(wv)}')
print(f' R = I ⇒ σ_v = 1。逐向比值 σ²_v/σ²_w = {1.0/wv}')
print(f' ⇒ [A1] 锚点落在 {"非平凡区（>10）" if np.all(1.0/wv > 10) else "混合区：至少一个方向 ≤10"}'
      f'，最小比值 {np.min(1.0/wv):.2f}')

# 候选传感器 O = 标准基输出，V = I（与 C 的 Table 1 的"标准基"候选同口径）
Vsum = np.zeros((n, n))
for i in range(n):
    Cu = np.zeros((1, n)); Cu[0, i] = 1.0
    Vsum += Cu.T @ Cu          # V=I


def alpha_hat(Cbase):
    """Π = 基线传感器集的稳态**预测**信息矩阵；α̂ = λmin(Π)/λmax(Π+Σ_u V_u)。"""
    r = Cbase.shape[0]
    Pm = sym(solve_discrete_are(A.T, Cbase.T, W, np.eye(r)))
    Pi = np.linalg.inv(Pm)
    num = float(np.linalg.eigvalsh(sym(Pi))[0])
    den = float(np.linalg.eigvalsh(sym(Pi + Vsum))[-1])
    return num / den, Pi, Pm


print('\n== (1) [A3] 两版 Π 口径下的 α̂ 与隐含证书 ==')
bases = [('全 4 根（信息最全 ⇒ Π 最大）', np.eye(n)),
         ('去掉第 1 行', np.eye(n)[1:]),
         ('去掉第 1、2 行', np.eye(n)[2:]),
         ('只留第 1 行', np.eye(n)[:1])]
rows = []
for lab, Cb in bases:
    try:
        a, Pi, Pm = alpha_hat(Cb)
    except Exception as e:
        print(f' {lab:<28} DARE 失败（不可检测）：{type(e).__name__}')
        continue
    for rs_lab, frac in [('r/s=1（预算刚好）', 1.0), ('r/s=0.5', 0.5)]:
        cert = 1.0 / (1.0 - np.exp(-a * frac)) if a * frac > 1e-12 else np.inf
        rows.append((lab, a, rs_lab, cert))
    cert1 = 1.0 / (1.0 - np.exp(-a)) if a > 1e-12 else np.inf
    print(f' {lab:<28} α̂={a:.5f}   证书比(r/s=1) {cert1:8.3f}×   '
          f'(r/s=0.5) {1.0/(1.0-np.exp(-a*0.5)):8.3f}×   '
          f'λmin(Π)={np.linalg.eigvalsh(sym(Pi))[0]:.4f} λmax(Π+ΣV)={np.linalg.eigvalsh(sym(Pi+Vsum))[-1]:.4f}')

print('\n== (2) 与实测贪心差距对照（C 的地板侧规则：距最优 1.08–1.20×）==')
worst = min(r[3] for r in rows)
print(f' 发表证书能给出的最好（最小）比 = {worst:.3f}×；实测 1.08–1.20×')
print(f' [A2] 判据门槛 5×：{"证书弱到不可用，正文只能按 A2 的措辞写" if worst > 5 else "证书本身已够紧"}')
print(f' 差距量级：worst 证书 / 实测 = {worst/1.20:.2f}×')

print('\n== (3) 忠实探针：把**所有**传感器的测量噪声 σ²_v 一起送向 0（基线与候选同尺度）==')
print('  前一版这里只缩放基线的 Π、不缩放候选的 $V_u$，那不是他们的极限 ⇒ 作废并改正（板 §66 自首）。')


def alpha_vs_sigma_v(sv2, Cbase=None):
    """基线传感器与候选传感器共用同一测量噪声 σ²_v：
    Π = Pm⁻¹，$V_u=C_u^t(σ²_v I)^{-1}C_u$。"""
    Cb = np.eye(n) if Cbase is None else Cbase
    r = Cb.shape[0]
    Pm = sym(solve_discrete_are(A.T, Cb.T, W, sv2 * np.eye(r)))
    Pi = np.linalg.inv(Pm)
    Vs = np.eye(n) / sv2            # Σ_u e_u e_u^t / σ²_v = I/σ²_v
    num = float(np.linalg.eigvalsh(sym(Pi))[0])
    den = float(np.linalg.eigvalsh(sym(Pi + Vs))[-1])
    return num / den


for sv2 in (1.0, 1e-1, 1e-2, 1e-4, 1e-6, 1e-8, 1e-12):
    a = alpha_vs_sigma_v(sv2)
    cert = 1.0 / (1.0 - np.exp(-a)) if a > 1e-12 else np.inf
    print(f'  σ²_v={sv2:9.0e}  α̂={a:.8f}   证书比(r/s=1)={cert:12.5f}×')
print(' ⇒ 判决看数字：若 α̂ 随 σ²_v→0 落到 0，则在本锚点上支持他们 Remark 3 的机制（硬实例住在噪声less 侧）；')
print('   若不落 0，则**不许**引用为"复现"，只报实测。Fig.2 的非平凡区是 σ²_v/σ²_w>10，锚点最大只有 3.37。')

