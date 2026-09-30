r"""OpenQ1 / 最小舞台（1 阶单步版）：非二次代价 + u!=0 下，线性 vs 非线性编码器的最小信息率。

问题（单步"决策—虚拟传感器"舞台，即 Tanaka 2018 (TAC, arXiv:1510.04214) Step 2
虚拟传感器子问题在非二次代价下的版本）：

  源 X ~ N(0, S)（S=1）。控制器在收到编码 z 后做最优决策 a(z)=E[g(X)|z]，
  stage cost（非二次，关于 X 为三次量；u=a(z) 恒参与）：
      E[(g(X) - a)^2]，  g(x) = x + eps * x**3。
  编码器 z = phi(X) + v，v ~ N(0,V)；率 R = I(X;Z)（nats 或 bits）。

  (L) 线性编码器：phi(x)=c x。在率失真 (c,V) 上扫，得 D_lin(R)。
      闭式（Stein 引理 + 高斯条件方差）：
        D = 1 + 3 eps + 15 eps**2
            - (1 + 6 eps + 15 eps**2)**2 * c**2*S / (c**2*S + V)
        R = 0.5 log(1 + c**2*S/V)。
      （c 的标度自由度与 V 合并；无噪极限 c->inf 等价 V->0。）
  (N) 非线性编码器：phi(x)=x + q * x**3（Borel；q>=0 严格单调，含最优线性点 q=0），
      任意加性高斯噪声。率 = h(Z)-h(v)（数值熵，Gauss-Hermite），
      MMSE D = E Var(g(X)|Z)（Gauss Hermite 条件矩，或方差差恒等式）。
      另用 Monte Carlo 独立实现交叉核对。

判决：在同一 D 水平下 R_nl(D) 是否显著低于 R_lin(D)（> 数值噪声底，带收敛证书）。

用法：
  python -X utf8 -m p0.exp_openq1_1d
  python -X utf8 -m p0.exp_openq1_1d --eps 0.3 --mc 2000000
"""
import argparse
import numpy as np
import sys
from numpy.polynomial.hermite_e import hermegauss

sys.stdout.reconfigure(encoding="utf-8")
LN2 = np.log(2.0)


# ---------- 闭式：线性编码器 ----------
def lin_constants(eps):
    A1 = 1.0 + 6.0 * eps + 15.0 * eps ** 2          # E[X g(X)] / S
    D0 = 1.0 + 3.0 * eps + 15.0 * eps ** 2          # E[g(X)^2]
    return A1, D0


def lin_D_R(eps, ngrid=400):
    """返回 (D, R bits)，从无信息点 D0 到 0。"""
    A1, D0 = lin_constants(eps)
    s = np.linspace(0, 1, ngrid)                     # s = c^2 S/(c^2 S+V)
    D = D0 - A1 ** 2 * s
    R = -0.5 * np.log(1 - s) / LN2                   # 0.5 log2(1/(1-s))
    return D, R


def lin_R_at_D(eps, D):
    A1, D0 = lin_constants(eps)
    s = (D0 - D) / A1 ** 2
    if s <= 0:
        return 0.0
    return -0.5 * np.log(1 - s) / LN2


# ---------- Gauss-Hermite 数值实现：非线性编码器 ----------
class NLChannel:
    """z = x + q x^3 + v, x~N(0,S), v~N(0,V). Gauss-Hermite 求 h(Z) 与 MMSE。"""

    def __init__(self, eps, q, V, S=1.0, n=400):
        self.eps, self.q, self.V, self.S = eps, q, V, S
        self.z_nodes, self.w = hermegauss(n)         # N(0,1) 节点/权重
        self.w /= np.sqrt(2 * np.pi)
        self.x = np.sqrt(S) * self.z_nodes

    def _density_z(self, z):
        # p(z) = E_x N(z; phi(x), V)
        phi = self.x + self.q * self.x ** 3
        d = z[:, None] - phi[None, :]
        dens = np.exp(-0.5 * d ** 2 / self.V) / np.sqrt(2 * np.pi * self.V)
        return dens @ self.w

    def rate_bits(self, m=1000, zmax=None):
        phi = self.x + self.q * self.x ** 3
        sd = np.sqrt(np.var(phi) + self.V)
        if zmax is None:
            zmax = 8.0 * sd
        z = np.linspace(-zmax, zmax, m)
        p = self._density_z(z)
        p = np.maximum(p, 1e-300)
        # h(Z) = -∫ p log p
        h = -np.trapezoid(p * np.log(p), z)
        h_noise = 0.5 * np.log(2 * np.pi * np.e * self.V)
        return (h - h_noise) / LN2

    def distortion(self, k=120):
        """MMSE = E[Var(g(X)|Z)]。对每个 x 节点 i，取 Z|x_i 的 GH 节点 z_ij，
        在 z_ij 上求 E[g|z]（对 x 积分），平均后再对 x 积分。向量化。"""
        zk, ww = hermegauss(k)
        ww /= np.sqrt(2 * np.pi)
        gx = self.x + self.eps * self.x ** 3
        phi = self.x + self.q * self.x ** 3
        D = 0.0
        norm = 1.0 / np.sqrt(2 * np.pi * self.V)
        for i in range(len(self.x)):
            z = phi[i] + np.sqrt(self.V) * zk                    # (k,)
            d2 = (z[None, :] - phi[:, None]) ** 2               # (n, k)
            dens = np.exp(-0.5 * d2 / self.V) * norm            # N(z_j; phi(x_h), V)
            pz = dens.T @ self.w                                # (k,)
            Eg = dens.T @ (self.w * gx) / np.maximum(pz, 1e-300)
            Eg2 = dens.T @ (self.w * gx ** 2) / np.maximum(pz, 1e-300)
            cvar = Eg2 - Eg ** 2
            D += self.w[i] * (ww @ cvar)
        return D


def nl_montecarlo(eps, q, V, N=1_000_000, seed=0, S=1.0):
        """独立路径：MC 样本 + KDE 估 h(Z)；MMSE 用粒子近似条件矩。"""
        rng = np.random.default_rng(seed)
        x = rng.normal(0, np.sqrt(S), N)
        v = rng.normal(0, np.sqrt(V), N)
        gx = x + eps * x ** 3
        z = x + q * x ** 3 + v
        # MMSE via kNN-style local mean: bin-free nearest-neighbor regression on sorted z
        order = np.argsort(z)
        zs, gs = z[order], gx[order]
        k = 200
        # 滚动窗口 E[g|z]
        cs = np.convolve(gs, np.ones(k) / k, mode="same")
        cs2 = np.convolve(gs ** 2, np.ones(k) / k, mode="same")
        D = np.mean(cs2 - cs ** 2)
        # h(Z)：Gaussian KDE，Scott 带宽
        h_bw = 1.06 * np.std(z) * N ** (-0.2)
        # 子采样估计 h 以控内存
        idx = rng.choice(N, size=min(N, 200_000), replace=False)
        nb = 5000
        ent = 0.0
        for s in range(0, len(idx), nb):
            zb = z[idx[s:s + nb]][:, None]
            d = (zb - z[None, :200_000]) / h_bw
            p = np.mean(np.exp(-0.5 * d ** 2) / np.sqrt(2 * np.pi) / h_bw, axis=1)
            ent -= np.mean(np.log(np.maximum(p, 1e-300)))
        R = (ent - 0.5 * np.log(2 * np.pi * np.e * V)) / LN2
        return D, R


def run(eps=0.3, mc=1_000_000):
    A1, D0 = lin_constants(eps)
    print(f"eps={eps}  E[g^2]=D0={D0:.8f}  (E[Xg])^2={A1**2:.8f}")
    D_lin, R_lin = lin_D_R(eps)

    # 非线性编码器：扫描 V，固定 q；再对少数点优化 q
    Vs = np.concatenate([np.geomspace(4.0, 0.02, 14), np.geomspace(0.02, 1e-3, 6)])
    rows = []
    print(f"\n{'q':>6} {'V':>9} {'D':>12} {'R_nl bit':>10} {'R_lin@D':>10} {'gap bit':>9}")
    for q in [0.0, eps, 1.5 * eps, 2 * eps]:
        for V in Vs:
            ch = NLChannel(eps, q, V, n=300)
            D = ch.distortion(k=80)
            R = ch.rate_bits(m=700)
            Rl = lin_R_at_D(eps, D)
            rows.append((q, V, D, R, Rl, Rl - R))
            print(f"{q:6.3f} {V:9.4f} {D:12.6f} {R:10.5f} {Rl:10.5f} {Rl-R:9.5f}")

    # 在固定 D 水平附近挑最好的非线性（最大 gap），给 MC 交叉核对
    rows_t = [r for r in rows if r[5] > 0]
    rows_t.sort(key=lambda r: -r[5])
    best = rows_t[0]
    q, V = best[0], best[1]
    print(f"\n最佳非线性点 q={q}, V={V}: GH  D={best[2]:.6f} R_nl={best[3]:.5f} gap={best[5]:.5f}")
    Dmc, Rmc = nl_montecarlo(eps, q, V, N=mc)
    Rl_mc = lin_R_at_D(eps, Dmc)
    print(f"MC 独立核对: D={Dmc:.6f} R_nl={Rmc:.5f}  R_lin@D={Rl_mc:.5f} gap={Rl_mc-Rmc:.5f} (N={mc})")

    # 收敛证书：GH 阶数加倍
    ch1 = NLChannel(eps, q, V, n=150)
    ch2 = NLChannel(eps, q, V, n=450)
    d1 = (ch1.distortion(k=60), ch1.rate_bits(m=400))
    d2 = (ch2.distortion(k=120), ch2.rate_bits(m=1200))
    print(f"收敛证书 nGH=150/450: |D_D|={abs(d1[0]-d2[0]):.2e} |R_D|={abs(d1[1]-d2[1]):.2e}")

    # 在几个共同 D 水平上汇总（用 GH 扫描里 gap 最大曲线 q=eps 的最佳包络）
    print("\n率差汇总（GH 扫描；每 D 取最优 q）：")
    levels = [D0 - 0.15 * A1 ** 2, D0 - 0.4 * A1 ** 2, D0 - 0.7 * A1 ** 2]
    for Dtarget in levels:
        cand = [r for r in rows if abs(r[2] - Dtarget) < 0.25]
        if not cand:
            continue
        b = min(cand, key=lambda r: r[3])
        print(f"D={b[2]:8.4f}  R_lin={b[4]:8.5f}  R_nl={b[3]:8.5f}  gap={b[5]:8.5f} bit")
    return rows, best, (Dmc, Rmc, Rl_mc)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--eps", type=float, default=0.3)
    ap.add_argument("--mc", type=int, default=1_000_000)
    a = ap.parse_args()
    run(a.eps, a.mc)
