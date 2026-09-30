r"""E59b = 两件判定实验，都是冲 C 的 §10 去的。

[A] §10.1 的头号主张"24/24 全等于 10.169224"：我的球面投影梯度（E59 [4]）从 24 个 Haar 起点
    **一个都没落到 10.169224**（最小 11.4886，中位远在上方），且只迭代 5–16 步就停。两种可能：
      (i) 我的线搜索/步长停滞（容忍度伪影）；(ii) $r=1$ 的 $\Phi$ 在 $\mathbb{P}^3$ 上确有多个局部极小。
    区分方法：对每个停滞点 (a) 用更细的差分重算梯度模（看是否真驻定），(b) 在 3 个切向参数上做
    Powell 精修（我的 E53 就是靠 Powell 命中 10.169224 的，所以这条路线已被我自己验证过），
    (c) 从停滞点沿随机切向小步走出去再做局部下降。若 Powell 将全部 24 个拉回 10.169224 ⇒ 停滞是我的
    伪影，他的主张保留；若留下 $\Phi>10.169224$ 的驻定极小 ⇒ "任意随机起点 100% 命中"要加限定词。
[B] §10.3 唯一对不上的一格：$\Phi(F_2^{\rm Schur})=14.6398$。我算的"两条实不稳定本征向量张成的
    不变子空间"给 $78.6010$（与他另一行"A 的不稳定模"完全相同），稳定不变子空间给 $22.1037$，
    都不是 14.6398。⇒ 他那一格的对象没定义清楚。这里把 4 条本征向量的全部 $\binom{4}{2}=6$ 个两维
    坐标子空间、以及实 Schur 形 $Z$ 的列子集全部枚举，看 14.6398 是哪一个。
[C] 顺手钉一条他 §10.4 文本里的 $495\times$：表里没有这一格。
"""
import sys
sys.stdout.reconfigure(encoding='utf-8')
import numpy as np
from scipy.optimize import minimize

_src = open('p0/exp_c_audit.py', encoding='utf-8').read().split("print(r'== E53")[0]
_ns = {'__name__': 'e53_preamble'}
exec(compile(_src, 'p0/exp_c_audit.py[preamble]', 'exec'), _ns)
A, W, TH, JC, n, sym = _ns['A'], _ns['W'], _ns['TH'], _ns['JC'], _ns['n'], _ns['sym']
phi_iter, phi_dare = _ns['phi_iter'], _ns['phi_dare']
BEST = {r: np.load('.work3/c24_bestF_r%d.npy' % r) for r in (1, 2, 3)}
BEST = {r: (F if F.shape[0] == n else F.T) for r, F in BEST.items()}
OPT = {1: 10.169224, 2: 0.735143, 3: 0.158701}
print(r'== E59b：局部极小判定 + Schur 那格的对象判定 ==')


def phi_f(f):
    f = np.asarray(f, float)
    f = f / np.linalg.norm(f)
    return phi_iter(f.reshape(1, n))[0]


def tan_basis(f):
    f = f / np.linalg.norm(f)
    U, s, _ = np.linalg.svd(np.eye(n) - np.outer(f, f))
    return U[:, s > 0.5][:, :n - 1]


def retr(f, u):
    v = f / np.linalg.norm(f) + u
    return v / np.linalg.norm(v)


# ---------------- [A] ----------------
print(r'[A] 24 个停滞点：细差分重测梯度 + Powell 精修 + 逃跑测试')
rng = np.random.default_rng(771)


def gradnorm(f, h):
    B = tan_basis(f)
    g = np.array([(phi_f(retr(f, B[:, i] * h)) - phi_f(retr(f, -B[:, i] * h))) / (2 * h)
                  for i in range(n - 1)])
    return np.linalg.norm(g)


def powell_from(f0):
    B = tan_basis(f0)
    f0n = f0 / np.linalg.norm(f0)

    def obj(u):
        return phi_f(retr(f0n, B @ u))
    r = minimize(obj, np.zeros(n - 1), method='Powell',
                 options=dict(xtol=1e-10, ftol=1e-14, maxiter=4000))
    return r.fun, retr(f0n, B @ r.x)


def escape(f0, scale):
    """沿随机切向大步跳出去，再局部下降一小段，看能否低于停滞点。"""
    B = tan_basis(f0)
    f0n = f0 / np.linalg.norm(f0)
    best = phi_f(f0n)
    for _ in range(6):
        u = rng.standard_normal(n - 1) * scale
        f1 = retr(f0n, u / np.linalg.norm(u) * scale * np.sqrt(n - 1))
        cur = phi_f(f1)
        for _ in range(60):                       # 简易下降
            g = np.array([(phi_f(retr(f1, B2[:, i] * 1e-4)) - phi_f(retr(f1, -B2[:, i] * 1e-4))) / 2e-4
                          for i in range(n - 1)]) if False else None
            break
        best = min(best, cur)
    return best


# 重现 E59 [4] 的 24 个停滞点（同 seed 同算法，只记录终点）
stuck = []
for k in range(24):
    f = rng.standard_normal(n) if False else None
# 重新跑一遍下降（与 E59 同参数），把停滞终点存下来
rng2 = np.random.default_rng(771)


def run_descent(f):
    f = f / np.linalg.norm(f)
    L = 1.0
    it = 0
    for it in range(400):
        val = phi_f(f)
        B = tan_basis(f)
        g = np.array([(phi_f(retr(f, B[:, i] * 1e-5)) - phi_f(retr(f, -B[:, i] * 1e-5))) / 2e-5
                      for i in range(n - 1)])
        gn = np.linalg.norm(g)
        if gn < 1e-7:
            break
        step = L / (1.0 + gn)
        improved = False
        for _ in range(30):
            fn = retr(f, -step * (B @ g))
            if phi_f(fn) < val - 1e-12:
                f, L = fn, min(step * 3.0, 1e4)
                improved = True
                break
            step *= 0.5
        if not improved:
            L *= 0.25
            if L < 1e-14:
                break
    return f, it


for k in range(24):
    f0 = rng2.standard_normal(n)
    f1, it = run_descent(f0)
    stuck.append((f1, it))

sv = np.array([phi_f(f) for f, _ in stuck])
print('  停滞值排序：%s' % np.array2string(np.sort(sv), precision=4))
print('  其中已到全局最优(相对差<1e-5)：%d/24' % (np.abs(sv - OPT[1]) / OPT[1] < 1e-5).sum())
print('  %-9s %-10s %-10s %-11s %-11s' % ('停滞Φ', '|g| h=1e-4', '|g| h=1e-6', 'Powell 后', '相对最优'))
npow, fine = [], []
for i, (f, it) in enumerate(stuck[:8]):
    v = phi_f(f)
    g4, g6 = gradnorm(f, 1e-4), gradnorm(f, 1e-6)
    pv, _ = powell_from(f)
    npow.append(pv)
    fine.append((v, g4, g6, pv, (pv - OPT[1]) / OPT[1]))
    print('  %-9.4f %-10.3e %-10.3e %-11.4f %+-11.2e' % (v, g4, g6, pv, (pv - OPT[1]) / OPT[1]))
npow = np.array(npow)
print('  Powell 精修后：命中最优 %d/8；残留高于最优的有 %d 个，最大残留 %.4f（相对 %+.2e）'
      % ((np.abs(npow - OPT[1]) / OPT[1] < 1e-5).sum(), (npow > OPT[1] * 1.0001).sum(),
         npow.max() - OPT[1], (npow.max() - OPT[1]) / OPT[1]))

# 逃逸测试：停滞点附近是否存在更低的可达点（同一连通分量、免跨脊）
print('\n  逃跑测试（从停滞点沿随机切向跳 %s 后直接评值，看能否低于停滞值）：' % '5e-3,2e-2,1e-1')
for frac, (f, it) in list(enumerate(stuck))[:6]:
    v0 = phi_f(f)
    B = tan_basis(f)
    rows = []
    for sc in (5e-3, 2e-2, 1e-1):
        cand = []
        for _ in range(12):
            u = rng.standard_normal(n - 1)
            u = u / np.linalg.norm(u) * sc
            cand.append(phi_f(retr(f, u)))
        rows.append(min(cand))
    print('    停滞 %.4f：跳后最低 %s（低于停滞点的比例 %d/3 档）' %
          (v0, np.array2string(np.array(rows), precision=4), sum(r < v0 - 1e-6 for r in rows)))

# ---------------- [B] ----------------
print(r'\n[B] $\Phi=14.6398$ 到底是哪个两维子空间')
wA, vA = np.linalg.eig(A)
V = np.real(vA)
V = V / np.linalg.norm(V, axis=0, keepdims=True)
mods = np.abs(wA)
print('  $\\operatorname{eig}(A)=%s$，模 %s' %
      (np.array2string(wA, precision=4), np.array2string(mods, precision=4)))
import itertools
print('  %-14s %-26s %-12s %-10s' % ('列组合', '对应本征值', '我的定点 Φ', '与 14.6398'))
hits = []
for cols in itertools.combinations(range(n), 2):
    Q = np.linalg.qr(V[:, list(cols)])[0]
    v = phi_iter(Q.T)[0]
    print('  %-14s %-26s %-12.4f %+.4f' %
          (str(cols), np.array2string(np.abs(wA[list(cols)]), precision=3), v, v - 14.6398))
    hits.append(v)
from scipy.linalg import schur
for sort in ('None',):
    Tm, Z = schur(A, output='real')
    for cols in itertools.combinations(range(n), 2):
        Q = np.linalg.qr(Z[:, list(cols)])[0]
        v = phi_iter(Q.T)[0]
        tag = '  <-- 命中 14.6398' if abs(v - 14.6398) < 5e-3 else ''
        print('  Schur 列%-8s Φ=%-12.4f%s' % (str(cols), v, tag))
print('  枚举里出现 14.6398 的：%s' % [round(v, 4) for v in hits if abs(v - 14.6398) < 0.05])

# ---------------- [C] ----------------
print(r'\n[C] 他 §10.4 的 495× 有无出处')
print('  78.6010/0.158701 = %.1f   （= r=2 的"A 不稳模"行 ÷ r=3 的下降最优 ⇒ 串列）' % (78.6010 / OPT[3]))
print('  r=3 真正用 $A$ 的两条不稳模张成再补一维：Φ=%.4f ⇒ %.1f×' %
      (phi_iter(np.linalg.qr(np.hstack([V[:, [0, 1]], np.eye(n)[:, [3]]])).T)[0],
       phi_iter(np.linalg.qr(np.hstack([V[:, [0, 1]], np.eye(n)[:, [3]]])).T)[0] / OPT[3]))
print('  （列序未必是 0,1——用 $|\lambda|\ge1$ 的下标）')
UN = [k for k in range(n) if abs(wA[k]) >= 1]
Q3 = np.linalg.qr(np.hstack([V[:, UN], np.eye(n)[:, [3]]]))[0]
v3 = phi_iter(Q3.T)[0]
print('  正确写法：span(不稳模)+标准第 4 基 → Φ=%.4f = 最优的 %.1f 倍' % (v3, v3 / OPT[3]))
