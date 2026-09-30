# -*- coding: utf-8 -*-
r"""e186：判决"跨支 prefactor 不同"这句话的**方向性**。

两车道都签过字的两句：
  * 我方 e179 [U4]：c_1 in [21.9,40.7] vs 我方 2/3 支 c in [116.7,333.5] ⇒ 判"平方律的 prefactor **不是**支无关常数"。
  * C 方补记二十四 86-A：Δ_min(55,r=1)=2.629e-3 ⇒ c_{1/2} ≤ 39.97 < 116.7，并写"**方向可靠**因为我的数是**上界**"；
    87-A 的修订联合句沿用"常数跨支不同"。
本脚本不跑优化，只做两件可复核的事：
  (1) 从落盘日志量出**同一格、两车道（或同车道两族）**的 Δ 上界读数之比 = 仪器松紧；
  (2) 把它与双方声称的**跨支效应量**并排，检查"上界带不相交"能不能判"真常数不同"。

口径（p0/e175_same_cell_delta.py:175--200）：dlt = I2 - SDP[t]，
  I2 = rank-r 设计在 25 起点搜索里取到的最小率；SDP[t] = **去掉秩约束**后的凸问题最优值。
  ⇒ I2 ≥ 真 rank-r 最优 ⇒ Δ 是上界；SDP[t] 是凸最优，且 e175 [X0] 已与 C 的 I_unc 对到 ≤1e-5。
  ⇒ 关键推论：**把秩约束直接丢弃的松弛恰好等于 SDP[t]，给出的 Δ ≥ 0 是空下界**；
    要判常数是否跨支相同，必须给**秩感知**的下界。[X4] 把这个需求换算成具体格的 bit 门槛。
"""
import re
import sys

OUT = []


def p(s=''):
    OUT.append(s)


def rd(f):
    return open(f, 'rb').read().decode('utf-8', errors='replace').replace('\r\n', '\n')


E178 = rd('p0/e178_out.txt')
E179 = rd('p0/e179_out.txt')
E175 = rd('p0/e175_same_cell_delta.py')
E179S = rd('p0/e179_rank1_branch.py')
C92 = rd('.work3/c92_out.txt')      # 只读对方车道日志（#51：引对方的数要去他自己的落盘核一行）

# ---------------------------------------------------------------- [X0] 口径来源
p('== [X0] Δ 为什么是上界（逐字取自本车道脚本，取不到就不出判决）==')
prov = []
for needle, tag in [('dlt = I2 - SDP[t]', '定义'),
                    ('措辞不得升格', '预注册措辞门')]:
    hit = [l.strip() for l in E175.split('\n') if needle in l]
    prov.append((tag, hit[0][:120] if hit else None))
for tag, line in prov:
    p('   %s：%s' % (tag, line if line else '**取不到 ⇒ 本节不成立**'))
hdr = [l.strip() for l in E179S.split('\n') if '不构成反证' in l]
p('   e179 文件头（我自己写下的方向纪律）：%s' % (hdr[0][:150] if hdr else '**缺**'))
n_prov = sum(1 for _, l in prov if l)
p('   有效来源行 %d/2；文件头纪律 %d 条' % (n_prov, 1 if hdr else 0))

# ---------------------------------------------------------------- [X1] 同格仪器松紧
p('')
p('== [X1] 同一格、两车道（或同车道两族）的 Δ 上界之比 = 仪器松紧（倍）==')
rows = []
in_v3 = False
for l in E178.split('\n'):
    if l.startswith('== [V3]'):
        in_v3 = True
        continue
    if in_v3 and l.startswith('=='):
        break
    if in_v3 and re.match(r'^\s+\d\d\.\d\d\s', l):
        m = re.search(u'（(c\\w+|e17\\d) 的 ([\\d.]+)', l)
        d = re.match(r'^\s*(\d\d\.\d\d)', l)
        if m and m.group(1).startswith('c9'):
            rows.append((d.group(1), m.group(1), float(m.group(2))))
for d, who, f in rows:
    p('   $D=%s$：C 方 %s 的读数是本格最小值的 %.1f$\\times$（e178 [V3] 末列）' % (d, who, f))
n_in = len(rows)
fac = [f for _, _, f in rows]
assert fac, '同格松紧读数 0 条 ⇒ 不许出判决（#46）'
p('   跨车道同格因子 %s ⇒ **最保守 %.1f$\\times$**、最大 %.1f$\\times$（有效 %d 格）'
  % (' '.join('%.1f' % f for f in fac), min(fac), max(fac), n_in))
c55 = [l for l in E178.split('\n') if l.strip().startswith('$D=55.00$') and '自身分歧' in l]
sd55 = None
if c55:
    m = re.search(r'自身分歧 \$([\d.]+)', c55[0])
    sd55 = float(m.group(1))
p('   C 自家两族（c92 vs c94）同格 $D{=}55.00$ 自我分歧 = %s$\\times$（e178 [V2]）'
  % (('%.1f' % sd55) if sd55 else '**未取到**'))

# ---------------------------------------------------------------- [X2] 效应量
p('')
def band(text, key):
    """紧跟 key 之后的 `[a,b]` 取两个浮点（纯字面定位，避开正则转义）。"""
    i = text.index(key) + len(key)
    j = text.index('[', i) + 1
    k = text.index(']', j)
    a, b = text[j:k].split(',')
    return float(a), float(b)


def fields(line, after, idx):
    """按行内字段出现顺序取第 idx 个空白分隔字段（#59：字段顺序而非正则）。"""
    i = line.index(after) + len(after)
    parts = [q for q in line[i:].split() if q]
    return parts[idx]


p('== [X2] 双方声称的"跨支常数不同"效应量（全部由落盘带端点相除）==')
u4line = [l for l in E179.split('\n') if '[U4]' in l][0]
lo1, hi1 = band(u4line, '$c_1')
lo2, hi2 = band(u4line, '支 $c')
row55l = [l for l in E179.split('\n') if l.strip().startswith('55.00')][0]
f55 = [q for q in row55l.split() if q]
# 列序（e179 表头）：D, delta_sw, r*, r, lambda1, lambda2, lambda2/lambda1, Delta1, seed, k/25, 残差
lam1_55, lam2_55, lam55, d55m = float(f55[4]), float(f55[5]), float(f55[6]), float(f55[7])
assert abs(lam2_55 / lam1_55 / lam55 - 1.0) < 2e-3, '列序解析错位：lambda2/lambda1 与前两列不自洽'
line_c92 = [l for l in C92.split('\n') if l.strip().startswith('55.00')][0]
fc92 = [q for q in line_c92.split() if q]
d55c = float(fc92[5].lstrip('+'))
cov55 = fc92[6].split('/')[0]
c55m, c55c = d55m / lam55 ** 2, d55c / lam55 ** 2
same55 = d55c / d55m
assert abs(c55c - 39.97) < 0.02, 'C 板 86-A 的 c=39.97 复算不上 ⇒ 列序或数据有变'
p('   [U4] 我方 $1/2$ 支 $c_1\\in[%.1f,%.1f]$，$2/3$ 支 $c\\in[%.1f,%.1f]$' % (lo1, hi1, lo2, hi2))
p('        带内极差（$2/3$ 支、同族工具）$=%.1f/%.1f=$ **%.3f$\\times$**' % (hi2, lo2, hi2 / lo2))
p('        两带**最近距离** $=%.1f/%.1f=$ **%.3f$\\times$**；极值比 $=%.1f/%.1f=$ %.3f$\\times$'
  % (lo2, hi1, lo2 / hi1, hi2, lo1, hi2 / lo1))
p('   同格 $D{=}55.00$（$\\lambda_1=%.4e$、$\\lambda_2=%.4e$、比值 $=%.4e$ 三列自洽）：'
  % (lam1_55, lam2_55, lam55))
p('        我方 $\\Delta=%.4e$ ⇒ $c=%.2f$；C $\\Delta=%.4e$（覆盖 %s/200）⇒ $c=%.2f$，'
  '与 C 板 86-A 的 $39.97$ **逐位对上**' % (d55m, c55m, d55c, cov55, c55c))
p('        ⇒ 同一格、同一分母，两车道的 $c$ 差 **%.3f$\\times$**（我方更紧）；'
  '分母相同 ⇒ 差全部来自 $\\Delta$ 上界的松紧，不含新物理' % same55)
p('   C 86-A 的"$\\ge2.9\\times$ 严格低于带下沿"复算：$116.7/%.2f=$ %.3f$\\times$' % (c55c, lo2 / c55c))

# ---------------------------------------------------------------- [X3] 判决
p('')
p('== [X3] 判决：上界带不相交能不能判真常数不同 ==')
instr, effect = min(fac), lo2 / hi1
fac_all = sorted(fac + [same55])
instr_min, instr_max = fac_all[0], fac_all[-1]
n_exceed = sum(1 for f in fac_all if f > effect)
p('   仪器尺度：同一格、跨车道（或跨族）的上界因子全体 %s（%d 个），其中 **%d 个大于效应尺度 %.2f$\\times$**；'
  'C 自家两族同格还差 %.1f$\\times$' % (' '.join('%.1f' % f for f in fac_all), len(fac_all), n_exceed, effect, sd55))
p('   ⇒ 所声称的跨支效应（%.2f$\\times$）落在仪器可造范围内（%.1f$\\times$ 到 %.3g$\\times$），'
  '且**上界互比本身无判决力**：普适常数只需 $c_{\\rm true}\\le\\min(%.1f,%.1f)$ 就与两带不相交完全兼容。'
  % (effect, instr_min, instr_max, lo1, lo2))
p('   逻辑：$\\Delta$ 是上界 ⇒ $c_{\\rm obs}=\\Delta_{\\rm ub}/\\lambda^{2}\\ge c_{\\rm true}$ 在两支都成立')
p('        ⇒ 普适常数只需满足 $c_{\\rm true}\\le\\min(%.1f,%.1f)$，与"两带不相交"**完全兼容**；'
  '两带上界互比只比谁的松弛更小（我方 [V2-更正] 早已对 C 用过这条，#53）。' % (lo1, lo2))
p('   [U6] **公开降级（措辞级；纪律账记为板账 #57）**：e179 [U4] 的判决句"平方律的 prefactor 不是支无关常数"'
  '**方向不成立**，改为"跨支常数是否相同：**未判**"。')
p('        C 86-A/87-A 的"常数跨支不同""方向可靠因为我的数是上界"是**同一处方向错误的对方版本**；'
  'C 的数作为**上界**没有错（$\\le39.97$ 我方复算通过），错在拿它去比我的上界。')
p('   仍然成立、我不撤回的部分：平方律形状（阈值端、两支持有，$N_1=1.00$、$\\mathrm{spr}=1.86$）；'
  '两带不相交这一**读数事实**；每支自己的合法上界 $c_{\\rm true}\\le%.1f$（$1/2$ 支）与 $\\le%.1f$（$2/3$ 支）。'
  % (lo1, lo2))

# ---------------------------------------------------------------- [X4] 可判条件
p('')
def kv(line, key):
    """行内 `key=<数>` 的字面值（去掉 $ 后按空白切字段，#59 字段出现顺序）。"""
    for tok in line.replace('$', '').split():
        if tok.startswith(key + '='):
            return float(tok.split('=', 1)[1])
    raise AssertionError('字段 %s 不在行内' % key)


p('== [X4] 要判掉这句话，下界证书需要紧到多少（换算到具体格，不是口号）==')
row34 = [l for l in E178.split('\n') if '$D=34.30$' in l and 'c=' in l][0]
lam34, d34 = kv(row34, '\\lambda_3/\\lambda_1'), kv(row34, '\\Delta_2')
row33 = [l for l in E178.split('\n') if '$D=33.00$' in l and 'c=' in l][0]
lam33, d33 = kv(row33, '\\lambda_3/\\lambda_1'), kv(row33, '\\Delta_2')
need_lo1 = lo1 * lam34 ** 2
p('   可判条件：取得 rank-$r$ 率的**下界** $L$ 后，若 $L_{2/3}/\\lambda^{2}>c_{1/2}$ 的最好上界（$%.1f$）'
  '⇒ 跨支不同得证。' % lo1)
p('   格 $D{=}34.30$（$2/3$ 支近阈，$\\lambda_3/\\lambda_1=%.4e$，现上界 $\\Delta=%.4e$，$c_{\\rm ub}=%.1f$）：'
  % (lam34, d34, d34 / lam34 ** 2))
p('        需 $\\Delta_{\\rm lower}\\ge$ **%.3e bit** $=$ 现上界的 **%.1f%%**；按 C 的 $c\\le%.2f$ 则只需 %.3e bit'
  % (need_lo1, 100.0 * need_lo1 / d34, c55c, c55c * lam34 ** 2))
p('   格 $D{=}33.00$（远端，$\\lambda_3/\\lambda_1=%.4e$，现上界 $\\Delta=%.4e$）：需 $\\Delta_{\\rm lower}\\ge$ %.3e bit'
  '（$=$ 现上界的 %.1f%%）' % (lam33, d33, c55c * lam33 ** 2, 100.0 * c55c * lam33 ** 2 / d33))
p('   为什么现在没有：$I_{\\rm unc}$ 一侧的 ${\\rm SDP}[t]$ 是凸最优（可信、且与 C 对到 $10^{-5}$），'
  '但"丢掉秩"这个最自然的松弛**正好回到 ${\\rm SDP}[t]$** ⇒ 只给 $\\Delta\\ge0$。')
p('   非空下界必须**知道秩**：对 $\\Gamma$ 谱冻结后的对角水填规划做拉格朗日对偶'
  '（$\\nu\\ge0$，$g(\\nu)=\\nu K+\\sum_i\\min_{0<r_i\\le1}(\\gamma_i r_i-\\nu\\ln r_i)$，'
  '因 feasible 集上 $\\log\\det R\\ge K$ 故 $g$ 合法）或对 rank-$r$ 的特征值互斥加约束。')
p('   ⇒ 预注册交付判据：证书在 $D{=}34.30$ 给 $\\Delta\\ge%.3e$ bit 即判"跨支不同"；'
  '给不到即判"未判"，两句都不许写"普适"。**不许**用 $\\Delta\\ge0$ 交差。' % need_lo1)

# ---------------------------------------------------------------- [X5] 条数
p('')
p('== [X5] 条数与退出（#46）==')
p('   [X0] 来源 %d/2+1；[X1] 同格松紧 %d 格（全部跨车道）+ C 自家分歧 %d 条；'
  '[X2] 同格 $c$ 复算 2 条；[X3] 判决依据因子 %d 个；[X4] 数值门槛 2 格'
  % (n_prov, n_in, 1 if sd55 else 0, n_in + 1))
assert n_prov == 2 and hdr and n_in >= 3 and sd55, '有效读数不足，不准出判决'
p('   本节不跑优化、不产新的 $\\Delta$ 读数；只做方向判定与门槛换算。')

txt = '\n'.join(OUT) + '\n'
open('p0/e186_out.txt', 'wb').write(txt.replace('\n', '\r\n').encode('utf-8'))
sys.stdout.write('e186 written %d lines\n' % len(OUT))
