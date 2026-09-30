# -*- coding: utf-8 -*-
# p0/claims_rc_lane.py —— R77-L 车道的主张清单，给 C 的 repro/check_claims.py **import** 用。
# 我不写 repro/（不是本车道）。用法：
#     from p0.claims_rc_lane import CLAIMS      # CLAIMS = [(id, 主张, 日志路径, 必现子串), ...]
# 直接跑本文件会自校验全部条目并打印条数（#46）。
# 纪律：#45 每条 needle 都必须能在落盘 stdout / 正文文件里当场 grep 到；改一个字要重跑。

CLAIMS = [
    # ---- e177（2/3 支：噪声底 + 三模型门）----
    ('RC-e177-N', '种子噪声底 N=1.53（公共 5 格里 4 格两种子比值=1.000）',
     'p0/e177_out.txt', ['$= $ 1.000', 'N=1.53']),
    ('RC-e177-M1', 'M1（Δ∝λ3）不过线：spr 657.39',
     'p0/e177_out.txt', ['spr $=657.39$']),
    ('RC-e177-M2', 'M2（∝λ3/λ1）不过线：spr 1343.91',
     'p0/e177_out.txt', ['spr $=1343.91$']),
    ('RC-e177-M3', 'M3（∝(λ3/λ1)^2）过线：spr 2.86 ≤ 2N=3.06',
     'p0/e177_out.txt', ['spr $=2.86$', '2N=3.06']),
    # ---- e178（带与有效指数）----
    ('RC-e178-band', '近阈值 4 格的带 c∈[116.7,138.7]（带宽 1.19×）',
     'p0/e178_out.txt', ['116.7', '138.7']),
    ('RC-e178-p', '有效指数：全程 2.40 → 末端 2.15（阈值端趋二次）',
     'p0/e178_out.txt', ['2.40', '2.15']),
    # ---- e179（1/2 支 = C 的支，我方工具）----
    ('RC-e179-N1', 'C 支的同族种子底 N1=1.00（公共 4 格比值 1.000）',
     'p0/e179_out.txt', ['公共格 4 个', '$=$ 1.000']),
    ('RC-e179-M3p', "M3' 在 1/2 支过线：spr 1.86（M1' 62.38、M2' 82.96 不过线）",
     'p0/e179_out.txt', ['spr $=1.86$', 'spr $=62.38$', 'spr $=82.96$']),
    ('RC-e179-c1', 'c1∈[21.9,40.7]：D=50 → 40.7、D=55 → 24.9',
     'p0/e179_out.txt', ['[21.9,40.7]']),
    ('RC-e179-tighter', '同格我更紧：Δ1(50)=4.387e-02、Δ1(55)=1.636e-03（C 报 5.734e-02 / 2.629e-03）',
     'p0/e179_out.txt', ['4.387e-02', '1.636e-03']),
    ('RC-e179-U4', '[U4] 跨支 prefactor 差 15.3× ⇒ 两支共享的只是形状、不是常数',
     'p0/e179_out.txt', ['15.3', '支无关常数']),
    # ---- e180（预注册外推检验）----
    ('RC-e180-grid', '新格 5 个与已跑格不相交、第三条种子流 20261001、r* 全 3',
     'p0/e180_out.txt', ['20261001', '不相交']),
    ('RC-e180-ladder', 'c 沿 δ_sw 降序严格下降：647.7 > 278.9 > 183.4 > 134.5 > 113.6',
     'p0/e180_out.txt', ['647.7', '278.9', '183.4', '134.5', '113.6']),
    ('RC-e180-hit', '带外预测命中 1 格：D=34.05 给 c=134.5 ∈ [116.7,138.7]',
     'p0/e180_out.txt', ['134.5']),
    # ---- e180b（功效复核：撤回两条判决）----
    ('RC-e180b-power', '有功效的否证需 c<116.7/1.53=76.3；实测缺口 1.027 ≤ N ⇒ 无功效',
     'p0/e180b_out.txt', ['c<116.7/1.53=76.3', '1.027']),
    ('RC-e180b-retract', '[P2-iii] 作废 + 两句撤回（含"方向是反的"）',
     'p0/e180b_out.txt', ['[P2-iii] **作废（无功效）**', '后者方向是反的', '撤回这两句']),
    ('RC-e180b-res', '可分辨相邻步 1／4（步比 2.32 / 1.52 / 1.36 / 1.18）',
     'p0/e180b_out.txt', ['可分辨相邻步 1／4', '2.32', '1.52', '1.36', '1.18']),
    ('RC-e180b-count', '判决行带条数：解析 5/应 5；带内 1、高于上沿 3、低于下沿 1；退出 0',
     'p0/e180b_out.txt', ['解析 5 行／应 5']),
    # ---- 正文（与 C 在同一份 preprint/ 上作业）----
    ('RC-tex-abst', '摘要 (iv)：closed-form value 只在 wall prior，那个常数属于被撤回的 frontier-wide 形式'
     '（2026-09-30 摘要压缩后措辞已更新，主张不变）',
     'preprint/preprint.tex',
     ['yielding a closed-form value at the \\emph{wall prior}',
      'would have fixed the constant of the wall-approach']),
    ('RC-tex-tabIV', '表 IV（tab:bound）题注与脚注都写 "not a frontier bound"',
     'preprint/frag_rc_waterfill.tex',
     ['at the wall prior at $I=3$ (a value at', 'not a frontier bound', '\\label{tab:bound}']),
    ('RC-tex-XI', '§XI：withdrawn, not ``unproved\'\'（115 in-cone 反例在内）',
     'preprint/sec_related.tex', ['is withdrawn, not', '$115$ reachable in-cone']),
    # ---- e181–e185（本轮：文献取数 / 书目解析 / 补 DOI / TAC-2022 适用域 / 渲染验句）----
    ('RC-e181-fetch', '三篇目标逐条双通道取数：crossref 3/3、openalex 3/3、摘要 3 条、'
     '开放全文只 1 篇落盘（另两篇无 OA 链接 ⇒ 明写"不下"，不编造全文）',
     'p0/e181_out.txt', ['目标 3 条；crossref 成功 3；openalex 成功 3；摘要可得 3；OA 全文落盘 1',
                        '无 OA 链接 ⇒ 不下（不编造全文）']),
    ('RC-e182-bib', '旧书目 24 条的locator全部当场解析：DOI 9/9、arXiv 15/15、'
     '标题等值 0 条不等；无定位符 1 条（nair2004stabilizability）；会议条目探测到候选期刊升版 3 条',
     'p0/e182_out.txt', ['DOI 解析 9／9', 'arXiv 解析 15／15', '候选升版 3 条',
                         "['nair2004stabilizability']"]),
    ('RC-e182b-doi', '补 DOI：双门槛（标题等值＋作者重叠）可补 14 条、未命中 0；'
     '手工升级项 1 条（e181 双通道核过的改标题期刊版）',
     'p0/e182b_out.txt', ['可补条目 14；未命中／失败 0', '待核 12；可补 DOI 14（去重后 14）']),
    ('RC-e183-tac03', 'TAC-2022（Stavrou–Skoglund）Prop.1 的三个结构情形在本锚点植物上命中 0/3：'
     'αI 残差 1.77、对称残差 1.87、A=W 残差 6.83 ⇒ 他们的式 (14) 不适用于本植物',
     'p0/e183_out.txt', ['三情形合计命中 **0／3**', '情形 1  $A=\\alpha I_p$', '6.8300e+00']),
    ('RC-e183-comm', '他们论证真正依赖的两件事在本锚点都不成立：$\\max|AW-WA|=7.7065e+00$（不能同时'
     '正交对角化）、eig($A^2$) 虚部最大 3.978e-01 ⇒ (14) 的 $\\mu_{A^2,i}$ 不是实数',
     'p0/e183_out.txt', ['7.7065e+00', '虚部最大 $3.978e-01$', '结构情形命中 0/3；交换残差 7.706e+00']),
    ('RC-e185-mode', '双栏抽取假象定量：15 条抽查短语在栏序模式下 15/15 命中、在 -layout 下 8 条 MISS；'
     '132 个句子的字母串原样连续出现数 = 栏序 101 vs 按行 27 ⇒ 先前"丢句"是 pdftotext 把左右栏'
     '同一行的碎片拼在一起，不是正文缺失（本行为摘要改写后 e185 重跑的新读数，改前 100/32）',
     'p0/e185_out.txt', ['layout=MISS', 'default=OK', '栏序=101 / 按行=27']),
    ('RC-e185-sent', '句子完整性（判决，摘要改写后已对新 PDF 重跑）：132 个 prose 句子的每个字母都在 '
     '12 页渲染里按序出现，字母级缺口 0（允许 \\ref 插入与浮动体跨块，共 2 次跨块跳转，改前 3 次）；'
     'caption 3/3、书目 28 条字母级完整',
     'p0/e185_out.txt', ['字母级缺口句=0', '跨块跳转=2', '书目 28 条完整',
                         'caption 3/3 完整']),
    ('RC-e185-ctrl', '门限有效性对照：从渲染文本删去一句的中段 12 个字母 ⇒ gaps_of 当场报缺失；'
     '同句对未改动渲染无缺口 ⇒ 0 缺口不是空转出来的',
     'p0/e185_out.txt', ['gaps_of 报缺失=True', "报出的缺口=['auseno']", '无缺口 True']),
    ('RC-e185-aux', '编号权威是 .aux 不是静态数数：eq:wf→(8)、eq:cond→(9)，'
     '正文引用 36 个 label 全在 aux 中有编号（渲染正文里 "(8) needs one disambiguation"、"(9) is conditional" 原地可见）',
     'p0/e185_out.txt', ['eq:wf     -> aux 记录 (8)', 'eq:cond   -> aux 记录 (9)',
                         '引用但 aux 无编号=none']),
    ('RC-tex-cite-tac', '§VII 里 TAC-2022 期刊版与其会议版同时被引（"closing a question"一句）',
     'preprint/frag_rc_waterfill.tex',
     ['\\cite{stavrou2022tac} (conference version \\cite{stavrou2018asymptotic})']),
    ('RC-tex-cite-sched', '§II 的调度线与 1/2 支对照引 shi2012scheduling（"其 rate 是占空比"），'
     '传感器选择线引 li2026oss',
     'preprint/sec_related.tex', ['shi2012scheduling', 'li2026oss']),
    ('RC-tex-bib28', 'sec_bib.tex 现有 \\bibitem{stavrou2022tac}（期刊版新条目）',
     'preprint/sec_bib.tex', ['\\bibitem{stavrou2022tac}']),
    ('RC-e186-upperdir', '公开降级（措辞级，板账 #57）：e179 [U4] 那句"平方律的 prefactor 不是支无关常数"'
     '方向不成立 —— $\\Delta$ 是上界，两条上界带不相交推不出真常数不同；改判为"跨支常数是否相同：未判"。'
     '对方 86-A/87-A 的"常数跨支不同"是同一处方向错误的对方版本，其数作为上界没有错',
     'p0/e186_out.txt', ['方向不成立', '跨支常数是否相同：**未判**', '同一处方向错误的对方版本']),
    ('RC-e186-instrument', '仪器尺度量出来了：同格跨车道（或同车道跨族）的上界因子 1.6/8.7/11.0/122.5/1821705.4，'
     '5 个里 4 个大于所声称的效应尺度 2.87$\\times$（两带最近距离 2.867$\\times$）；普适常数只需 '
     '$c_{\\rm true}\\le\\min(21.9,116.7)$ 就与两带不相交完全兼容 ⇒ 上界互比本身无判决力',
     'p0/e186_out.txt', ['1821705.4', '2.867$\\times$', '4 个大于效应尺度 2.87$\\times$',
                         '与"两带不相交"**完全兼容**']),
    ('RC-e186-lowerneed', '要判掉那句话需要 rank-aware 下界证书，门槛换算到具体格：$D{=}34.30$ 需 '
     '$\\Delta_{\\rm lower}\\ge4.927e$-05 bit（现上界的 18.8%）；去秩约束恰好吃回 SDP[t] ⇒ 只有 $\\Delta\\ge0$ '
     '的空下界，候选路线是谱冻结对角水填规划的拉格朗日对偶；判据预注册：给不到就判"未判"，两句都不许写"普适"',
     'p0/e186_out.txt', ['4.927e-05 bit', '拉格朗日对偶', '**不许**用 $\\Delta\\ge0$ 交差',
                         '给不到即判"未判"']),
    ('RC-jaigp-front', '投稿格式（非预印本）：封面署名 Xiangrui Meng + Qwen3.8-Flash + DeepSeek-V4.1-Flash，'
     '四条 \\thanks 交代人类作者职责、两条车道各自的脚本与日志位置、代码可得性，以及本文对已刊 '
     'LCSYS-2026 那条是 comment-and-extension 而非其版本；未用的草稿宏 \\Pv/\\TODO 已删',
     'preprint/preprint.tex', ['Xiangrui~Meng', 'Qwen3.8-Flash', 'DeepSeek-V4.1-Flash',
                               'is a comment on, and reproduction study of']),
    ('RC-jaigp-flat', '上传包对账：单文件源编译 0 错误、页数与多文件版一致（12 页）；label 44 定义/36 引用/缺失 0，'
     'bibitem 28/28/缺失 0；正文 preprint/to be settled/\\TODO/\\Pv/\\today 计 0（"draft" 3 处是对自家早期版本的'
     '自述，逐条列出后放行）；摘要单段、散文 250 词；PDF 与源包均在表单限额内',
     'p0/e187_out.txt', ['引用但未定义 0 个', '散文 250 词', '页数 一致',
                         "['fig_wall.pdf', 'preprint.tex']"]),
]


def verify():
    import os
    bad, cache = [], {}
    for cid, claim, path, needles in CLAIMS:
        if not os.path.exists(path):
            bad.append((cid, '文件不存在 ' + path))
            continue
        if path not in cache:
            with open(path, 'rb') as f:
                cache[path] = f.read().decode('utf-8', 'replace')
        hay = cache[path]
        for nd in needles:
            if nd not in hay:
                bad.append((cid, '缺 needle %r in %s' % (nd[:40], path)))
    return bad


if __name__ == '__main__':
    import io, sys
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')
    bad = verify()
    files = sorted(set(c[2] for c in CLAIMS))
    print('== p0/claims_rc_lane.py 自校验 ==')
    print('   主张 %d 条，覆盖文件 %d 个：%s' % (len(CLAIMS), len(files), ', '.join(files)))
    print('   失败 %d 条 $-$ $-$ %s' % (len(bad), '全部通过' if not bad else ''))
    for cid, why in bad:
        print('   ✗ %s：%s' % (cid, why))
    print('   退出码 %d（0 $=$ 全过）' % (1 if bad else 0))
    sys.exit(1 if bad else 0)
