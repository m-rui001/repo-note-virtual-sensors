# JAIGP 投稿材料清单（本车道 p0 产物）

期刊：JAIGP — Journal for AI Generated Papers（jaigp.org）。全部字段取自
`preprint/preprint.tex` 与 `p0/e187_out.txt`，粘贴文本见 `p0/jaigp_title.txt`、
`p0/jaigp_abstract.txt`。核对脚本：`p0/e187_jaigp_package.py`（退出码 0）。

## 1. 表单字段

| 字段 | 内容 | 来源 |
|---|---|---|
| Title | On Task-Restricted Virtual Sensors: Reproduction, Infeasibility Interval, Non-Factorisation, and the Rate-Side Ladder | `e187_out.txt` [P4] |
| Abstract | 单段，散文 250 词（不计行内数学式；含式 277 词）→ `p0/jaigp_abstract.txt` | [P4] |
| Keywords（正文 IEEEkeywords，表单无此格） | Remote state estimation, sensor scheduling, information form Kalman filter, non-convex optimisation, directed information, reverse water-filling | `preprint.tex` |
| Research Fields | 按 OpenAlex 建议输入，候选：networked control systems、estimation and filtering、information theory、rate-distortion theory、task-oriented communication、convex and semidefinite optimization、semantic communication、AI-generated research | 需现场匹配下拉项 |
| Academic Categories（1–5，必填） | 建议 4 条：Control Engineering / Systems & Control；Information & Coding Theory；Applied Mathematics（optimization）；Artificial Intelligence。⚠ 本站首页未公开候选值，须在表单现场按下拉列表逐项对齐，不要凭记忆填写 | [P4] + 首页核查 |
| PDF | `preprint/preprint.pdf`，12 页，351,560 B（0.34 MB，上限 20 MB），`pdflatex` 0 错误、0 未定义引用 | [P1][P2][P5] |
| LaTeX 源 | `preprint/jaigp_submission.zip`，0.04 MB（上限 30 MB），条目 `preprint.tex` + `fig_wall.pdf`；包内只有一个 `.tex` 且不 `\input` 任何文件，转换管线不会选错主文件；扁平版编译 0 错误、页数与多文件版一致（12 页） | [P1][P5] |
| Cover image | 可选，未生成。若要，建议用本文图 `fig_wall.pdf` 渲染 1200×800 JPG，不用装饰性生成图 | 待用户定 |
| 验证邮件 | 提交后 7 天内点验证链接，文章才公开 | 站点规则 |

## 2. Human Prompters

- 姓名：Xiangrui Meng（孟祥瑞），上海师范大学 数理学院（College of Mathematics
  and Sciences, Shanghai Normal University）
- 勾选 Lead Author：是
- ORCID：0009-0006-8902-0519 ⚠ **待本人确认**：这个号码没有出现在本工作区的任何
  文件里，我无法独立核实；填表前请从自己的 ORCID 记录复制。
- Google Scholar / rankless 链接：可选，未提供。

## 3. AI Co-Authors（名字/版本/角色）

| 名称 | 角色（论文里怎么写的事实） | 车道 |
|---|---|---|
| Qwen3.8-Flash | 驱动一条独立验证车道：写出 `p0/` 下的脚本与运行日志，提出并否决命题，逐条复核对方的数字 | lane D（`p0/`） |
| DeepSeek-V4.1-Flash | 驱动另一条车道：写出 `.work3/` 下的脚本与日志，负责地板侧与速率侧的推导、以及本文表格的逐数 QA | lane C（`.work3/`） |

⚠ **两个模型名与 C/D 车道的对应关系待本人确认**。我在正文 `\thanks` 里按"lane D = `p0/`、
lane C = `.work3/`"写了，这是从附录 `sec_appendix.tex` 的车道定义直接对上的；如果
两个模型名对应反了，改 `preprint.tex` 的 `\thanks` 两处即可，不涉及任何结论。

## 4. 版权与"是否重复发表"

- 本文不是任何已发表论文的版本。它是一篇 comment-and-extension note：沿用
  \cite{refpaper}（Moirangthem, Natarajan & Branicky, *Minimum directed information
  in LQG control under task-restricted sensing*, IEEE Control Systems Letters,
  vol. 10, pp. 1621–1626, 2026, doi: 10.1109/LCSYS.2026.3710208）的模型、符号与
  四态算例，只复用其公开算例与其 Fig. 1 的数字化读数，并在封面上明确写出这层关系
  （`preprint.tex` 第四条 `\thanks`）。该条 doi/卷/页已由本车道两条通道复算通过
  （`p0/e182_out.txt`）。
- 站点首页未声明 prior-publication 政策，也没有公开许可证条款（页脚只有
  "© 2026 JAIGP. All rights reserved."）。⚠ 若 JAIGP 后台要求转让版权，须先确认
  它与 IEEE 版权不冲突——这一点我不能从首页判定，需查作者协议原文。

## 5. 本轮已完成的格式改动（逐条可核）

1. 封面 `\author` 从 "Preprint draft — collaborative two-lane corpus, \today"
   改为真实署名 + 4 条 `\thanks`（人类作者/两个模型车道职责与可核对的日志位置/
   代码可得性/与已刊论文的关系）。
2. 删掉未使用的草稿宏 `\Pv`、`\TODO`（全文 0 处使用，`e187_out.txt` [P3]），
   `\today` 归零。
3. 摘要 480 词（含式口径）压到散文 250 词；被删掉的细节都还在引言的 8 条贡献
   清单里（`sec_intro.tex`），没有主张丢失。
4. 修好一句断句：`frag_D_repro_infeasible_nofactor.tex:14` 原为
   "The reference example of is a four-state plant with"（缺 `\cite{refpaper}`）。
5. 注释里的 "merged preprint" 字样清理；正文 `draft` 只剩 3 处对自家早期版本的
   自述句（撤回账的一部分，逐条列在 [P3]），不删。
6. 引用完整性：label 定义 44 / 引用 36 / 缺失 0；bibitem 定义 28 / 引用 28 /
   缺失 0 / 未引用 0。

## 6. 还没做完的事

- MIT public 仓库未建：`git push` 不能自动建仓，本机没有 `gh`，需要网页端或 API
  token 建一次空仓库，然后我按清单裁剪后推（`papers/` 52 MB、`p0/_cache_*.pkl`
  等不进仓库；`.work3/` 与 `E:\pdf\community.md` 未经你明确同意不公开）。
- 仓库地址要在填表时补回 `preprint.tex` 第三条 `\thanks`（现在写的是"地址随投稿
  表提供"，没有编造 URL）。
- `p0/claims_rc_lane.py` 需要为本轮追加条目。
