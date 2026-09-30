---

## 62 [2026-09-30 20:5x | R47-D] 对 C 的 `run_all.sh` 复现审计：**范围问题已修，落盘问题没修**（`printf` 打到 stdout，不进 `$LOG`）；另外我们两车道现在各有**一份注记**

### A. 承认（C 这一轮做对的两件事）

1. **§53-A 刺 1 的"范围"部分被解决了**：`.work3/run_all.sh`（00:26 写）现在跑 **19 个脚本**，
   把我点名的 `c53_criterion / c54_alloc / c55_n4robust / c59_foldshape` 全放进循环，
   而且每个脚本后面 `printf '%-16s EXIT=%-3s %4ss\n'` —— 这正是我要求的格式。
   我上一轮说"§51(C) 的冻结跑只覆盖 13/19"，这条刺 C 接了。
2. **`c67_qa.py` 是我没见过的好设计**：24 个正文数字逐个回 `repro_all.log` 找，缺一个就
   `sys.exit(1)`。这比"脚本跑绿"强，因为它核的是**正文引用的量**而不是脚本的退出码。
   （它自己还没有输出件，见 C-3。）
3. **C 的 `\bibitem{ts}` 元数据比我的对**：1503.01848 的正式版是
   *Proc. 54th IEEE CDC, 2015, pp. 4486–4491*（Crossref 命中 `10.1109/cdc.2015.7402920`）。
   我原来只写 "arXiv:1503.01848, 2015"。已按 C 的口径升级我的 `tanaka2015joint` 条目。
   **这是一次"对方对我错"的记账**，不是客套。

### B. 挑刺（一处，具体到行）

`run_all.sh` 的退出码**不会进日志**。第 10 行
`timeout 1800 python -u ".work3/$s.py" >> "$LOG" 2>&1; ec=$?` 只把 python 的输出重定向了；
紧接着的 `printf ... EXIT=...` 和第 13 行的 `python -u .work3/c67_qa.py | tail -3` 都写到
**stdout**。⇒ 直接 `bash .work3/run_all.sh` 跑完之后，`repro_all.log` 里**仍然没有 `EXIT=` 行**，
和我 §53-A 刺 1 描述的现象一模一样；差别只是这次"证据"在终端上、跑完即失。
一行修法（任选）：
```bash
printf '%-16s EXIT=%-3s %4ss\n' "$s" "$ec" "$(( $(date +%s) - t0 ))" >> "$LOG" 2>&1
# 或整条：bash .work3/run_all.sh 2>&1 | tee -a .work3/repro_all.log
```
所以 §51(C) 的"13/13 脚本 EXIT=0"目前**仍然不可复核**，而且这次我能指出为什么修了脚本还是不可复核。
按纪律 20 的"证据要落在被引用的产物里"，这条判**未结**；补跑一次并把 stdout 并进 `$LOG` 即结。
（C-3 同一处：`c67_qa.py` 的 24 条对账也还没有落盘件，跑完请 `| tee .work3/qa_out.txt`。）

### C. 结构性问题：两份注记、两套引用

- 我在 `p0/note/note.tex`（11 节，§XI Related work + 17 条 `\bibitem`，6 页 / Overfull 0 / undefined 0）；
  C 在 `.work3/note_clean.tex`（9 节，含 `\section{Related-work positioning}` 的三条 itemize + 3 条 `\bibitem`，
  以及 `\section{Reproduction}`）。**两份都活着**，时间戳分别是 01:33（我）与 00:15（C）。
- 风险不是重复劳动，是**同一批数字以两种口径进入两份正文**：C 的 `note_clean.tex` 里没有
  `refpaper` 这个 key（它引 `\bibitem{ts}` 打头），而 `frag_D` 用的是 `\cite[Fig.~1]{refpaper}`；
  合并时编号会整体错位。
- 我的建议（请 C 或用户定）：**以 `p0/note/note.tex` 为唯一正文**，把 C 的
  `Related-work positioning` 三条 itemize 并进我的 §XI（它们讲的是"我们不与线性类 SDP 竞争设计"
  与"AM–GM 下界的地位"，与我的三条线不重叠，是补充），把 C 的 `\section{Reproduction}` 与
  `\section{What we do not claim}` 原样搬成 §XII/§XIII；`.work3/note_clean.tex` 转为**只读的历史件**。
  我不在 C 的车道里动手，等 C 自己搬或用户点头。
- 顺带一条硬事实：`.work3/note_clean.tex` 里 `refpaper`/`Moirangthem`/`task-restricted` **各 0 次命中**，
  三条 `\bibitem` 全是别人的论文 ⇒ 那份稿子作为"注记（comment）"投稿的话，**没有引用它要评论的那篇**。
  我的 `note.tex` 有 `refpaper`（DOI 齐全）并在 §III 用 `\cite[Fig.~1]{refpaper}` 指被复绘的图。
  合并时这一条必须按我的来，否则整篇的文体不成立。
- 顺带一条**引用防撞**：C 的 `\bibitem{gap}`（arXiv:1604.01227, Tanaka–Johansson–Oechtering–Sandberg–Skoglund,
  "Rate of prefix-free codes in LQG control systems"）与我引的 `\bibitem{cuvelier2021side}`
  （arXiv:2101.09329, Cuvelier–Tanaka, "Rate of prefix-free codes in LQG control systems **with side information**"）
  **标题只差四个词**、作者链重叠、结论不同（前者给"前缀码可达率的下界 $=$ 某个 DI 下确界、SDP 可算"，
  后者给"解码器免费持有状态子集"的率–代价前沿）。审稿人极易当成同一篇。
  我已把 1604.01227 收进我的引用集（key `tanaka2016prefix`），并在 §XI 里写明
  **"本文率轴上画的都是有向信息，它是前缀码比特率的下界 \cite{tanaka2016prefix}，不是码率本身"**——
  这句话本来就该在正文里（我的所有率数字都是 DI），现在它顺带把两篇分开。

### D. 我这轮的净交付

`p0/note/note.tex`：§XI 三段 Related work（计算线 / 限制坐在哪一侧 / 率侧）+ **17 条 `\bibitem`**
（§61 写的是 16，随后把 C 的 1604.01227 收进来并加了 DI≠码率那句），元数据逐条当场核：
13 条走 arXiv API（含 `<journal_ref>`，据此把 1606.01946/47 写成 CDC 2016、2101.09329 写成 CISS 2021）、
4 条走 Crossref（TAC 49(9)/TIT 52(8)/SICOMP 43(2)/CDC 2015 的卷期页与 DOI）、
L-CSS 本体走出版社抽取件的作者行 + DOI、Nair–Evans 另由 Sabag 原文的 [33] 交叉确认。
`frag_D` 新增 `rmk:dual`（§61-B 的前占切割）；`papers/notes/1503.01848.md` 式 (7) 的正负序改对。
编译：6 页、Overfull 0、undefined 0，渲染核对过 §XI 与参考文献栏。
**页数已经到顶**：C 的 §VIII/§IX `\TODO` 再进正文就会超 6 页，必须先做减法（我倾向合并 §X 与 §IX 的两份限定语清单）。
