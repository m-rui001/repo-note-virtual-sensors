
### I. 本轮正文**有**改动（与 §63-F 的"正文零改动"相反，因为这次修的是我自己挖出的引用缺陷）

`p0/note/note.tex` 两处，共 8 行：
1. `\bibitem{kostina2019rate}`：V.~Kostina and B.~Hassibi, "Rate-cost tradeoffs in control,"
   *IEEE Trans. Autom. Control*, vol. 64, pp. 4525--4540, 2019, doi: 10.1109/TAC.2019.2912256。
   **卷号与页码按 Crossref 当场核；issue 号没核到所以不写**（对比：同一清单里 `tanaka2018di` 有 no.~1，
   那条是核过的）。引用表从 **17 条变 18 条**，§61 那一格里的"17"作为历史记录保留不改。
2. §XI "Rate side" 段插入 6 行定位句：KH 的率–代价函数是**最小化率**那一侧；其部分观测形式
   **只在观测矩阵满秩时非平凡**，受限传感下右端塌成一个与代价无关的常数，
   所以"本注记研究的预算侧渐近线读不出任何已发表界"。
   措辞上我**刻意没有**写 `rank C=n` 这种符号（正文的 `C` 还没定义成观测矩阵，避免自造记号冲突），
   也**没有**写"算不出"（KH Remark 3 说 Tanaka 的 SDP 是精确的，见 B 段）。

渲染核对（`pdflatex` ×2）：**exit=0、6 页、Overfull 0、undefined 0**；
从 PDF 抽文验证新句与文献号对得上：正文印作 "nearest published converse is the rate-cost function of **[9]**"，
表内 `[8] Sabag–Tian–Kostina–Hassibi` / `[9] Kostina and Hassibi` / `[17] Atay–Chandrasekaran–Kostina`，
**引子不引母的缺口已闭合**，页数仍在 L-CSS 上限内（§VIII/§IX 的 `\TODO` 若要落，还是得先做 §61-D 说的 §X/§IX 合并）。
