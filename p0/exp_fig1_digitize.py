"""E25：把 L-CSS 2026 Fig.1 逐像素读出来，与 Prop D/E 的预言同台对表。

Fig.1（`p0/fig/lcss_p5_x327.png`，从 PDF 页 5 抽出的 2043×1374 位图）画四条曲线：
  黑 "Tanaka unrestricted"、红 "2 unstable"、蓝 "2 unstable + 1 stable"、品红 "2 unstable + 2 stable"，
  外加一条灰虚线 "Unstable-mode lower bound"。

与本仓库的对应关系（这是关键，不是猜）：
  红 = 我的 F2 = 排序 Schur 的两个不稳定列（rank 2，动态封闭）→ Prop D 地板 46.1231；
  蓝 = 我的 F3 = 尾部3（1 稳定 + 2 不稳定，rank 3）→ 地板 36.6032；
  品红 = rank 4 = 满秩 = Prop 3 的 Δ≡0 情形 → 地板 D_min^unc = 31.4833；
  黑 = 无约束，三锚点已由 E24 复现到 +0.040/+0.006/+0.001 bit（D=33/40/80）。
#10 的起表预言：只有红曲线在图窗内起竖（46.12 ↔ 图上 ≈47），蓝与品红的地板 36.60/31.48 被横轴起点 40 裁掉。

标定：图框左边 = D 40、右边 = D 90、底边 = 1 bit、顶边 = 5 bit（与刻度标签同位）。
内部校验：灰虚线的读出值应等于 $\\sum_{|\\lambda|>1}\\log_2|\\lambda|=1.169$ bit（原文 Theorem 2 的下界）。
"""
import sys
import numpy as np
from PIL import Image
sys.stdout.reconfigure(encoding='utf-8')

IMG = 'p0/fig/lcss_p5_x327.png'
D_LEFT, D_RIGHT, I_BOT, I_TOP = 40.0, 90.0, 1.0, 5.0
PRED_FLOOR = {'red': 46.1231, 'blue': 36.6032, 'magenta': 31.4833}
MINE = {  # E23 (B) 真最优链：D = 地板+x → 我的 I_opt（bits/sample）
    48.1231: 3.8746, 47.5231: 4.3208, 46.9231: 5.0569, 46.5001: 6.0901,
    46.3231: 6.9823, 46.2231: 7.9695, 46.1731: 8.9631, 46.1531: 9.6975,
}


def frame(a):
    """找图框：最长的水平/垂直暗线。"""
    dark = (a[..., 0] < 120) & (a[..., 1] < 120) & (a[..., 2] < 120)
    H, W = dark.shape
    rows = np.argsort(-dark.sum(1))[:6]
    cols = np.argsort(-dark.sum(0))[:6]
    rh = [r for r in rows if dark[r].sum() > .7 * W]
    cl = [c for c in cols if dark[:, c].sum() > .5 * H]
    return min(rh), max(rh), min(cl), max(cl)


def series(a, key, box, legend=None):
    r0, r1, c0, c1 = box
    R, G, B = a[..., 0].astype(int), a[..., 1].astype(int), a[..., 2].astype(int)
    mx = np.maximum(np.maximum(R, G), B)
    mn = np.minimum(np.minimum(R, G), B)
    if key == 'red':
        m = (R > 150) & (G < 110) & (B < 110)
    elif key == 'blue':
        m = (B > 150) & (R < 110) & (G < 160)
    elif key == 'magenta':
        m = (R > 150) & (B > 150) & (G < 120)
    elif key == 'black':
        m = (mx < 90)
    elif key == 'dash':
        m = (mx > 110) & (mx < 215) & ((mx - mn) < 28)      # 中性灰，排除彩色曲线
    mg = 30 if key == 'black' else 3                        # MATLAB 刻度朝内，黑曲线要留足边距
    m[:r0 + mg, :] = m[r1 - mg:, :] = False
    m[:, :c0 + mg] = m[:, c1 - mg:] = False
    if legend is not None and key != 'dash':
        lr0, lr1, lc0, lc1 = legend
        m[lr0 - 4:lr1 + 5, lc0 - 4:lc1 + 5] = False
    out = {}
    for c in range(c0 + mg, c1 - mg):
        rr = np.where(m[r0 + mg:r1 - mg, c])[0]
        if len(rr):
            # 一列可能同时穿过曲线与标记；取最长连续段的中点
            split = np.where(np.diff(rr) > 3)[0]
            segs = np.split(rr, split + 1)
            sg = max(segs, key=len)
            out[c] = (sg.mean() + r0 + mg, len(sg))
    return out


def _runs(idx):
    if not len(idx):
        return np.array([], int), np.array([], int)
    brk = np.where(np.diff(idx) > 1)[0]
    segs = np.split(idx, brk + 1)
    return np.array([s[0] for s in segs]), np.array([len(s) for s in segs])


def _longest(idx, L):
    st, ln = _runs(idx)
    return (np.max(ln) - 1) if len(ln) else 0          # 间隙 = 段长-1，这里要的是最长连续段


def legend_box(a, box):
    """图例框 = 图框内部的**连续**长黑直线（用最长连续段判据，虚线不会被误判）。"""
    r0, r1, c0, c1 = box
    blk = (a.max(2).astype(int) < 90)
    sub = blk[r0 + 6:r1 - 6, c0 + 6:c1 - 6]
    H, W = sub.shape
    rows = [i + r0 + 6 for i in range(H)
            if len(_runs(np.where(sub[i])[0])[1]) and _longest(np.where(sub[i])[0], W) > .2 * W]
    cols = [j + c0 + 6 for j in range(W)
            if len(_runs(np.where(sub[:, j])[0])[1]) and _longest(np.where(sub[:, j])[0], H) > .15 * H]
    if len(rows) < 2 or len(cols) < 2:
        return None
    return min(rows), max(rows), min(cols), max(cols)


def dash_row(a, box):
    """Unstable-mode lower bound 那条虚线：短段数最多的行。"""
    r0, r1, c0, c1 = box
    mx = a.max(2).astype(int)
    mn = a.min(2).astype(int)
    g = (mx < 215) & ((mx - mn) < 28)
    best, bn = None, 0
    for i in range(r0 + 6, r1 - 6):
        st, ln = _runs(np.where(g[i, c0 + 6:c1 - 6])[0])
        k = int(np.sum((ln >= 4) & (ln <= 60)))
        if k > bn:
            best, bn = i, k
    return best, bn


def to_data(px, py, box):
    r0, r1, c0, c1 = box
    return (D_LEFT + (px - c0) * (D_RIGHT - D_LEFT) / (c1 - c0),
            I_BOT + (r1 - py) * (I_TOP - I_BOT) / (r1 - r0))


if __name__ == '__main__':
    im = Image.open(IMG).convert('RGB')
    a = np.asarray(im)
    box = frame(a)
    leg = legend_box(a, box)
    print('图框像素 (r0,r1,c0,c1)=%s  图例框=%s  图像 %dx%d'
          % (box, leg, a.shape[1], a.shape[0]))
    # 标定自检 1：灰虚线（Unstable-mode lower bound）应读成 1.169 bit
    rd, nrun = dash_row(a, box)
    print('标定自检① 虚线行=%s（%d 个短段）读出 = %.4f bit（Theorem 2 下界 1.169）→ 误差 %+.4f'
          % (rd, nrun, to_data(0, rd, box)[1], to_data(0, rd, box)[1] - 1.169))

    print('\n%-9s %8s %8s %10s %10s   %s' % ('曲线', '左端D', '顶端I', '起竖D@I=4', '地板预言', '与预言差'))
    curves = {}
    for key in ('red', 'blue', 'magenta', 'black'):
        s = series(a, key, box, leg)
        pts = []
        for c, (r, n) in sorted(s.items()):
            D, I = to_data(c, r, box)
            pts.append((D, I))
        pts = np.array(pts)
        if not len(pts):
            print('%-9s  未采到' % key)
            continue
        curves[key] = pts
        i4 = pts[(pts[:, 1] > 3.8) & (pts[:, 1] < 4.2)]
        print('%-9s %8.3f %8.3f %10s %10.4f   %s'
              % (key, pts[:, 0].min(), pts[:, 1].max(),
                 '%.3f' % i4[:, 0].min() if len(i4) else '—',
                 PRED_FLOOR.get(key, np.nan),
                 '%+.3f' % (i4[:, 0].min() - PRED_FLOOR[key]) if (len(i4) and key in PRED_FLOOR) else '—'))

    # 标定自检 2：黑曲线（Tanaka 无约束）在 D=40 应读成 3.266（原文 §V 已发表值，E24 复现 3.2719）
    Bk = curves['black']
    for tgt, ref in ((40.0, 3.266), (80.0, 1.602)):
        sel = Bk[np.abs(Bk[:, 0] - tgt) < .15]
        if len(sel):
            print('标定自检② 黑曲线 D=%.0f 读出 %.4f bit（原文 %.3f / E24 复现 %.4f）→ 差 %+0.4f'
                  % (tgt, sel[:, 1].mean(), ref, ref + .006, sel[:, 1].mean() - ref))

    print('\n红曲线（=我的 F2）逐点读出 vs 我的真最优链 I_opt：')
    print('%9s %9s %9s %8s' % ('图上D', '图上I', '我的I', '差'))
    R = curves['red']
    for Dm, Im in sorted(MINE.items()):
        sel = R[np.abs(R[:, 0] - Dm) < .12]
        if len(sel):
            print('%9.4f %9.4f %9.4f %+8.4f' % (Dm, sel[:, 1].mean(), Im, Im - sel[:, 1].mean()))
    print('\n红曲线在 I=5（图顶）处的 D = %.3f；我的地板 %.4f；#10 预言起表 47.0±0.5'
          % (R[np.argmax(R[:, 1])][0], PRED_FLOOR['red']))
    np.save('p0/fig/fig1_red.npy', curves['red'])
    for k in ('blue', 'magenta', 'black'):
        np.save('p0/fig/fig1_%s.npy' % k, curves[k])
    print('四条曲线已存 p0/fig/fig1_*.npy（列：D, I_dir）')
