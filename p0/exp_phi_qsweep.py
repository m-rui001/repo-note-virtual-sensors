"""E10（2026-09-28 第 3 轮，B 侧）：Φ(F;Q) 对 Q 的依赖 —— 架构排序会不会随 Q 翻转？
   顺带把"锥左端点预测"做成可被原文图表直接打脸的形式：D_min(K_F)=tr(WPc)+Φ(F;Q)。

方法全部用 E9 的无噪声定点（Γ→∞ 代理已对账到 <1e-3）。
"""
import numpy as np
from p0.p0_replicate_letter import (A, B, W, n, sym, ctrl, ss_filter)
from p0.exp_trv_penalty import schur_task
from p0.exp_phi_exact import noiseless_posterior

np.set_printoptions(precision=4, suppress=True)
F2, _ = schur_task(2)
F3, _ = schur_task(3)
Fs = {'k=2 不稳定Schur(封闭)': F2, 'k=3 +1稳定(不封闭)': F3,
      'e1,e4': np.array([[1., 0, 0, 0], [0, 0, 0, 1.]]),
      'e1,e2': np.array([[1., 0, 0, 0], [0, 1., 0, 0]]),
      '随机秩2#1': np.array([[0.8, -0.3, 0.5, 0.2], [0.1, 0.9, -0.4, 0.3]])}

QS = {'I': np.eye(n),
      'F2ᵀF2': sym(F2.T @ F2),
      'F3ᵀF3': sym(F3.T @ F3),
      'diag(1,1,10,10)': np.diag([1., 1, 10, 10]),
      'diag(10,10,1,1)': np.diag([10., 10, 1, 1]),
      'e1e2权重': np.diag([1., 1, 0, 0]),
      'Θ-like(KᵀK)': None}


def floor_and_phi(F, Qmat):
    try:
        Pc, K, Theta = ctrl(Qmat)
    except Exception as e:
        return None, None, None, 'LQR不可解'
    unc = float(np.trace(W @ Pc))
    P, res, t = noiseless_posterior(F)
    if not np.isfinite(P).all() or res > 1e-8:
        return None, unc, None, '定点未收敛(res=%.1e)' % res
    phi = float(np.trace(Theta @ P))
    return phi, unc, unc + phi, ''


if __name__ == '__main__':
    Pc0, K0, Th0 = ctrl(np.eye(n))
    QS['Θ-like(KᵀK)'] = sym(K0.T @ K0)

    print('%-22s' % 'F \\ Q', ''.join('%14s' % k for k in QS))
    tab = {}
    for tag, F in Fs.items():
        row = []
        for qk, Qm in QS.items():
            phi, unc, dm, err = floor_and_phi(F, Qm)
            row.append((phi, unc, dm, err))
        tab[tag] = row
        print('%-22s' % tag, ''.join('%14s' % (('%.2f' % r[0]) if r[0] is not None else r[3][:13]) for r in row))

    print('\n---- 同 Q 下的 D_min^unc（无约束地板）与相对超出 %')
    for j, qk in enumerate(QS):
        ucs = [tab[t][j][1] for t in Fs]
        print('  Q=%-14s D_min^unc=%9.4f' % (qk, np.mean([u for u in ucs if u])), end='')
        for tag in Fs:
            phi, unc, dm, err = tab[tag][j]
            if phi is None:
                print(' | %s:%s' % (tag.split()[0], err[:12]), end='')
            else:
                print(' | %s:%+.1f%%' % (tag.split()[0], 100 * phi / unc), end='')
        print()

    print('\n---- 排序是否随 Q 翻转（同一 Q 下按 Φ 升序给出架构名次）')
    for j, qk in enumerate(QS):
        vals = [(tab[t][j][0], t) for t in Fs if tab[t][j][0] is not None]
        vals.sort()
        print('  Q=%-14s 由便宜到贵: %s' % (qk, ' < '.join('%s(%.2f)' % (t.split()[0], v) for v, t in vals)))

    print('\n---- 可判决的图表预测：D_min(K_F) 就是 Δ(D) 曲线的左端点')
    for tag in ('k=2 不稳定Schur(封闭)', 'k=3 +1稳定(不封闭)'):
        for qk in ('I', 'F2ᵀF2'):
            phi, unc, dm, err = tab[tag][list(QS).index(qk)]
            if dm:
                print('  %-22s Q=%-8s D_min(K_F)=%9.4f  → D=50 %s ; D=80 %s' %
                      (tag, qk, dm, '可行' if 50 > dm else '不可行', '可行' if 80 > dm else '不可行'))
