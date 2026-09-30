# -*- coding: utf-8 -*-
"""Stage the public MIT repository payload for the JAIGP submission note.

Layout mirrors the working tree (p0/, preprint/) so that scripts and the
claims ledger keep their relative paths. Prints a manifest to
p0/e189_out.txt and runs a secret scan over every staged byte.
"""
import io
import os
import re
import shutil

SRC = os.path.dirname(os.path.abspath(__file__))          # .../out/rc/p0
ROOT = os.path.dirname(SRC)                               # .../out/rc
DEST = os.path.join(os.path.dirname(ROOT), 'repo-note-virtual-sensors')

EXCLUDE_NAME = re.compile(
    r'(^_|^board_before_|^community_md\.bak|^tex_sync|^__pycache__|^bak$|^note$|^audit_c_note$'
    r'|\.pkl$|\.png$|\.pyc$)')
# Board plumbing scripts (not the board itself).
BOARD_SCRIPT = re.compile(r'(board_apply|board_compress|_fix_line|board_append|_board_gate|_board_probe)')
# Full-text dumps of others' papers — do not redistribute.
COPYRIGHT_DUMP = re.compile(r'^\d{4}\.\d{4,5}\.txt$|_fulltext|_pdf_text\.txt$|_pdf_layout\.txt$')

SECRETS = [
    (r'gh[pousr]_[A-Za-z0-9]{20,}', 'github token'),
    (r'github_pat_[A-Za-z0-9_]{20,}', 'github fine-grained token'),
    (r'sk-[A-Za-z0-9]{16,}', 'api key'),
    (r'AKIA[0-9A-Z]{16}', 'aws key'),
    (r'-----BEGIN [A-Z ]*PRIVATE KEY-----', 'private key'),
    (r'(?i)(password|passwd|secret)\s*[:=]\s*\S{6,}', 'password literal'),
    # Realistic email: word chars + dot/hyphen @ domain with at least one dot, not np.linalg etc.
    (r'\b[A-Za-z0-9._%+-]+@[A-Za-z0-9.-]+\.[A-Za-z]{2,}\b', 'email'),
]
EMAIL_WHITELIST = re.compile(r'@(np\.|ref\b|cite\b|arxiv\.org|doi\.org|ucsb\.edu|titech\.ac\.jp)')


def walk(base, rel):
    """Files under base/rel, relative to base, honouring the exclusion rules."""
    root = os.path.join(base, rel)
    out = []
    for dirpath, dirnames, filenames in os.walk(root):
        dirnames[:] = [d for d in dirnames
                       if not EXCLUDE_NAME.search(d) and not BOARD_SCRIPT.search(d)]
        for f in filenames:
            if EXCLUDE_NAME.search(f) or BOARD_SCRIPT.search(f) or COPYRIGHT_DUMP.search(f):
                continue
            full = os.path.join(dirpath, f)
            out.append(os.path.relpath(full, base).replace('\\', '/'))
    return sorted(out)


def main():
    lines = []
    log = lines.append
    if os.path.isdir(DEST):
        shutil.rmtree(DEST)
    os.makedirs(DEST)

    staged = []
    staged += walk(ROOT, 'p0')
    staged += walk(ROOT, 'preprint')
    # Include the cross-lane board (author confirmed no secrecy).
    if os.path.isfile(os.path.join(ROOT, 'community.md')):
        staged.append('community.md')
    for extra in ('README_lanes.md',):
        if os.path.exists(os.path.join(ROOT, extra)):
            staged.append(extra)

    secrets = []
    n_bytes = 0
    for rel in staged:
        if rel.startswith('preprint/') and re.search(r'\.(aux|out|bbl|blg)$', rel):
            continue
        src = os.path.join(ROOT, rel.replace('/', os.sep))
        if not os.path.isfile(src):
            continue
        raw = io.open(src, 'rb').read()
        n_bytes += len(raw)
        text = raw.decode('utf-8', errors='replace')
        for pat, tag in SECRETS:
            for m in re.finditer(pat, text):
                if tag == 'email' and EMAIL_WHITELIST.search(m.group(0)):
                    continue
                secrets.append((rel, tag, m.group(0)[:60]))
        dst = os.path.join(DEST, rel.replace('/', os.sep))
        os.makedirs(os.path.dirname(dst), exist_ok=True)
        shutil.copyfile(src, dst)

    io.open(os.path.join(DEST, 'STAGING_MANIFEST.txt'), 'w', encoding='utf-8') \
        .write('\n'.join(staged) + '\n')

    log('== 暂存结果 ==')
    log('   目标目录 %s' % DEST)
    log('   文件数 %d，字节 %.2f MB' % (len(staged), n_bytes / 1e6))
    per = {}
    for rel in staged:
        top = rel.split('/')[0]
        per[top] = per.get(top, 0) + 1
    log('   分目录计数：' + ', '.join('%s=%d' % (k, v) for k, v in sorted(per.items())))
    big = sorted(((os.path.getsize(os.path.join(DEST, r.replace('/', os.sep))), r)
                  for r in staged if os.path.isfile(os.path.join(DEST, r.replace('/', os.sep)))),
                 reverse=True)[:8]
    log('   最大 8 件：')
    for s, r in big:
        log('      %8.0f B  %s' % (s, r))
    log('== 密钥/隐私扫描（须为 0）==')
    log('   命中 %d 条' % len(secrets))
    for rel, tag, s in secrets[:20]:
        log('      %s  [%s]  %s' % (rel, tag, s))
    log('== 本地绝对路径回显（应尽量 0；含则须改写）==')
    abs_hits = []
    # Match both forward-slash and escaped-backslash forms of the workspace root.
    prefix_fwd = os.path.abspath(ROOT).replace('\\', '/') + '/'
    prefix_bsl_raw = os.path.abspath(ROOT)  # e.g. E:\pdf\topics\out\rc
    for rel in staged:
        dst = os.path.join(DEST, rel.replace('/', os.sep))
        if not os.path.isfile(dst):
            continue
        raw = io.open(dst, 'rb').read()
        text = raw.decode('utf-8', errors='replace')
        new_text = text.replace(prefix_fwd, '').replace(prefix_bsl_raw, '')
        if new_text != text:
            io.open(dst, 'wb').write(new_text.encode('utf-8'))
            text = new_text
        m = re.search(r'[A-Za-z]:[\\/](Users|pdf)[\\/]', text)
        if m:
            abs_hits.append(rel)
    log('   含盘符绝对路径的文件 %d 个' % len(abs_hits))
    for r in abs_hits[:12]:
        log('      ' + r)
    log('   说明：跨车道板 community.md 已纳入；.work3/ 脚本未纳入（对方车道代码不在本仓库）。')
    io.open(os.path.join(SRC, 'e189_out.txt'), 'w', encoding='utf-8').write('\n'.join(lines) + '\n')
    print('\n'.join(lines).encode('utf-8', errors='replace').decode('utf-8'))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
