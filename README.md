# On Task-Restricted Virtual Sensors

Reproduction package for the JAIGP submission note:
**"On Task-Restricted Virtual Sensors: Reproduction, Infeasibility Interval, Non-Factorisation, and the Rate-Side Ladder"**
by Xiangrui Meng, Qwen3.8-Flash, and DeepSeek-V4.1-Flash.

## What this repo contains

| Directory | Contents |
|-----------|----------|
| `preprint/` | LaTeX source (single self-contained `.tex`), compiled PDF, wall-figure PDF |
| `p0/`     | Lane-D verification scripts (`e*.py`), evidence receipts (`*_out.txt`), data files (`fig/*.npy`, `results/`) |
| `agent-logs/` | Cross-lane collaboration board (`community.md`), work log, claims ledger |

## Quick start

```bash
pip install -r requirements.txt
cd p0
python e187_jaigp_package.py    # flatten + QA the LaTeX source
python claims_rc_lane.py verify  # self-check the claims ledger (38 claims, 0 failures)
```

## Two-lane design

This note was produced by two independent verification lanes driven by the same human prompter:

- **Lane D** (Qwen3.8-Flash): scripts in `p0/`, evidence logs `p0/*_out.txt`.
- **Lane C** (DeepSeek-V4.1-Flash): scripts in a separate workspace; results are quoted in the cross-lane board `community.md`.

The full cross-lane board is included here because all artefacts belong to the same principal (the human author). Lane-C's *code* is not redistributed — only its published numbers and conclusions appear in the board and in lane-D's evidence logs.

## Key results (from the paper)

1. The published four-state example's rate-cost frontiers are optima over the whole admissible cone, to within 0.07 cost units inside the digitisation band.
2. No isotropic-noise family reaches them (1.11–16.4 units above).
3. The admissibility wall is a horizontal asymptote of every frontier; the interval `(Δ_wall, Δ_min(R)]` is a complete infeasible region.
4. The wall does **not** order the frontiers (111/120 random rank-one designs invert the ordering at fixed rate).
5. At a frozen prior the trade-off reduces to one convex program on the information form, solved by reverse water-filling.

See `preprint/preprint.pdf` for the full note.

## Licence

MIT — see `LICENSE`.

## Author

Xiangrui Meng, College of Mathematics and Sciences, Shanghai Normal University, Shanghai, China.
ORCID iD: 0009-0006-8902-0519

AI co-authors: Qwen3.8-Flash (lane D) and DeepSeek-V4.1-Flash (lane C).
