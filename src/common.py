"""Shared paths and small helpers for the PXD072052 NC-EV atlas pipeline."""

from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
DATA_RAW = ROOT / "data" / "raw"
DATA_REF = ROOT / "data" / "ref"
RESULTS = ROOT / "results"
QC = RESULTS / "qc"
REPORT = ROOT / "report"

for _d in (DATA_RAW, DATA_REF, RESULTS, QC, REPORT):
    _d.mkdir(parents=True, exist_ok=True)

EV_MARKERS = ["CD9", "CD63", "CD81", "CD82", "PDCD6IP", "TSG101", "SDCBP", "FLOT1", "FLOT2"]
SOLUBLE_MARKERS = ["ALB", "HBB", "HBA1", "FGA", "FGB", "TF", "APOA1"]


def bh_fdr(pvals):
    """Benjamini-Hochberg adjusted p-values; NaN in, NaN out."""
    p = np.asarray(pvals, dtype=float)
    q = np.full(p.shape, np.nan)
    ok = ~np.isnan(p)
    n = int(ok.sum())
    if n == 0:
        return q
    idx = np.flatnonzero(ok)[np.argsort(p[ok], kind="stable")]
    ranked = p[idx] * n / np.arange(1, n + 1)
    q[idx] = np.minimum.accumulate(ranked[::-1])[::-1].clip(max=1.0)
    return q
