"""Reproducibility filter and the EV-vs-SOL corona axis.

Retained EV proteome = detected with a confident ID in >=4 of 7 EV preps, matching
the source publication. That set is also the enrichment background for step 06.

The EV/SOL contrast is paired within prep. Proteins absent from every SOL run are
not a missing-data problem to impute away - absence is the measurement - so they
are tested on detection counts and flagged EV-exclusive rather than given a ratio.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import DATA_RAW, RESULTS, bh_fdr  # noqa: E402

MIN_EV_PREPS = 4


def main():
    samples = pd.read_csv(DATA_RAW / "samples.tsv", sep="\t")
    log2 = pd.read_csv(RESULTS / "abundance_log2_norm.tsv", sep="\t", index_col=0)
    found = pd.read_csv(RESULTS / "found_status.tsv", sep="\t", index_col=0)
    ev = pd.read_csv(RESULTS / "protein_evidence.tsv", sep="\t", index_col=0)
    high = found.eq("High")

    pairs = (
        samples.pivot(index="prep", columns="fraction", values="f_number")
        .loc[:, ["EV", "SOL"]]
    )
    ev_f, sol_f = pairs.EV.tolist(), pairs.SOL.tolist()
    prep_batch = samples.drop_duplicates("prep").set_index("prep").batch

    out = pd.DataFrame(index=log2.index)
    out["n_ev_high"] = high[ev_f].sum(axis=1)
    out["n_sol_high"] = high[sol_f].sum(axis=1)
    out["n_ev_quantified"] = log2[ev_f].notna().sum(axis=1)
    out["n_sol_quantified"] = log2[sol_f].notna().sum(axis=1)
    out["donor_reproducibility"] = (out.n_ev_high / len(ev_f)).round(3)
    out["mean_ev_log2"] = log2[ev_f].mean(axis=1).round(3)
    out["mean_sol_log2"] = log2[sol_f].mean(axis=1).round(3)

    linear_ev = np.power(2.0, log2[ev_f])
    out["cv_ev_pct"] = (linear_ev.std(axis=1, ddof=1) / linear_ev.mean(axis=1) * 100).round(1)

    batches = sorted(prep_batch.unique())
    for b in batches:
        preps = prep_batch.index[prep_batch == b]
        out[f"n_ev_high_{b}"] = high[pairs.loc[preps, "EV"]].sum(axis=1)
    out["n_batches_detected"] = sum(
        (out[f"n_ev_high_{b}"] > 0).astype(int) for b in batches
    )

    # --- paired EV/SOL contrast within prep ---
    # SOL runs carry 10-40x more total protein than EV runs, so the ratio is only
    # meaningful on composition-normalised values ("what share of each fraction is
    # this protein"). The raw ratio is kept alongside: normalisation shifts every
    # log2FC by a fixed offset and that shift should be visible, not buried.
    diff = pd.DataFrame(
        {p: log2[pairs.loc[p, "EV"]] - log2[pairs.loc[p, "SOL"]] for p in pairs.index}
    )
    raw_log2 = np.log2(pd.read_csv(RESULTS / "abundance_raw.tsv", sep="\t", index_col=0).where(lambda d: d > 0))
    diff_raw = pd.DataFrame(
        {p: raw_log2[pairs.loc[p, "EV"]] - raw_log2[pairs.loc[p, "SOL"]] for p in pairs.index}
    )
    out["n_pairs_usable"] = diff.notna().sum(axis=1)
    out["log2FC_EV_SOL"] = diff.mean(axis=1).round(3)
    out["log2FC_EV_SOL_unnormalised"] = diff_raw.mean(axis=1).round(3)

    usable = diff.notna().sum(axis=1) >= 3
    t_p = pd.Series(np.nan, index=diff.index)
    vals = diff[usable]
    if len(vals):
        res = stats.ttest_1samp(vals.values, 0.0, axis=1, nan_policy="omit")
        t_p[usable] = np.asarray(res.pvalue, dtype=float)
    out["p_paired_ttest"] = t_p
    out["q_paired_ttest"] = bh_fdr(t_p.values)

    # per-batch sign agreement on the ratio
    signs = []
    for b in batches:
        preps = [p for p in pairs.index if prep_batch[p] == b]
        signs.append(np.sign(diff[preps].mean(axis=1)))
    out["batch_consistent"] = (signs[0] == signs[1]) & signs[0].notna() & signs[1].notna()

    # detection-based test, the only valid one for on/off proteins
    n = len(ev_f)
    fisher_p = [
        stats.fisher_exact([[a, n - a], [b, n - b]])[1]
        for a, b in zip(out.n_ev_high, out.n_sol_high)
    ]
    out["p_detection_fisher"] = fisher_p
    out["q_detection_fisher"] = bh_fdr(fisher_p)

    out["ev_exclusive"] = (out.n_sol_high == 0) & (out.n_ev_high >= MIN_EV_PREPS)
    out["sol_exclusive"] = (out.n_ev_high == 0) & (out.n_sol_high >= MIN_EV_PREPS)
    out["retained"] = out.n_ev_high >= MIN_EV_PREPS

    master = ev.drop(columns=[c for c in out.columns if c in ev.columns]).join(out)
    master.to_csv(RESULTS / "ev_sol_contrast.tsv", sep="\t")
    retained = master[master.retained]
    retained.to_csv(RESULTS / "retained_ev_proteome.tsv", sep="\t")

    sig = (out.q_paired_ttest < 0.05) & out.retained
    print(f"retained EV proteome: {len(retained)}/{len(out)} "
          f"(>= {MIN_EV_PREPS}/{len(ev_f)} EV preps with confident ID)")
    print(f"  detected in all 7 EV preps : {(out.retained & (out.n_ev_high == 7)).sum()}")
    print(f"  EV-exclusive (0 SOL runs)  : {out.ev_exclusive.sum()}")
    print(f"  median CV across EV preps  : {retained.cv_ev_pct.median():.1f}%")
    print(f"  batch-consistent direction : {(retained.batch_consistent).sum()}/{len(retained)}")
    print(f"\nEV/SOL contrast (retained, q<0.05): {sig.sum()} significant")
    print(f"  EV-enriched  : {(sig & (out.log2FC_EV_SOL > 1)).sum()} (log2FC > 1)")
    print(f"  SOL-enriched : {(sig & (out.log2FC_EV_SOL < -1)).sum()} (log2FC < -1)")


if __name__ == "__main__":
    main()
