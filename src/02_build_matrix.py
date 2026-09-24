"""Parse the PD Proteins.txt export into abundance / detection matrices + QC.

Guards the F-number -> sample assignment: PD numbers files in consensus-workflow
input order, which the pdStudy does not explicitly confirm, so we test it against
marker biology before anything downstream trusts it.
"""

import re
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import DATA_RAW, EV_MARKERS, QC, RESULTS, SOLUBLE_MARKERS  # noqa: E402

FRUNS = [f"F{i}" for i in range(1, 15)]


def load():
    df = pd.read_csv(DATA_RAW / "Proteins.txt", sep="\t", dtype=str, encoding="utf-8-sig")
    df.columns = [c.strip().strip('"') for c in df.columns]
    for c in df.columns:
        df[c] = df[c].str.strip().str.strip('"').str.strip()
    return df


def main():
    df = load()
    samples = pd.read_csv(DATA_RAW / "samples.tsv", sep="\t")
    print(f"{len(df)} protein groups x {len(samples)} runs")

    acc = df["Accession"]
    abundance = pd.DataFrame(
        {f: pd.to_numeric(df[f"Abundance: {f}: Sample"], errors="coerce") for f in FRUNS}
    ).set_index(acc)
    found = pd.DataFrame(
        {f: df[f"Found in Sample: {f}: Sample"] for f in FRUNS}
    ).set_index(acc)

    # "High" = confident identification; "Peak Found" = quantified via matching only.
    high = found.eq("High")
    print(f"detection sanity: {(high.sum(axis=1) == 14).sum()} proteins High in all 14 runs, "
          f"{(high.sum(axis=1) == 1).sum()} in exactly one")
    if (high.sum(axis=1) == 14).sum() == 0:
        raise SystemExit("FAIL: no proteins High anywhere - value stripping is broken")

    ev_f = samples.loc[samples.fraction == "EV", "f_number"].tolist()
    sol_f = samples.loc[samples.fraction == "SOL", "f_number"].tolist()

    # --- median-ratio normalisation on log2, anchored on complete-case proteins ---
    log2 = np.log2(abundance.where(abundance > 0))
    complete = log2.dropna()
    offsets = (complete.sub(complete.mean(axis=1), axis=0)).median()
    log2_norm = log2.sub(offsets, axis=1)
    print(f"normalisation anchored on {len(complete)} complete-case proteins; "
          f"offsets {offsets.min():.2f}..{offsets.max():.2f}")

    # --- F-mapping gate: EV markers must favour EV runs, soluble markers SOL runs ---
    symbol = df["Description"].str.extract(r"GN=(\S+)", expand=False)
    sym_idx = pd.Series(symbol.values, index=acc.values)
    verdict = []
    for name, markers, expect in [
        ("EV", EV_MARKERS, "EV"),
        ("soluble", SOLUBLE_MARKERS, "SOL"),
    ]:
        hits = sym_idx[sym_idx.isin(markers)].index
        d = log2_norm.loc[hits]
        delta = (d[ev_f].mean(axis=1) - d[sol_f].mean(axis=1)).mean()
        got = "EV" if delta > 0 else "SOL"
        verdict.append(got == expect)
        print(f"F-mapping check: {len(hits)} {name} markers "
              f"{sorted(sym_idx[hits].unique())} mean log2(EV/SOL)={delta:+.2f} -> favours {got}")
    if not all(verdict):
        raise SystemExit("FAIL: F-number to sample mapping contradicts marker biology")
    print("F-mapping gate PASSED")

    # --- per-protein evidence and protein-group ambiguity ---
    num = lambda c: pd.to_numeric(df[c], errors="coerce")  # noqa: E731
    ev = pd.DataFrame(
        {
            "accession": acc,
            "description": df["Description"],
            "pig_symbol_from_export": symbol,
            "master_status": df["Master"],
            "fdr_confidence": df["Protein FDR Confidence: Combined"],
            "exp_q_value": num("Exp. q-value: Combined"),
            "sum_pep_score": num("Sum PEP Score"),
            "coverage_pct": num("Coverage [%]"),
            "n_peptides": num("# Peptides"),
            "n_unique_peptides": num("# Unique Peptides"),
            "n_razor_peptides": num("# Razor Peptides"),
            "n_psms": num("# PSMs"),
            "n_protein_groups": num("# Protein Groups"),
            "mw_kda": num("MW [kDa]"),
        }
    ).set_index("accession")

    ev["unique_peptide_fraction"] = (ev.n_unique_peptides / ev.n_peptides).round(4)
    ev["razor_only"] = ev.n_unique_peptides == 0
    ev["low_unique_evidence"] = ev.n_unique_peptides < 2
    dup = ev.pig_symbol_from_export.dropna().duplicated(keep=False)
    dup_syms = set(ev.pig_symbol_from_export.dropna()[dup])
    ev["split_group_symbol"] = ev.pig_symbol_from_export.isin(dup_syms)
    print(f"ambiguity: {ev.razor_only.sum()} razor-only, {ev.low_unique_evidence.sum()} <2 unique "
          f"peptides, {len(dup_syms)} symbols split across >1 group")

    ev["n_ev_high"] = high[ev_f].sum(axis=1)
    ev["n_sol_high"] = high[sol_f].sum(axis=1)
    ev["n_ev_peak_found"] = found[ev_f].eq("Peak Found").sum(axis=1)
    ev["n_sol_peak_found"] = found[sol_f].eq("Peak Found").sum(axis=1)

    abundance.to_csv(RESULTS / "abundance_raw.tsv", sep="\t")
    log2_norm.round(4).to_csv(RESULTS / "abundance_log2_norm.tsv", sep="\t")
    found.to_csv(RESULTS / "found_status.tsv", sep="\t")
    ev.to_csv(RESULTS / "protein_evidence.tsv", sep="\t")

    qc_plots(log2_norm, high, samples, ev_f, sol_f)
    print(f"\nwrote matrices to {RESULTS}")


def qc_plots(log2_norm, high, samples, ev_f, sol_f):
    colors = ["#2563eb" if f in ev_f else "#ea580c" for f in FRUNS]
    labels = [
        f"{r.f_number} {r.prep} {r.fraction}" for r in samples.itertuples()
    ]

    fig, axes = plt.subplots(1, 3, figsize=(17, 5))
    axes[0].bar(labels, high.sum().values, color=colors)
    axes[0].set_title("Confident IDs per run")
    axes[0].set_ylabel("proteins with High ID")
    axes[0].tick_params(axis="x", rotation=90)

    axes[1].boxplot(
        [log2_norm[f].dropna().values for f in FRUNS], tick_labels=labels, showfliers=False
    )
    axes[1].set_title("Normalised log2 abundance")
    axes[1].tick_params(axis="x", rotation=90)

    mat = log2_norm.dropna()
    centered = mat.sub(mat.mean(axis=1), axis=0).values.T
    u, s, _ = np.linalg.svd(centered, full_matrices=False)
    pcs = u[:, :2] * s[:2]
    var = (s**2 / (s**2).sum())[:2] * 100
    for i, f in enumerate(FRUNS):
        axes[2].scatter(pcs[i, 0], pcs[i, 1], color=colors[i], s=60)
        axes[2].annotate(labels[i].split()[1], (pcs[i, 0], pcs[i, 1]), fontsize=7)
    axes[2].set_title(f"PCA of runs (blue=EV, orange=SOL)\n{len(mat)} complete-case proteins")
    axes[2].set_xlabel(f"PC1 {var[0]:.0f}%")
    axes[2].set_ylabel(f"PC2 {var[1]:.0f}%")

    fig.tight_layout()
    fig.savefig(QC / "run_qc.png", dpi=140)
    plt.close(fig)


if __name__ == "__main__":
    main()
