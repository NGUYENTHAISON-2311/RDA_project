"""Summary figures and the compact data payload used to build the HTML report."""

import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import DATA_RAW, QC, REPORT, RESULTS  # noqa: E402

EV_PANEL = ["CD63", "CD9", "CD81", "TSG101", "PDCD6IP", "FLOT1", "FLOT2",
            "CHMP4B", "VPS4A", "RAB27A", "CAV1", "AHNAK"]


PATHWAY_SHORT = {
    "GAG biosynthesis - chondroitin / dermatan sulfate": "Chondroitin / dermatan sulfate",
    "GAG biosynthesis - heparan sulfate / heparin": "Heparan sulfate / heparin",
    "GAG biosynthesis - keratan sulfate": "Keratan sulfate",
    "GAG degradation": "GAG degradation",
    "Amino sugar & nucleotide sugar metabolism": "UDP-sugar precursor supply",
    "Sulfur metabolism (PAPS donor)": "Sulfate donor (PAPS)",
}


def gag_payload(gag, steps, enr, df):
    det = gag[gag.detected]
    coverage = []
    for name, short in PATHWAYS_ORDER:
        members = gag[gag.kegg_pathway_names.fillna("").str.contains(name, regex=False)]
        coverage.append({
            "pathway": short,
            "total": int(len(members)),
            "detected": int(members.detected.sum()),
            "retained": int(members.retained.sum()),
        })

    # every detected GAG-specific gene, ordered soluble -> vesicle for the figure
    plot = det[(det.role == "gag_specific") & ~det.shared_metabolism]
    plot = plot.dropna(subset=["log2FC_EV_SOL"]).sort_values("log2FC_EV_SOL")
    core = df[df.matrisome_category == "Proteoglycans"].dropna(subset=["log2FC_EV_SOL"])
    core = core.sort_values("log2FC_EV_SOL", ascending=False).drop_duplicates("human_symbol")

    return {
        "coverage": coverage,
        "steps": [
            {"step": t.functional_step, "detected": int(t.detected), "retained": int(t.retained),
             "median_log2fc": float(t.median_log2FC), "sol": int(t.sol_enriched),
             "ev": int(t.ev_enriched), "in_step": int(t.genes_in_step)}
            for _, t in steps.iterrows()
        ],
        "genes": [
            {"symbol": t.symbol, "log2fc": round(float(t.log2FC_EV_SOL), 2),
             "step": t.functional_step, "retained": bool(t.retained),
             "ev": int(t.n_ev_high), "sol": int(t.n_sol_high),
             "pep": int(t.n_unique_peptides),
             "q": None if pd.isna(t.q_paired_ttest) else float(t.q_paired_ttest)}
            for _, t in plot.iterrows()
        ],
        "enrichment": [
            {"query": t.query_set, "term": t.term, "observed": int(t.observed),
             "in_bg": int(t.term_size_in_background), "fold": t.fold_enrichment,
             "q": float(t.q_value), "underpowered": bool(t.underpowered)}
            for _, t in enr[enr.q_value < 0.05].sort_values("p_value").iterrows()
        ],
        "absent": sorted(gag[~gag.detected & (gag.role == "gag_specific")].symbol),
        "core_proteins": [
            {"symbol": t.human_symbol, "log2fc": round(float(t.log2FC_EV_SOL), 2),
             "retained": bool(t.retained), "ev": int(t.n_ev_high), "sol": int(t.n_sol_high)}
            for _, t in core.iterrows()
        ],
        "summary": {
            "retained_non_shared": int(len(
                det[det.retained & ~det.shared_metabolism].dropna(subset=["log2FC_EV_SOL"]))),
            "median_log2fc": round(float(
                det[det.retained & ~det.shared_metabolism].log2FC_EV_SOL.median()), 2),
            "sol_enriched": int((det[det.retained & ~det.shared_metabolism]
                                 .class_empirical == "SOL-enriched").sum()),
            "ev_enriched": int(det[det.retained & ~det.shared_metabolism]
                               .class_empirical.isin(["EV-enriched", "EV-exclusive"]).sum()),
        },
    }


PATHWAYS_ORDER = [
    ("GAG biosynthesis - chondroitin / dermatan sulfate", "Chondroitin / dermatan sulfate"),
    ("GAG biosynthesis - heparan sulfate / heparin", "Heparan sulfate / heparin"),
    ("GAG biosynthesis - keratan sulfate", "Keratan sulfate"),
    ("GAG degradation", "GAG degradation"),
    ("Amino sugar & nucleotide sugar metabolism", "UDP-sugar precursor supply"),
    ("Sulfur metabolism (PAPS donor)", "Sulfate donor (PAPS)"),
]


def main():
    df = pd.read_csv(RESULTS / "ev_atlas_master.tsv", sep="\t", index_col=0)
    mat = pd.read_csv(RESULTS / "enrichment_matrisome.tsv", sep="\t")
    go = pd.read_csv(RESULTS / "enrichment_go.tsv", sep="\t")
    reac = pd.read_csv(RESULTS / "enrichment_reactome.tsv", sep="\t")
    gag = pd.read_csv(RESULTS / "gag_pathway_genes.tsv", sep="\t")
    gag_steps = pd.read_csv(RESULTS / "gag_functional_steps.tsv", sep="\t")
    gag_enr = pd.read_csv(RESULTS / "gag_enrichment.tsv", sep="\t")
    samples = pd.read_csv(DATA_RAW / "samples.tsv", sep="\t")
    omap = pd.read_csv(RESULTS / "ortholog_map.tsv", sep="\t")
    r = df[df.retained]

    fig, ax = plt.subplots(2, 2, figsize=(14, 10))

    v = r.dropna(subset=["log2FC_EV_SOL", "q_paired_ttest"])
    ax[0, 0].scatter(v.log2FC_EV_SOL, -np.log10(v.q_paired_ttest.clip(lower=1e-12)),
                     s=6, c="#cbd5e1", label="other")
    for mask, color, lab in [
        (v.is_ev_associated, "#2563eb", "curated EV marker"),
        (v.is_corona_candidate, "#ea580c", "corona candidate"),
    ]:
        s = v[mask]
        ax[0, 0].scatter(s.log2FC_EV_SOL, -np.log10(s.q_paired_ttest.clip(lower=1e-12)),
                         s=18, c=color, label=lab)
    ax[0, 0].axvline(0, color="#94a3b8", lw=0.8)
    ax[0, 0].axhline(-np.log10(0.05), color="#94a3b8", lw=0.8, ls="--")
    ax[0, 0].set_xlabel("log2 EV / SOL (composition-normalised)")
    ax[0, 0].set_ylabel("-log10 q (paired t-test)")
    ax[0, 0].set_title("EV vs soluble fraction, paired across 7 preps")
    ax[0, 0].legend(fontsize=8)

    panel = (r[r.human_symbol.isin(EV_PANEL)]
             .groupby("human_symbol").log2FC_EV_SOL.max().sort_values())
    ax[0, 1].barh(panel.index, panel.values, color="#2563eb")
    ax[0, 1].axvline(0, color="#475569", lw=0.8)
    ax[0, 1].set_xlabel("log2 EV / SOL")
    ax[0, 1].set_title("Canonical EV markers (validation)")

    counts = r.n_ev_high.value_counts().sort_index()
    ax[1, 0].bar(counts.index, counts.values, color="#0d9488")
    ax[1, 0].set_xlabel("EV preps with confident ID (of 7)")
    ax[1, 0].set_ylabel("protein groups")
    ax[1, 0].set_title(f"Donor reproducibility, retained set (n={len(r)})")

    sig = mat[(mat.q_value < 0.05) & ~mat.circular & (mat.query_set == "empirically_SOL_enriched")]
    sig = sig.sort_values("fold_enrichment")
    ax[1, 1].barh(sig.term, sig.fold_enrichment, color="#ea580c")
    ax[1, 1].set_xlabel("fold enrichment vs detected EV proteome")
    ax[1, 1].set_title("Matrisome enrichment of the SOL-enriched set")

    fig.tight_layout()
    fig.savefig(QC / "atlas_summary.png", dpi=140)
    plt.close(fig)

    def top(frame, qs, src=None, n=8):
        s = frame[frame.query_set == qs]
        if src:
            s = s[s.source == src]
        s = s.sort_values("p_value").head(n)
        return [
            {"name": t["name"], "fold": float(t.fold_enrichment),
             "p": float(t.p_value), "size": int(t.intersection_size)}
            for _, t in s.iterrows()
        ]

    payload = {
        "dataset": {
            "accession": "PXD072052",
            "species": "Sus scrofa domesticus",
            "runs": len(samples),
            "preps": int(samples.prep.nunique()),
            "batches": sorted(samples.batch.unique()),
            "protein_groups": int(len(df)),
            "retained": int(len(r)),
            "retained_symbols": int(r.human_symbol.nunique()),
        },
        "mapping": {
            "total": int(len(omap)),
            "mapped": int(omap.human_symbol.notna().sum()),
            "deleted_from_uniprot": int(
                len((Path("data/ref/deleted_accessions.txt").read_text().split()))
            ),
            "methods": omap.mapping_method.value_counts().to_dict(),
        },
        "layers": {
            lab: {"n": int(r[col].sum()),
                  "ev_side": int((r[col] & r.class_empirical.isin(["EV-enriched", "EV-exclusive"])).sum())}
            for lab, col in [
                ("EV-associated (curated)", "is_ev_associated"),
                ("ECM & proteoglycan", "is_ecm_proteoglycan"),
                ("Matricellular", "is_matricellular"),
                ("Growth factor / cytokine regulator", "is_growth_factor_regulator"),
                ("Surface / corona-associated", "is_surface_corona"),
                ("Intracellular / soluble co-isolate", "is_intracellular_soluble"),
                ("Corona candidate", "is_corona_candidate"),
            ]
        },
        "reproducibility": {int(k): int(v) for k, v in r.n_ev_high.value_counts().sort_index().items()},
        "run_totals": {
            f: round(float(v) / 1e9, 1)
            for f, v in pd.read_csv(RESULTS / "abundance_raw.tsv", sep="\t", index_col=0).sum().items()
        },
        "batch_consistent": int(r.batch_consistent.sum()),
        "normalisation_shift": round(float(
            r.log2FC_EV_SOL.median() - r.log2FC_EV_SOL_unnormalised.median()
        ), 2),
        "empirical": r.class_empirical.value_counts().to_dict(),
        "agreement": r.class_agreement.value_counts().to_dict(),
        "corona_interpretation": df.corona_interpretation.value_counts().to_dict(),
        "ev_markers": [
            {"symbol": s, "log2fc": float(v)}
            for s, v in r[r.is_ev_associated].groupby("human_symbol").log2FC_EV_SOL.max()
            .sort_values(ascending=False).items()
        ],
        "top_corona": [
            {"symbol": t.human_symbol, "log2fc": float(t.log2FC_EV_SOL),
             "category": (t.matrisome_category if isinstance(t.matrisome_category, str) else None),
             "n_ev": int(t.n_ev_high), "n_sol": int(t.n_sol_high)}
            for _, t in r[r.corona_interpretation == "EV-enriched corona"]
            .sort_values("log2FC_EV_SOL", ascending=False).head(15).iterrows()
        ],
        "soluble_dominant": [
            {"symbol": t.human_symbol, "log2fc": float(t.log2FC_EV_SOL),
             "category": (t.matrisome_category if isinstance(t.matrisome_category, str) else None),
             "n_ev": int(t.n_ev_high), "n_sol": int(t.n_sol_high)}
            for _, t in r[r.corona_interpretation == "EV-bound, soluble-dominant"]
            .sort_values("log2FC_EV_SOL").head(15).iterrows()
        ],
        "volcano": [
            {"x": round(float(t.log2FC_EV_SOL), 2),
             "y": round(float(-np.log10(max(t.q_paired_ttest, 1e-12))), 2),
             "s": t.human_symbol if isinstance(t.human_symbol, str) else "",
             "k": ("ev" if t.is_ev_associated else "corona" if t.is_corona_candidate else "other")}
            for _, t in v.iterrows()
        ],
        "enrichment": {
            "reactome_ev_enriched": top(reac, "empirically_EV_enriched"),
            "reactome_corona": top(reac, "corona_candidate"),
            "go_ev_enriched": top(go, "empirically_EV_enriched", "GO:CC"),
            "go_sol_enriched": top(go, "empirically_SOL_enriched", "GO:MF"),
            "matrisome_noncircular": [
                {"query": t.query_set, "term": t.term, "fold": float(t.fold_enrichment),
                 "q": float(t.q_value), "observed": int(t.observed), "query_size": int(t.query_size)}
                for _, t in mat[(mat.q_value < 0.05) & ~mat.circular]
                .sort_values("p_value").head(12).iterrows()
            ],
        },
        "background_check": {
            # Compartment terms only. The BP term "extracellular exosome biogenesis"
            # is a small, specific ESCRT term and is legitimately enriched in the
            # curated marker set - it is a result, not a background failure.
            "max_exosome_cc_fold": float(
                go[(go.source == "GO:CC")
                   & go.name.str.contains("exosome|extracellular vesicle", case=False)]
                .fold_enrichment.max()
            ),
            "biogenesis": [
                {"name": t["name"], "fold": float(t.fold_enrichment),
                 "p": float(t.p_value), "size": int(t.intersection_size),
                 "term_size": int(t.term_size), "query": t.query_set}
                for _, t in go[go.name.str.contains("vesicle biogenesis|exosome biogenesis",
                                                    case=False)].iterrows()
            ],
        },
        "samples": samples.to_dict("records"),
        "gag": gag_payload(gag, gag_steps, gag_enr, df),
    }
    (REPORT / "atlas_data.json").write_text(json.dumps(payload, indent=1))
    # the published page loads this as a script: artifact CSP blocks fetch/XHR
    (REPORT / "atlas_data.js").write_text(
        "window.ATLAS = " + json.dumps(payload, separators=(",", ":")) + ";"
    )
    print(f"wrote {REPORT / 'atlas_data.json'} and {QC / 'atlas_summary.png'}")
    print(f"  volcano points: {len(payload['volcano'])}, "
          f"EV markers: {len(payload['ev_markers'])}, "
          f"top corona: {len(payload['top_corona'])}")


if __name__ == "__main__":
    main()
