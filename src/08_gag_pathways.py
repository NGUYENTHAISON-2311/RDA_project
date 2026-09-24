"""Glycosaminoglycan pathway readout from KEGG.

Notochordal cells sit in a proteoglycan-rich nucleus pulposus and GAG content is
the functional readout of disc regeneration, so the question is whether these
vesicles carry the machinery that builds GAG chains - not just the proteoglycan
core proteins the atlas already annotates.

The primary output is a membership-and-evidence table per gene, not a p-value.
These KEGG pathways hold 10-38 genes each; tested against a 2,931-symbol
background under g:Profiler's g_SCS correction none of them survives, which would
report "nothing here" about pathways whose genes are demonstrably present and
strongly directional. Fisher with BH is the sensitive complement, and genes that
are absent get explicit rows - the missing hyaluronan synthases are part of the
answer, and a table of hits alone cannot show them.
"""

import sys
from pathlib import Path

import pandas as pd
import requests
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import DATA_REF, RESULTS, bh_fdr  # noqa: E402

KEGG = "https://rest.kegg.jp"
PATHWAYS = {
    "hsa00532": "GAG biosynthesis - chondroitin / dermatan sulfate",
    "hsa00534": "GAG biosynthesis - heparan sulfate / heparin",
    "hsa00533": "GAG biosynthesis - keratan sulfate",
    "hsa00531": "GAG degradation",
    "hsa00520": "Amino sugar & nucleotide sugar metabolism",
    "hsa00920": "Sulfur metabolism (PAPS donor)",
}
CORE_GAG = {"hsa00532", "hsa00534", "hsa00533", "hsa00531"}
GLYCOLYSIS = "hsa00010"

# Curated grouping by biosynthetic step. This is our curation, not a KEGG
# classification - KEGG groups by pathway map, which splits the shared linker
# enzymes across the CS and HS maps and says nothing about where in the chain an
# enzyme acts.
FUNCTIONAL_STEP = {
    "linker initiation": ["XYLT1", "XYLT2", "B4GALT7", "B3GALT6", "B3GAT3", "FAM20B"],
    "CS/DS polymerisation": ["CSGALNACT1", "CSGALNACT2", "CHSY1", "CHSY3", "CHPF", "CHPF2", "DSE", "DSEL"],
    "CS/DS sulfation": ["CHST3", "CHST7", "CHST11", "CHST12", "CHST13", "CHST14", "CHST15", "UST"],
    "HS polymerisation": ["EXT1", "EXT2", "EXTL1", "EXTL2", "EXTL3"],
    "HS modification": ["NDST1", "NDST2", "NDST3", "NDST4", "GLCE", "HS2ST1",
                        "HS6ST1", "HS6ST2", "HS6ST3", "HS3ST1", "HS3ST2", "HS3ST3A1",
                        "HS3ST3B1", "HS3ST4", "HS3ST5", "HS3ST6"],
    "keratan sulfate": ["B3GNT1", "B3GNT2", "B3GNT7", "B4GALT1", "B4GALT4", "CHST1",
                        "CHST2", "CHST4", "CHST5", "CHST6", "ST3GAL1", "ST3GAL2", "FUT8"],
    "hyaluronan synthesis": ["HAS1", "HAS2", "HAS3"],
    "sulfate donor (PAPS)": ["PAPSS1", "PAPSS2", "BPNT1", "BPNT2", "SLC35B2", "SLC35B3"],
}
# Abundant central-carbon enzymes that appear in the amino-sugar map. Their
# presence says nothing about GAG synthesis specifically and they must not be
# counted as GAG machinery.
SHARED_METABOLISM = {"HK1", "HK2", "HK3", "GPI", "PGM1", "PGM2", "PGM3", "GCK",
                     "GALK1", "GALK2", "GALT", "GALE", "GNPDA1", "GNPDA2", "AMDHD2"}


def kegg_gene_sets():
    """pathway -> {symbol}, cached. KEGG returns 'hsa:<entrez>' identifiers."""
    cache = DATA_REF / "kegg_gag.tsv"
    if cache.exists():
        df = pd.read_csv(cache, sep="\t", dtype=str)
    else:
        gi = pd.read_csv(
            DATA_REF / "Homo_sapiens.gene_info.gz", sep="\t", dtype=str,
            compression="gzip", header=0, usecols=[1, 2], names=["GeneID", "Symbol"],
        )
        entrez2symbol = dict(zip(gi.GeneID, gi.Symbol))
        rows = []
        for pid in list(PATHWAYS) + [GLYCOLYSIS]:
            r = requests.get(f"{KEGG}/link/hsa/{pid}", timeout=120)
            r.raise_for_status()
            for line in r.text.strip().splitlines():
                if not line:
                    continue
                entrez = line.split("\t")[1].replace("hsa:", "")
                sym = entrez2symbol.get(entrez)
                if sym:
                    rows.append({"pathway": pid, "entrez": entrez, "symbol": sym})
        df = pd.DataFrame(rows).drop_duplicates()
        df.to_csv(cache, sep="\t", index=False)
    sets = {p: set(g.symbol) for p, g in df.groupby("pathway")}
    for pid, name in PATHWAYS.items():
        print(f"  {pid} {name:46} {len(sets.get(pid, ())):3} genes")
    return sets


def default_step(in_paths):
    """Step for genes the curated map does not name, taken from their pathway.

    Without this the lysosomal glycosidases and sulfatases (GUSB, HEXA/B, NAGLU,
    SGSH, GNS, GLB1) fall through to 'precursor supply' and the degradation arm
    disappears from the rollup entirely.
    """
    if "hsa00531" in in_paths:
        return "degradation"
    if "hsa00920" in in_paths:
        return "sulfate donor (PAPS)"
    return "precursor supply"


def fisher(query, group, background, label, term):
    """One-sided Fisher for group membership in query, vs the EV proteome."""
    gb = group & background
    a = len(gb & query)
    b = len(query) - a
    c = len(gb) - a
    d = len(background) - len(query) - c
    expected = len(query) * len(gb) / len(background) if background else 0
    return {
        "query_set": label, "term": term,
        "observed": a, "expected": round(expected, 2),
        "fold_enrichment": round(a / expected, 2) if expected else None,
        "term_size_in_background": len(gb),
        "query_size": len(query), "background_size": len(background),
        "p_value": stats.fisher_exact([[a, b], [c, d]], alternative="greater")[1],
        "underpowered": len(gb) < 3,
    }


def main():
    print("KEGG gene sets:")
    sets = kegg_gene_sets()
    glyc = sets.get(GLYCOLYSIS, set())

    df = pd.read_csv(RESULTS / "ev_atlas_master.tsv", sep="\t", index_col=0)
    detected = df.dropna(subset=["human_symbol"]).copy()
    retained = detected[detected.retained]
    background = set(retained.human_symbol)
    print(f"\nbackground = retained EV proteome: {len(background)} human symbols")

    step_of = {g: step for step, genes in FUNCTIONAL_STEP.items() for g in genes}

    # One row per gene, including genes never detected, and including curated
    # members that KEGG does not file under these six maps - the hyaluronan
    # synthases sit outside them entirely, and their absence is part of the answer.
    kegg_genes = set().union(*(sets[p] for p in PATHWAYS if p in sets))
    curated_genes = {g for members in FUNCTIONAL_STEP.values() for g in members}
    all_genes = sorted(kegg_genes | curated_genes)
    best = (
        detected.sort_values("n_unique_peptides", ascending=False)
        .drop_duplicates("human_symbol").set_index("human_symbol")
    )

    rows = []
    for g in all_genes:
        in_paths = [p for p in PATHWAYS if g in sets.get(p, ())]
        rec = {
            "symbol": g,
            "source": "kegg" if in_paths else "curated_only",
            "kegg_pathways": ";".join(in_paths),
            "kegg_pathway_names": ";".join(PATHWAYS[p] for p in in_paths),
            "role": ("gag_specific" if (set(in_paths) & CORE_GAG) or not in_paths
                     else "support"),
            "functional_step": step_of.get(g) or default_step(in_paths),
            "shared_metabolism": g in SHARED_METABOLISM or g in glyc,
            "detected": g in best.index,
        }
        if g in best.index:
            t = best.loc[g]
            rec.update({
                "accession": t.name if isinstance(t.name, str) else None,
                "retained": bool(t.retained),
                "n_ev_high": int(t.n_ev_high), "n_sol_high": int(t.n_sol_high),
                "log2FC_EV_SOL": t.log2FC_EV_SOL,
                "log2FC_EV_SOL_unnormalised": t.log2FC_EV_SOL_unnormalised,
                "q_paired_ttest": t.q_paired_ttest,
                "q_detection_fisher": t.q_detection_fisher,
                "class_empirical": t.class_empirical,
                "corona_interpretation": t.corona_interpretation,
                "n_unique_peptides": int(t.n_unique_peptides),
                "n_psms": int(t.n_psms), "coverage_pct": t.coverage_pct,
                "mapping_method": t.mapping_method,
                "mapping_confidence": t.mapping_confidence,
            })
        else:
            rec["retained"] = False
        rows.append(rec)

    genes = pd.DataFrame(rows)
    genes.insert(1, "accession", genes.pop("accession") if "accession" in genes else None)
    genes.to_csv(RESULTS / "gag_pathway_genes.tsv", sep="\t", index=False)

    print("\npathway coverage:")
    print(f"  {'pathway':46} {'total':>5} {'detected':>8} {'retained':>8}")
    for p, name in PATHWAYS.items():
        g = sets.get(p, set())
        print(f"  {name:46} {len(g):5} {len(g & set(best.index)):8} {len(g & background):8}")

    # --- enrichment: pathways and functional steps, against the EV proteome ----
    queries = {
        "retained": background,
        "SOL_enriched": set(retained[retained.class_empirical == "SOL-enriched"].human_symbol),
        "EV_enriched_or_exclusive": set(
            retained[retained.class_empirical.isin(["EV-enriched", "EV-exclusive"])].human_symbol
        ),
    }
    tests = []
    for label, q in queries.items():
        if label == "retained":
            continue
        for p, name in PATHWAYS.items():
            tests.append(fisher(q, sets.get(p, set()), background, label, name))
        for step, members in FUNCTIONAL_STEP.items():
            tests.append(fisher(q, set(members), background, label, f"step: {step}"))
    enr = pd.DataFrame(tests)
    enr["q_value"] = bh_fdr(enr.p_value.values)
    enr.sort_values("p_value").to_csv(RESULTS / "gag_enrichment.tsv", sep="\t", index=False)

    print("\nenrichment vs detected EV proteome (q < 0.05):")
    sig = enr[enr.q_value < 0.05].sort_values("p_value")
    for _, t in sig.iterrows():
        flag = "  [underpowered]" if t.underpowered else ""
        print(f"  {t.query_set:26} {t.term:48} {t.observed:3}/{t.term_size_in_background:<3} "
              f"fold={t.fold_enrichment:6} q={t.q_value:.2e}{flag}")
    if sig.empty:
        print("  none")

    # --- the headline direction ------------------------------------------------
    gag = genes[genes.detected & genes.retained & ~genes.shared_metabolism]
    gag = gag.dropna(subset=["log2FC_EV_SOL"])
    print(f"\nGAG-relevant proteins retained on vesicles (excluding shared metabolism): {len(gag)}")
    print(f"  median log2FC EV/SOL : {gag.log2FC_EV_SOL.median():.2f}")
    print(f"  SOL-enriched         : {(gag.class_empirical == 'SOL-enriched').sum()}")
    print(f"  EV-enriched          : {gag.class_empirical.isin(['EV-enriched', 'EV-exclusive']).sum()}")
    print("\n  strongest evidence (by unique peptides):")
    for _, t in gag.sort_values("n_unique_peptides", ascending=False).head(10).iterrows():
        print(f"    {t.symbol:12} {t.functional_step:22} log2FC={t.log2FC_EV_SOL:7.2f} "
              f"EV {t.n_ev_high}/7  SOL {t.n_sol_high}/7  {int(t.n_unique_peptides):3} unique pep")

    # Rollup by biosynthetic step. This is where the partition splits: the
    # cytosolic sugar-donor and sulfate-donor enzymes sit near zero and are
    # present in every EV prep, while the Golgi-lumenal chain-building and
    # sulfating enzymes sit far to the soluble side. The split is read off the
    # data, not asserted from compartment annotation.
    det = genes[genes.detected]
    roll = det.groupby("functional_step").apply(
        lambda s: pd.Series({
            "genes_in_step": int(len(genes[genes.functional_step == s.name])),
            "detected": int(len(s)),
            "retained": int(s.retained.sum()),
            "median_log2FC": round(float(s.log2FC_EV_SOL.median()), 2),
            "sol_enriched": int((s.class_empirical == "SOL-enriched").sum()),
            "ev_enriched": int(s.class_empirical.isin(["EV-enriched", "EV-exclusive"]).sum()),
            "max_unique_peptides": int(s.n_unique_peptides.max()),
        }), include_groups=False
    ).sort_values("median_log2FC")
    roll.to_csv(RESULTS / "gag_functional_steps.tsv", sep="\t")
    print("\nby biosynthetic step (detected genes, sorted soluble -> vesicle):")
    print(f"  {'step':24} {'det':>4} {'ret':>4} {'medFC':>7} {'SOL':>4} {'EV':>3}")
    for step, t in roll.iterrows():
        print(f"  {step:24} {int(t.detected):4} {int(t.retained):4} "
              f"{t.median_log2FC:7.2f} {int(t.sol_enriched):4} {int(t.ev_enriched):3}")

    absent = genes[~genes.detected & (genes.role == "gag_specific")]
    print(f"\nGAG-specific genes not detected anywhere in the deposit: {len(absent)}")
    print("   ", ", ".join(absent.symbol.tolist()))

    # --- core proteins the machinery acts on -----------------------------------
    pg = detected[detected.matrisome_category.isin(["Proteoglycans"])].copy()
    cols = ["human_symbol", "matrisome_category", "retained", "n_ev_high", "n_sol_high",
            "log2FC_EV_SOL", "q_paired_ttest", "class_empirical", "corona_interpretation",
            "n_unique_peptides"]
    pg[cols].sort_values("log2FC_EV_SOL", ascending=False).to_csv(
        RESULTS / "gag_core_proteins.tsv", sep="\t", index=False)
    print(f"\nproteoglycan core proteins cross-referenced: {len(pg)} "
          f"({int(pg.retained.sum())} retained)")


if __name__ == "__main__":
    main()
