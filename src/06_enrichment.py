"""Enrichment against the detected EV proteome, never the whole human proteome.

Using all ~20k human genes as background would report "extracellular exosome" and
"extracellular matrix" as enriched in any EV prep, which is circular. The
background here is the retained EV proteome itself, so a term is only enriched if
it stands out among proteins this preparation actually contains.

GO + Reactome go through g:Profiler with a custom background; the matrisome is not
a g:Profiler source, so it is tested directly by Fisher exact against that same
background.
"""

import sys
from pathlib import Path

import pandas as pd
import requests
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import RESULTS, bh_fdr  # noqa: E402

GPROFILER = "https://biit.cs.ut.ee/gprofiler/api/gost/profile/"
# KEGG is included for completeness, but g_SCS is conservative on small pathways
# against a background this size - the glycosaminoglycan maps, for instance, do
# not survive it despite their genes being clearly present and directional.
# Step 08 runs the sensitive BH-corrected Fisher complement for those.
SOURCES = ["GO:BP", "GO:CC", "GO:MF", "REAC", "KEGG"]


def gprofiler(query, background, label):
    r = requests.post(
        GPROFILER,
        json={
            "organism": "hsapiens",
            "query": sorted(query),
            "sources": SOURCES,
            "user_threshold": 0.05,
            "significance_threshold_method": "g_SCS",
            "domain_scope": "custom_annotated",
            "background": sorted(background),
            "no_evidences": True,
            "no_iea": False,
        },
        headers={"User-Agent": "NC-EV-atlas/PXD072052"},
        timeout=600,
    )
    r.raise_for_status()
    res = r.json().get("result", [])
    if not res:
        return pd.DataFrame()
    df = pd.DataFrame(res)
    keep = ["source", "native", "name", "p_value", "term_size", "query_size",
            "intersection_size", "effective_domain_size"]
    df = df[[c for c in keep if c in df.columns]].copy()
    df["fold_enrichment"] = (
        (df.intersection_size / df.query_size)
        / (df.term_size / df.effective_domain_size)
    ).round(2)
    df["query_set"] = label
    return df.sort_values("p_value")


# Query sets whose definition already involves matrisome membership. Testing them
# for matrisome enrichment is circular, so results are kept but marked, never
# headlined - the informative tests are the ones defined by the EV/SOL data alone.
MATRISOME_CIRCULAR = {
    "corona_candidate", "ecm_proteoglycan", "matricellular",
    "growth_factor_cytokine_regulator", "corona_and_EV_enriched",
}


def matrisome_enrichment(df, query_mask, background_mask, label):
    rows = []
    bg = df[background_mask]
    q = df[query_mask]
    n_bg, n_q = len(bg), len(q)
    for col in ["matrisome_division", "matrisome_category"]:
        for term in bg[col].dropna().unique():
            a = int((q[col] == term).sum())
            b = n_q - a
            c = int((bg[col] == term).sum()) - a
            d = n_bg - n_q - c
            if a == 0:
                continue
            p = stats.fisher_exact([[a, b], [c, d]], alternative="greater")[1]
            expected = n_q * (a + c) / n_bg
            rows.append(
                {
                    "query_set": label, "source": col, "term": term,
                    "observed": a, "expected": round(expected, 1),
                    "fold_enrichment": round(a / expected, 2) if expected else None,
                    "term_size_in_background": a + c,
                    "query_size": n_q, "background_size": n_bg, "p_value": p,
                }
            )
    out = pd.DataFrame(rows)
    if len(out):
        out["q_value"] = bh_fdr(out.p_value.values)
        out["circular"] = label in MATRISOME_CIRCULAR
    return out


def main():
    df = pd.read_csv(RESULTS / "ev_atlas_master.tsv", sep="\t", index_col=0)
    retained = df[df.retained & df.human_symbol.notna()]
    background = set(retained.human_symbol)
    print(f"background = retained EV proteome: {len(background)} human symbols "
          f"(from {len(retained)} protein groups)")

    ev_side = retained.class_empirical.isin(["EV-enriched", "EV-exclusive"])
    query_sets = {
        "corona_candidate": retained.is_corona_candidate,
        "ecm_proteoglycan": retained.is_ecm_proteoglycan,
        "matricellular": retained.is_matricellular,
        "growth_factor_cytokine_regulator": retained.is_growth_factor_regulator,
        "surface_corona_associated": retained.is_surface_corona,
        "ev_associated": retained.is_ev_associated,
        "intracellular_soluble": retained.is_intracellular_soluble,
        "empirically_EV_enriched": ev_side,
        "empirically_SOL_enriched": retained.class_empirical.eq("SOL-enriched"),
        "corona_and_EV_enriched": retained.is_corona_candidate & ev_side,
    }

    go_frames, mat_frames = [], []
    for label, mask in query_sets.items():
        genes = set(retained.loc[mask, "human_symbol"].dropna())
        if len(genes) < 5:
            print(f"  {label:36} skipped ({len(genes)} genes)")
            continue
        try:
            res = gprofiler(genes, background, label)
        except Exception as e:
            print(f"  {label:36} g:Profiler failed: {e}")
            res = pd.DataFrame()
        go_frames.append(res)
        mat_frames.append(
            matrisome_enrichment(retained, mask, retained.retained.astype(bool), label)
        )
        top = res.iloc[0]["name"] if len(res) else "-"
        print(f"  {label:36} {len(genes):5} genes -> {len(res):4} terms   top: {top[:52]}")

    go = pd.concat([f for f in go_frames if len(f)], ignore_index=True) if any(len(f) for f in go_frames) else pd.DataFrame()
    if len(go):
        go[go.source.str.startswith("GO")].to_csv(RESULTS / "enrichment_go.tsv", sep="\t", index=False)
        go[go.source == "REAC"].to_csv(RESULTS / "enrichment_reactome.tsv", sep="\t", index=False)
        go[go.source == "KEGG"].to_csv(RESULTS / "enrichment_kegg.tsv", sep="\t", index=False)
    mat = pd.concat([f for f in mat_frames if len(f)], ignore_index=True)
    mat.to_csv(RESULTS / "enrichment_matrisome.tsv", sep="\t", index=False)

    print(f"\nGO terms: {len(go[go.source.str.startswith('GO')]) if len(go) else 0}, "
          f"Reactome: {len(go[go.source == 'REAC']) if len(go) else 0}, "
          f"matrisome tests: {len(mat)}")

    if len(go):
        print("\ntop Reactome (corona_and_EV_enriched):")
        sub = go[(go.source == "REAC") & (go.query_set == "corona_and_EV_enriched")].head(10)
        for _, r in sub.iterrows():
            print(f"  {r['name'][:62]:62} fold={r.fold_enrichment:5} p={r.p_value:.2e}")

    print("\nmatrisome enrichment, non-circular query sets only (q<0.05):")
    sig = mat[(mat.q_value < 0.05) & ~mat.circular].sort_values("p_value")
    for _, r in sig.head(12).iterrows():
        print(f"  {r.query_set:34} {r.term:24} {r.observed:4}/{r.query_size:<5} "
              f"fold={r.fold_enrichment:5} q={r.q_value:.2e}")
    print(f"  ({(mat.circular & (mat.q_value < 0.05)).sum()} further hits suppressed as circular)")

    # Sanity: with an EV-proteome background, generic exosome terms must collapse in
    # effect size. Against the whole human proteome they would run fold >5; anything
    # like that here would mean the custom background silently failed to apply.
    if len(go):
        exo = go[(go.source == "GO:CC")
                 & go.name.str.contains("exosome|extracellular vesicle", case=False)]
        print(f"\nbackground sanity check - generic EV/exosome CC terms: {len(exo)} significant, "
              f"max fold enrichment {exo.fold_enrichment.max() if len(exo) else 0}")
        if len(exo) and exo.fold_enrichment.max() > 4:
            print("  WARNING: exosome terms strongly enriched - custom background may not have applied")


if __name__ == "__main__":
    main()
