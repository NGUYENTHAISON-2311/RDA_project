"""Annotate the atlas into six biological layers, twice over.

The layers are deliberately non-exclusive - a proteoglycan can also be a growth
factor reservoir - so each protein carries independent flags plus one primary
class resolved by precedence.

Two opinions are produced side by side:
  class_knowledge  from curated lists and GO (what this protein is)
  class_empirical  from the paired EV/SOL contrast (where it actually partitions)
Disagreements are the point, not noise: an "intracellular" protein that is
strongly EV-enriched is either real cargo or a corona passenger, and this study's
thesis says not to discard it by reflex.
"""

import gzip
import io
import re
import sys
import tarfile
from pathlib import Path

import pandas as pd
import requests

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import DATA_REF, RESULTS  # noqa: E402

MATRISOME_TGZ = (
    "https://raw.githubusercontent.com/Matrisome/matrisome-annotator/HEAD/"
    "MatrisomeAnnotatorLinuxV1.tgz"
)
GOA = "http://current.geneontology.org/annotations/goa_human.gaf.gz"
OBO = "http://current.geneontology.org/ontology/go-basic.obo"

GO_EV = {"GO:0070062": "extracellular exosome", "GO:1903561": "extracellular vesicle"}
GO_SURFACE = {"GO:0009986": "cell surface", "GO:0005886": "plasma membrane",
              "GO:0005576": "extracellular region", "GO:0031012": "extracellular matrix"}
GO_GROWTH = {"GO:0008083": "growth factor activity", "GO:0005125": "cytokine activity",
             "GO:0019838": "growth factor binding", "GO:0050431": "TGF-beta binding",
             "GO:0036122": "BMP binding", "GO:0005160": "TGF-beta receptor binding"}
# Deliberately specific organelles. Broad terms (cytosol, ER) propagate to most of
# the proteome and would make "intracellular contaminant" the default verdict.
GO_INTRACELLULAR = {"GO:0005634": "nucleus", "GO:0005840": "ribosome",
                    "GO:0005739": "mitochondrion", "GO:0000502": "proteasome complex",
                    "GO:0005681": "spliceosomal complex", "GO:0005730": "nucleolus",
                    "GO:0005694": "chromosome"}

EV_MARKERS = {
    "CD9", "CD63", "CD81", "CD82", "CD37", "CD53", "TSPAN8", "TSPAN14",
    "PDCD6IP", "TSG101", "VPS4A", "VPS4B", "CHMP1A", "CHMP2A", "CHMP4B", "VTA1",
    "SDCBP", "FLOT1", "FLOT2", "CAV1", "ANXA11", "RAB7A", "RAB5A", "RAB11A",
    "RAB27A", "RAB27B", "ARF6", "EHD1", "EHD4", "SDC1", "SDC4", "AHNAK",
    "HSPA8", "HSP90AA1", "HSP90AB1", "ACTB", "GAPDH", "ITGB1", "ITGA6", "MFGE8",
}
MATRICELLULAR = {
    "THBS1", "THBS2", "THBS3", "THBS4", "SPARC", "SPARCL1", "TNC", "TNR", "TNXB",
    "TNN", "POSTN", "CCN1", "CCN2", "CCN3", "CCN4", "CCN5", "CCN6",
    "CYR61", "CTGF", "NOV", "WISP1", "WISP2", "WISP3",
    "SPP1", "FBLN1", "FBLN2", "FBLN5", "LGALS1", "LGALS3", "EDIL3", "TGFBI",
}
BLOOD = {
    "ALB", "HBB", "HBA1", "HBA2", "HBD", "FGA", "FGB", "FGG", "TF", "HP", "HPX",
    "APOA1", "APOA2", "APOB", "APOE", "C3", "C4A", "C4B", "C1QA", "C1QB", "C1QC",
    "SERPINA1", "A2M", "AHSG", "TTR", "ORM1", "PLG", "F2",
}


def download(url, dest, binary=True):
    if not dest.exists():
        print(f"  downloading {dest.name} ...")
        r = requests.get(url, timeout=1800)
        r.raise_for_status()
        dest.write_bytes(r.content if binary else r.text.encode())
    return dest


def load_matrisome():
    dest = DATA_REF / "HsMatrisomeAnnotations.csv"
    if not dest.exists():
        r = requests.get(MATRISOME_TGZ, timeout=300)
        r.raise_for_status()
        t = tarfile.open(fileobj=io.BytesIO(r.content))
        # the archive is Mac-made: skip the AppleDouble "._" resource forks
        member = next(
            m for m in t.getmembers()
            if m.name.endswith("HsMatrisomeAnnotations.csv") and "/._" not in m.name
        )
        dest.write_bytes(t.extractfile(member).read())
    df = pd.read_csv(dest, encoding="latin-1")
    print(f"  matrisome: {len(df)} genes, {df.Division.nunique()} divisions, "
          f"{df.Category.nunique()} categories")
    return df


def symbol_normalizer():
    """Old symbol -> current NCBI symbol (the matrisome list predates CCN renaming)."""
    gi = pd.read_csv(
        DATA_REF / "Homo_sapiens.gene_info.gz", sep="\t", dtype=str, compression="gzip",
        header=0, usecols=[2, 4], names=["Symbol", "Synonyms"],
    )
    current = set(gi.Symbol)
    syn2cur = {}
    for sym, syns in zip(gi.Symbol, gi.Synonyms.fillna("-")):
        for s in syns.split("|"):
            if s != "-" and s not in current:
                syn2cur.setdefault(s, sym)
    return lambda s: s if s in current else syn2cur.get(s, s)


def load_go():
    """symbol -> propagated GO ids (is_a / part_of closure)."""
    gaf = download(GOA, DATA_REF / "goa_human.gaf.gz")
    obo = download(OBO, DATA_REF / "go-basic.obo")

    parents, term, aspect = {}, None, {}
    for line in obo.read_text(encoding="utf-8", errors="replace").splitlines():
        if line == "[Term]":
            term = None
        elif line.startswith("id: GO:"):
            term = line[4:].strip()
            parents.setdefault(term, set())
        elif term and line.startswith("is_a: GO:"):
            parents[term].add(line[6:16])
        elif term and line.startswith("relationship: part_of GO:"):
            parents[term].add(line.split()[2])
        elif term and line.startswith("namespace:"):
            aspect[term] = line.split(": ", 1)[1]

    closure_cache = {}

    def closure(go_id):
        if go_id in closure_cache:
            return closure_cache[go_id]
        seen, stack = set(), [go_id]
        while stack:
            g = stack.pop()
            if g in seen:
                continue
            seen.add(g)
            stack.extend(parents.get(g, ()))
        closure_cache[go_id] = seen
        return seen

    sym2go = {}
    with gzip.open(gaf, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("!"):
                continue
            f = line.split("\t")
            if len(f) < 10 or f[3].startswith("NOT"):
                continue
            sym2go.setdefault(f[2], set()).update(closure(f[4]))
    print(f"  GO: {len(sym2go)} human symbols annotated, {len(parents)} terms")
    return sym2go


def main():
    m = pd.read_csv(RESULTS / "ortholog_map.tsv", sep="\t", index_col=0)
    contrast = pd.read_csv(RESULTS / "ev_sol_contrast.tsv", sep="\t", index_col=0)
    matrisome = load_matrisome()
    sym2go = load_go()

    df = contrast.join(m[[c for c in m.columns if c not in contrast.columns]])
    sym = df.human_symbol

    normalize = symbol_normalizer()
    matrisome["symbol"] = matrisome.EntrezGeneSymbol.map(normalize)
    renamed = (matrisome.symbol != matrisome.EntrezGeneSymbol).sum()
    print(f"  matrisome symbols updated to current nomenclature: {renamed}")
    mat = matrisome.drop_duplicates("symbol").set_index("symbol")
    df["matrisome_division"] = sym.map(mat.Division)
    df["matrisome_category"] = sym.map(mat.Category)

    go_sets = sym.map(lambda s: sym2go.get(s, set()) if isinstance(s, str) else set())
    has = lambda terms: go_sets.map(lambda g: bool(g & set(terms)))  # noqa: E731

    text = (df.subcellular_location.fillna("") + " " + df.keywords.fillna(""))
    df["has_signal_peptide"] = df.signal_peptide.notna()
    df["has_transmembrane"] = df.transmembrane.notna()
    df["is_gpi_anchored"] = df.lipidation.fillna("").str.contains("GPI", case=False)
    df["is_secreted"] = text.str.contains("Secreted", case=False)

    # "extracellular exosome" is annotated to thousands of human proteins, largely
    # from EV proteomics itself. In an EV prep it is near-universal and circular, so
    # it is kept as a weak flag while the category proper uses curated MISEV markers.
    df["in_exosome_go"] = has(GO_EV)
    df["is_ev_associated"] = sym.isin(EV_MARKERS)
    df["is_ecm_proteoglycan"] = df.matrisome_division.eq("Core matrisome") | has(["GO:0031012"])
    df["is_proteoglycan"] = df.matrisome_category.eq("Proteoglycans")
    df["is_matricellular"] = sym.isin(MATRICELLULAR)
    df["is_growth_factor_regulator"] = (
        df.matrisome_category.isin(["Secreted Factors", "ECM Regulators"]) | has(GO_GROWTH)
    )
    df["is_surface_corona"] = (
        has(GO_SURFACE) | df.has_signal_peptide | df.is_gpi_anchored | df.is_secreted
    )
    df["is_intracellular_soluble"] = (
        (has(GO_INTRACELLULAR) & ~df.has_signal_peptide & ~df.has_transmembrane
         & ~df.is_ecm_proteoglycan)
        | sym.isin(BLOOD)
    )

    # corona proper: secreted/ECM riding on the vesicle, not integral to its membrane
    df["is_corona_candidate"] = (
        (df.is_ecm_proteoglycan | df.is_matricellular | df.is_growth_factor_regulator
         | df.is_secreted)
        & ~df.has_transmembrane
    )

    precedence = [
        ("EV-associated", "is_ev_associated"),
        ("Matricellular", "is_matricellular"),
        ("ECM & proteoglycan", "is_ecm_proteoglycan"),
        ("Growth factor / cytokine regulator", "is_growth_factor_regulator"),
        ("Surface / corona-associated", "is_surface_corona"),
        ("Intracellular / soluble co-isolate", "is_intracellular_soluble"),
    ]
    df["class_knowledge"] = "Unclassified"
    for label, col in reversed(precedence):
        df.loc[df[col], "class_knowledge"] = label

    # --- empirical partitioning from the paired EV/SOL contrast ---
    sig = (df.q_paired_ttest < 0.05) | (df.q_detection_fisher < 0.05)
    df["class_empirical"] = "Ambiguous"
    df.loc[sig & (df.log2FC_EV_SOL > 1), "class_empirical"] = "EV-enriched"
    df.loc[sig & (df.log2FC_EV_SOL < -1), "class_empirical"] = "SOL-enriched"
    df.loc[sig & df.log2FC_EV_SOL.between(-1, 1), "class_empirical"] = "Co-present"
    df.loc[df.ev_exclusive, "class_empirical"] = "EV-exclusive"

    ev_side = df.class_empirical.isin(["EV-enriched", "EV-exclusive"])
    sol_side = df.class_empirical.eq("SOL-enriched")
    knowledge_ev = df.class_knowledge.ne("Intracellular / soluble co-isolate") & df.class_knowledge.ne("Unclassified")
    knowledge_sol = df.class_knowledge.eq("Intracellular / soluble co-isolate")

    df["class_agreement"] = "Undetermined"
    df.loc[(knowledge_ev & ev_side) | (knowledge_sol & sol_side), "class_agreement"] = "Agree"
    df.loc[(knowledge_ev & sol_side) | (knowledge_sol & ev_side), "class_agreement"] = "Disagree"

    # How a corona/ECM protein actually partitions. A negative EV/SOL ratio is NOT
    # evidence of contamination here: the SOL fraction is the soluble secretome by
    # construction, so matrix proteins are expected to be more concentrated there.
    # What marks a corona protein is reproducible presence ON the EVs regardless.
    df["corona_interpretation"] = pd.NA
    corona = df.is_corona_candidate
    df.loc[corona & df.retained & df.class_empirical.isin(["EV-enriched", "EV-exclusive"]),
           "corona_interpretation"] = "EV-enriched corona"
    df.loc[corona & df.retained & df.class_empirical.eq("SOL-enriched"),
           "corona_interpretation"] = "EV-bound, soluble-dominant"
    df.loc[corona & df.retained & df.class_empirical.isin(["Co-present", "Ambiguous"]),
           "corona_interpretation"] = "EV-associated, undifferentiated"
    df.loc[corona & ~df.retained, "corona_interpretation"] = "predominantly soluble"

    df.to_csv(RESULTS / "ev_atlas_master.tsv", sep="\t")

    retained = df[df.retained]
    print(f"\nretained EV proteome: {len(retained)} protein groups, "
          f"{retained.human_symbol.nunique()} unique human symbols")
    print("\nannotation layers (retained set):")
    for label, col in precedence + [("Corona candidate", "is_corona_candidate"),
                                    ("Proteoglycan", "is_proteoglycan"),
                                    ("in exosome GO (weak flag)", "in_exosome_go")]:
        n = int(retained[col].sum())
        ev_n = int((retained[col] & ev_side).sum())
        print(f"  {label:38} {n:5}   EV-enriched/exclusive: {ev_n}")

    print("\nprimary class (retained):")
    print(retained.class_knowledge.value_counts().to_string())
    print("\nempirical class (retained):")
    print(retained.class_empirical.value_counts().to_string())
    print("\nagreement (retained):")
    print(retained.class_agreement.value_counts().to_string())

    print("\ncorona interpretation (all corona candidates):")
    print(df.corona_interpretation.value_counts().to_string())

    dis = retained[retained.class_agreement == "Disagree"].sort_values("log2FC_EV_SOL", ascending=False)
    cols = ["human_symbol", "protein_name", "class_knowledge", "class_empirical",
            "log2FC_EV_SOL", "q_paired_ttest", "n_ev_high", "n_sol_high", "mean_ev_log2"]
    dis[cols].to_csv(RESULTS / "classification_disagreements.tsv", sep="\t")
    print(f"\nwrote {len(dis)} disagreements for review")

    for label, col in precedence + [("corona_candidate", "is_corona_candidate")]:
        slug = re.sub(r"[^a-z0-9]+", "_", label.lower()).strip("_")
        retained[retained[col]].to_csv(RESULTS / f"category_{slug}.tsv", sep="\t")


if __name__ == "__main__":
    main()
