"""Map porcine accessions to human gene symbols.

Two problems make this the load-bearing step:

  * 24% of deposited rows carry no GN= symbol, and those hide core EV/ECM markers.
  * 57% of the accessions have since been DELETED from UniProt - the Sus scrofa
    TrEMBL proteome was pruned after the study's Dec-2023 search database. A
    current-UniProt-only lookup silently loses more than half the atlas.

So pig identity is resolved from four sources (current UniProt, UniParc archive,
the deposited description, UniSave history), then carried to human by NCBI
orthology, symbol passthrough, or Ensembl Compara. Every protein keeps an audit
trail of which route produced its symbol.
"""

import io
import re
import sys
import time
from pathlib import Path

import pandas as pd
import requests

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import DATA_REF, RESULTS  # noqa: E402

UNIPROT_FIELDS = (
    "accession,protein_name,gene_names,gene_primary,xref_geneid,xref_ensembl,"
    "cc_subcellular_location,ft_signal,ft_transmem,ft_lipid,keyword,go_id,length"
)
NCBI = "https://ftp.ncbi.nlm.nih.gov/gene/DATA/"
# Markers that must survive mapping. FN1 is deliberately absent: no fibronectin
# entry exists anywhere in this deposit (only FN3-domain proteins), so expecting
# it would mask a real finding as a mapping bug.
CHECKPOINT = ["CD9", "CD81", "TSG101", "DCN", "BGN", "THBS1", "AHNAK", "HSPG2"]
SESSION = requests.Session()


def cached(name, fn):
    path = DATA_REF / name
    if path.exists():
        return pd.read_csv(path, sep="\t", dtype=str)
    df = fn()
    df.to_csv(path, sep="\t", index=False)
    return df


def fetch_uniprot(accessions):
    def go():
        frames = []
        for i in range(0, len(accessions), 100):
            batch = accessions[i : i + 100]
            r = SESSION.get(
                "https://rest.uniprot.org/uniprotkb/accessions",
                params={"accessions": ",".join(batch), "fields": UNIPROT_FIELDS, "format": "tsv"},
                timeout=180,
            )
            r.raise_for_status()
            if r.text.strip():
                frames.append(pd.read_csv(io.StringIO(r.text), sep="\t", dtype=str))
            print(f"  UniProt {i + len(batch)}/{len(accessions)}", end="\r")
            time.sleep(0.15)
        return pd.concat(frames, ignore_index=True)

    df = cached("uniprot_pig.tsv", go)
    print(f"  UniProt active entries: {len(df)}")
    return df.rename(columns={"Entry": "accession"}).set_index("accession")


def fetch_uniparc(accessions):
    """UniParc keeps deleted UniProtKB entries, with gene and protein names."""

    def go():
        rows = []
        for i in range(0, len(accessions), 50):
            batch = accessions[i : i + 50]
            r = SESSION.get(
                "https://rest.uniprot.org/uniparc/search",
                params={
                    "query": " OR ".join(batch),
                    "fields": "upi,gene,protein,accession",
                    "format": "tsv",
                    "size": "500",
                },
                timeout=180,
            )
            r.raise_for_status()
            df = pd.read_csv(io.StringIO(r.text), sep="\t", dtype=str)
            wanted = set(batch)
            for _, row in df.iterrows():
                members = {
                    t.split(".")[0]
                    for t in re.split(r"[;,\s]+", str(row.get("UniProtKB") or ""))
                    if t
                }
                for acc in members & wanted:
                    rows.append(
                        {
                            "accession": acc,
                            "uniparc_gene": (str(row["Gene names"]).split(";")[0].strip()
                                             if pd.notna(row.get("Gene names")) else None),
                            "uniparc_protein": (str(row["Protein names"]).split(";")[0].strip()
                                                if pd.notna(row.get("Protein names")) else None),
                        }
                    )
            print(f"  UniParc {i + len(batch)}/{len(accessions)}", end="\r")
            time.sleep(0.15)
        return pd.DataFrame(rows).drop_duplicates("accession")

    df = cached("uniparc_pig.tsv", go)
    print(f"  UniParc recovered: {len(df)} of {len(accessions)} deleted accessions")
    return df.set_index("accession")


def fetch_unisave(accessions):
    """Last archived UniProt version: GN name and Ensembl gene for stubborn cases."""

    def go():
        rows = []
        for n, acc in enumerate(accessions, 1):
            try:
                r = SESSION.get(
                    f"https://rest.uniprot.org/unisave/{acc}",
                    params={"format": "txt"}, timeout=120,
                )
                text = r.text if r.status_code == 200 else ""
            except Exception:
                text = ""
            gn = re.findall(r"^GN   Name=([^;{\s]+)", text, re.M)
            eg = re.findall(r"ENSSSCG\d+", text)
            rows.append(
                {
                    "accession": acc,
                    "unisave_gene": gn[-1] if gn else None,
                    "unisave_ensembl": eg[-1] if eg else None,
                }
            )
            if n % 25 == 0:
                print(f"  UniSave {n}/{len(accessions)}", end="\r")
            time.sleep(0.05)
        return pd.DataFrame(rows)

    df = cached("unisave_pig.tsv", go)
    print(f"  UniSave: {df.unisave_gene.notna().sum()} genes, "
          f"{df.unisave_ensembl.notna().sum()} Ensembl ids from {len(df)} queried")
    return df.set_index("accession")


def human_swissprot():
    def go():
        r = SESSION.get(
            "https://rest.uniprot.org/uniprotkb/stream",
            params={
                "query": "reviewed:true AND organism_id:9606",
                "fields": "accession,gene_primary,protein_name",
                "format": "tsv",
            },
            timeout=900,
        )
        r.raise_for_status()
        return pd.read_csv(io.StringIO(r.text), sep="\t", dtype=str)

    return cached("human_swissprot.tsv", go)


def _split_names(s):
    """UniProt packs alternative names in parentheses: 'Foo (Bar) (EC 1.2.3.4)'."""
    out, depth, cur = [], 0, ""
    for ch in s:
        if ch == "(":
            depth += 1
            if depth == 1:
                if cur.strip():
                    out.append(cur.strip())
                cur = ""
                continue
        elif ch == ")":
            depth -= 1
            if depth == 0:
                if cur.strip():
                    out.append(cur.strip())
                cur = ""
                continue
        cur += ch
    if cur.strip():
        out.append(cur.strip())
    return out


def _norm_name(n):
    n = re.sub(r"\b(isoform|variant)\s*X?[0-9]+\b", "", n, flags=re.I)
    n = re.sub(r"-like\b", "", n, flags=re.I)
    n = re.sub(r"\bEC [0-9.\-]+", "", n)
    n = re.sub(r"\b(Fragment|precursor|preproprotein)\b", "", n, flags=re.I)
    return re.sub(r"[^a-z0-9]", "", n.lower())


def protein_name_resolver():
    """Map a descriptive protein name to a human symbol.

    Pig RefSeq-derived entries keep human-style names ('Biglycan', 'Myosin-9'),
    so this recovers identities that have no gene name or cross-reference at all.
    The full name is matched before any parenthetical fragment - splitting first
    would resolve 'Collagen alpha-1(VI) chain' to COL1A1 instead of COL6A1.
    """
    sp = human_swissprot()
    full, part = {}, {}
    for sym, pn in zip(sp["Gene Names (primary)"], sp["Protein names"]):
        if not isinstance(sym, str) or not isinstance(pn, str):
            continue
        full.setdefault(_norm_name(pn), sym)
        names = _split_names(pn)
        if names:
            full.setdefault(_norm_name(names[0]), sym)
        for p in names[1:]:
            k = _norm_name(p)
            if len(k) > 8:
                part.setdefault(k, sym)

    def resolve(pn):
        if not isinstance(pn, str):
            return None
        if _norm_name(pn) in full:
            return full[_norm_name(pn)]
        names = _split_names(pn)
        if names and _norm_name(names[0]) in full:
            return full[_norm_name(names[0])]
        for p in sorted(names, key=len, reverse=True):
            k = _norm_name(p)
            if len(k) > 8 and k in part:
                return part[k]
        return None

    print(f"  human SwissProt name index: {len(full)} names")
    return resolve


def download(url, dest):
    if not dest.exists():
        print(f"  downloading {dest.name} ...")
        with SESSION.get(url, stream=True, timeout=1800) as r:
            r.raise_for_status()
            with open(dest, "wb") as fh:
                for chunk in r.iter_content(1 << 20):
                    fh.write(chunk)
    return dest


def ncbi_tables():
    orth = download(NCBI + "gene_orthologs.gz", DATA_REF / "gene_orthologs.gz")
    info = download(
        NCBI + "GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz",
        DATA_REF / "Homo_sapiens.gene_info.gz",
    )
    o = pd.read_csv(orth, sep="\t", dtype=str, compression="gzip")
    o.columns = [c.lstrip("#") for c in o.columns]
    fwd = o[(o.tax_id == "9606") & (o.Other_tax_id == "9823")]
    rev = o[(o.tax_id == "9823") & (o.Other_tax_id == "9606")]
    pig2human = dict(zip(fwd.Other_GeneID, fwd.GeneID))
    pig2human.update(dict(zip(rev.GeneID, rev.Other_GeneID)))

    gi = pd.read_csv(
        info, sep="\t", dtype=str, compression="gzip", header=0,
        usecols=[1, 2, 4, 5], names=["GeneID", "Symbol", "Synonyms", "dbXrefs"],
    )
    entrez2symbol = dict(zip(gi.GeneID, gi.Symbol))
    approved = set(gi.Symbol)
    synonym2symbol = {}
    for sym, syns in zip(gi.Symbol, gi.Synonyms.fillna("-")):
        for s in syns.split("|"):
            if s != "-" and s not in approved:
                synonym2symbol.setdefault(s, sym)
    ensg2symbol = {}
    for sym, xr in zip(gi.Symbol, gi.dbXrefs.fillna("")):
        for m in re.findall(r"Ensembl:(ENSG\d+)", xr):
            ensg2symbol[m] = sym
    print(f"  NCBI: {len(pig2human)} pig-human ortholog pairs, {len(approved)} human symbols, "
          f"{len(ensg2symbol)} Ensembl-mapped")
    return pig2human, entrez2symbol, approved, synonym2symbol, ensg2symbol


def ensembl_gene_of(transcript_ids):
    """Batch transcript/protein -> parent gene via Ensembl POST lookup."""
    out = {}
    ids = [i for i in transcript_ids if i]
    for i in range(0, len(ids), 500):
        batch = ids[i : i + 500]
        try:
            r = SESSION.post(
                "https://rest.ensembl.org/lookup/id",
                json={"ids": batch},
                headers={"Content-Type": "application/json", "Accept": "application/json"},
                timeout=180,
            )
            if r.status_code != 200:
                continue
            for k, v in r.json().items():
                if not v:
                    continue
                out[k] = v.get("Parent") if v.get("object_type") == "Transcript" else v.get("id")
        except Exception:
            pass
        print(f"  Ensembl lookup {i + len(batch)}/{len(ids)}", end="\r")
    return out


def compara(ensembl_genes, ensg2symbol):
    def go():
        rows = []
        for n, g in enumerate(ensembl_genes, 1):
            sym = htype = None
            try:
                r = SESSION.get(
                    f"https://rest.ensembl.org/homology/id/sus_scrofa/{g}",
                    params={"target_species": "homo_sapiens", "type": "orthologues",
                            "format": "condensed"},
                    headers={"Content-Type": "application/json"},
                    timeout=90,
                )
                if r.status_code == 200:
                    homs = r.json()["data"][0]["homologies"]
                    if homs:
                        best = min(homs, key=lambda h: 0 if h.get("type") == "ortholog_one2one" else 1)
                        sym = ensg2symbol.get(best.get("id"))
                        htype = best.get("type")
            except Exception:
                pass
            rows.append({"pig_ensembl": g, "compara_symbol": sym, "compara_type": htype})
            if n % 25 == 0:
                print(f"  Compara {n}/{len(ensembl_genes)}", end="\r")
        return pd.DataFrame(rows)

    df = cached("compara.tsv", go)
    print(f"  Compara: {df.compara_symbol.notna().sum()}/{len(df)} resolved")
    return df.set_index("pig_ensembl")


def main():
    ev = pd.read_csv(RESULTS / "protein_evidence.tsv", sep="\t", dtype=str).set_index("accession")
    accessions = ev.index.tolist()
    print(f"mapping {len(accessions)} porcine accessions\n")

    up = fetch_uniprot(accessions)
    deleted = [a for a in accessions if a not in up.index]
    print(f"  deleted from current UniProt: {len(deleted)} ({len(deleted) / len(accessions):.1%})")
    uparc = fetch_uniparc(deleted)

    m = pd.DataFrame(index=pd.Index(accessions, name="accession"))
    m["symbol_deposited"] = ev.pig_symbol_from_export
    m["symbol_uniprot"] = up.get("Gene Names (primary)")
    m["symbol_uniparc"] = uparc.get("uniparc_gene")
    m["protein_name"] = up.get("Protein names").combine_first(
        uparc.get("uniparc_protein")).combine_first(ev.description.str.replace(r" OS=.*", "", regex=True))
    for col, src in [("subcellular_location", "Subcellular location [CC]"), ("signal_peptide", "Signal peptide"),
                     ("transmembrane", "Transmembrane"), ("lipidation", "Lipidation"),
                     ("keywords", "Keywords"), ("go_ids", "Gene Ontology IDs")]:
        m[col] = up.get(src)
    m["pig_entrez"] = up.get("GeneID").str.split(";").str[0] if "GeneID" in up else None
    m["pig_ensembl_xref"] = up.get("Ensembl")

    m["pig_symbol"] = m.symbol_uniprot.combine_first(m.symbol_uniparc).combine_first(m.symbol_deposited)
    print(f"  pig symbol resolved for {m.pig_symbol.notna().sum()}/{len(m)}")

    pig2human, entrez2symbol, approved, synonym2symbol, ensg2symbol = ncbi_tables()

    m["human_symbol"] = None
    m["mapping_method"] = None
    m["ortholog_type"] = None

    # layer A: pig Entrez -> human Entrez (authoritative orthology)
    he = m.pig_entrez.map(lambda g: pig2human.get(g) if isinstance(g, str) else None)
    hs = he.map(lambda g: entrez2symbol.get(g) if isinstance(g, str) else None)
    m["human_entrez"] = he
    take = hs.notna()
    m.loc[take, ["human_symbol", "mapping_method", "ortholog_type"]] = pd.DataFrame(
        {"human_symbol": hs[take], "mapping_method": "ncbi_ortholog", "ortholog_type": "1:1"}
    )
    print(f"  layer A NCBI orthology: {take.sum()}")

    # layer B: pig symbol already an approved human symbol, or a known synonym
    need = m.human_symbol.isna() & m.pig_symbol.notna()
    direct = need & m.pig_symbol.isin(approved)
    m.loc[direct, "human_symbol"] = m.loc[direct, "pig_symbol"]
    m.loc[direct, ["mapping_method", "ortholog_type"]] = ["symbol_exact", "1:1"]
    need = m.human_symbol.isna() & m.pig_symbol.notna()
    syn = need & m.pig_symbol.isin(synonym2symbol)
    m.loc[syn, "human_symbol"] = m.loc[syn, "pig_symbol"].map(synonym2symbol)
    m.loc[syn, ["mapping_method", "ortholog_type"]] = ["symbol_synonym", "1:1"]
    print(f"  layer B symbol passthrough: {direct.sum()} exact + {syn.sum()} synonym")

    # layer C: Ensembl Compara for the remainder
    residual = m.index[m.human_symbol.isna()].tolist()
    print(f"  residual after A+B: {len(residual)}")
    if residual:
        us = fetch_unisave(residual)
        m.loc[residual, "unisave_gene"] = us.unisave_gene
        m.loc[residual, "unisave_ensembl"] = us.unisave_ensembl
        m["pig_symbol"] = m.pig_symbol.combine_first(m.get("unisave_gene"))

        # retry layer B with UniSave-recovered symbols
        need = m.human_symbol.isna() & m.pig_symbol.notna()
        direct = need & m.pig_symbol.isin(approved)
        m.loc[direct, "human_symbol"] = m.loc[direct, "pig_symbol"]
        m.loc[direct, ["mapping_method", "ortholog_type"]] = ["symbol_exact", "1:1"]
        syn = m.human_symbol.isna() & m.pig_symbol.isin(synonym2symbol)
        m.loc[syn, "human_symbol"] = m.loc[syn, "pig_symbol"].map(synonym2symbol)
        m.loc[syn, ["mapping_method", "ortholog_type"]] = ["symbol_synonym", "1:1"]
        print(f"  layer B retry after UniSave: +{direct.sum() + syn.sum()}")

        still = m.index[m.human_symbol.isna()].tolist()
        tx = {}
        for acc in still:
            xr = m.at[acc, "pig_ensembl_xref"]
            if isinstance(xr, str):
                ids = re.findall(r"ENSSSC[TG]\d+", xr)
                if ids:
                    tx[acc] = ids[0]
        parents = ensembl_gene_of(sorted(set(tx.values())))
        m["pig_ensembl"] = pd.Series(
            {acc: parents.get(t, t if t.startswith("ENSSSCG") else None) for acc, t in tx.items()}
        )
        if "unisave_ensembl" in m:
            m["pig_ensembl"] = m.pig_ensembl.combine_first(m.unisave_ensembl)

        genes = sorted(m.loc[m.human_symbol.isna(), "pig_ensembl"].dropna().unique())
        print(f"  Compara candidates: {len(genes)} pig genes")
        if genes:
            cp = compara(genes, ensg2symbol)
            hit = m.human_symbol.isna() & m.pig_ensembl.isin(cp.index[cp.compara_symbol.notna()])
            m.loc[hit, "human_symbol"] = m.loc[hit, "pig_ensembl"].map(cp.compara_symbol)
            m.loc[hit, "mapping_method"] = "ensembl_compara"
            m.loc[hit, "ortholog_type"] = m.loc[hit, "pig_ensembl"].map(
                cp.compara_type.map(lambda t: "1:1" if t == "ortholog_one2one" else "1:many")
            )
            print(f"  layer C Compara: {hit.sum()}")

    # layer D: descriptive protein name -> human symbol, for entries with no
    # gene name and no usable cross-reference left
    resolve = protein_name_resolver()
    still = m.human_symbol.isna()
    names = m.loc[still, "protein_name"].map(resolve)
    hit = names.notna()
    m.loc[names[hit].index, "human_symbol"] = names[hit]
    m.loc[names[hit].index, "mapping_method"] = "protein_name"
    m.loc[names[hit].index, "ortholog_type"] = "inferred"
    print(f"  layer D protein-name match: {hit.sum()}")

    m["mapping_method"] = m.mapping_method.fillna("unmapped")
    m["mapping_confidence"] = m.mapping_method.map(
        {"ncbi_ortholog": "high", "ensembl_compara": "high", "symbol_exact": "high",
         "symbol_synonym": "medium", "protein_name": "medium", "unmapped": "none"}
    )
    m.reset_index().to_csv(RESULTS / "ortholog_map.tsv", sep="\t", index=False)

    mapped = m.human_symbol.notna().sum()
    print(f"\nmapped {mapped}/{len(m)} ({mapped / len(m):.1%}); "
          f"{m.human_symbol.nunique()} unique human symbols")
    print(m.mapping_method.value_counts().to_string())
    found = [c for c in CHECKPOINT if c in set(m.human_symbol.dropna())]
    print(f"checkpoint recovered {len(found)}/{len(CHECKPOINT)}: {found}")
    if len(found) < len(CHECKPOINT):
        print(f"WARNING: missing {[c for c in CHECKPOINT if c not in found]}")
    if mapped / len(m) < 0.90:
        print("WARNING: mapping rate below 90%")


if __name__ == "__main__":
    main()
