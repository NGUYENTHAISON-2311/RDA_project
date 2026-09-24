# NC-EV Proteomic Atlas — PXD072052

A corona-aware reanalysis of [PXD072052](https://www.ebi.ac.uk/pride/archive/projects/PXD072052),
the proteomic atlas behind *"Proteomic atlas of notochordal cell-derived extracellular vesicles
uncovers functional protein corona affecting intracellular pathways in chondrocyte-like cells"*
(porcine NC-EVs, Orbitrap Eclipse, Proteome Discoverer 2.5).

The pipeline maps porcine protein groups to human orthologs, filters for donor reproducibility,
annotates six biological layers, and runs enrichment **against the detected EV proteome** rather than
the full human proteome — so "ECM is enriched" means enriched relative to what an EV preparation
actually contains, not relative to all of biology.

## Design

14 runs = 7 preparations × paired **EV** and **SOL** (soluble) fractions, across 2 batches
(TRY_0463: P13–P14; TRY_0500: P20–P24). The pairing is what makes an empirical
corona-versus-co-isolate call possible instead of an argument about annotation lists.

## Getting the data without downloading 223 GB

The processed results are deposited as a single ZIP64 archive of 223,172,625,425 bytes, mostly
Proteome Discoverer result databases of 56–77 GB each. The protein table inside it is 2.4 MB.

`src/01_fetch_pride.py` reads the ZIP64 central directory from the tail of the remote file, locates
the two members it needs by name, and inflates only those byte ranges — about 700 KB transferred.
No raw files are downloaded and no spectra are re-searched.

## Pipeline

| Step | Does |
|---|---|
| `01_fetch_pride.py` | Range-request the protein table and study metadata out of the 223 GB archive; build the sample sheet |
| `02_build_matrix.py` | Parse the PD export into abundance / detection matrices; **gate** the F-number→sample mapping against marker biology |
| `03_map_orthologs.py` | Porcine → human symbols through four identity sources and four mapping layers |
| `04_filter_reproducible.py` | Retain ≥4 of 7 EV preparations; paired EV/SOL contrast within preparation |
| `05_annotate_categories.py` | Six annotation layers, plus knowledge-based and measurement-based classifications side by side |
| `06_enrichment.py` | GO / Reactome / KEGG via g:Profiler with a custom background; matrisome by Fisher exact |
| `07_report.py` | Summary figures and the report's data payload |
| `08_gag_pathways.py` | KEGG glycosaminoglycan pathway membership, evidence and enrichment |

Run in order:

```bash
python src/01_fetch_pride.py
python src/02_build_matrix.py
# ... through 08
```

Everything downloaded is cached under `data/`; delete it to force a clean re-fetch. Both `data/` and
`results/` are generated and therefore untracked.

## Two problems that shaped the analysis

**57% of the accessions no longer exist.** The study searched a December 2023 *Sus scrofa* database;
UniProt has since pruned the porcine TrEMBL proteome, so 2,364 of 4,116 accessions are now *deleted*
entries. A current-UniProt lookup alone silently loses more than half the atlas, and a quarter of the
deposited rows carry no gene symbol at all. Identity is resolved from current UniProt → the UniParc
archive → the deposited description → the UniSave version history, then carried to human by NCBI
orthology, symbol passthrough, Ensembl Compara, or matching descriptive protein names against human
SwissProt. Final mapping rate: **97.2%**. The most abundant protein in the dataset is deposited as
"Uncharacterized protein" with a deleted accession — it is AHNAK.

**The soluble fraction carries 10–40× more protein than the vesicle fraction.** Ratios are therefore
computed on composition-normalised values; normalisation shifts the median log2 ratio by ~2.4 units
against the raw numbers, so both are kept in the output. A negative EV/SOL ratio is *not* a verdict
of contamination — the SOL fraction is the soluble secretome by construction, so matrix proteins
belong there. What marks a corona protein is reproducible presence on the vesicles regardless.

## Validation built into the pipeline

- `02` refuses to continue unless EV markers favour the EV runs and plasma proteins favour the SOL
  runs — the F-number→sample assignment is tested, not trusted.
- `03` fails loudly if canonical markers (CD9, CD81, TSG101, DCN, BGN, THBS1, AHNAK, HSPG2) do not
  survive mapping, or if the mapping rate drops below 90%.
- `06` checks that generic "extracellular exosome" terms collapse in effect size (max 2.9× fold);
  anything higher means the custom background silently failed to apply.

## Results

4,116 protein groups → **3,048 retained**, 2,931 unique human symbols.

Every canonical EV marker is EV-enriched (CD63 +4.1, CD9 +3.4, CD81 +2.3; ALIX and RAB27A
EV-exclusive). The EV-enriched corona is an annexin and S100 layer. Structural matrix (thrombospondins,
tenascin-X, periostin, perlecan, versican) is soluble-dominant, while a defined subset — aggrecan,
COMP, fibromodulin, CILP, collagen II, decorin, biglycan, SPARC — is detected on vesicles in all
seven preparations.

GAG-synthesis machinery is present and confidently identified (XYLT1 at 78 unique peptides), but the
Golgi-lumenal chain builders and sulfotransferases partition to the soluble fraction, while the
cytosolic nucleotide-sugar and PAPS-donor enzymes sit near zero and are retained almost completely.
No chain-building or sulfating enzyme is EV-enriched.

## Requirements

Python 3.12+ with `pandas`, `numpy`, `scipy`, `requests`, `matplotlib`. No API keys — PRIDE, UniProt,
NCBI, Ensembl, KEGG, g:Profiler and the Gene Ontology are all queried anonymously.

## Data sources

PRIDE · UniProt / UniParc / UniSave · NCBI Gene orthologs · Ensembl Compara · Gene Ontology ·
Reactome and KEGG via g:Profiler · the human matrisome masterlist (Naba lab).
