# Execution, coordinates, and output schema

## Requirements and inputs

Python 3.10+ standard library, NCBI BLAST+ (`blastn`), Infernal (`cmsearch`, `cmalign`), and SeqKit for independent extraction QA. The development checks use Python 3.13, BLAST 2.16.0+, Infernal 1.1.5, and SeqKit 2.13.0; another version requires recording its version and reviewing changed results. Scripts do not install programs, download references, or modify source FASTA/GFF/index files.

First select the assembly and decide its nuclear scope. Supply uncompressed FASTA; unpack compressed downloads into a new work directory and preserve their original hashes. GFF is optional: missing GFF means context unassessed. If a full-assembly GFF accompanies a nuclear FASTA, provide only explicitly justified non-nuclear excluded GFF seqids. An unexplained GFF/FASTA ID mismatch is an error, not a reason to silently drop annotation.

Minimal JSON config (replace the input paths and reference identity; these are examples, not runnable placeholder data):

```json
{
  "key": "species",
  "prefix": "Sp",
  "species": "Genus species",
  "assembly": "assembly accession.version and cultivar",
  "genome": "inputs/selected_nuclear.fa",
  "gff": "inputs/matched.gff3",
  "excluded_gff_seqids": [],
  "query_panel": "plant-u6-identification/resources/queries/U6_all_queries.fasta",
  "historical_query": "plant-u6-identification/resources/queries/AtU6_legacy.fasta",
  "u6_model": "plant-u6-identification/resources/models/RF00026_U6.cm",
  "atac_model": "plant-u6-identification/resources/models/RF00619_U6atac_exclusion.cm",
  "threads": 2,
  "blastn": "blastn",
  "cmsearch": "cmsearch",
  "cmalign": "cmalign"
}
```

Relative input paths resolve against the config file's directory, not the current shell directory. Executables are PATH names or explicit executable paths. Use `"gff": null` when unavailable. The example assumes the skill folder is beside the config; use correct paths if it is installed elsewhere. For historical-At-only BLAST, set `query_panel` to the same `AtU6_legacy.fasta`; this does not disable RF00026 or RF00619.

`key` is a filename-safe identifier beginning with a letter; use only letters, digits, underscores and hyphens. `prefix` is the species abbreviation (e.g. `Sp`); locus IDs append `U6L0001` etc. Names are run-local, not stable cross-assembly identifiers. `reference_scope` can explicitly describe nuclear selection or a diagnostic subset. A GFF `##sequence-region` can legitimately describe only part of an included contig; its interval must be in bounds, and partial coverage is recorded rather than treated as complete annotation.

Optional `activity_evidence` is a list of externally reviewed records with `seqid`, `strand`, integer `body_start1`/`body_end1`, `source` (citation), and `evidence_type`. Allowed types are `engineered_promoter_activity`, `native_transcription`, or `genome_editing_with_promoter`. Additional notes may describe genotype, construct and assay. Only exact body-coordinate matches are linked; all supplied and matched records remain in `qa/activity_evidence_mapping.json`. The script does not verify papers or prove those claims. Review the source and mapping before supplying them; a predicted motif is not an activity record. Geometry stays unchanged even when experimental-source metadata is linked.

Run with real paths substituted:

```bash
python3 /path/to/plant-u6-identification/scripts/run_u6.py --config /path/to/config.json --out /path/to/new_run
python3 /path/to/plant-u6-identification/scripts/validate_outputs.py --run /path/to/new_run --seqkit /path/to/seqkit --out /path/to/new_independent_qa
```

Do not launch a full run for a reference-selection-only request. Preview `--help` when adapting command usage. The normal pipeline intentionally refuses an existing output directory. To reuse a completed expensive search in a fresh run, add `--reuse-search /path/to/old_run/01_search`. The runner validates completion, raw-output hashes, input hashes, and the frozen scoring contract. If any gate fails, retain the failed records, diagnose the mismatch, and run the affected work into a fresh version; do not edit a ledger to force acceptance.

## Output layers

`README_where_to_find.md` is the single entrypoint. The numbered folders keep reference audit, raw searches, curation, promoter evidence and final delivery separate. `04_delivery/` contains the user-facing catalogue:

- `<key>_U6_candidates.tsv`: every retained canonical-family candidate, not merely a shortlist.
- `<key>_U6_loci.bed` and `<key>_U6_body.fa`: corresponding body intervals and transcript-oriented sequences.
- `<key>_U6_promoter_300.fa`, `_500.fa`, `_1000.fa`: primary-boundary upstream sequences.
- Separate strict historical-At subset, all search loci, U6atac/review evidence, and motif/context/boundary evidence tables.
- `U6_candidate_report.md`, `Methods_and_Results_current_stage.md`, command/version records and input manifests.

The generated entrypoint gives the exact paths for the current version. Empty candidate sets remain valid searched outcomes, not a reason to invent candidates; any empty/truncated upstream has an explicit status. FASTA with an empty record alone is not adequate documentation: use the corresponding promoter table to distinguish zero upstream sequence from failed extraction.

## Coordinate and missingness contract

| Item | Convention |
|---|---|
| GFF, locus and promoter table intervals | 1-based closed, start≤end on either strand |
| BED | 0-based start, end-exclusive; strand in column 6 |
| `predicted_tss1` | genomic coordinate of hypothetical +1, not experimentally measured TSS |
| Plus-strand upstream length L | `[t−L, t−1]`, clipped only at reference edge |
| Minus-strand upstream length L | reverse-complement `[t+1, t+L]`, clipped at reference edge |
| Transcript-relative positions | upstream −L..−1, first RNA base +1; no zero |
| Source promoter/vector block | preserve original length and mapping; may contain +1, never assume |
| Missing scalar | `NA` in TSV, `null` in JSON |
| Structured TSV field | JSON object/list inside the cell; distinguish nested `null` from 0 |
| Negative finding | zero only within the actually assessed annotation/search scope |

Minimum evidence carried per candidate: assembly and locus identity; seqid/strand/body interval; RF00026/RF00619 evidence and classification; strict historical-query status; alignment-derived body integrity; hypothetical +1 and boundary origin; promoter coordinates and extraction status; detected motif sequences/positions/mismatches and profile-specific geometry; termination context; source-specific annotation overlaps; activity evidence, interpretation, and review reasons. PCR uniqueness and homeolog grouping remain unassessed unless an additional authorized analysis supplies them.

Keep local-body BLAST used for boundary mapping separate from the original genome-wide search. Its E-values are not genome-wide discovery statistics. Alternative boundary exports are not additional genomic loci.

## QA and reporting

Primary QA checks pipeline bookkeeping and input preservation. Independent QA reconstructs genomic extraction requests from result tables, checks the body/BED/promoter relationships, and compares SeqKit-extracted sequences without using the primary extraction implementation. Corrupted-sequence and coordinate controls test that the validator fails when it should.

Neither QA layer certifies promoter activity, native TSS, expression rank, exhaustive functional sensitivity, assembly correctness, or PCR success. Missing SeqKit is an unrun independent check, not PASS. The final narrative must state:

1. Which reference and exact search/profile versions were used.
2. Full family count, strict-At subset count, separate U6atac/ambiguous/weak sets.
3. Integrity, boundary, motif and annotation limitations relevant to experimental choices.
4. Actual execution and QA status, with any stage not run.

Use deterministic cohort counts and factual methods text; do not describe all computational family candidates as expressed genes or all motif-positive windows as active promoters.
