---
name: plant-u6-identification
description: Identify plant nuclear U6 snRNA loci with nucleotide BLAST and Rfam covariance models, distinguish U6atac, and extract and assess strand-correct upstream promoters. Use for new-species U6 discovery, candidate review, or reproducible promoter identification; not protein-coding families or RNA-expression quantification.
metadata:
  version: "1.0.1"
  repository: "https://github.com/yz3394/Genome-wide-plant-u6-identification"
---

# Plant U6 identification

Deliver a reference-bounded U6 candidate catalogue and auditable promoter evidence for plant genome-editing research. Separate RNA-family support, body integrity, boundary confidence, promoter geometry, annotation context, and experimental activity. Computational support does not establish transcription or promoter strength.

## Select the requested stage

- **Reference selection or planning only:** audit species/ploidy/cultivar, assembly scope and matched annotation; do not begin searches before selection is accepted. Nuclear unplaced scaffolds matter. Freeze source accessions, versions, checksums, and deliberate exclusions.
- **Discovery and promoter extraction:** read [the evidence rules](references/evidence_and_profiles.md) and [execution/schema guide](references/schema_and_coordinates.md), then use the packaged scripts. Do not rewrite coordinate or clustering logic ad hoc.
- **Review existing results:** inspect raw evidence and provenance first; read [validation cases](references/validation_cases.md) when evaluating method limitations or rerunning a regression. A review request alone does not authorize overwriting or rerunning canonical results.

## Execute

1. Confirm the user-selected reference and search scope. Audit an uncompressed nuclear FASTA and optional matching GFF3. Record missing annotation categories as unassessed, not clean. Do not silently rename contigs or exclude unplaced records.
2. Configure `scripts/run_u6.py` with portable paths and a fresh output directory, following the schema guide. Bundled resources contain the frozen four-query panel, historical AtU6 query, and Rfam 15.1 RF00026/RF00619 models; their provenance and hashes are in `resources/manifest.json`. Respect requests for an At-only BLAST panel. RF00026 and RF00619 still scan the entire selected reference independently.
3. Run the pipeline. Require actual successful search logs, matching hashes, and complete expected outputs. For an interrupted project, reuse only a validated completed search with `--reuse-search`; retain failed/incomplete directories and restart into a new version. Never treat a directory or a lone completed flag as proof of completion.
4. Run `scripts/validate_outputs.py` in a separate fresh QA directory using SeqKit. Report execution success, independent extraction QA, and biological uncertainty separately. Missing dependencies or unrun checks must remain explicit; do not install tools without the necessary authorization.
5. Deliver `README_where_to_find.md`, `04_delivery/`, and QA links. The candidates table is the **entire RF00026-supported canonical catalogue**, not just strict-At hits or a motif shortlist. Preserve all search loci, U6atac, ambiguous and weak evidence in separate tables. Explain proposed experimental priorities with reasons; do not imply measured activity.

## Scientific invariants

- U6 is noncoding: BLASTN on the transcribed body, not BLASTP or promoter-sequence queries. Multiple query agreement is optional corroboration, not a two-query gate.
- Classify RF00026 versus RF00619 at the same strand/locus with real overlap and each frozen model's GA. Weak RF00619 or a nearby distinct U6atac cannot exclude a U6 candidate. Keep dual-GA ambiguity for review.
- Single-HSP identity, exact integer-derived query coverage, and score must support any strict BLAST label together. Do not combine maxima or transfer bit-score thresholds across scoring profiles.
- Use actual aligned query position 1 if supported; otherwise name the body-edge hypothesis. Neither is a measured target-species TSS. Keep alternative boundaries and incomplete submitted bodies explicit.
- Upstream FASTA excludes +1 and is transcript-oriented. Source/vector blocks may include +1; inspect their boundaries rather than assuming a name defines them.
- USE/TATA detection and geometry are separate from family and activity. Retain noncanonical candidates; do not widen thresholds simply to pass a known example. The current generic profile is provisional and Arabidopsis-derived.
- Self-U6 annotation and foreign exon/CDS/gene-span/repeat overlaps are different. A conflict is a review flag, not automatic disproof of family membership. Missing GFF/repeat/UTR evidence, PCR specificity and homeolog assignments remain NA unless independently assessed.

## Boundaries

Do not automatically launch RNA-seq analysis, PCR/construct design, experimental genome editing, or infer high-expression promoters. Respect stage-limited requests. Preserve originals and version outputs. This skill is a research-grade candidate/evidence workflow, not a guarantee of exhaustive functional U6 discovery in every assembly.

For a requested skill fix or upgrade, read [release maintenance](references/release_maintenance.md). Keep source publication, installed copies and research outputs separate; analysis itself never uploads data.
