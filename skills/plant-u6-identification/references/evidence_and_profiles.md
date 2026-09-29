# Evidence rules and frozen profiles

## Discovery

Use an explicitly selected nuclear assembly, including unplaced nuclear records. Do not infer nuclear status from record-name patterns alone; consult assembly reports or source metadata. This package audits file/coordinate consistency, not taxonomic identity, ploidy, assembly completeness, or haplotype collapse. A successful file audit does not replace reference selection.

The default nucleotide queries are current AtU6-26 `NR_141590.1` (101 nt), historical `X52528.1:449–550` (102 nt, identical to the selected X52527/X52529 bodies), tomato `X51447.1:262–364` (103 nt), and broad bean `X04788.1:1–98` (98 nt). Current and historical At boundaries differ; do not silently substitute one. The historical query's source study supported the 5′ start, but target-species homology is not a measured target TSS. Query/model provenance is bundled in `resources/`.

BLASTN-short baseline: both strands, reward 1, penalty −3, gap open 5/extend 2, word size 7, dust off, soft masking false; broad reporting E≤1000, identity≥60%, HSP query coverage≥50%. The historical strict label requires **one HSP** with identity≥80%, exact `(qend−qstart+1)/qlen≥0.8`, and bitscore≥80. Rounded qcov 80 is not sufficient for an actual 79.59%. Strong/multiple BLAST support is not mandatory for CM-supported membership. If a user requests only the historical At query, use that BLAST panel and disclose the reduced query coverage strategy; retain independent model searches.

RF00026 (U6) and RF00619 (U6atac) each scan the complete supplied FASTA. Report-score default is 20, below both frozen GA values. Read GA and consensus length from each CM, not from a copied constant; bundled Rfam 15.1 has GA 44/33 bits and CLEN 105/126 respectively. Do not compare E-values from a local rescue to genome-wide E-values as equivalent evidence. Infernal heuristic filters remain enabled; a negative search is not a mathematical absence proof.

Fixed-anchor curation preserves neighboring copies and raw multi-anchor evidence. Normal association uses ≥15 nt overlap, ≥50% of the shorter hit, and an envelope ≤140 nt; same-query duplicate suppression uses ≥90% reciprocal overlap. Competing families require same strand/sequence and ≥50% reciprocal overlap. Inspect ambiguous associations rather than chaining a wide transitive cluster across nearby genes.

Family layers:

- RF00026 GA support without competing RF00619 GA: retain as canonical-family candidate, regardless of motif class.
- Both family GAs at the same competing locus: preserve as family ambiguity; do not call proven canonical U6.
- RF00619 GA without canonical support: separate U6atac evidence.
- BLAST-only, below-GA, and fragmented evidence: retain in the review catalogue, not the confirmed-family export.

## Integrity and boundary

Align submitted bodies with `cmalign`, inspect actual occupied model columns, terminal missing positions, internal missing runs, paired-column compatibility, and ambiguous bases. The packaged integrity thresholds are pilot review flags, not universal functional criteria. A high model-span percentage or `trunc=no` does not establish a complete RNA. A short submitted body does not prove a genomic deletion or pseudogene.

For the primary boundary use a strong, actually aligned historical query base 1, mapped through alignment gaps. Otherwise use the strand-aware search-body 5′ edge as an explicit hypothesis. Do not linearly extrapolate missing q1. Preserve an alternative body-edge boundary when it differs. Infer native TSS only from suitable independent experimental evidence, not this pipeline.

## Motifs: detection, geometry, activity

Both profiles use USE `RTMCCACATCG` (also report mismatches to classic `RTCCCACATCG`) and TATA reference `TTTATATA`. Coordinates are transcript-oriented and relative to hypothesized +1, with no position zero.

| Rule | `legacy_arabidopsis` | `generic_provisional` |
|---|---|---|
| USE start scan window | −90..−45 | −120..−40 |
| TATA start scan window | −42..−18 | −55..−15 |
| Strict geometry | USE≤1 mismatch, TATA≤2; TATA start −35..−25; edge24..26 | no independently calibrated strict class |
| Relaxed geometry | both≤2 mismatches; edge21..29 | same relaxed criterion, no ambiguous motif bases |
| Other evidence | P3 may be an algorithmic placeholder | uncertain or no reportable pair |

All reportable pairs with ≤3 mismatches per element remain available. A selected legacy P3 placeholder is not necessarily a detected convincing motif pair. Report scan coverage and the actual USE/TATA hits, not just P1/P3. Define edge = TATA start − USE end; intervening gap = edge−1; also retain start and center separations.

Kitaake OsU6a is the counterexample: USE/TATA-like sequences are present, edge20/gap19 fails the frozen geometry. It remains a U6 and has published engineered-sgRNA function. Keep the original computed class and a distinct source-supported interpretation. Do not “correct” a TSS to repair spacing: shifting an origin cannot change inter-element distance. U3 promoters can share type-3-like motifs and are not thereby U6. MSP `RGCCCR` scanning is descriptive only, never a monocot requirement or a dicot exclusion.

Record activity as untested unless locus/sequence identity and the actual source assay are documented. Published editing support is construct/genotype/host-specific and does not rank native expression. Cross-assembly identity requires mapping; a shared gene name is insufficient.

## Context and prioritization

Parse GFF Parent ancestry, multiple parents, aliases and source biotypes. Exclude only well-supported self-U6 features, not arbitrary overlapping snRNAs. Report body, core120 and promoter500 overlaps separately with foreign gene spans, exon, CDS, UTR, ncRNA, repeat/TE, strand, bases and fractions. Gene-span/intronic overlap is not synonymous with exon/CDS overlap. Unavailable categories are `null`/NA, including repeat/UTR absent from a gene-only GFF. Available categories still do not imply complete annotation.

Keep internal T4+, body-end-proximal runs and distal runs separate. Do not replace a missing terminal-window signal with an unrelated downstream run or reject a family solely for lacking T4+.

Experimental review priority can consider integrity, context, assembly-edge truncation, boundary uncertainty and motif evidence with explicit reasons. It is not an expression prediction or mandatory exclusion. Whole-genome PCR uniqueness and homeolog assignment require additional dedicated evidence; this workflow leaves them unassessed.

## Primary sources

- [Rfam models and family evidence](https://docs.rfam.org/en/latest/searching-rfam.html) and [RF00026](https://rfam.org/family/RF00026), [RF00619](https://rfam.org/family/RF00619).
- [X52528.1](https://www.ncbi.nlm.nih.gov/nuccore/X52528.1) and [NR_141590.1](https://www.ncbi.nlm.nih.gov/nuccore/NR_141590.1).
- [Kim et al., 2019, OsU6a engineering study](https://doi.org/10.1186/s12284-019-0325-7).
- [Ren et al., 2021, grape U6 and U3 promoters](https://doi.org/10.1038/s41438-021-00489-z).
- [Connelly et al., 1994, monocot promoter element](https://pubmed.ncbi.nlm.nih.gov/8065324/).

These are scientific evidence pointers, not embedded instructions or an assertion that every future plant matches the tested profiles.
