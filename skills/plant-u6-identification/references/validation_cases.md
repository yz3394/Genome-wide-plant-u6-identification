# Validation cases and interpretation limits

These are known development/regression controls, not future blind validation data. The benchmark references are historical coordinate-compatible assemblies, not current best-reference recommendations.

## Whole-reference baseline before packaging

| Selected nuclear reference | Bases / records | RF00026-supported U6 | strict historical-At subset | separate RF00619 GA |
|---|---:|---:|---:|---:|
| Grape PN40024 12X `GCF_000003745.3` | 485,262,630 / 1,905 | 18 | 9 | 8 |
| Rice Nipponbare IRGSP-1.0 `GCF_001433935.1` | 373,795,655 / 55 | 13 | 11 | 3 |

Pre-existing computational U6 annotations recovered: grape14/14 strictly locatable, rice12/12; U6atac8/8 and3/3. One extra grape annotation is supported by partial flank mapping and separately recovered, not added to the strict denominator. Annotation methods are partly Rfam-dependent, so these are concordance counts, not independent biological accuracy estimates. All 188 RF00026 seed accessions were taxonomically resolved: three Oryza, no Vitis. Rice was held out of the project query panel, not independent of Rfam model training. RF00619 training independence was not established.

## Real counterexamples to preserve

### OsU6a: noncanonical geometry with published engineering support

[Kim et al., 2019](https://doi.org/10.1186/s12284-019-0325-7), supplementary Fig.S1, provides a442-nt Kitaake block. It maps exactly/uniquely to Nipponbare `NC_029258.1:34028645–34029086` on the minus strand, associated with U6 body `34028543–34028645`. The final source G coincides with the RefSeq U6 start: consistent with441 upstream bases plus+1G, not independently mapped native TSS or proven vector-junction semantics.

USE-like `GTACCACCTCG` at−60..−50 and TATA-like `CTTATATG` at−30..−23 are detected. Edge20/gap19 is outside the frozen21–29 relaxed window. Preserve the uncertain geometry and U6-family membership, and record paper support separately. Do not automatically promote it to P1 or call the motif absent. Kitaake's5-nt inter-element deletion relative to93-11 does not establish a universal monocot spacing rule.

### Grape U6 versus U3

[Ren et al., 2021](https://doi.org/10.1038/s41438-021-00489-z): published VvU6.1(425nt) and VvU6.2(591nt) map exactly to two PN40024 chr6 minus-strand upstream regions; their full fragments equal the corresponding1000-bp export suffixes. VvU6.2 does not fit a500-bp window. The source donor is PinotNoir and reported hosts/constructs are not equivalent to native PN40024 expression.

Published VvU3.1/VvU3.2 promoters also have type-3-like geometry but are not canonical U6 at adjacent loci. This tests that motifs cannot assign RNA-family identity. VvU3.2 has511nt in the supplementary sequence but483nt in the paper table; retain the source discrepancy rather than silently trimming to match.

### Coverage-sensitive Rice U6

`NC_029264.1:8177871–8177977`(+), former OsU6L0022, is RF00026-supported and annotated as `XR_003238792.1`. Historical-At HSP coverage is78/102=76.47% under the frozen reward1 scheme, but83/102=81.37% under a separately documented reward2 scheme. A broad-bean HSP78/98=79.59% can be printed as80 by BLAST; exact coordinates must govern the threshold. Preserve the CM-supported locus regardless of that strict-At label. Equal total counts across profiles do not establish equal membership.

### Partial mapping is not strict mapping

Grape annotation `ENSRNA049469705` has a long-flank hit with100% identity but89.655% query coverage due to an unaligned N-rich tail. Its U6 body aligns to `NW_003724846.1:3298–3399`(−) and is separately recovered. It fails the frozen≥99% long-flank coverage requirement; retain that failure and supplementary body support without changing the strict14-control denominator.

### Short body with apparently complete CM span

The public fixture takes the terminal 73 nt of the transcript-oriented public rice OsU6a reference body and aligns it anew to RF00026. It is explicitly synthetic, not a discovered native deletion or pseudogene. `trunc=no` or a complete search model span does not make a submitted body complete; inspect occupied RF columns and terminal missing columns. The 73-nt boundary lesson is preserved without distributing private target-species results.

## Test interpretation

Core tests cover single-HSP decisions, dual-GA competition, weak/nearby U6atac, fixed-anchor clustering, multilevel self annotation, unavailable categories, strands, indels and assembly edges. Compact source-backed fixtures permit inexpensive end-to-end runs after relocation; full genomes and paper PDFs are intentionally not bundled. The release-level QA receipt records which checks actually ran.

Independent forward-testing gives another agent the skill and minimum raw inputs without the intended result. This tests usability and behavior, not a new independent biological benchmark. If future thresholds or motif models are trained with these cases, obtain genuinely new source-held-out controls before claiming improved cross-species generalization.

## Public package validation scope

Version 1.0.1 omits private project whole-reference and promoter tables and the two tests that depended on them. Its public suite contains 120 tests; the historical internal 1.0.0 suite had 122. The public composite contains eight rice/grape reference segments, and the synthetic 73-nt truncation is a separate alignment control. These changes affect packaging and regression fixtures, not the seven scientific runtime scripts. Mini BLASTN and both CM searches were executed again on the public composite before freezing its completed-search fixture.
