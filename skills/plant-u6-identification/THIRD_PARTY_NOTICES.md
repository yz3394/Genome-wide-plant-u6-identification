# Third-party resources and public fixtures

Original workflow code/instructions and external data are different materials.
No license for this skill's original code transfers rights over third-party data.

- `resources/models/`: frozen Rfam 15.1 RF00026 and RF00619 covariance models.
  Rfam states that its data are available under CC0; retain Rfam attribution and
  accession/version/hash provenance. [Rfam license](https://docs.rfam.org/en/latest/#license).
  The files are unchanged, including provider-generated headers/build paths.
- `resources/queries/`: small public nucleotide records, with accession,
  coordinate and source provenance in the resource manifests. NCBI places no
  restrictions on use/distribution of molecular data but warns that submitters
  may assert rights; this is not a blanket transfer of third-party rights.
  [NCBI data policies](https://www.ncbi.nlm.nih.gov/home/about/policies/#data).
- `tests/fixtures/`: short public rice/grape reference extracts, published
  sequence controls, synthetic cases and derived tool outputs. Exact source
  mappings/citations are in `source_manifest.json`, `known_cases.json`
  and the fixture provenance records. Full genomes and paper PDFs are excluded.
  The synthetic 73-nt truncated RNA is a software test, not a natural locus.
- Published sequence controls cite [Kim et al. (2019)](https://doi.org/10.1186/s12284-019-0325-7)
  and [Ren et al. (2021)](https://doi.org/10.1038/s41438-021-00489-z).
  Preserve their source attribution; no ownership of those sources is claimed.
- External BLAST+, Infernal and SeqKit executables are not distributed. Install
  them from their providers and comply with their respective terms. Version/help
  text in fixture logs identifies the tool used, not a bundled software license.

Personal paths in public derivative diagnostic logs are replaced by symbolic
`/source/` paths and hashes recalculated, as recorded in public fixture provenance.
These are sanitized derivatives of actual executions, not original unmodified
raw logs. The original private project archive remains outside this repository.
