# Plant U6 identification skill

`plant-u6-identification` 1.0.1 is a research-grade, reference-bounded workflow for plant nuclear U6 snRNA discovery and upstream-promoter assessment. Its [entrypoint](skills/plant-u6-identification/SKILL.md) separates RNA-family support, submitted-body integrity, boundary confidence, promoter geometry, annotation overlap and experimental-source evidence.

中文概要：使用 U6 转录区 BLASTN、独立全参考 RF00026/RF00619 检索，输出完整候选家族、单独的严格 AtU6 子集和链方向正确的 300/500/1000 bp 上游序列。USE/TATA 是结构线索，不是启动子表达强度或活性的证明。不能保证鉴定出所有功能性 U6，也不能替代实验验证。

## Install on another computer

Ask Codex:

> 使用 $skill-installer，从 yz3394/Genome-wide-plant-u6-identification 安装 skills/plant-u6-identification，指定版本 plant-u6-identification-v1.0.1。先检查是否已有同名技能，保留已有修改，再检查文件完整性和依赖。

Alternatively clone the repository at that tag, then copy the entire `skills/plant-u6-identification/` directory into your user-level `.agents/skills/` directory. Do not copy only `SKILL.md`, and do not copy an absolute symlink from another computer. Preserve an existing installation before replacing it; avoid multiple conflicting copies. Codex's local skill discovery and installer are described in the [official guide](https://learn.chatgpt.com/docs/build-skills).

This repository is a source distribution, not a plugin-marketplace installation. It does not install bioinformatics software or automatically update other computers.

## Dependencies and execution

- Python 3.10+; runtime code uses the standard library.
- NCBI BLAST+ `blastn`; tested baseline 2.16.0+.
- Infernal `cmsearch` and `cmalign`; tested baseline 1.1.5.
- SeqKit for independent extraction QA; tested baseline 2.13.0.

Install suitable builds for your operating system and record their versions. macOS was used for release testing; Linux/Windows execution is not independently certified by this release. Other versions may change search results and need review. No executable binaries are bundled.

Follow the [configuration and coordinate guide](skills/plant-u6-identification/references/schema_and_coordinates.md). Input paths are relative to the config file; executable paths can be PATH names or explicit locations. First select and freeze the genome and annotation; do not assume a bundled benchmark assembly is the best reference for a new study.

```bash
python3 skills/plant-u6-identification/scripts/run_u6.py --config /path/to/config.json --out /path/to/new_run
python3 skills/plant-u6-identification/scripts/validate_outputs.py --run /path/to/new_run --seqkit seqkit --out /path/to/new_qa
```

Existing output directories are refused. Inputs are preserved. Completed-search reuse requires input/output hashes and command-contract validation, not just a `completed` flag. The normal result entrypoint is `README_where_to_find.md`; `04_delivery/` holds candidate tables, BED, RNA bodies, promoters, methods/results text and provenance.

## Verification and public-package boundary

From the repository root, with dependencies on PATH:

```bash
python3 -m unittest discover -s skills/plant-u6-identification/tests -p 'test_*.py' -v
```

Release-specific executed checks are recorded in [the QA summary](releases/plant-u6-identification-v1.0.1-validation.json). A missing-tool skip is not a pass of that tool-dependent check.

For 1.0.1: 120 unit tests passed with zero skips; a fresh eight-segment, 16,882-bp public-reference mini search recovered five U6 candidates and passed 40 case checks. Independent QA passed 100 checks and 21 SeqKit sequence extractions. These are technical checks on selected controls, not genome-wide sensitivity or activity measurements.

The project-local 1.0.0 archive is retained privately. Version 1.0.1 is its public-sharing derivative: the seven runtime scripts and scientific thresholds are unchanged. Private coffee/Catharanthus whole-reference search tables, selected promoter sequences and two private-data regression tests are excluded. Public reference/published rice and grape controls and synthetic edge cases remain. Personal paths in diagnostic provenance are replaced with explicit symbolic source paths; derivative hashes are recomputed and the provenance records this transformation. This must not be described as byte-identical original raw logs.

Whole-reference rice/grape counts in [validation cases](skills/plant-u6-identification/references/validation_cases.md) summarize earlier local technical validation, not searches rerun in GitHub or universal biological accuracy estimates. No promoter-expression rankings, native TSS measurements or wet-lab activity are claimed.

## Updates, sharing and rollback

The GitHub source is canonical for future skill maintenance. The owner's standing request is to publish each verified U6 upgrade here with its version, changelog, tests and immutable tag; see [MAINTENANCE.md](MAINTENANCE.md). This is an upgrade workflow, not a background synchronization service. Runtime analysis never uploads study files.

Other computers explicitly install a chosen tag. To update or roll back, preserve local changes, replace the complete skill with the intended tagged version and rerun a small check. Each research project should record the exact tag/commit, dependency versions and input/output hashes. Never overwrite an old study solely because the skill was upgraded.

Third-party models and sequence resources retain their own provenance and terms; see [resource manifest](skills/plant-u6-identification/resources/manifest.json) and [third-party notices](skills/plant-u6-identification/THIRD_PARTY_NOTICES.md). The public manifest is versioned 1.0.1; the underlying query/model bytes remain frozen.

Original code and documentation are released under the [MIT License](LICENSE), as authorized by the repository owner. The complete skill folder also includes the license so it remains available after installation. MIT does not replace the separate terms/provenance of third-party models, sequences or source publications.
