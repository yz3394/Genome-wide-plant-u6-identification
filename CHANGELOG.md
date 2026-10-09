# Changelog

## Documentation update — 2026-10-09

- Translate the remaining Chinese README summary, installation prompt and repository maintenance instructions into English for sharing with colleagues.
- Keep the skill at version 1.0.1: the skill instructions, runtime scripts, scientific thresholds, resources, fixtures and historical validation records are unchanged.

## 1.0.1 — 2026-09-28 — first public distribution

- Publish `plant-u6-identification` in the owner's dedicated repository, with installation, dependency and future-release instructions.
- Release original code/documentation under the owner-approved MIT License; preserve separate third-party terms and scientific attribution.
- Preserve all seven runtime scripts and frozen Rfam models from the private 1.0.0 baseline byte-for-byte. No discovery, curation, coordinate, motif or annotation-conflict threshold changed.
- Exclude private coffee/Catharanthus whole-reference raw search tables, selected promoter sequences and two associated private-data tests. Keep the original local release untouched.
- Use eight public rice/grape mini-reference segments; replace the private short-body case with an explicitly synthetic 73-nt truncation of public OsU6a. Recompute its actual `cmalign` evidence instead of inheriting the old case's metrics.
- Rerun all three mini-reference searches, sanitize personal paths in derivative provenance, and recompute derivative hashes. Provider-supplied Rfam CM files retain original bytes.
- Public checks: 120 unit tests passed, zero skips; fresh mini run has five U6 candidates and 40 passing case checks; independent QA has 100 passing checks and 21 SeqKit extractions. Official skill-format validation passed.
- Document that public code is not a promoter-activity, expression or exhaustive functional-gene classifier. Cross-platform execution and new-species biological accuracy remain unverified.

## 1.0.0 — 2026-09-28 — private project baseline

The project-local package and detailed research QA remain preserved outside this repository. It included 122 tests, including two regressions that require private study tables. It is not distributed here. Earlier full-reference rice/grape validation is described with its limitations in the packaged references, not claimed as a new GitHub execution.
