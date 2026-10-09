# Upgrades and GitHub synchronization

This repository is the source for ongoing maintenance of `plant-u6-identification`. On 2026-09-28, the owner explicitly requested that future verified skill upgrades be synchronized to this GitHub repository. This agreement does not provide a background synchronization service or authorize uploading research data.

1. Synchronize with the remote repository and inspect the working tree, installed version and the user's uncommitted changes. Preserve existing content.
2. Record the issue, dependency versions, a minimal reproduction with documented provenance, and observed versus expected behavior. Distinguish software defects, methodological rules and species-specific exceptions.
3. Make a scoped revision and add a regression test. Keep U6/U6atac classification, single-HSP thresholds, integrity assessed from actual model alignments, coordinates, motifs and annotation conflicts as separate evidence layers.
4. Run the full public unit-test suite, any affected small searches using real reference sequences, and independent SeqKit QA. When membership decisions or coordinates change, compare old and new members, boundaries and promoter sequences. Substantial method changes also require independent validation cases that were not used to tune parameters.
5. Check the public scope and provenance: exclude research-project tables, whole genomes, PDFs, tokens and personal paths; label synthetic controls explicitly. Preserve third-party sources and licenses. Do not automatically treat a new model version as an equivalent replacement.
6. Update `metadata.version` in `SKILL.md`, `CHANGELOG.md`, the README and the version-specific QA summary. Generate the corresponding ZIP and SHA-256 file manifest. Check the local package, remote source and actual installed copy separately.
7. Commit and push to `yz3394/Genome-wide-plant-u6-identification`, and create the `plant-u6-identification-vX.Y.Z` tag/Release. Do not move existing tags or overwrite historical research results. Browser uploads are also acceptable; afterward, synchronize the commits created through the browser to the local repository.
8. Read back the remote file set, hashes, commit and tag before reporting the synchronization result and installation status. If network access, authentication, conflicts or validation failures block completion, preserve local changes and explicitly report `not synchronized`; do not skip checks.

Version semantics: use a patch version for compatible fixes or distribution revisions, and a minor version for compatible capability extensions. Use a major version for incompatible changes to default rules or outputs, and explain how to compare earlier analyses. A higher version number is not evidence of greater biological accuracy.

Suggested issue-record fields: version/commit, species and reference version, shareable inputs or hashes, reproduction command, expected/observed results, scope of impact, fix, QA/unverified checks, membership/sequence differences, and GitHub and installation status.

Other computers explicitly select a tag to install or update; existing installations are not overwritten automatically. Preserve local changes before restoring an older tag. For each analysis, record the skill tag, Git commit, software versions and actual file hashes.
