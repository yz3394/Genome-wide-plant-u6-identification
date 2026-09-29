# Maintaining the skill source

Source repository: https://github.com/yz3394/Genome-wide-plant-u6-identification

For a requested implementation fix or upgrade, use the repository source rather
than silently changing only an installed copy. Follow repository `AGENTS.md` and
`MAINTENANCE.md`: reproduce the problem, preserve evidence boundaries, add a
targeted regression, run relevant checks, update version/changelog, and verify
the remote commit and new immutable tag after publication. Preserve study inputs,
outputs, old releases and unrelated user edits.

The repository owner requested synchronization of future verified upgrades.
This applies to that owner's authorized maintenance work, not automatic upload
of another user's changes/data. Other users should follow their own repository
and publishing authorization. Ordinary analysis and read-only review do not
trigger publication; runtime scripts have no GitHub upload behavior.

If validation, authentication, connectivity or merge conflicts block release,
report `not synchronized` and retain the work. A local commit, pushed source,
Release and installed version are separate states. Updating GitHub does not
automatically update installations on other computers.
