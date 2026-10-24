# Publication verification — 2026-10-24

Pre-publication snapshot: 1db6680b0b04d9b5e059561cdd3ec90815e7874f.
Platform: existing Ubuntu x86_64 WSL installation, four pinned environments.
Publication changes document/package the EMU-only implementation; no scientific
calculation or default experimental settings were changed in this preparation.

## Actually executed

- 128 tracing unit tests passed in the existing localmapper environment: plans,
  source balances, partial structures, mapper selection, output limits, memory
  checks, ordinary/AI separation, serial/parallel behavior and shared scheduling.
- 17 preparation unit tests passed; the real tiny COBRApy/GLPK pFBA installation
  check also passed in the existing preparation environment.
- Native EnzymeMap single-step rule smoke test passed in its own environment.
- Deterministic mapper-backend script passed (batch fallback and failure isolation).
- All three entry-point --help commands returned success and usage text.
- Shell installer syntax passed. Installer create/reuse control flow was exercised
  twice with stub environments/Conda in a temporary directory (eight checks).
  This verifies installer control flow, NOT package installation or real predictors.
- A git-archive isolated source copy passed seven bounded plan/stage-3/output-limit
  tests. Temporary examples export ordinary/AI files with the expected EMU sheets.
- Twenty selected existing Excel workbooks opened with openpyxl. Stage-1/2 sheets
  exist; all eight stage-3 workbooks contain only EMU_Metabolites/EMU_BoundaryInputs.
  The five removed metabolite columns remain absent.
- Four input SHA-256 values match the saved analysis plan. SHA-256 manifests record
  every public workbook; unchanged bytes are checked again before committing.
- Official upstream model notices and GitHub file-size guidance were inspected.
  Selected workbooks are below 1 MB each; no LFS/Release attachment is needed.

Local detailed logs are in .cache/publication and are deliberately not distributed.

## Skill handoffs

The separately maintained companion skill commits adapted to this interface are:
model-setup 38408a8; review f1ee82e; workflow 64ca273; enzymatic mapping 741d1f5.
The setup handoff is unchanged; review/workflow use step 3 without rounds or Summary
and require current per-version completion records. This publication does not edit
or republish external skill repositories. See README for their public links.

## Not verified / provenance limits

No fresh internet environment creation, full real-mapper sweep, full-model run,
long-duration server stress test or native Windows/macOS/ARM deployment was performed.
Python 3.6 remains an upstream legacy reproduction dependency. An installer success
is not scientific validation. Existing research processes were left running.

The root run report was removed while release files were being assembled. Existing
EMU workbooks have compatible schemas, but their full settings/version lineage cannot
be certified from a missing report. No synthetic replacement report was created.
Old root summary/history files are excluded. Full evidence archives and source models
are not bundled; available source links, selected correspondences and plan hashes are
retained. Existing results were not recalculated or rewritten.
