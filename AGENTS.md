# Project handoff for local agents

Read README.md, INSTALL.md and docs/interfaces.md before running stages. Existing
workbooks are indexed in results/README.md. Preserve user data and running jobs.
Do not launch study computations merely to install or inspect the project.

Choose the starting point from available material: research question/model -> step 1;
prepared flux Excel + plan -> step 2; mapped Excel + compatible plan -> step 3.
Use explicit --plan for custom cases under results/results_by_AI/<case>/.
Never replace a user's plan/results as an installation check.

Step 3 is size-1 EMU-inspired fixed-flux source balance only, not a full isotopomer
simulator. No enumeration, rounds or Summary collection. Keep compartments and Unknown.
See docs/interfaces.md for workbook contracts, mapping hashes and completion checks.
Read current carbon_flow_run.json and process exit status; require emu_completed,
not merely emu_solved or the existence of an old workbook. Resource failures are not
convergence. Ordinary and _AIGenerated outputs are separate.

Ordinary scripts require no AI. Companion skills are optional: model-setup prepares
conditions, review produces reviewed mapping copies and runs step 3, enzymatic mapping
establishes reaction correspondences, workflow orchestrates stages. Repository links
are in README. Follow installed skills only when compatible with the current interface.
Retain original mapper sheets/Mapping_Selection during review; edits go to Final_Mapping
in a separate copy. Save evidence, hashes, plans, versions and actual run status.

Use bounded temporary-directory tests. Never run all predictors or replace published
results for a smoke test. Do not upload environments, caches, locks, credentials or
host paths. Results are computational outputs, not fermentation measurements.
