## Batched carbon-flow output checks (2026-10-16)

The current implementation stages rows in a temporary SQLite database and streams
candidate XLSX workbooks. It discards an incomplete round when the candidate
reaches the configured size limit or exceeds one worksheet's data-row limit.
There is no worksheet pagination. Ordinary and AI-reviewed propagation share this
implementation.

Validated in Ubuntu_Lab using the existing LocalMapper tracing environment:

- `python -m unittest scripts.test_batched_output scripts.test_output_limits scripts.test_memory_guard`: 21 tests passed.
- `python -m scripts.test_carbon_tracing`: passed.
- `python -m scripts.test_stages`: passed, including ordinary/AI files and parallel execution.
- `python -m scripts.test_analysis_plan`: passed with isolated case plans and fake mapper adapters.
- Workflow, model-setup and AI-review skill metadata validation passed for installed and maintained copies.
- `git diff --check`: passed.

Checks compare staged records, ordering, ProductionIDs and convergence with the
existing in-memory algorithm, including mixed sources and CO2/bicarbonate
termination. Exact 50,000,000-byte and overshoot decisions use mocked file sizes;
a real saved XLSX is also tested against a small limit. Size rejection, row
rejection, memory pressure and simulated disk failure preserve the preceding
complete output; rejection in the first round reports no output for the run.

These checks did not rerun real mapper predictions or full research analyses.
The 17 existing three-carbon-source case workbooks match their saved SHA-256
hashes. Six root carbon-flow workbooks differ from the older migration snapshot,
but their modification times precede this batching change; they were not rewritten
by these fixture-only tests. No staging directories remained in results.

The size check occurs at batch/layer/round checkpoints, not after each byte.
Rebuilding candidate workbooks uses additional disk space and I/O. Accumulated
label states and indexes still require RAM, so the independent memory guard can
stop a run. This is not resumable computation and no full-scale performance
benchmark was run for the new writer.

# Validation of the COBRApy development workflow

This file contains historical validation records. Earlier post-save-size and
pagination policies are superseded by the batch-checkpoint policy in the
[current README](../README.md). `scripts/test_batched_output.py` and
`scripts/test_output_limits.py` exercise the current behavior.

The current step 1 uses COBRApy 0.30.0 with GLPK through swiglpk 5.0.12 on
Linux/WSL x86_64. It preserves the reference model conditions, six named blocks,
input schema and version-3 analysis-plan interface. Source MAT files are read
as arrays by COBRApy/SciPy. No alternative implementation is maintained here.

The environment check solves a real small pFBA problem. Model preparation
checks optimal status, finite fluxes, positive biological objective, steady-state
residual below 1e-6 and reaction-bound compliance. The primary production flux
is read from the selected reaction, not pFBA's total-flux objective value.

`scripts/test_model_preparation.py` exercises model loading, pFBA, medium filtering,
export/provenance and failure handling. The existing mapping, propagation, source
numbering, inorganic-carbon endpoint and reviewed-workbook tests remain separate.
These checks do not certify every biological reconstruction or atom mapping.

Real-model comparison results are recorded below after execution. Matching
objective values alone does not imply an identical nonzero network: equally
optimal pFBA representatives may differ. Previous results are comparison data,
not constraints added to force the new solver's answer.

## Observed checks (2026-10-15)

A fresh isolated preparation environment was installed from its specification;
the real COBRApy/GLPK smoke test and `pip check` passed. The default
`python3 step1_prepare_carbon_tracing_inputs.py` command selected that interpreter,
optimized all six real models and generated six inputs plus the version-3 plan.
Seven focused preparation tests and all eight existing deterministic regression
programs passed. The custom iML1515 growth example also exported in isolation.
The existing compatible tracing interpreter was used for downstream regressions;
a clean reinstallation of the three unchanged mapper environments was not run.

The actual six-model checks verified the non-IUPD/IUPD shared fixed growth, doubled
non-IUPD glucose bounds, 2.09-fold IUPD isoprenol bounds, nonfunctional knocked-out
genes and the expected blocked reactions. The unchanged step-2 reader accepted
all six inputs with offline structure completion; input hashes were unchanged.
Shared metabolite names, formulas, charges and available InChI/InChIKey values
matched the reference exports. Calibri fonts and plain tables were retained.

Reference comparison uses rtol=1e-6 and atol=1e-8 for fixed growth, production
flux and total absolute flux. All six comparisons passed those numerical
criteria. Reaction sets are **not identical**:

| Model | Reference reactions | COBRApy reactions | Removed / added |
|---|---:|---:|---:|
| Saccharomyces_cerevisiae_non_IUPD | 493 | 505 | 32 / 44 |
| Saccharomyces_cerevisiae_IUPD | 486 | 492 | 23 / 29 |
| Escherichia_coli_non_IUPD | 434 | 434 | 2 / 2 |
| Escherichia_coli_IUPD | 436 | 436 | 5 / 5 |
| Yarrowia_lipolytica_non_IUPD | 390 | 389 | 4 / 3 |
| Yarrowia_lipolytica_IUPD | 386 | 385 | 4 / 3 |

E. coli non-IUPD, for example, uses PYK2/PYK3 in place of NDPK1/NDPK2 while retaining
the same numerical objectives. The larger yeast differences include alternative
pathways and small fluxes near the export cutoff. Such differences are consistent
with alternate optimal representatives and numerical solver behavior; objective
agreement alone does not certify every position-specific downstream conclusion.
No constraints were added to force an old reaction set. These development
exports must not be called an exact reproduction of previous carbon-flow data.

Detailed per-reaction comparisons, input hashes, versions and the environment
freeze are local runtime records under `data/validation/`, excluded from version
control. The original reference inputs were read without modification. Source
model MAT files are data files; their Python reader needs no external runtime.

Full-network mapper predictions and carbon propagation on the newly selected
networks were not rerun. No model data, environments, credentials or private
research results were added to tracked files. Validation was on Linux/WSL x86_64;
native Windows, macOS and ARM installations remain unverified. No changes were
pushed to GitHub.

## Readability-only step-1 revision (2026-10-15)

The six strain blocks now show COBRApy knockout, reaction stoichiometry, separate
biomass bounds and `pfba()` calls directly. The pFBA wrapper and nested
stoichiometry comprehensions were removed; loading, validation and export remain
shared. The custom-model example follows the same direct pFBA convention.

Before and after this revision, all six actual models were serialized and their
complete pFBA flux vectors compared using the same preparation environment.
Model data, inactive genes and every reaction flux were identical (maximum
absolute flux difference: 0 for each model). These checks compare the two Python
revisions, not the older reference implementation. All seven preparation tests
and static checks passed. Existing input exports and the analysis plan were not
rewritten; full-network mapping and carbon propagation were not rerun.

## Automatic structure discovery (2026-10-15)

Step 1 now records source-model paths, source SHA256 and MAT variable names
automatically. Step 2 reads supported native structures and compatible same-folder
JSON/YAML annotations without an analysis-name-to-YAML configuration table.
The optional E. coli companion setting was removed from step 1; step 2 loads its
InChIKeys automatically instead. Structure completion does not change pFBA bounds,
reaction stoichiometry, objectives, or label-propagation rules.

Eight preparation tests and nine discovery tests passed. The discovery checks
cover native MAT/JSON SMILES, differently named YAML companions, identity matching,
directory isolation, source metadata records, source hash changes, manual/input
priority, equivalent SMILES, conflicting structures and identifier-based lookup.
All eight existing deterministic regression programs also passed, including the
six-model plan-to-mapping-to-propagation integration check and step-4 settings.

For the six real reference networks, complete pFBA flux vectors were identical to
the preceding Python revision (maximum absolute difference 0). Model data were
identical except that E. coli companion InChIKey enrichment moved to step 2.
Offline automatic structure selection matched the previous explicit-YAML selection
for every exported metabolite in all six networks; E. coli InChIKeys were also
restored exactly when the old exported identifier column was removed in a test.
Input Excel bytes were unchanged. The existing generated plan was refreshed with
source provenance without rewriting its workbooks. Detailed local-only evidence
is under `data/validation/auto_structure_*`.

No new dependencies, online structure retrieval, real mapper predictions, or
full-network carbon-flow runs were needed for this revision. Native SBML structure
notes and arbitrary YAML layouts are not parsed by the automatic annotation reader.
Unsupported structures can still use the existing database/manual workflow.
No original test-project files were changed and no changes were published.

## Per-strain carbon-source settings (2026-10-15)

Each step-1 strain block now defines its own ordered `carbon_sources` list and
exports it together with that model's solution. The separate global source table
was removed. All nine preparation tests passed, including two independent source
lists with different introduction orders. The six real models retained exactly
the same source lists and complete flux vectors (maximum absolute difference 0).
No existing input, mapping or carbon-flow workbook was rewritten.

## Uniform model blocks and validation scope (2026-10-15)

All six strain blocks use the same reading order: model loading, ordered carbon
sources, culture/model modifications, growth constraints, production pFBA and
registration for export. Fixed-version enforcement was removed. non-IUPD biomass
objectives are assigned explicitly instead of requiring a matching source objective.
The Yarrowia glucose setting is local to each block rather than a shared constant.

The generic requirement that a gene knockout must change reaction bounds was
removed; redundant or already inactive genes can be valid knockouts. COBRApy still
raises on missing gene/reaction IDs. Maximum-growth solutions are validated before
use. Duplicate production checks were removed from the model blocks; the existing
exporter still verifies optimal status, finite fluxes, steady state and bounds
before creating output.

All ten preparation tests passed, including invalid-solution rejection before any
output was created. A fresh comparison of all six real models found identical
serialized model data, inactive genes, ordered source lists and complete pFBA flux
vectors (maximum absolute difference 0). Local evidence is in
`data/validation/uniform_blocks_*`. Existing workbooks and the generated analysis
plan were not rewritten. Mapping and propagation were not rerun; no changes were
published.

## non-IUPD terminology (2026-10-15)

The comparison group now uses `non-IUPD` in descriptions and `non_IUPD` in code
identifiers, analysis names and filenames. The three existing comparison input
workbooks were renamed without rewriting them; all input SHA256 values remained
unchanged. Manual structure scopes and the generated analysis plan were updated
together. Chemical structures, source order and source-model provenance were
preserved.

Reversing the name substitutions produces the identical parsed step-1 computation.
All ten preparation tests and eight existing regression programs passed. Offline
checks on the six actual inputs verified identical selected structures and review
flags before/after renaming, unchanged source lists, and the expected normal/AI
filenames. Details are recorded locally in `data/validation/non_iupd_*`.
Real-model optimization, mapper prediction and full-network propagation were not
rerun for this naming-only revision. No reference-project files were changed or
published.

## Runtime memory protection (2026-10-15)

Run the focused safety tests with:

    .envs/localmapper/bin/python -m scripts.test_memory_guard

All seven tests passed: reserve/temporary-allocation checks, throttling and
unreadable-memory failure, nested cgroup limits, stopping inside combination
enumeration, unchanged propagation results, preservation of saved round 1 when
round 2 stops, and guarded Excel writes/step-4 setting inheritance.

Existing carbon-tracing and script-settings regressions passed. The six-model
analysis-plan fixture passed with mocked safe-memory readings. The stage/AI-copy
regression passed with mocked safe-memory readings and serial execution forced;
the parallel-process path was not rerun while a live analysis occupied about
33 GiB. A real low-memory reading (about 0.9 GiB free) also stopped a fixture
before loading a mapping workbook, as intended.

The committed pre-change propagation implementation was compared against this
revision on a small mixed-source cycle with repeated stoichiometric copies and
both source introduction orders. All table values, row order, ProductionIDs,
new-state counts and convergence results matched exactly (10, 11 and 12 rows
in rounds 1, 2 and 3 for each source order).

These checks used synthetic fixtures, not six full real-model propagation runs.
No real mapper was rerun, existing research workbooks were not overwritten,
and the active analysis was not interrupted. No changes were published.

## Strict output limit and Summary collection (2026-10-15, historical)

The strict rollback policy tested in this section was replaced on 2026-10-16
by the 50 MB post-save stopping policy described in the current README. The
comparisons below document the earlier implementation and preserved outputs.

Run the focused checks with:

```bash
.envs/localmapper/bin/python -m scripts.test_output_limits
```

Five tests cover exact byte limits and rejected-save preservation, complete-round
rollback, first-round rejection without reusing stale outputs, Calibri,
separate Original/AI Summary copies, source-workbook preservation, step-4 setting
inheritance, and collecting saved rounds before re-raising a memory-stop error.
The eight related regression scripts also passed: output_limits, memory_guard,
carbon_tracing, script_settings, analysis_plan, stages, parallel_carbon and
consensus_mapping. They include real concurrent processes and synthetic AI
mapping changes; no external mapper was called by propagation.

All six real step-2 workbooks were propagated for one round in an isolated
validation directory. Every first-round carbon-flow value and original Summary
count matched the existing output exactly. Input and mapping file hashes were
unchanged. This was not a full rerun of later rounds or real AI-reviewed models.

Existing ordinary outputs were preserved under
`results/archive/before_size_limit_fit_20261015_213122/`. Their already-computed
rounds were used to produce bounded workbooks without another mapper
or multi-round propagation run. The Saccharomyces_cerevisiae_non_IUPD output
retains rounds 1-8 (about 16.65 MB) rather than rounds 1-9 (about 273 MB);
round 9 was removed in full, including all its pages. The other five outputs
retain all ten rounds. Retained flow-sheet uncompressed sizes/CRCs and all
original Summary observations were verified against the archive. The copied
Summary sheets match their per-model sources, and all new outputs are below
100,000,000 bytes.

Automatic curve fitting was subsequently removed at the user's request.
Current Summary sheets contain only calculated observations. The temporary
function/statistics columns were removed from the six current workbooks, and
the aggregate was regenerated. Their saved rounds, original Summary values,
Calibri formatting and all flow-sheet contents remain unchanged. The output
limit no longer reserves space for fitting metadata.

Internal historical records in `data/validation/flow_limit_fit/` retain the earlier
one-round comparisons. Removal checks and current regression logs are under
`data/validation/remove_auto_fitting/`; both directories are ignored by Git.
Historical archives are preserved. The old test project was not changed and no
GitHub push was performed.

## Separate analysis plans and automation interfaces (2026-10-15)

The default six-model script and no-argument stage commands remain available.
An optional copied COBRApy recipe can generate a separate analysis plan and
results under `results/results_by_AI/<name>/`; steps 2 and 3 select it with `--plan`. Each analysis
shares the installed environments and original model data. Generated plans add
preparation provenance, exchange conditions and actual source uptake. Propagation
adds one batch status report, without changing Excel tables or tracing rules.
Stage scripts contain no LLM/skill calls. The documented future orchestrator
includes external mapping review by default, followed by step 3 for normal and reviewed mappings.

Thirteen preparation tests passed, including multiple allowed carbon exchanges,
an unused supplied source, isolated plan invalidation and preparation-script
hashes. The analysis-plan suite exercised two separate analyses through mocked
mapping and two propagation rounds, checked a zero-uptake source without changing
the mapped-network denominator, and verified interpreter relaunch arguments.
Under the then-current strict size policy, it also verified that a rejected
first round reports zero current saves and preserves an older workbook. The
five output-limit tests passed with the new
size-stop and memory-stop report checks.

Nine other regression programs passed: memory_guard, carbon_tracing,
script_settings, stages, parallel_carbon, consensus_mapping, model_structures,
smiles_completion and mapper_backends. These include concurrent propagation and
AI-copy fixtures; real external mapping predictors were not rerun.

The new glucose/xylose example was solved with real COBRApy/GLPK in isolation.
The single-source model used 10 mmol/gDW/h glucose; the simultaneously supplied
model used 5 glucose and 6 xylose. Biological growth objectives were approximately
0.876997214 and 0.871354700, respectively. This validates input preparation,
not a new full-network mapping or carbon-propagation result.

All six study inputs were independently regenerated using a selected validation
plan. `Nonzero` and `Metabolites` tables matched the existing exports, including
row order and numerical values within 1e-9 tolerance. The default plan and all
12 existing root-level result workbooks retained identical SHA256 hashes.
Preparation comparison and fixture logs are kept locally under the ignored
`data/validation/case_interfaces/` directory. The original test checkout was
not modified; nothing was pushed to GitHub. No new skill was created.

## Post-save 50 MB stopping threshold (2026-10-16)

Steps 3 and 4 now save each complete round before checking its file size. The
default threshold is 50,000,000 bytes. A round reaching or exceeding it is kept
in full; no next round starts for that model. The run report records
`file_size_limit` ahead of a simultaneous convergence condition. This replaces
the historical strict rollback behavior without changing propagation, mapping
selection, Summary counting or memory safeguards.

All seven output-limit tests passed, covering exact-threshold, overshoot and
first-round size stops, below-threshold continuation, step-4 setting inheritance,
failed-save preservation, Summary copying, and memory-stop reporting. Boundary
tests mock observed sizes while saving and reading actual small workbooks; no
large workbook was created. Separate generic capped-save behavior remains tested.
The analysis-plan fixture also used a real small threshold to verify that round
1 is retained, its unchanged flow/Summary values are copied, and the report
records one saved round and its output path.

Eight regression programs passed: output_limits, analysis_plan, memory_guard,
carbon_tracing, script_settings, stages, parallel_carbon and consensus_mapping.
This includes concurrent propagation and reviewed-copy fixtures. Original model
inputs, the default analysis plan, the reference preparation script and existing
root-level workbooks were not rewritten. Logs are kept locally under
`data/validation/post_save_size_50mb/`. No real mapper or full research-model
propagation was rerun, and no GitHub push was performed.
