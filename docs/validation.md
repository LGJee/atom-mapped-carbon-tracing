# Release validation

Validation performed on 2026-10-13. Results below separate actual software
execution from deterministic tests using mock adapters.

## MATLAB reference export

Executed the published reference script with MATLAB R2025b, the installed
COBRA Toolbox and GLPK, using the unmodified official Yeast9 v9.1.1 model.
The upstream YAML was also retrieved and compared byte-for-byte with the study
copy; it contains 1,787 SMILES annotations. Its SHA256 is
`2ddfaaa73ead41ab243127d1b34fa56db2c86052f6eabcd12bb4eea93947dc17`.

Both new workbook tables matched the original study inputs, including fluxes
within numerical comparison tolerances, equations, metabolite identities and
SMILES overrides:

| Strain | Nonzero reactions | Participating metabolite rows |
|---|---:|---:|
| WT | 493 | 522 |
| IUPD | 486 | 516 |

The original model, code, inputs and result files were not overwritten.

## Scientific behavior before and after refactoring

The original and release implementations read the same full WT/IUPD inputs
and the same previously generated three-tool candidates. No new predictions
were substituted during this comparison. Dataframes were compared, not merely
their sizes. Selected correspondence strings, mapper choice, support groups,
review/AI flags, fluxes, both rounds' state IDs, labels, labelled SMILES and
Summary values matched.

| Strain | Usable mappings | Round 1 rows | Round 2 rows | Mapped-network carbon metabolites |
|---|---:|---:|---:|---:|
| WT | 335 | 252 | 334 | 334 |
| IUPD | 335 | 266 | 342 | 337 |

This verifies preservation of the implementation's scientific behavior. It is
not biochemical validation of all candidate mappings.

## Real three-tool example

Executed the supplied synthetic acetaldehyde example with actual LocalMapper,
RXNMapper and EnzymeMap calls, followed by one round of carbon propagation.
EnzymeMap used the downloaded official v2 CSV and SHA256-verified rule file.
All three produced usable chemical-reaction mappings. The transport was handled
separately by identity correspondence. Consensus selected EnzymeMap for the
reduction; the labelled products were `acetaldehyde[c]` and `ethanol[c]`, each
carrying C1 and C2. `AIGenerated` remained false.

The example is synthetic and does not assert pFBA-derived fluxes or experimental
isotope measurements. The full biological networks were not freshly predicted
again solely for the refactoring comparison.
The separate `tests.test_enzymemap_rules` integration check forced an empty
database lookup and successfully executed the official single-step rule,
restoring the supplied structures with complete carbon correspondence.

## Deterministic regression checks

The five deterministic test programs (separate from the native-rule integration) cover:

- Exact substrate-state pairing, cleavage, multiple product copies and repeated
  carbon-origin labels.
- Structural alignment, symmetric equivalents, complete carbon coverage,
  selection priority and preservation of raw mapper predictions.
- Batch failure and individual retry, isolated errors and validated transport.
- Six concurrent mock prediction processes with a completion barrier.
- Separate strain propagation processes and equivalence to serial execution.
- Prediction cache validation, flux direction and input-overwrite protection.
- Equivalent combined and separate map/trace execution; trace makes no mapper
  calls and removes stale higher-round sheets.
- Rejection of duplicate IDs, zero/nonfinite flux and conflicting structures.
- Source-specific Summary counts and mapped-network denominator rules.

All checks passed. Python files passed syntax, unused-import and unused-variable
checks; the setup shell script passed `bash -n`. Example input/output sheets were
rendered and inspected. Copies retain plain styles, readable column widths,
wrapped equations and no added colour fills, filters or frozen panes.

## Environment and platform scope

The checked execution platform was Linux x86_64 under Ubuntu WSL on Windows.
The three tool environments are isolated. Critical tracing versions are pinned
to the original study stack; exact tool/model/data references are in the
environment files and main configuration.

The first clean RXNMapper installation exposed an incompatibility with
setuptools 84, which removed `pkg_resources`. Both the RXNMapper and EnzymeMap
specifications now pin setuptools 80.9.0. This is required by the selected
legacy tool stack and should not be removed without rechecking imports.
The clean EnzymeMap installation also required an explicit `requests=2.32.5`
dependency, which was implicit in the original workstation. This dependency
is now included in the environment specification.
All three isolated environments were newly created for release validation,
their functional checks passed, and the real three-tool example was rerun
using these environments. Existing workstation environments remained intact.
Rerunning the setup script reused all three environments and passed without
creating replacements. New-environment raw prediction sheets, selections,
final mappings, carbon flow and Summary matched the original-environment
example results.

Native Linux on a separate machine, native Windows Python, macOS, ARM and GPU
execution were not independently tested. Environment specifications support
Linux/WSL x86_64 only. The research tracing environment uses end-of-life Python
3.6 to preserve the checked serialization behavior.

## Code cleanup

Removed machine-specific Windows-to-WSL launch switching, personal interpreter
paths, the legacy nonzero-sheet mapping writer, old split-tool output lookup,
obsolete IUPD split-source sheet cleanup, unused root/output constants and unused
imports. Consolidated atomic workbook replacement and plain column-width
handling. Retained adapters, input/structure validation, failure propagation,
resource limits and the original label-transfer algorithm.

No AI review stage, trained mapper, web service, additional scientific correction
rules or skill package was added.

## Local reference-entry check (2026-10-14)

The test checkout renamed `pFBA_network.m` to
`prepare_carbon_tracing_inputs.m` without changing its file contents. The
scientific defaults in `config.json` and all Python implementation files
remained unchanged.

Executed the renamed script using Windows MATLAB R2025b, COBRA and GLPK,
accessing the Ubuntu_Lab project through the WSL shared folder. The official
Yeast9 v9.1.1 MAT and YAML were downloaded into the test project's `model/`.
Both output workbooks matched the original study inputs: reaction IDs,
equations, fluxes within 1e-10 comparison tolerances, metabolite identities
and SMILES. WT has 493 reactions and 522 metabolite rows; IUPD has 486
reactions and 516 metabolite rows.

The default Python configuration resolved both exported files, accepted
their sheets and retained the reference source order and 10-round limit.
Documentation links and section anchors passed checks. Full-network atom
mapping and ten-round propagation were not repeated for this rename and
documentation-only change; no computational algorithm was modified.
These changes remain local pending publication approval.

## Local script-settings check (2026-10-14)

The two Python entry files now own their stage-specific settings; external JSON
configuration files were removed. Step 2 does not require carbon sources or
propagation rounds. Step 3 independently reads the original reaction/structure
tables and Final_Mapping from step 2, validates the correspondences, and writes
only Summary and Carbon_Flow sheets without changing the mapping workbook.
Its model-process/memory limits now govern the propagation process pool.

The label-transfer and mapping-selection algorithms were not changed. Stage
regressions compare combined versus separate computation, protect the mapping
workbook byte-for-byte, and trace successfully while the original step-1 input
is temporarily unavailable. A settings test verifies stage ownership and
command-line precedence. The synthetic real-tool example exercises the new
entries using --example, with no configuration file.

Observed checks for this local change: all six deterministic programs and the
native EnzymeMap rule integration passed; the actual three-tool synthetic
example passed through both new entries. Direct two-model propagation ran in
separate processes and matched serial results. The two-round example's full
Carbon_Flow tables and Summary matched the pre-change workbook exactly.
Python 3.6 syntax and unused-code checks passed. Active documentation links
and anchors passed; the changes have not been published to GitHub.

## Local model-first filename check (2026-10-14)

Result filenames now follow `<model>_<step>_<contents>.xlsx`, keeping one
model's three stages together when sorted by name. MATLAB export paths,
Python readers/writers, example fixtures, tests and documentation were updated.
Existing example workbooks were renamed with SHA256 checks confirming their
contents were unchanged. All seven test programs passed, followed by the
renamed two-round example. No scientific algorithm was changed or full
biological-network run repeated; these changes remain unpublished.

## Local automatic-input check (2026-10-14)

Step 2 now discovers `<model>_1_*.xlsx` in `INPUT_DIR` and infers the two input
sheets from their required headers. Per-model filenames and sheet settings were
removed from its entry script. Step 3 reuses the same sheet-detection helper.
Duplicate model files, ambiguous sheets and missing exports fail explicitly;
lock files, temporary workbooks and later-stage outputs are ignored.

All seven test programs passed, including custom model names and sheet names,
model selection, duplicate/empty input checks and the independent example path.
The real-tool synthetic example reused cached predictions and passed mapping
selection and two propagation rounds. Syntax/unused-code checks passed. No
scientific algorithm changed; full WT/IUPD mapping and propagation were not
rerun for this input-discovery change. Changes remain local and unpublished.

## Local stage-3 discovery and settings check (2026-10-14)

Step 3 now detects `<model>_2_atom_mapping.xlsx` and defaults to step 2's output
folder and tracing interpreter. Per-model mapping filenames were removed.
User-selected source lists remain in CARBON_SOURCES at the top of step 3;
missing source settings raise an error before propagation. Input/output,
analysis, resources and installation settings are grouped above execution code.
MATLAB file names, exchange/objective reaction IDs, model names, reference gene
settings and culture conditions are likewise centralized without changing defaults.

All seven regression programs and the cached synthetic mapping/two-round
propagation example passed. Native Python 3.6 syntax and unused-code checks
passed. An additional MATLAB export attempt was blocked during COBRA's GitHub
submodule initialization by a connection failure; the revised MATLAB export was
therefore not runtime-verified in this check. Full WT/IUPD mapping and propagation
were not repeated. No scientific selection or propagation algorithm changed.
These changes remain in the test checkout and have not been published.


## Local structure-completion check (2026-10-14)

Step 2 now applies model-scoped verified overrides from data/manual_smiles.csv
and optionally queries PubChem for missing structures. Exact PubChemCID/InChIKey
lookup is preferred; name-selected candidates always require review. Source,
CID, compatible candidates and lookup status are retained in Metabolites, and
affected reaction review flags reach both selection and final mapping tables.
Raw input workbooks are protected and the structure cache is excluded from Git.

Eight test programs passed, including identifier/name selection, structure
priority, model scoping, cache reuse, transient failure recovery, reaction-review
propagation and two-round tracing of enriched structures. An actual PubChem
request returned isoprenol CID 12988 with CC(=C)CCO. The real three-tool synthetic
mapping example and its two-round carbon-flow/Summary regression passed.

In manual-only mode, moving the MATLAB overrides preserved all 522 WT and 516
IUPD effective structures and all 493 WT/486 IUPD prepared reaction SMILES and
statuses against the original study inputs. No full WT/IUPD online-completion
or mapping/propagation run was repeated. MATLAB runtime export remains untested
for these edits because the earlier COBRA initialization failed at its submodule
download. Python 3.6 syntax and unused-code checks passed. Other databases and
cross-model chemical identity accuracy are not claimed to be validated.
Changes remain local pending publication instruction.


## Real iML1515 no-SMILES adaptation check (2026-10-14)

Ran examples/step1_prepare_custom_model_inputs.m with Windows MATLAB R2025b
and the already-installed COBRA/GLPK solver. The local test harness validated
GLPK directly without retrying the network-dependent COBRA submodule setup.
Model conditions: glucose uptake limit 10, oxygen limit 1000 mmol/gDW/h, other
carbon-containing boundary uptake closed, original inorganic medium/maintenance
bounds retained, maximum biomass followed by minimum total absolute flux.
Growth was 0.876997214/h; input contained 432 active reactions and 439 metabolite
identities, formulae, charges and available InChIKeys, with no SMILES columns.

The initial CLI rejected the input because it required legacy SMILES fields.
Made structure columns optional, added generic SMILES_Model support and retained
legacy SMILES_Yeast9 precedence compatibility. Metabolite ID/name remain required.
Real retrieval then exposed PubChem's charged formula syntax (C3H5O6P-2 versus
model formula C3H5O6P and separate charge -2). Comparing elemental formula without
the trailing charge, and retaining separate charge validation, raised completion
from 66/439 to 365/439. No retrieved structure or mapping rule was altered.

All three real mappers completed. Final_Mapping includes 300/432 reactions;
197 selected reactions require review. Selected sources: EnzymeMap 196,
RXNMapper 61, Identity transport 38, LocalMapper 5. One propagation round produced
251 records, 187 labelled states and 178 distinct labelled carbon metabolites.
Summary reports 178/293=60.75%, with the denominator restricted to carbon
metabolites in the usable mapped network, not the whole downloaded model.
The remaining 74 metabolites have no compatible PubChem structure under the
current identity/formula/charge constraints; mapping and biochemical correctness
are not claimed for every model reaction.

Verified the source input SHA256 remained
b4bc462d34415fbca852bb45ffc6f5d494445bac76e348c14390573552dac308.
Verified product labels occur only on carbon atoms, labels agree with the
LabeledSMILES atom-map fields and use the glucose C1-C6 range, and each
ProductionID identifies a consistent metabolite/labelled structure.
Eight regression programs, Python 3.6 syntax and unused-code checks passed.
Yarrowia was not analyzed in this check. Changes and outputs remain local;
no GitHub publication or experimental isotope-abundance validation was done.

## Single-entry configuration and six WT/IUPD inputs (2026-10-14)

All per-analysis settings now live at the top of step1_prepare_carbon_tracing_inputs.m.
MATLAB generates data/analysis_plan.json; the Python entries load it without
parsing MATLAB source. Portable project-relative paths and exact model lists
are used. A SHA256 links each exported input to stage-2 Mapping_Info; stage 3
rejects a same-name mapping from a different export. Neither Python stage alters
MATLAB inputs. Stage 3 does not invoke or initialize mapper tools.

Executed step 1 with actual MATLAB R2025b, installed COBRA and validated GLPK:

| Model | Objective | Nonzero reactions | Participating metabolites |
|---|---:|---:|---:|
| WT | Squalene, 0.163903667 | 493 | 522 |
| IUPD | Squalene, 0.348333333 | 486 | 516 |
| Ecoli_WT | Growth, 0.876997214/h | 432 | 439 |
| Ecoli_IUPD | Growth, 0.877225831/h | 429 | 436 |
| Yarrowia_WT | Growth, 1.05691644/h | 388 | 413 |
| Yarrowia_IUPD | Growth, 1.07420521/h | 386 | 412 |

The Yeast9 pFBA tables and effective structures match the pre-refactor inputs.
Legacy input workbooks still had some manual SMILES cells; those now reside in
the model-scoped manual CSV and are applied by step 2. The old input files were
preserved outside the test checkout for comparison. No scientific carbon-label
propagation or consensus rules were changed by this handoff refactor.

Confirmed MECDPDH5 is inactive after E. coli ispG deletion (b2515), and iYali
y000559/y000560 are inactive after ERG13 deletion (YALI0F30481g). All IUPD inputs
include active isoprenol exchange/transport and both introduced phosphorylation
reactions. These checks validate the configured GEM simulations, not actual
strain viability or every alternative pathway in biology.

The new deterministic six-model test executes mapping and two-round propagation
with an existing test adapter, verifies plain output sheets and unchanged input
and mapping files, checks subset/round overrides, excludes unrelated old files,
and rejects stale same-name mappings. The installation example was also freshly
mapped with all three real tools. Nine regression programs and static Python
checks pass. Full three-tool mapping and multi-round propagation of all six
genome-scale networks were not repeated in this configuration refactor.
Only the local test checkout was changed; no GitHub publication occurred.

## Scientific species names (2026-10-14)

Renamed the S. cerevisiae analysis defaults from WT/IUPD to
Saccharomyces_cerevisiae_WT/Saccharomyces_cerevisiae_IUPD, including
the manual-structure scopes and generated plan. Historical validation logs
above retain their original names. Input reaction worksheets now use Nonzero,
because full species names plus a suffix exceed Excel's 31-character limit.
The six-model stage test also exercises the full species names. Model conditions,
mapping selection and carbon-label propagation rules are unchanged.

## Squalene objective for all six models (2026-10-14)

The earlier growth-only defaults for E. coli and Y. lipolytica were replaced:
all six comparison models now fix growth at 50% of each strain's maximum,
maximize squalene formation and then minimize total absolute flux. Yeast9 uses
r_1012 and iYali uses y001012. Since iML1515 lacks SQS, both E. coli comparison
models receive functional heterologous SQS activity (RHEA:32295) and a squalene
output; the WT label does not mean an unmodified native squalene-producing strain.
The added SQS reaction is elemental- and charge-balanced for the native iML1515
substrate forms. No specific enzyme construct or experimental validation is implied.

Real MATLAB R2025b/COBRA/GLPK execution in an independent validation copy passed
for all six models. Assertions checked each actual growth flux against 0.5 times
its separately optimized maximum, checked the pFBA objective against an independent
squalene-maximization solve, and checked steady-state residuals below 1e-6.
Observed squalene synthase fluxes (mmol/gDW/h): WT 0.163903667, IUPD 0.348333333,
Ecoli_WT 0.625040214, Ecoli_IUPD 0.348158692, Yarrowia_WT 0.646577207,
Yarrowia_IUPD 0.348333333. Native synthase flux includes downstream consumption;
it must not be described as measured squalene accumulation or extracellular titre.

The additional species retain glucose availability 10 mmol/gDW/h in both growth
and production phases, avoiding an infeasible reduction to 2 after fixing growth.
Yeast9 retains its original study-specific change from 1 to 2 for WT. Thus the
optimization sequence is shared but medium conditions are not identical across
species. Existing main-thread results/analysis_plan were not replaced by validation.
The Python stages were unchanged; full mapping and carbon propagation were not
rerun for these new production inputs, and nothing was published to GitHub.
## Consistent scientific species names (2026-10-14)

All six default model identifiers now use the full scientific species name
followed by WT or IUPD: Saccharomyces_cerevisiae, Escherichia_coli and
Yarrowia_lipolytica. MATLAB config variables use the lowercase genus names
saccharomyces, escherichia and yarrowia. Model/data download paths are unchanged.
The scoped manual-structure table, generated plan, current usage documentation
and six-model stage regression test follow these identifiers. Historical logs
retain their original names. This is a naming-only change; numerical settings,
objectives and propagation rules are unchanged.

## Six explicit COBRA model blocks (2026-10-14)

Replaced the generic species loop and the loadModel, addIU, productionModels and
solvePair model-building helpers with six independent, sequential model blocks.
Each loads its original MAT model and directly calls changeRxnBounds,
changeObjective, deleteModelGenes, addMetabolite, addReaction and optimizeCbModel.
IUPD never inherits an already-fixed WT growth bound. Common Excel/structure
annotation export and SHA256/analysis-plan bookkeeping remain shared.

Ran the original and refactored versions with real MATLAB R2025b, installed
COBRA and GLPK in separate ignored test-output folders. All six reaction tables,
fluxes and metabolite tables matched within 1e-9 relative/1e-10 absolute
floating-point tolerances. Effective structures, prepared reaction SMILES and
mapping statuses also matched for all six networks with online lookup disabled.
Current sequential model calculations match those numerically tested blocks.
A two-reaction glucose example using the same real LocalMapper correspondences
produced identical one-round Carbon_Flow and Summary. The six-model Python stage
regression passed. The full genome-scale mapper/propagation analyses were not
rerun; neither original research outputs nor the live analysis plan were changed
by these isolated tests. Changes remain local and unpublished.
An isolated run containing only the current E. coli IUPD block, without any WT
block, matched its six-model export and generated a one-model analysis plan.

## Short working variables and separate export section (2026-10-14)

Simplified the six explicit model blocks to use modelWT/modelIUPD and direct
COBRA calls with visible numerical bounds. Each saves its model and solution
under species-specific names before the next block starts. Moved medium-filter
loops, model/knockout/solution checks and output bookkeeping to local helpers;
the six export calls now follow the model calculations. Existing validation is
retained, including finite solutions and steady-state residuals below 1e-6.

Executed all six current blocks with MATLAB R2025b, COBRA and GLPK in a separate
ignored output folder. Compared Nonzero and Metabolites sheets against the
pre-refactor baseline: all tables match within 1e-9 relative/1e-10 absolute
floating-point tolerances. Generated input SHA256 hashes validate; model lists,
source order and Python-stage settings match apart from isolated output paths.
No Python algorithm changed. Full genome-scale mapping and carbon propagation
were not repeated in this readability change. Original research results and the
live analysis plan were untouched by validation; no GitHub publication occurred.
The existing six-model analysis-plan regression also passed, using its deterministic test adapter for mapping and two-round propagation; this is distinct from real mapper validation.

## Stage-specific execution settings (2026-10-14)

MATLAB step 1 now contains only input/model settings, ordered carbon sources
and the propagation-round limit. Its version-2 generated plan records models,
input paths/hashes, the shared output directory and rounds. Mapper tools,
batch/atom limits, mapping resources, structure retrieval and mapper/data paths
belong to step 2; propagation resources and the tracing interpreter belong to
step 3. Python execution settings are no longer overwritten by the MATLAB plan.
Older plans are rejected with a request to rerun step 1; no automatic migration
of their execution settings is performed.

Real MATLAB R2025b/COBRA/GLPK execution passed for all six model blocks in an
isolated test directory. Nonzero and Metabolites tables match the preceding
version within 1e-9 relative/1e-10 absolute tolerances. Input hashes, carbon-source
order and effective default Python-stage settings match; the actual generated
version-2 plan loads for both Python stages.

Plan, script-settings and independent-stage regressions passed. These cover
Python parameter precedence, six-model isolation, two-round propagation,
command-line overrides, stale-mapping rejection, invalid-round/legacy-plan
rejection and unchanged input/mapping files. Static Python checks passed.
The small installation example was freshly predicted by real LocalMapper,
RXNMapper and EnzymeMap, selected mappings successfully and completed two
propagation rounds (convergence at round 2). Full six-network mapper and carbon
propagation runs were not repeated. Tests did not replace the live plan or
original research results. All changes remain local and unpublished.

## Rounds in step 3; optional external-input example (2026-10-14)

Removed the remaining cfg.propagation.rounds and cfg.externalInputs blocks,
including their consumers, from the main MATLAB exporter. The version-3 plan
contains only version, models (input paths, hashes and ordered sources) and
output_dir. Step 3's ROUNDS setting now controls normal propagation; changing it
requires no MATLAB export. Earlier version-1/2 plans need regeneration.
The optional examples/step1_use_existing_inputs.py registers existing workbooks
instead of running the MATLAB exporter. It reuses the shared plan validation,
records input hashes and replaces the plan only after validation succeeds.
It neither modifies the Excel inputs nor requires manual JSON editing.

Real MATLAB R2025b/COBRA/GLPK export passed for all six models in an isolated
output folder. Reaction, flux and metabolite tables matched the preceding
version within 1e-9 relative/1e-10 absolute tolerances, and input hashes validated.
The actual reduced plan loaded in both Python stages; changing step 3's ROUNDS
from 10 to 15 changed the resolved limit with the same plan.
Plan and script-setting regressions passed, including six-model two-round
propagation with the deterministic test adapter, command-line round overrides,
legacy-plan/invalid-round rejection, stale-mapping protection, unchanged inputs,
and preservation of the previous external-input plan after invalid sources.
Static Python checks passed. No mapping or carbon-propagation algorithm changed;
full genome-scale mapping/propagation was not repeated. Live outputs/plan and
original research files were not replaced by validation. No GitHub publication.

## WT native-glucose growth and doubled production availability (2026-10-14)

Scientific-condition change requested for the three WT models: retain the source
model's original glucose lower and upper bounds while maximizing growth, fix
growth to 50% of that maximum, double the native glucose uptake limit, maximize
squalene and compute L1 pFBA. Native bounds are Yeast9 r_1714 [-1,1000], iML1515
EX_glc__D_e [-10,1000] and iYali y001714 [-1,-1]. Their production uptake limits
are now 2, 20 and 2 mmol/gDW/h, respectively. The fixed iYali upper bound is
released to 1000 only in the production phase. Other medium settings and all
three IUPD model blocks are unchanged; IUPD glucose limits remain 1, 10 and 10.
The limits are permitted uptake capacities, not forced consumption rates, apart
from the retained fixed native iYali glucose bound during growth optimization.

Executed all six models using MATLAB R2025b, installed COBRA and GLPK in a separate
ignored output folder. Extra assertions verified both initial WT glucose bounds
against the source MAT model, biomass lower/upper bounds and actual biomass flux
at 0.5 times maximum growth, doubled production glucose limits, the pFBA objective
against independent squalene maximization and steady-state residuals below 1e-6.
Observed WT maximum growth / fixed growth / squalene synthase flux:
S. cerevisiae 0.0811002406 / 0.0405501203 / 0.163903667;
E. coli 0.876997214 / 0.438498607 / 1.91055422;
Y. lipolytica 0.105691644 / 0.0528458218 / 0.185189123.
All exported input hashes validated. The three IUPD exports and S. cerevisiae WT
tables match the preceding version within 1e-9 relative/1e-10 absolute tolerances.
E. coli and Y. lipolytica WT predictions intentionally change under the new
conditions. These fluxes are model predictions, not experimental titres.

The live analysis plan and original research results were not replaced by the
isolated test. Rerun normal step 1, then step 2 and step 3 to analyze the new WT
networks. Full mapping and carbon propagation were not repeated in this change;
the Python algorithms were untouched. No GitHub publication occurred.

## IUPD native glucose, proportional isoprenol and paired WT growth (2026-10-14)

Scientific-condition change: each IUPD now retains both native glucose bounds,
sets the isoprenol uptake limit to 2.09 times the native glucose uptake limit,
fixes biomass at 50% of the matching WT maximum, then maximizes squalene and
computes L1 pFBA. IUPD no longer computes or uses its own maximum growth rate.
The WT maximum is saved before WT glucose availability is doubled, using a
species-specific reference variable. Each IUPD block must follow its WT block;
models still load their original MAT files independently. The introduced IU
reactions and species-specific knockouts are unchanged.

IUPD glucose / isoprenol uptake limits (mmol/gDW/h) are S. cerevisiae 1 / 2.09,
E. coli 10 / 20.9 and Y. lipolytica 1 / 2.09. Native iYali IUPD glucose remains
fixed at -1, including its upper bound; it is not relaxed as in WT production.
The 2.09 factor constrains allowed uptake capacity, not the actual uptake-flux
ratio or isotope abundance.

Ran all six models with MATLAB R2025b, COBRA and GLPK in an isolated output folder.
Assertions verified both IUPD glucose bounds against original MAT data, the
2.09-fold isoprenol lower bound, equal paired WT/IUPD biomass bounds and fluxes,
the half-WT growth reference, steady-state residuals below 1e-6 and the pFBA
objective against independent squalene maximization. IUPD growth / squalene
synthase fluxes were S. cerevisiae 0.0405501203 / 0.348333333,
E. coli 0.438498607 / 3.48315874 and Y. lipolytica 0.0528458218 / 0.348333333.
All six input hashes validated; the three WT reaction, flux and metabolite
exports match the preceding version within 1e-9 relative/1e-10 absolute tolerance.
These values are model predictions, not experimental strain validation or titres.

The scientific input settings changed deliberately; Python mapping selection and
carbon-propagation algorithms were not modified. Full mapping/propagation was
not repeated. Validation did not replace the live plan, original source models
or research results. Rerun normal step 1 followed by steps 2 and 3 for the new
IUPD analysis. No GitHub publication occurred.

## Remove repeated source-model settings (2026-10-14)

Checked all three original MAT files: their oxygen lower bounds are already
-1000 and each native objective consists solely of its reference biomass reaction
with coefficient 1. Removed six repeated oxygen lower-bound calls and six initial
biomass-objective resets. WT uses its loaded biomass objective, with compatibility
assertions preventing accidental maximization of another source objective. IUPD
uses the paired WT growth reference and sets its squalene objective directly.
Removed three unused WT generateRules calls; IUPD gene-rule preparation before
knockout remains. The shared carbon-medium helper now skips uptake boundaries
that are already closed instead of reapplying the same bound.

Kept genuine model changes: half-WT growth constraints, WT glucose-limit doubling,
iYali WT's production-only glucose upper-bound release, proportional IUPD
isoprenol uptake, knockouts, added IU/SQS reactions and final squalene objectives.
All six models ran with real MATLAB R2025b, COBRA and GLPK in isolated outputs.
Extra assertions checked inherited oxygen bounds and equal paired growth bounds;
existing solution/steady-state checks passed. Every Nonzero/Metabolites table
matches the preceding version within 1e-9 relative/1e-10 absolute tolerance.
Input hashes and ordered carbon sources also validated.
This cleanup did not change scientific conditions or Python algorithms.
Full mapping and carbon propagation were not repeated. Original source models,
research results and the live analysis plan were not replaced by validation.
All changes remain local; no GitHub publication occurred.

## Yeast9 YAML structures prepared in step 2 (2026-10-14)

Removed the Yeast9-specific MATLAB readYAMLSmiles function, YAML path setting,
SMILES export branch and the extra exportModel argument from step 1. Step 1 now
exports identities, formulas, charges and available identifiers for all models;
it does not read or retrieve SMILES. Yeast9 consequently gains the same basic
formula/charge metadata already exported for the other models.
Step 2 configures YEAST9_YAML_FILES by analysis name. The existing
carbon_flow_tools/smiles_completion.py module reads the first SMILES entry per
metabolite ID from the pinned Yeast9 9.1.1 YAML layout, with no new dependency.
The reader is not a generic model-YAML importer. Unconfigured analyses never
receive Yeast9 annotations. Manual overrides and existing input structures are
preserved, local annotation loading works offline and PubChem remains the
fallback for missing structures. StructureSource/StructureReference record YAML
origin in the enriched step-2 workbook; raw step-1 inputs remain unchanged.

Real MATLAB R2025b/COBRA/GLPK exported all six inputs without the YAML reader.
All reaction tables match the preceding version within 1e-9 relative/1e-10
absolute tolerances. The new Python reader exactly matches every original
SMILES_Yeast9 entry in both prior S. cerevisiae exports. With online lookup
disabled, all six effective SMILES lists, prepared reaction SMILES and statuses
match the preceding version, including scoped manual overrides. Other species
do not receive a SMILES_Yeast9 column. Input hashes validated.
A real LocalMapper run on a two-reaction glucose input exported without SMILES
loaded the structures in step 2 and produced unchanged one-round Carbon_Flow
and Summary compared with the preceding fixed-mapping small example. Raw input
and mapping file hashes were preserved during propagation.
Structure-completion, script-settings and six-model plan regressions passed,
covering first-entry parsing, unsupported YAML layout rejection, model scoping,
offline annotations, manual/input priority, PubChem mocks, caching, stage isolation
and immutable inputs. Static Python checks passed. Full genome-scale three-tool
mapping and live PubChem completion were not repeated; extra formula/charge
metadata may refine fallback candidate filtering where Yeast9 lacks structures.
Scientific model conditions, selection rules and carbon-propagation algorithms
are unchanged. Original research files and the live plan were not replaced by
validation. Changes remain local and unpublished.

## Merge test helpers into scripts (2026-10-14)

Moved the nine test programs and their package initializer from tests/ into
scripts/, alongside environment setup and data-download helpers. Updated all
active cross-test imports, spawned adapter module names and documented commands
to scripts.test_*. The installation example now uses scripts/results/example/.
Existing ignored validation outputs moved to scripts/results/, with portable
paths in saved analysis_plan.json files updated. SHA256 checks before and after
the move and after regression testing confirmed all 68 Excel files unchanged.
The empty tests/ directory was removed. Historical entries above retain their
original test module names; current commands are in the user guide.

All nine relocated regression programs passed: carbon tracing, mapper backends,
consensus selection, parallel propagation, independent stages, SMILES completion,
script settings, analysis-plan handoff and native EnzymeMap rule fallback.
The first eight ran in the existing LocalMapper environment; the last ran in
the existing EnzymeMap environment. Python syntax/unused-code checks also passed.
No scientific algorithm or model settings changed. Full genome-scale model or
mapper runs were not repeated. No GitHub publication occurred. Original edited
files, test logs and the workbook hash manifest are kept outside the project.

## Profile-guided propagation optimization (2026-10-14)

Measured one real S. cerevisiae WT network: 493 input reactions, 335 selected
mapping rows, 522 metabolite rows, eight propagation rounds and 16,073 rows in
round 8. Existing selected mappings were reused; this benchmark excludes mapper
prediction and structure lookup. Results were written to a private validation
folder, leaving current analysis outputs and the input workbooks untouched.

Two sequential fresh-interpreter wall-clock runs per configuration, without
cProfile, gave the following medians on the same Ubuntu WSL x86_64 machine:

| Phase | Before (s) | After (s) |
|---|---:|---:|
| Total, including workbook input and validation | 11.87 | 4.47 |
| Propagation computation | 5.66 | 0.84 |
| Eight workbook saves | 4.67 | 2.09 |

The measured total speedup was 2.66x. Exact base SMILES and positional labels
now key two bounded LRU caches (8,192 entries each per process), reusing only
immutable conversion results. These are not mapper prediction caches. The
35,164 conversion requests produced about 92% hits in this sample. lxml 4.9.4
enables openpyxl's native XML backend without replacing the workbook writer.
Scientific rules, state enumeration, ID assignment, output order, per-round
saves, atomic file replacement and existing model-level concurrency are unchanged.

All eight deterministic regression programs passed with the final backend.
Checks include iterable/generator assignments, position-sensitive labels,
symmetric canonicalization, repeated origin labels and cache capacity bounds.
Every output cell, type, style and dimension in all eight Carbon_Flow sheets
and Summary matched baseline. Original research input and private mapping
snapshot SHA256 values were unchanged. The real LocalMapper environment smoke
test also passed with the native backend installed.

This is one eight-round workload, not a guarantee for twenty rounds, all six
models or other hardware. Combinatorial state growth and cumulative Excel
serialization remain; bounded helper reuse does not eliminate those costs.
No within-round reaction process pool or new user setting was added. Profiling
scripts, research data and timing artifacts remain outside the project. Changes
remain local and unpublished.


## iYali5 model and paired structure annotations (2026-10-14)

The local test project now uses ordinary iYali5 for Yarrowia pFBA and the paired
eciYali5 YAML solely for SMILES. Upstream commit
3b97d38d92327a099f3b0f6f9cd92e8aeab78c10, CC BY 4.0; file hashes and sources are
recorded in data/yarrowia/iYali5_source.json. Earlier iYali files/results remain
untouched. This replacement is not a reconstruction of the published study's
original Yarrowia results.

The new MAT model has 1,986 reactions and 1,687 metabolites. Its glucose exchange
1714 is closed in the release; step 1 explicitly selects a glucose medium with
limit 1 (WT production 2), closes other carbon uptake, and retains the paired-WT
half-growth / IUPD isoprenol-limit strategy. Native oxygen bounds are retained.
SQS is 1012; ATP/ADP/IPP are m141/m143/m491; ERG13 YALI0F30481g constrains 559/560.

The existing reader was extended to the supplied COBRA/RAVEN YAML layouts, with
model-name scoping and ID/name/compartment/available-formula checks. Manual and
existing input structures precede native annotations, then online retrieval.
The existing manual table supplies explicit [H+] for iYali5's H+_p+1 entries;
no generic polymer/lipid structure was invented. No new dependency or change to
carbon-transfer rules was introduced. A separate import compatibility fix reads
numeric ReactionID values as strings in prediction reuse and Final_Mapping.

MATLAB R2025b/COBRA/GLPK successfully optimized/exported all six model blocks to
isolated results/iyali5_validation outputs. Yarrowia WT/IUPD exported 390/386
nonzero reactions and 407/404 metabolites. Both fixes use growth 0.0640660659/h;
both SQS optima were 0.00138260977 mmol/gDW/h under these model constraints.
This result was recorded rather than changing bounds to force an IUPD increase.

With online lookup disabled, WT used 352 paired YAML structures plus 8 proton
overrides; IUPD used 346 YAML structures plus 11 overrides (including introduced
IU metabolites). Each retained 47 unresolved structures. Reaction preparation
produced 223/219 mapper-ready chemical reactions and 99/98 identity mappings;
these counts are eligibility, not verified mapping accuracy. Yeast9 WT/IUPD
selected SMILES matched the preceding reader on all 522/516 exported metabolites.

A real three-reaction glucose-pathway smoke test called LocalMapper, RXNMapper
and EnzymeMap and produced 3/3 complete carbon maps (including identity
transport); one result retained ReviewRequired. EnzymeMap used rules after no
database matches. One-round tracing generated glucose[c], glucose 6-phosphate[c]
and fructose 6-phosphate[c]. Numeric-ID reuse was checked for all three tools.
Structure, settings, analysis-plan, stage/numeric-ID, carbon-transfer and consensus
regressions passed. Current input/result files were not replaced; rerun step 1
before using the replacement for a full analysis. Full-network mapping/tracing,
biochemical accuracy of every native SMILES, and GitHub publication were not done.

## Optional AI-reviewed workbook propagation (2026-10-14)

Added step4_carbon_flow_ai.py as a small entry reusing step 3's settings and
the existing propagation/validation implementation. It reads only available
<model>_2_atom_mapping_AIGenerated.xlsx copies from the current analysis plan
and writes <model>_3_carbon_flow_AIGenerated.xlsx. Normal step-3 outputs remain
separate. Source-input hashes, participant structures, active direction and
carbon correspondences undergo the same checks as ordinary mapping workbooks.
The entry does not invoke AI, mappers or mapping selection. The external skill
and reviewed workbook updates are a separate upstream operation.

All eight deterministic regression programs passed. Extended stage tests use
synthetic reviewed mappings to prove changed correspondences change carbon-label
positions, one/two reviewed models work, unavailable copies do not fall back to
ordinary files, stale hashes are rejected, and interpreter relaunch retains
the step-4 entry. Input/mapping/ordinary-output hashes remain unchanged.
The actual `python3 step4_carbon_flow_ai.py --example --rounds 1` command passed
through its interpreter relaunch using a temporary example copy; this checks
routing, not scientific AI prediction. The fixture copy and generated result
were removed from the project after checks. Syntax/unused-code checks passed.
No real AI mapping review or full six-model propagation was performed for this
entry. A skill-creator prompt was prepared separately; no skill was installed
and no changes were published to GitHub.

## Explicit Calibri and Yarrowia source-count audit (2026-10-14)

All 19 current result/input workbook font tables named Calibri, but seven
carbon-flow workbooks (including a user-created copy) inherited the minor theme
font. Their themes contain a SimSun/Hans definition. The shared save helper now
sets explicit Calibri and removes font scheme/charset inheritance, rebuilding
the font index after normalization. Repeated-save tests preserve values, wrapping,
bold/italic attributes and valid font indexing. Plain output styling is retained.
Current seven flow files received only a styles.xml font patch. SHA256 checks
confirmed every other ZIP entry unchanged. The twelve input/mapping files were
already explicit Calibri and remain byte-identical. All six input hashes still
match the analysis plan. Representative Summary/flow crops were rendered and
inspected; native Microsoft Excel's UI was not inspected.

For Yarrowia_lipolytica_IUPD, independent counts from all ten saved flow sheets
exactly matched Summary: isoprenol reaches 86 metabolite identities in round 1
and 129/315 (40.95%) in round 10. Every saved CarbonLabels multiset agreed with
LabeledSMILES map labels. Across 180603 flow records, ProductionID identity, substrate
state names, stoichiometric source-ID counts and source-label inventories passed.
This found no counting mismatch or invented source labels.

The saved pathway carries isoprenol carbon through sterols to formate (reaction
317), CO2 (445), bicarbonate (1667), oxaloacetate (958) and carbamoyl phosphate
(250). Thus the current qualitative reachability calculation includes metabolites
outside terpene synthesis. Important mappings in this chain (462, 698, 317) still
have ReviewRequired=TRUE. Numerical consistency does not establish biochemical
correctness of those atom correspondences; no mapping, source rule, denominator
or count was silently corrected or excluded. Evidence-based mapping review and
subsequent propagation remain necessary to resolve that scientific uncertainty.

## Inorganic-carbon label endpoints (2026-10-15)

Step 3 now defaults to `ALLOW_INORGANIC_CARBON_RECYCLING = False`; step 4
inherits it. The common propagation function identifies CO2, bicarbonate,
carbonic acid and carbonate by structure in all compartments. Labelled products
remain in the output/catalog and Summary, but are not added to available
labelled substrate states. Carboxylation reactions continue transferring labels
from their other substrates. Terminal substrates therefore appear as SourceID 0
downstream. Mapped-network Summary denominators and all mapper outputs are
unchanged. Explicit inorganic-carbon sources require the setting to be True.
This is a qualitative analysis-scope assumption, not a claim that all CO2
escapes or that biological refixation is negligible.

All eight regression programs passed, including tests for both source orders,
CO2/bicarbonate termination, compartment transport, subsequent rounds, retention
of pyruvate labels in oxaloacetate, settings shared with step 4, an inorganic
source, and restoring recycling. Syntax/unused-code checks passed.

The actual default step-3 entry completed ten rounds for all six current models
using existing mappings. Independent record counts matched every Summary row;
no terminal ProductionID was consumed downstream. Input/mapping workbook and
analysis-plan SHA256 values were unchanged. Previous carbon-flow workbooks and
the code diff were retained in a private audit directory outside the repository.
Explicit Calibri fonts were checked, and representative Summary/flow crops were
rendered and visually inspected; native Excel UI behavior was not tested.

Round-10 labelled metabolite counts (recycling enabled before / disabled after):

| Model | Mapped carbon denominator | Glucose before / after | Isoprenol before / after |
| --- | ---: | ---: | ---: |
| Saccharomyces_cerevisiae_WT | 334 | 239 / 233 | -- |
| Saccharomyces_cerevisiae_IUPD | 336 | 231 / 226 | 15 / 15 |
| Escherichia_coli_WT | 292 | 226 / 221 | -- |
| Escherichia_coli_IUPD | 292 | 225 / 220 | 14 / 14 |
| Yarrowia_lipolytica_WT | 318 | 189 / 184 | -- |
| Yarrowia_lipolytica_IUPD | 315 | 167 / 162 | 129 / 20 |

Yarrowia IUPD C11 destinations independently decreased from 129 to 20 at round
10; isoprenol-labelled coverage changed from 40.95% to 6.35%. A separate ten-round
run with recycling re-enabled reproduced every sheet name, cell value and table
dimension in all six pre-change carbon-flow workbooks exactly. No new mapper
predictions, MATLAB optimization, biochemical correction or GitHub publication
was performed. These checks validate the rule and compatibility, not the
biochemical correctness of the selected mappings or quantitative isotope flux.

All eight deterministic regression programs and syntax/unused-code checks passed.
No new mapper predictions or AI corrections were run. Original results and audit
evidence are backed up privately outside the project. Changes remain local and
were not published to GitHub.
