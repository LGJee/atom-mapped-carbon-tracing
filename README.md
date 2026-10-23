# Atom-Mapped Carbon Tracing

Prepare pFBA-active metabolic networks with **COBRApy**, compare three atom-mapping
tools and trace possible carbon-source label propagation. In
`step1_prepare_carbon_tracing_inputs.py`, each of the six named strain blocks follows
the same sequence: load model, list carbon sources, set medium and modifications,
fix growth, run production pFBA, and register output. Shared solver/output settings
remain at the top. Steps 2 and 3 read the generated analysis plan. All computation
runs in Python on Linux/WSL.

```text
step1_prepare_carbon_tracing_inputs.py
    -> <model>_1_pFBA_nonzero_reactions.xlsx + data/analysis_plan.json
step2_atom_mapping.py
    -> <model>_2_atom_mapping.xlsx
step3_carbon_flow.py
    -> <model>_3_carbon_flow.xlsx
```

## 1. Install Python environments

In Linux/WSL x86_64, install [Miniforge](https://github.com/conda-forge/miniforge)
and Git, then run:

```bash
git clone https://github.com/LGJee/atom-mapped-carbon-tracing.git
cd atom-mapped-carbon-tracing
bash scripts/setup_environments.sh
python3 scripts/download_enzymemap.py
```

Setup creates four isolated environments: `preparation` for COBRApy/GLPK,
`localmapper` for LocalMapper and tracing, `rxnmapper`, and `enzymemap`.
Compatible existing environments are checked and reused; none are deleted.
The preparation check solves a real small pFBA problem; mapper checks call their
actual interfaces. No commercial solver or additional numerical application is
required. See the [installation details](docs/user-guide.md#1-install-python-environments-and-enzymemap-data).

## 2. Obtain the reference models

Download these external files into `data/`, preserving their upstream licenses:

| Species | Source | Local files |
|---|---|---|
| S. cerevisiae | [Yeast9 v9.1.1](https://github.com/SysBioChalmers/yeast-GEM/tree/v9.1.1/model) | `data/Saccharomyces_cerevisiae/yeast-GEM.mat`, `yeast-GEM.yml` |
| E. coli | [BiGG iML1515](https://bigg.ucsd.edu/models/iML1515) | `data/Escherichia_coli/iML1515.mat`, `iML1515.json` |

Keep each source model in its own subfolder under `data/`, with its companion
JSON/YAML files in the same folder. Configure the computational model path only
once in step 1. The generated plan records its source path and SHA256; step 2
automatically reads supported native structures and scans that folder for matching
annotations. Companion filenames need not match the computational model.

MAT is an upstream model **data format**, read directly by COBRApy and SciPy.
Model files, large datasets and environments are not bundled.

## 3. Configure and run step 1

Edit the user settings and relevant named model blocks in
[step1_prepare_carbon_tracing_inputs.py](step1_prepare_carbon_tracing_inputs.py):

```bash
python3 step1_prepare_carbon_tracing_inputs.py
```

Step 1 reports preparation by condition, the actual condition count, numbered
exports, each target reaction and its flux, and the generated workbook path.
Target fluxes use the source model's units. The pFBA total absolute flux is a
separate diagnostic, available with `--verbose`; it is not the target-product flux.
Only known legacy MAT field/compartment notices are hidden by default; other
warnings remain visible:

```bash
python3 step1_prepare_carbon_tracing_inputs.py --verbose
```

Step 1 uses the configured `.envs/preparation/bin/python` when it is installed,
even if the launching interpreter also has COBRApy. You can invoke it explicitly.
Paths are resolved from the project directory, independently of the working folder.

The four numbered blocks use standard COBRApy operations directly:
`gene.knock_out()`, `Metabolite`, `Reaction.add_metabolites()`, explicit
`lower_bound` / `upper_bound` assignments, `model.optimize()` and `pfba()`.
Each strain's modifications stay together, even where operations repeat between
strains. Model versions are not enforced; source paths and hashes are recorded.
The non-IUPD growth objective
is set explicitly. Maximum-growth solutions are checked before use; production
solutions are validated by the shared exporter before writing any result.
See the [COBRApy model-building tutorial](https://cobrapy.readthedocs.io/en/latest/building_model.html)
for these operations.

The default blocks prepare non-IUPD/IUPD pairs for Saccharomyces_cerevisiae
and Escherichia_coli. Each IUPD block independently reloads
its original model and uses 50% of its paired non-IUPD maximum growth rate. Both
comparators maximize squalene, then minimize the sum of absolute reaction fluxes
using `cobra.flux_analysis.pfba`. Its `objective_value` is the minimized flux sum;
the biological production value is read from the squalene reaction flux.

Here, **non-IUPD** denotes the comparison model without IU-specific engineering;
it does not imply an unmodified wild type. Shared production modifications can
still be present, such as introduced squalene synthase in both E. coli models.
Code identifiers and filenames use `non_IUPD`; descriptions use `non-IUPD`.

| Species | Growth-phase glucose limit | non-IUPD production limit | IUPD isoprenol limit |
|---|---:|---:|---:|
| S. cerevisiae | 1 | 2 | 2.09 |
| E. coli | 10 | 20 | 20.9 |

Limits are mmol/gDW/h. Both species inherit native glucose bounds.
IUPD retains the growth-phase glucose bounds. Native oxygen settings are
inherited. Other carbon uptake is closed in the E. coli comparison.
E. coli non-IUPD/IUPD both receive functional SQS activity and
a squalene output. These are model-specific conditions, not universal media
or experimentally validated strains. See [model notes](docs/yeast9-example.md).

Each workbook contains `Nonzero` and `Metabolites`. Step 1 exports reaction IDs,
names, signed fluxes and equations, plus metabolite IDs, compartment-qualified
names, formulas, charges and available identifiers. It does not retrieve SMILES.
Step 2 handles structures from the source model directory, manual overrides and PubChem.
No analysis-name-to-YAML table needs editing in step 2.

Step 1 writes `data/analysis_plan.json` only after every export succeeds. The
plan records exact input hashes, model provenance, relative paths, ordered carbon
sources and output folder. Do not edit it manually. Rerunning step 1 replaces its input exports;
rerun step 2 before tracing those inputs. Step 2 controls mapping execution;
step 3 controls propagation rounds and resources. Their defaults need no edits.

pFBA selects one optimal representative, not all possible pathways. Equivalent
optima can differ between solver versions. Record versions and compare nonzero
reaction sets as well as objective values when validating a new environment.

## 4. Predict and select atom mappings

```bash
python3 step2_atom_mapping.py
```

The script reads step 1's plan, runs LocalMapper, RXNMapper and EnzymeMap with
CPU/memory-limited concurrency, then selects final correspondences. Each model
has one `<model>_2_atom_mapping.xlsx` containing the three original predictions,
`Mapping_Selection`, `Final_Mapping` and the enriched metabolite table.

Step 2 applies species/name overrides from `data/manual_smiles.csv`, preserves
existing input structures, then uses compatible source/companion annotations before
online structure completion. It checks metabolite identity rather than choosing
a file by species name. Conflicting structures and name-selected database results
mark affected reactions for review. Structure sources and candidates are recorded
in the enriched metabolite table. Inputs remain unchanged. Inspect `ReviewRequired=TRUE` and, if necessary,
correct `Final_Mapping` before step 3. Rerunning step 2 can replace those edits.
`AIGenerated` is stored only in `Final_Mapping` and defaults to FALSE.
The three-tool mapping stage does not perform AI review.

Formula/SMILES carbon-count conflicts require structure or repeat-unit review
before accepting a whole-reaction map. Independent consistent cofactor pairs may
still be recovered. Selection then uses complete carbon agreement; ties/conflicts use
EnzymeMap > RXNMapper > LocalMapper (scores are not comparable). When a whole
prediction fails validation, step 2 tries an exactly structure-aligned carbon-only
view of that prediction. It keeps only unique carbon-to-carbon correspondences;
partial views remain usable and require review. Raw predictions are unchanged.
Only if no mapper supplies usable carbon correspondence does
`carbon_flow_tools/carbon_mapping_rules.py` try conservative skeleton rules:
unique identity retention, A/G/C/U nucleoside mono/di/triphosphates, NAD(H) or
NADP(H) within each family, and uniquely positioned CoA-core carbons (not acyl carbon).
Full structure connectivity, stereochemical compatibility, equal transfer amounts
and unique donor/product pairing are required. Unsupported carbon remains unresolved.
See [scientific basis and limitations](docs/carbon-mapping-rules.md).

Every reaction appears in Final_Mapping. MappingCoverage is Complete, Partial,
Unresolved or Boundary. Legacy NoCarbon records remain readable. Partial mapper views use sparse MappedReactionSMILES;
rule-only partials show a reduced unit subreaction in the mapped-SMILES columns, with actual transfers in Partial_Mapping (zero-based atom
positions in the exact saved structures, original molar coefficient and rule
version). They are not a fabricated complete biomass map. All rule fallbacks and
partial mapper views have ReviewRequired=TRUE and AIGenerated=FALSE. SelectionReason
states the selected tool or rule and its limitation. Step 3 uses partial transfers
for label reachability without rounding fractional coefficients; it never treats
the unassigned carbon positions as known. Summary counts carbon-containing
metabolite IDs with usable atom correspondences, now including partial donors and
products, rather than all participants in an unresolved pseudo-reaction.

AI review changes only Final_Mapping in a separate copy. Confirmed partial rules
stay Partial; rejected or scientifically unresolved rules become Unresolved and
are not propagated even though their original audit records remain. A validated
complete AI replacement becomes Complete and supersedes partial records.


## 5. Analyze carbon origins

```bash
python3 step3_carbon_flow.py
```

Step 3 runs two independent methods from the same selected mapping:
Enumeration gives possible carbon-source label states; the restored size-1
EMU-inspired solver gives fixed-flux mean carbon-origin fractions and retains
Unknown. It is not a full isotopomer simulator or an experimental measurement.
Set `ANALYSIS_METHODS` to `["Enumeration"]`, `["EMU"]`, or both at the top of
`step3_carbon_flow.py`. Method-specific settings are grouped separately.

Enumeration reads both normal and optional AI-reviewed step-2 mappings for models
in the selected plan, using the same ordered carbon sources for both versions.
It never calls mappers or changes input/mapping workbooks. Each
`<model>_3_carbon_flow.xlsx` combines `Summary`, a `Carbon_Flow` sheet with
`Carbon_Flow_Round` as its first column,
`EMU_Metabolites` and `EMU_BoundaryInputs`. AI copies use `_AIGenerated.xlsx` and
both methods read that copy's Final_Mapping. No separate `_EMU` file is created.
The default maximum is 10 Enumeration rounds, stopping earlier on convergence.
`EXPORT_DETAIL_ROUNDS = 3` exports only R1-R3 details; later completed rounds
still update Summary. Set this to the requested maximum to export all rounds,
or 0 to keep Summary only. EMU does not use Enumeration rounds or ProductionID.
See [EMU assumptions and limitations](docs/source-contributions.md).

`MAX_EXCEL_SIZE_MB = 50` applies separately to each output. Generated rows are staged in temporary disk storage rather than retained
as a whole-round DataFrame or an in-memory Excel workbook. `OUTPUT_BATCH_ROWS =
50000` is an upper batch bound; byte-based checks can flush earlier. Size checks
also occur at layer and round boundaries.
At a checkpoint, the program streams a real candidate XLSX and measures its
compressed size. If it reaches 50,000,000 bytes, enumeration stops and the
unfinished round is discarded. Only complete rounds below the limit atomically
replace the output; the preceding complete result is retained unchanged.

There is no worksheet pagination. An exported round exceeding 1,048,575 data rows cannot fit
in one sheet and is discarded with `row_limit`. Size stops use `file_size_limit`,
never convergence. If round one cannot be saved, the report records zero completed
Enumeration rounds. Combined output is published only with a complete Enumeration round; if none
can be saved, prior output is retained. Select EMU alone to export it independently.

Checks occur at batch boundaries, not at every byte: temporary candidates may
exceed the limit before rejection. Rebuilding sorted candidate files costs disk
space and I/O. Label states and their indexes are still needed in memory; this is
not full external-memory propagation, and independent memory safety checks remain.

Each model/version keeps its own Summary sheet. Step 3 does not generate a
separate summary-collection workbook. The run manifest is `carbon_flow_run.json`.
Inspect its current status and saved rounds rather than trusting an older output file.

Propagation checks available Linux/WSL memory during enumeration, state merging
and streaming Excel writing. MIN_AVAILABLE_MEMORY_GB = 2.0 at the top of step 3 reserves
at least 2 GiB or 10% of detected memory, whichever is larger. Tighter cgroup-v2
limits are respected. Normal and AI tasks use the same protection.
If the reserve is crossed, the task exits with an error; completed rounds already
saved remain available, and the incomplete round does not replace them. This is
a safety stop, not convergence or a resumable checkpoint. Rerunning starts from
round 1. Free memory, select fewer models, or set MAX_CARBON_PROCESSES = 1.
No states or combinations are silently discarded.

By default, labels reaching CO2 or bicarbonate are recorded as endpoints but
are not reused in downstream reactions. Carbonic acid and carbonate are treated
the same way in every compartment. Carboxylation reactions still transfer labels
from their other substrates. This excludes inorganic-carbon recycling from the
reachability analysis; it does not assume that all CO2 escapes or that refixation
is biologically negligible. To include recycling or trace an inorganic-carbon
source, set `ALLOW_INORGANIC_CARBON_RECYCLING = True` at the top of
`step3_carbon_flow.py`. Both mapping versions use this setting. Summary still counts labelled
inorganic-carbon endpoints, and its mapped-network denominator is unchanged.

For an optional subset or one-run override:

```bash
python3 step2_atom_mapping.py Escherichia_coli_non_IUPD Escherichia_coli_IUPD
python3 step3_carbon_flow.py Escherichia_coli_non_IUPD Escherichia_coli_IUPD --rounds 1
```



Later rounds do not build or stage complete reaction-provenance rows. They retain
one minimum-step entry per unique product state for that round's Summary.
The global state catalogue, positional assignments and earliest steps are stored
in temporary SQLite tables, with bounded caches and ordered pages of 256 states.
No whole-state copy is needed for convergence checks. Layer/round snapshot
boundaries, ProductionIDs and exported ordering are unchanged. The maximum round
count remains 10; reducing output does not eliminate combinatorial computation.

Output buffers flush at the configured row bound or an estimated 8 MiB, whichever
comes first. Candidate chunks also account for molecule size, not just combination
count. Working-disk checks reserve at least 1 GiB or 1% of the filesystem capacity.
Low disk space is a resource failure, not convergence. Last complete workbooks
remain available; temporary working databases are not resumable checkpoints.
Normal cleanup removes only the current run's temporary directory. Files left by
an older interrupted run are not deleted automatically.

Memory checks also cover model/input loading, mapper batches and rule workers,
state-page decoding, Summary scans, EMU arrays/solver preparation and Excel export.
Workbook-load reservations use uncompressed archive sizes. These checks and
allocation estimates are safeguards, not hard OS limits: an external process,
a large native-library allocation or host pressure hidden from WSL can still
exhaust resources. In particular, sparse solver fill-in is estimated conservatively.

### Console output and estimated time

Steps 2 and 3 show short progress messages and result locations by default.
For example: `Example / RXNMapper: 50% (174/348 reactions).`
Mapper dependency warnings, library messages and per-reaction errors are saved in
`logs/<model>_<mapper>.log` under the analysis result directory (latest run per tool).
The shared EnzymeMap run writes `logs/EnzymeMap.log`.
A short notice identifies reactions that need checking; final workbook summaries
still report the review count. These messages do not change calculation or selection.
Use `--verbose` to show library output and technical progress directly in the terminal. Routine batch, layer and resource details are available with
`--verbose`; expected input/configuration errors otherwise show their context
and an action instead of a traceback:

```bash
python3 step2_atom_mapping.py --verbose
python3 step3_carbon_flow.py --verbose
```

Progress refreshes at computation checkpoints about every 15 seconds. Approximate
remaining time uses recent completed work, including samples between printed updates,
and resets when the phase or reaction changes. It waits for initial samples and
withdraws the estimate after a prolonged gap without completions. Variable reaction
costs can still move the estimate up or down; times are rounded rather than exact.
Step 2 counts unique mapper requests after duplicate reuse. EnzymeMap database
matching and rule prediction have separate counters; only rule prediction is timed.
Step 3 counts completed, ordered enumeration chunks, never prefetched combinations.
Its estimates cover the current reaction or Excel row preparation, excluding final
workbook compression. They do not predict growth in later rounds or full-run time.
Structure queries, model/database loading, solver work and scheduling waits have no
reliable ETA. Default output shows counts or elapsed time in those phases; verbose
output may show `Estimating...` or `Unavailable`. A blocking library call may delay
updates. Parallel tasks identify their model/tool or mapping version on each line.
These are changing execution estimates, not convergence or biological rates.

### Optional: propagate using AI-reviewed mapping copies

After external review, place `<model>_2_atom_mapping_AIGenerated.xlsx` copies
in the configured results directory, keeping the original model names and input
hashes. Successful replacements belong in `Final_Mapping.MappedReactionSMILES`,
with `AIGenerated=TRUE`; leave `Mapping_Selection` unchanged as the original tool-selection record. Preserve raw tool predictions and all unchanged reactions.

```bash
python3 step3_carbon_flow.py
```

Step 3 automatically processes both `<model>_2_atom_mapping.xlsx` and any
`<model>_2_atom_mapping_AIGenerated.xlsx` in the selected plan's output folder.
Each version produces independent Enumeration outputs, retaining the
`_AIGenerated` suffix for reviewed copies. Either mapping version can be present
alone. Excel lock
files, temporary files and unrelated names are ignored. Only models in the plan
are included; model selection applies to both available versions of that model.
Both versions must pass input-hash, direction and structure checks. The suffix
selects a file; it does not certify correctness. Step 3 performs no AI review or
mapper prediction. The separate step-4 entry has been removed.


### Automatic task scheduling

Steps 2 and 3 share `carbon_flow_tools/task_scheduler.py`; users do not run it
separately. Existing `MAX_MAPPER_PROCESSES`, `MAPPER_MEMORY_GB`,
`MAX_CARBON_PROCESSES` and `CARBON_MEMORY_GB` settings remain at the top of the
entry scripts. No external scheduler or extra dependency is required.

- Independent model/mapper jobs and mapping-workbook analyses use a bounded queue.
  Larger estimated jobs start first; returned results keep their original order.
  Original precedes AI for a model when only one slot is available.
- CPU affinity, cgroup CPU quotas and available memory bound the capacity. Free
  memory is checked again before every job starts. Active per-task estimates are
  reserved conservatively, so memory pressure can reduce concurrency further.
  These estimates are not measured peak guarantees; runtime MemoryGuard stays active.
- LocalMapper and RXNMapper jobs run concurrently as before. EnzymeMap then uses
  a separate phase, sharing one database scan across pending models with the same
  adapter, interpreter and options. Its independent rule searches use a process pool
  sized to available CPUs, memory and reaction count; `MAX_MAPPER_PROCESSES` can cap
  this pool too. It is not fixed to 32 threads. Each rule worker loads the pinned
  rules once, uses one native-library thread, and does not scan the database.
- Exact duplicate reaction SMILES are predicted once per mapper call and restored
  to every original row. EnzymeMap shares this reuse across the merged models.
  Direction, component multiplicity, charges and stereochemistry remain part of the
  exact input. Different configurations never share predictions. Existing workbook
  reuse still applies on later runs; no new permanent prediction cache is created.
- EnzymeMap reserves 1 GiB for its parent and estimates 0.5 GiB per rule worker,
  with the shared scheduler checking free memory before each reaction starts.
  These are estimates, not hard memory caps; difficult reactions can require more.
  One reaction's internal official rule order is unchanged. Rule workers and spawned
  analysis processes inherit thread limits before numerical-library imports.
- A task failure stops new admissions; already running jobs finish normally before
  the original error is raised. An idle queue with insufficient memory fails clearly
  rather than waiting forever. Unknown memory readings permit only serial admission.
- Rounds and reaction traversal remain ordered. Large Enumeration reactions can
  compute independent candidate chunks in parallel; the parent consumes chunks in
  original order, assigns ProductionIDs, merges states and writes SQLite/Excel.
  Workers never update the shared catalogue or start the next round. Enumeration
  and EMU remain separate sequential calculations within the same model job.
- Step 3's `MAX_ENUMERATION_PROCESSES = None` uses the CPU/memory share assigned
  to that model, not every CPU for every model. Set it to `1` for serial candidates.
  Small reactions stay serial to avoid startup overhead. Molecule templates and
  fixed atom transfers are reused without changing canonical labelled structures.

Scheduling improves use of independent work, not the combinatorial size of one
reaction. Small jobs may be dominated by process startup and Excel writing; no
fixed speedup or full CPU utilization is guaranteed. Use `--verbose` for queue and
resource details. Ordinary progress remains concise.


## 6. Read the results

Every model has three numbered workbooks: input reactions, selected atom mapping,
and carbon flow. `ProductionID` identifies a compartment-qualified metabolite
and its labeled structure, rather than a count of molecules. `CarbonLabels`
lists source labels; `LabeledSMILES` records their positions. Summary percentages
use carbon-containing metabolites participating in usable final mapped
reactions as the denominator, not the whole genome-scale model.

Results describe **possible carbon-source label propagation and reachability**,
not experimental isotope abundance or quantitative carbon allocation.
Mapping agreement and confidence scores do not establish biochemical correctness.
Compare non-IUPD/IUPD within each species under explicitly stated conditions; differences
in model coverage and medium prevent interpreting cross-species percentages as
direct biological rankings.

## Custom inputs and installation check

The default commands reproduce the four-model study. A different analysis can
keep its own preparation script, generated plan and results under
`results/results_by_AI/<analysis-name>/`, sharing this project's environments and model data.
The case folder directly contains its script, plan, reports and Excel workbooks;
there is no nested results folder. The name identifies an AI-organized analysis,
not AI-generated mappings; `_AIGenerated` still identifies reviewed mapping copies.
This leaves the study script and its outputs unchanged. A future orchestration
skill can manage these locations; users do not need to write JSON.

The [glucose/xylose example](examples/step1_prepare_custom_model_inputs.py)
contains two explicit COBRApy model blocks for single-source and simultaneous
mixed-source growth pFBA. See the [custom analysis guide](docs/user-guide.md#3-configure-a-custom-network-and-carbon-sources)
for commands and the [calling contract](docs/interfaces.md) for automation.

Steps 1-3 run fixed code and do not call an LLM or a skill. The intended future
AI-orchestrated workflow includes mapping review after step 2 by default:
`carbon-tracing-review-skill` delegates flagged reactions to
`enzymatic-atom-mapping`, saves reviewed copies, then step 3 propagates both normal
and reviewed copies for comparison. This repository provides the
handoff files; automatic skill orchestration is not implemented here.

For a different model, adapt a named COBRApy block and its source list in step 1.
For Excel inputs prepared elsewhere, use
[examples/step1_use_existing_inputs.py](examples/step1_use_existing_inputs.py)
instead of the model exporter. Set propagation rounds with `ROUNDS` at the top
of step3_carbon_flow.py; the default is 10. See the
[detailed guide](docs/user-guide.md#3-configure-a-custom-network-and-carbon-sources)
for the input schema and mixed-source configuration.

The small installation example does not require source models or an analysis plan:

```bash
python3 step2_atom_mapping.py --example
python3 step3_carbon_flow.py --example
```

It writes to `scripts/results/example/`. This is an installation check, not the
reference study. The implementation lives in `carbon_flow_tools/`; normal users
edit step 1 for model/input preparation and step 3 to change rounds. Third-party sources and licenses are recorded in
[THIRD_PARTY.md](THIRD_PARTY.md); validation is in [docs/validation.md](docs/validation.md).

The `scripts/` folder contains environment setup, data-download utilities and
optional `test_*.py` regression checks. Run tests as modules from the project root,
for example `.envs/localmapper/bin/python -m scripts.test_carbon_tracing`.

The original equation and all raw predictions remain in step 2. A nonempty
prediction that fails validation is retained for review with Unresolved coverage;
it is not silently used for propagation. Mapping_Selection shows only reaction
information, selected tool, review flag, reason and mapped SMILES. Every selection
row is retained in Final_Mapping. AI review also targets missing mapped strings,
even if ReviewRequired is FALSE, except explicit NoCarbon/Boundary records.


Manual structures in data/manual_smiles.csv take precedence over retrieved
structures. Verified partial carbon correspondences in data/manual_atom_mappings.csv
take precedence over mapper selections when exact inputs match; the current file
contains a user-designated IspG carbon-backbone record. Missing-SMILES participants
may be omitted from mapper-only inputs, retaining original equations/fluxes and
marking partial coverage for optional AI review. See
[manual data and reduced inputs](docs/carbon-mapping-rules.md#user-maintained-manual-data-has-priority).


Manual SMILES use exactly three columns: Species, MetaboliteName, SMILES. Species
is the full underscore-separated binomial used by result prefixes, without the
non_IUPD/IUPD condition suffix. For example Escherichia_coli shares one isoprenol
record across conditions and compartments. Terminal compartment suffixes are
ignored for name matching; names are otherwise exact. Local model metabolite IDs
are not manual-table keys. Conflicting structures for one species/name, or an
explicit model formula/charge incompatible with an override, raise an error.
The original detailed workbook IDs and review archives remain available for audit.


AI-assisted discoveries are recorded separately in data/ai_smiles.csv and
data/ai_atom_mappings.csv, each with a Source column for evidence URLs and archive
references. Manual structures/mappings keep priority. Saved AI structures only fill
missing input structures and require review; the AI mapping CSV is a reviewable
index, not an automatic override. See [manual and AI records](docs/carbon-mapping-rules.md#user-maintained-manual-data-has-priority).
