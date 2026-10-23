# Detailed workflow and custom analyses

Follow the [quick start](../README.md) for the four reference models. All stages
run in Python; COBRApy and GLPK prepare the pFBA networks.

```text
Metabolic model + COBRApy -> pFBA Excel
                              |
            three independent atom-mapping tools
                              |
                  consensus / Final_Mapping
                              |
              carbon propagation -> Excel results
```

## 1. Install Python environments and EnzymeMap data

Run these commands in a **Linux or WSL terminal** with Conda and Git available:

```bash
bash scripts/setup_environments.sh
python3 scripts/download_enzymemap.py
```

The specifications target Linux/WSL **x86_64**; native Windows, macOS and ARM
environment installation is not supported by them. In WSL, installing under
the Linux filesystem is faster than extracting environments under `/mnt/c`
or `/mnt/d`.

If Conda is not available, install [Miniforge](https://github.com/conda-forge/miniforge)
and initialize the current shell before setup. On WSL, a Windows localhost
proxy may not be reachable from Linux; use a proxy address reachable from WSL
when required by your network. This is a local network setting, not a project
configuration value.

If a network failure interrupts the pip part of environment creation, complete
the affected environment before rerunning its checks. For example:

```bash
conda env update --prefix .envs/enzymemap --file environments/enzymemap.yml
bash scripts/setup_environments.sh
```

This retains the existing environment. Substitute the affected mapper's name
for `enzymemap` if necessary; do not delete unrelated environments.

Four isolated environments are created in `.envs/`. Existing environments at
those paths are checked and reused, never deleted. The `preparation` environment pins COBRApy 0.30.0, swiglpk 5.0.12, NumPy 1.26.4, SciPy 1.13.1, pandas 2.2.3 and openpyxl 3.1.5 on Python 3.11. Its smoke test solves a real GLPK pFBA model. Checks also include actual
LocalMapper/RXNMapper calls and actual rdchiral rule execution for EnzymeMap's
dependencies. The full EnzymeMap adapter is exercised by the example after
downloading its data. An incompatible environment fails the check.

The LocalMapper/tracing stack pins Python 3.6.13, RDKit 2021.9.4, pandas 1.1.5
and openpyxl 3.0.10 to preserve study serialization. lxml 4.9.4 enables the
native Excel XML backend. To update an existing tracing environment, run
`.envs/localmapper/bin/python -m pip install lxml==4.9.4`, then rerun the setup
checks. Python 3.6 is end-of-life;
this is a research reproduction environment, not an internet-facing service.
All specifications are in [environments/](../environments).

For a custom prefix, run `bash scripts/setup_environments.sh /path/to/envs`
and set `ENVIRONMENT_DIR` in step2_atom_mapping.py and `TRACING_PYTHON` in
step3_carbon_flow.py. Existing compatible mapper environments can be referenced
by step 2's `PYTHON_OVERRIDES` dictionary containing Python paths for LocalMapper,
RXNMapper and EnzymeMap. Step 3 requires the compatible tracing environment;
it calls no mapper.

The downloader obtains official **EnzymeMap v2** CSV and rules from commit
`c91578856f11312d7941b290f8b773eb41b54973`. The compressed CSV is about 13 MB.
Before unpickling rules, SHA256 must match
`21cdebccdda5ea4308555b6dd8dfb27413950725ace66e7b8cb7ed3dba2c1e1a`.
Only load trusted official rule files. Data and model weights are not bundled.
On a timeout, rerun the downloader or obtain files from the
[official pinned data directory](https://github.com/hesther/enzymemap/tree/c91578856f11312d7941b290f8b773eb41b54973/data)
and configure their paths. Rebuilding the original BRENDA-derived database is
outside this workflow.

## 2. Prepare input workbooks

Step 1 keeps six numbered strain blocks. Edit the relevant block to change its
medium, knockout, introduced reactions or objective. Reaction stoichiometry is
written explicitly as `{metabolite: coefficient}`: negative coefficients consume
substrates and positive coefficients produce products. Setting both biomass
bounds to the same value fixes growth. The direct `pfba()` call then retains
optimal production while minimizing total absolute flux. Numerical checks and
export helpers do not construct or alter these biological modifications.

Provide one Excel workbook per metabolic network or culture condition, with
nonzero reaction fluxes and metabolite identities. Structures are optional at export. The Python workflow accepts
any model that supplies the following schema; non-IUPD/IUPD names and Yeast9 reaction
IDs are not required.

| Sheet | Required columns |
|---|---|
| Reaction sheet (any name) | `ReactionID`, `ReactionName`, `Flux`, `Equation` |
| Metabolite sheet (any name) | `MetaboliteID`, `MetaboliteName` |

No SMILES columns are required. If the model already has structures, optionally
export `SMILES_Model`; the legacy `SMILES_Yeast9` column remains supported.
`SMILES_Manual` takes priority when supplied. Optional `MolecularFormula`,
`Charge`, `InChIKey` and `PubChemCID` aid step-2 structure retrieval.
Step 1 records the source model automatically when `read_model()` loads it;
the generated analysis plan carries this provenance alongside the Excel input. Step 1 does not
retrieve SMILES or require a second annotation-file setting. Put each original
model and its companion files in a separate `data/<Full_species_name>/` directory. Do not mix
unrelated model versions in that folder.

Step 2 verifies the recorded source SHA256 and examines that file plus JSON/YAML
files in the same directory. It does not recursively scan other species. JSON
source/license records without model metabolites are ignored. Annotations can
come from native MAT fields (`metSmiles`, `metSMILES`, `smiles`, `SMILES`), COBRA
JSON metabolite `smiles`/`annotation.smiles`, or the supplied COBRA/RAVEN ordered-map
YAML layouts. A scalar SMILES or the first list entry is used. Unsupported layouts
do not establish a structure; existing input structures and database lookup remain
available. SBML model preparation is supported, but native SBML structure-note
extraction is not implemented.

Priority is manual table, existing input structure, original model, then compatible
companions. Same-stem companions precede other filenames, followed by alphabetical
order. IDs, names, compartments and available formulas must agree; matching IDs
alone do not identify companion compounds. Different canonical structures, including
stereochemical or charge differences, set `StructureReviewRequired`. The preferred
candidate remains available for mapping; affected reactions inherit
`ReviewRequired=TRUE`. `StructureReference`, `StructureCandidates` and
`StructureStatus` record the choice and alternatives. This is not biochemical
validation. Conflicting database identifiers are not automatically selected.

Yeast9 supplies S. cerevisiae structures. Compatible companion InChIKeys are
loaded automatically for exact database lookup, including iML1515 JSON annotations.
No generic lipid/polymer structure is invented.
PubChem supplies remaining missing structures. Local sources still work when
`COMPLETE_MISSING_SMILES=False` disables online lookup.

Rerun step 1 once to add source provenance to a plan created by an older version.
Externally prepared Excel inputs without source provenance can still use existing
structures, the manual table and PubChem; no original model is guessed from their names.
Names must include compartments, such as `D-glucose[c]`, and match the names in
reaction equations. Equations use `=>`, `->` or `<=>`, terms separated by ` + `,
and positive integer stoichiometry, e.g. `A[c] + 2 B[c] => C[c]`.
Flux sign determines the active direction. Fractional reactions and missing
structures cannot form complete mapper inputs. After trying usable mapper carbon
correspondences, step 2 may recover conservative partial skeleton transfers;
other carbon remains unresolved. Generic biomass/macromolecule pools are not
filled by PubChem name search with an arbitrary small-molecule representative.
See [carbon mapping rules and evidence](carbon-mapping-rules.md).

[step1_prepare_carbon_tracing_inputs.py](../step1_prepare_carbon_tracing_inputs.py)
contains the model, medium, knockout, objective and carbon-source settings.
The default setup exports the non-IUPD/IUPD pairs for S. cerevisiae and E. coli.
Both models of each species fix growth at 50% of the non-IUPD maximum and maximize
squalene formation before minimizing total absolute flux. The E. coli comparison
models both include introduced SQS activity and a squalene output.
Both species preserve native glucose bounds during non-IUPD growth optimization.
After fixing half of the maximum growth rate, non-IUPD glucose uptake limits
double: S. cerevisiae 1 -> 2 and E. coli 10 -> 20 mmol/gDW/h.
IUPD retains the growth-phase glucose limits (1 and 10), permits isoprenol uptake
up to 2.09 times those limits (2.09 and 20.9), and uses the same fixed growth rate
as its paired non-IUPD, not half its own maximum. Uptake limits permit uptake;
they do not force it. Native oxygen bounds are inherited (lower bound -1000 in
both reference models). Each non-IUPD block sets its biomass objective explicitly. Set a
different medium bound or objective explicitly when adapting another model.
IUPD uses its paired non-IUPD growth reference before selecting the squalene objective.
pFBA supplies one representative solution, not all possible pathways.

The COBRApy exporter also writes `data/analysis_plan.json` with portable paths,
ordered sources, source-model provenance, the output folder and input hashes. Do not edit this generated file; change
step 1 and rerun it. Copy it with the Excel inputs and the referenced data subfolders when moving between computers.

## 3. Configure a custom network and carbon sources

For model/input preparation, edit step 1. Shared solver/output settings are at the top.
Every named strain block follows the same six sections:

1. Set `model_file` and load an independent model.
2. Set the ordered `carbon_sources` list.
3. Set culture conditions and introduce model modifications.
4. Set growth constraints: calculate the non-IUPD maximum or use the paired non-IUPD reference.
5. Set the production objective and run pFBA.
6. Register the model, solution and sources for export.

Source-list order controls label introduction; it does not change the pFBA medium.
Source paths and hashes are recorded, without requiring a specific model version. Required
reaction and gene IDs still need to exist. Growth solutions are checked before
setting biomass bounds. The exporter checks production status, finite fluxes,
steady-state residuals and flux bounds before creating a workbook. Each block saves its
model, solution and sources together for automatic export. There is no separate
model-to-carbon-source table to maintain. Validation and export helpers are in
`carbon_flow_tools/model_preparation.py`. Set maximum propagation rounds with `ROUNDS` at
the top of step3_carbon_flow.py. Execution settings belong to the Python steps:
step 2 configures tools, batch/atom limits, mapping resources, structure completion
and mapper/data paths; step 3 configures rounds, propagation resources and its interpreter.
The generated plan does not override those execution settings. The atom limit
counts both reaction sides. Optional command-line overrides take precedence for
that run. Older generated plans are rejected with a request to rerun step 1.

For a new analysis, keep a separate preparation script rather than editing the
four-model study script. Copy the supplied readable example once:

```bash
mkdir -p results/results_by_AI/my_analysis
cp examples/step1_prepare_custom_model_inputs.py results/results_by_AI/my_analysis/step1_prepare_inputs.py
# Edit this copy's COBRApy model blocks and carbon sources, then run:
python3 results/results_by_AI/my_analysis/step1_prepare_inputs.py
python3 step2_atom_mapping.py --plan results/results_by_AI/my_analysis/analysis_plan.json
python3 step3_carbon_flow.py --plan results/results_by_AI/my_analysis/analysis_plan.json
```

The copied script places its generated plan and three numbered output workbooks
per model under `results/results_by_AI/my_analysis/`. Environments, implementation modules and
raw models remain shared. No project copy or hand-written JSON is needed.
The example compares glucose alone with simultaneously supplied glucose and
xylose in iML1515, maximizing growth. Its limits are illustrative capacities,
not forced uptake or the study's squalene conditions. The generated plan records
actual exchange fluxes and warns about declared sources with no net uptake.
Such a source remains in the metabolite table and source list; unused sources
do not enlarge the mapped-network Summary denominator. Valid source structures
are still needed in step 2, including for unused declared sources.

`--plan` also works with the default step-1 script and with step 3. If step 1
receives a nondefault plan path, its outputs default to a `results/` folder next
to that plan; `--output-dir` overrides this location. Paths inside a plan are
relative to the project root. Only the selected plan is replaced or invalidated;
the default `data/analysis_plan.json` is untouched by a separate analysis.
See [interfaces](interfaces.md) for the small calling contract used by future skills.

To run fewer reference models, omit their complete model blocks. Keep the matching non-IUPD block
when retaining an IUPD block: its non-IUPD maximum supplies the shared growth reference.
Each model still loads its original MAT file; IUPD does not inherit non-IUPD production
bounds.

For Excel inputs prepared elsewhere, edit `INPUTS` and `OUTPUT_DIR` at the top
of [examples/step1_use_existing_inputs.py](../examples/step1_use_existing_inputs.py),
then run it instead of the model exporter:

```bash
python3 examples/step1_use_existing_inputs.py
python3 step2_atom_mapping.py
python3 step3_carbon_flow.py
```

Each entry supplies an analysis name, an Excel path and ordered carbon-source
names. The example generates the same portable plan and hashes without changing
the input workbooks. It replaces the plan's analysis list; include all networks
you want to process in `INPUTS`. No manual JSON editing is needed.

The input schema above remains unchanged. Reaction and metabolite sheets are
identified by their headers. Steps 2 and 3 use the exact model list in the plan;
unrelated old files are ignored. Each analysis name must be unique and contain
letters, digits or underscores, starting with a letter. Paths are relative to
the project, not the terminal working directory.

Source names must exactly match compartment-qualified metabolites. Source order
controls label introduction, not physical uptake order or diauxic growth.
Medium/uptake constraints belong to the upstream model. Mixed carbon sources
are represented by multiple entries in the ordered source list.

Run step 1 again after changing its inputs or sources. Changes to
Python execution settings do not require a new input export. If export fails, no usable plan
remains: fix the error and rerun step 1 before proceeding. `--example` bypasses
the plan for the supplied installation test.




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

Step 1 reports prepared condition counts, numbered exports, target reaction IDs
and actual target fluxes in model units. Its pFBA total absolute flux is a separate
diagnostic shown with `python3 step1_prepare_carbon_tracing_inputs.py --verbose`.
The same option restores known legacy MAT field/compartment notices; unrelated
warnings and errors remain visible.

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


## 4. Run carbon propagation

All three stage scripts run fixed code without LLM or skill calls. A future
AI orchestrator is expected to include reaction review by default after step 2:
the review skill processes flagged and blank-mapping rows (except explicit NoCarbon/Boundary), delegates biochemical mapping
to Enzymatic Atom Mapping Skill, and writes `_AIGenerated` copies. Step 3 then propagates both normal and reviewed copies independently. Failed
reviews retain their original mappings and review flags. Neither deploying the
environments nor running step 2 alone executes this review.

Run step 3 after the mapping workbook has been generated:

```bash
python3 step3_carbon_flow.py --rounds 10
```

This command overrides step 3's `ROUNDS` setting to run up to 10 propagation rounds. Step 3 needs
the mapping workbook, generated plan and tracing interpreter; it does not need
the original step-1 workbook or the prediction tools. Input SHA256 recorded in
Mapping_Info must match the plan, preventing accidental use of a same-name old
mapping after step 1 is rerun. Run step 2 again if this check fails.
Later rounds reuse the selected
atom correspondences without calling the prediction tools again.

`ALLOW_INORGANIC_CARBON_RECYCLING = False` in step 3 is the default analysis
assumption. CO2, bicarbonate, carbonic acid and carbonate are identified by
structure, independently of model names and compartments. Their labelled
production records remain in Carbon_Flow and Summary, but those states cannot
supply labels to later reactions or rounds. Their substrate `SourceIDs` are
therefore `0` in downstream reactions, even if an earlier endpoint record exists.
Other substrates of the same reaction still transfer their labels normally;
reaction fluxes and the selected mappings are unchanged. For example,
labelled pyruvate plus bicarbonate produces oxaloacetate retaining the three
pyruvate labels, with no inherited bicarbonate label. The mapped-network Summary
denominator is unchanged.

Set `ALLOW_INORGANIC_CARBON_RECYCLING = True` to restore the original recycling
behavior, compare the two assumptions, or study CO2/bicarbonate as a carbon
source. An inorganic source with recycling disabled raises an explanatory error
rather than silently giving no propagation. Normal and AI tasks use the same setting.
The default deliberately limits the scope of qualitative reachability; it is
not evidence that CO2 all escapes or that its biological contribution is small.

Each rerun starts from configured sources and replaces old carbon-flow and
Summary sheets. It does not resume a previous unfinished enumeration.

Memory protection is active for normal and AI outputs in step 3, independently of the startup
concurrency estimate. The common implementation checks /proc/meminfo and
available cgroup-v2 ancestor limits before work, at intervals within substrate
combination enumeration and state merging, and during Excel row creation.
The minimum free reserve is max(MIN_AVAILABLE_MEMORY_GB, 10% of detected total)
with a default of 2 GiB. Batch flushing, disk-backed state merging and streaming workbook writing
also check the memory reserve. These estimates and periodic checks
are safeguards, not hard allocation bounds: a sudden allocation, another process,
or Windows-host memory pressure not reflected in WSL can still exhaust memory.

A memory safety stop raises an error, preserves the last successfully saved
workbook, and reports how many rounds were saved in this run. If none were saved,
any existing workbook is from a previous run. An unfinished round is never
reported as converged or complete, and no states/records are truncated.
The incomplete round is not saved and there is no automatic resume; a rerun
starts at round 1. Free memory, run fewer models together, or set
MAX_CARBON_PROCESSES = 1 before rerunning. Changing a script does not retrofit
memory protection into an already-running process.

Each propagation batch writes `carbon_flow_run.json` in its output folder.
It records each mapping workbook's model name, mapping path/hash, requested/saved
rounds, output path and stop reason, including both normal and AI versions. A file-size stop
does not imply convergence; zero saved rounds means any existing workbook is
from an earlier run in a combined task. EMU-only tasks save independently.
Always check the command's exit status as well: an input
validation error can occur before a new report is written. See
[report fields](interfaces.md#propagation-status).

- Round one introduces sources in configured order, retaining the states
  reached from earlier sources. Each reaction is processed at most
  once per source phase; current-layer products become usable next layer.
- Later rounds enumerate from the full accumulated state snapshot and merge
  results at the round boundary. Stop at the requested round or when neither
  states nor earliest steps change.
- Enumerate all known labelled substrate states. Use an unlabelled state only
  when no labelled state exists. This is a reachability rule, not a physical
  abundance or substrate-competition model.
- Map repeated stoichiometric copies separately and consider positional
  permutations to retain asymmetric products. Later-round repeated-term IDs
  are sorted; identical output rows are removed.
- Transfer labels through fixed, structurally aligned correspondences. Include
  all labelled products, including cofactors. Reuse `ProductionID` only when
  compartment-qualified name and canonical labelled SMILES both match.

Source labels occupy consecutive, nonoverlapping ranges in list order. A
six-carbon source followed by a two-carbon source receives C1-C6 and C7-C8.
Supported glucose and isoprenol structures have dedicated numbering rules;
other sources use carbon order in the supplied single-component SMILES, which
must not be interpreted as standard biochemical numbering. Glucose requires
the supported six-membered ring form. See the
[reference example](yeast9-example.md) for its numbering details.
Reaction-local map IDs are bridges, not carbon-origin labels.

### Optional AI-reviewed mappings

External review may supply `<model>_2_atom_mapping_AIGenerated.xlsx` copies.
Keep the originals, raw tool predictions, copied inputs and hidden Mapping_Info
with its SourceInputSHA256. Replace validated correspondences in Final_Mapping,
set AIGenerated to TRUE only for successful replacements, and leave
Mapping_Selection unchanged as the original tool-selection record. Every reaction
now has a final row; set MappingCoverage=Complete for full accepted replacements.
On older files, append successful formerly missing final rows. Confirmed rules
remain Partial and AIGenerated=FALSE. Scientifically rejected/unresolved rules
must become Unresolved and must not propagate; keep saved partial records in Mapping_Info as audit context. AI provenance
alone does not clear ReviewRequired or establish biochemical correctness.

```bash
python3 step3_carbon_flow.py
python3 step3_carbon_flow.py MyCulture --rounds 1
```

Place these copies in the plan's output directory (normally results/).
Step 3 processes `<model>_2_atom_mapping.xlsx` and its optional
`_AIGenerated` copy separately. Either version may be present alone.
Both Enumeration versions use the same source lists, rounds and resource settings and retain
input-hash/direction/structure checks. They produce independent
`<model>_3_carbon_flow.xlsx` and `<model>_3_carbon_flow_AIGenerated.xlsx` files.
The mapping inputs are unchanged. Lock files, temporary files, backups and models
outside the selected plan are not included. Step 3 runs no mappers or AI skill;
review and workbook updates happen separately before propagation. The separate
step-4 entry has been removed.

## 5. Inspect Excel results

Step 2 writes `<OUTPUT_DIR>/<analysis>_2_atom_mapping.xlsx`; step 3 reads it
and writes `<OUTPUT_DIR>/<analysis>_3_carbon_flow.xlsx`
(or `_AIGenerated.xlsx` for a reviewed copy). The default output
directory is `results`. Raw step-1 inputs use `<model>_1_` filenames.
Step-1 inputs remain unchanged. Close results in Excel before replacement;
for step 2, a denied rename leaves the completed `.tmp.xlsx`. Step 3 retains
the preceding complete output if publication fails and cleans its temporary
staging directory. No permanent intermediate result directories are created.

The step 2 workbook contains copied input tables, three raw tool sheets,
`Mapping_Selection`, `Final_Mapping` and hidden `Mapping_Info` diagnostics, including partial carbon correspondences.
It is unchanged by step 3. The combined step-3 workbook contains `Summary`,
one `Carbon_Flow` sheet containing R1-R3 details by default (the first column,
`Carbon_Flow_Round`, identifies each row's round), `EMU_Metabolites`, and `EMU_BoundaryInputs`.
`EXPORT_DETAIL_ROUNDS=3` limits detail export, not computation: later complete
rounds still update Summary. Set it to 0 for no detail sheets, or to ROUNDS for
all requested rounds. There is no pagination. Select methods with
`ANALYSIS_METHODS=["Enumeration", "EMU"]`; each can also run alone.
EMU uses the same ordinary/AI mapping version independently of Enumeration.
See [source-balance assumptions](source-contributions.md).
Input structures, mapper provenance and review flags remain in step 2.
`AIGenerated` appears only in `Final_Mapping` and defaults to FALSE.

Raw prediction columns are `ReactionID`, `ReactionName`, `Equation`,
`ReactionSMILES`, `MappedReactionSMILES`, `MappingDirection`, `MappingConfident`.
Long SMILES remain complete in their cells; use Excel's formula bar to inspect
or copy them without expanding row heights.

| Carbon-flow column | Meaning |
|---|---|
| `Step` | Earliest possible reaction depth from selected states, not elapsed time |
| `ReactionID` | Reaction supplying this transfer |
| `Source` | All substrates in the active flux direction |
| `SourceIDs` | Selected state IDs in substrate order; 0 means no known labelled state |
| `Production` | Product including compartment |
| `ProductionID` | Identifier of name plus canonical labelled SMILES |
| `CarbonLabels` | Origin labels present; repeats are retained |
| `LabeledSMILES` | Product structure with origin labels in carbon atom-map fields |

For `2 B`, `12&14` denotes two copies using states 12 and 14; ` + ` separates
equation terms. Seed sources also use ID 0 but are labelled internally, so their
transport row can have `SourceIDs=0` and a labelled product. IDs are not temporal
order and need not all appear as product IDs in every round.

Summary's denominator counts distinct carbon-containing **MetaboliteID** values
in usable final-mapping substrates and products, including identity transport.
Unmapped-only metabolites are excluded; mapped but unreachable ones remain.
Compartments are distinct. Counts use products in that round's records; mixed
products count for each applicable source. `All sources` is their union and
appears only with multiple sources. Increments compare adjacent rows and can be
negative. Zero denominator gives a blank percentage.

For multiple sources, `MixedSources_LabeledCarbonMetabolites` counts distinct
MetaboliteIDs with at least one recorded labeled state containing carbon labels
from two or more configured sources in that round. Separate single-source states
of the same metabolite do not qualify. Repeated states count once, compartments
remain distinct, and the mapped-network denominator is unchanged. This count
indicates possible mixed-origin label states, not measured co-assimilation or
carbon-exchange efficiency. It is omitted for single-source analyses.

The last saved round's Summary row also reports
`<source>_ExclusiveSourceMetabolites` for each configured source; earlier rows
leave these columns blank. A metabolite qualifies only if the union of labels
across **all its recorded states in that round** is nonempty and belongs to that
one source. A subset of that source's labels qualifies. A metabolite with separate
pure-source states from different sources does not qualify, even if no single
state mixes them. Counts use distinct mapped MetaboliteIDs, with compartments
kept separate, and do not change the denominator. This describes observed labels,
not proof that every carbon in the molecule comes from that source; unlabeled
positions and uncomputed states are not evidence of exclusivity.

Set `EXPORT_EXCLUSIVE_SOURCE_METABOLITES = True` at the top of
`step3_carbon_flow.py` to include an `Exclusive_Source_Metabolites` sheet with
`Round`, `CarbonSource`, `MetaboliteID` and `MetaboliteName` for the last complete
saved round. The default is `False`; counts are always included. The optional
sheet is subject to the same workbook size limit. Normal and AI-reviewed mapping
versions follow the same rules, and the Summary collector copies these columns.

Results describe **qualitative carbon-origin reachability under the chosen pFBA
network and mapping assumptions**, not measured 13C abundance, quantitative
carbon allocation, isotope kinetics or one molecule's physical history.
Multi-round combinations can grow rapidly even without rerunning mappers.
Uncertain mappings, incomplete reactions and polymer-as-monomer structures can
limit position-specific interpretation. Coverage concerns the mapped subnetwork.

## Workbook size limit

Normal and AI results use `MAX_EXCEL_SIZE_MB` at the top of step 3 (default
50 decimal MB, or 50,000,000 bytes). `OUTPUT_BATCH_ROWS` defaults to 50,000
generated records, with an additional estimated 8 MiB buffer bound. Exported rows are batched into temporary SQLite storage using the
Python standard library; no additional dependency or user-managed database is
required. Deduplication and output sorting occur on disk. Label states merge at
the same layer/round boundaries as before.

At each batch checkpoint, and at layer and round boundaries, a candidate XLSX is
streamed from staged rows and its actual compressed size is checked. Reaching
the limit stops further enumeration and discards the candidate round. A complete
round is published atomically only below the limit. The saved output and its
Summary therefore contain complete rounds only. Other models continue normally.

There is no pagination or truncation. If the combined Carbon_Flow worksheet exceeds 1,048,575
data rows, stop with `row_limit` and retain the preceding complete round. File
size uses `file_size_limit`; neither is convergence. If the first round is
rejected, zero Enumeration rounds are reported. For a combined run no new output
is published until a complete Enumeration round is saved. EMU can be exported
independently by selecting only EMU; emu_solved/emu_completed distinguish solving
from delivery. Existing historical files are not
deleted or relabeled.

Size checks are bounded by batches, not instantaneous per-byte checks. A temporary
candidate may exceed the limit before rejection; the final complete output stays
below it. Smaller batches detect growth sooner but rebuild the sorted workbook
more often. Temporary storage requires disk space and is removed when the task
ends. This removes full output-row tables from RAM; accumulated label states,
catalogs and other working data still require memory. The existing memory guard
can still stop a run safely. There is no automatic resume. Failed saves never
replace the preceding complete workbook.

Summary contains computed observations only. Step 3 does not fit
curves, adds predicted rounds or changes the denominator. Use these values for
manual plotting or fitting as needed.

Each model/version keeps its Summary in its own result workbook. Step 3 does
not generate a separate summary-collection workbook. Use `carbon_flow_run.json`
to distinguish current outputs from older files. A memory-stop error preserves the rounds successfully saved in that run;
if no rounds were saved, any earlier output remains unchanged.

## Small example and tests

After setup and data download:

```bash
.envs/localmapper/bin/python examples/make_example.py
python3 step2_atom_mapping.py --example
python3 step3_carbon_flow.py --example --rounds 1
```

The supplied example uses synthetic nonzero fluxes for acetaldehyde transport
and NADH-dependent reduction. It is not a Yeast9 pFBA result. The reduction is
actually submitted to all three tools. Expected labelled products are cytosolic
acetaldehyde and ethanol; identity transport is not three independent predictions.

```bash
.envs/preparation/bin/python -m scripts.test_model_preparation
.envs/localmapper/bin/python -m scripts.test_carbon_tracing
.envs/localmapper/bin/python -m scripts.test_mapper_backends
.envs/localmapper/bin/python -m scripts.test_consensus_mapping
.envs/localmapper/bin/python -m scripts.test_parallel_carbon
.envs/localmapper/bin/python -m scripts.test_stages
.envs/localmapper/bin/python -m unittest scripts.test_model_structures
.envs/localmapper/bin/python -m scripts.test_smiles_completion
.envs/localmapper/bin/python -m scripts.test_script_settings
.envs/localmapper/bin/python -m scripts.test_analysis_plan
.envs/enzymemap/bin/python -m scripts.test_enzymemap_rules
```

Scheduling/failure tests use explicit mock adapters; they do not replace the
real-tool example. See [docs/validation.md](validation.md) for observed
release checks and platform limits. Record the repository commit, configuration,
input hashes, tool versions and external data versions with manuscript analyses.

## Attribution

Authored code is under [MIT](../LICENSE); dependencies retain their own terms in
[THIRD_PARTY.md](../THIRD_PARTY.md). LocalMapper's pinned LICENSE is CC BY-NC-SA 4.0;
this repository's license does not remove its noncommercial restriction.
No upstream implementations, weights or full EnzymeMap dataset are
redistributed. Follow citation instructions for
[COBRApy](https://github.com/opencobra/cobrapy),
[Yeast9](https://github.com/SysBioChalmers/yeast-GEM),
[LocalMapper](https://github.com/snu-micc/LocalMapper),
[RXNMapper](https://github.com/rxn4chemistry/rxnmapper) and
[EnzymeMap](https://github.com/hesther/enzymemap).

The reference MAT files store the model as `model` (Yeast9) or
`iML1515` (E. coli). COBRApy reads the arrays directly. The shared reader also
accepts COBRApy JSON, SBML and YAML; adapt the model-specific IDs, validation and
medium when using another model. Native RAVEN YAML is not necessarily a COBRApy
network format; the step-2 structure reader handles its separate annotation role.
The custom growth example is
`examples/step1_prepare_custom_model_inputs.py`. Run it with
`.envs/preparation/bin/python`; it produces an input workbook and the same plan.

For an isolated output destination, use `--output-dir <folder>`. Step 3 performs
Enumeration and/or size-1 EMU-inspired balances, as selected at the script top.
Enumeration reports possible states and coverage. EMU reports model-dependent
mean source fractions, retains Unknown and does not simulate full isotopomer
abundances. It uses original reaction turnover and saved partial-transfer amounts;
formula-only structures never create invented positional links.


### Maintained manual structures and carbon correspondences

Use data/manual_smiles.csv for verified metabolite structures (highest structural
priority), and data/manual_atom_mappings.csv for verified partial carbon
correspondences (highest mapping priority). A Manual result is still a partial
map of the original reaction; its remaining participants are not invented.
Keep source evidence in the mapping CSV and review atom indices if SMILES change.
See [the rule contract](carbon-mapping-rules.md#user-maintained-manual-data-has-priority).
Missing structures may be omitted from mapper-only inputs; the original equations
and fluxes remain intact. Optional AI review assesses these provisional partial
results; ordinary script execution does not require AI.


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
index, not an automatic override. See [manual and AI records](carbon-mapping-rules.md#user-maintained-manual-data-has-priority).


### Inspecting and clearing AI cache tables

The two AI CSV files are the active, user-editable cache lists. Keep their four-column
headers; the first three columns match the corresponding manual table. Rows are
sorted and saved as UTF-8 with BOM for Excel. Delete data rows to disable those old
AI candidates, or remove a CSV to disable that AI cache. Archived evidence alone
must not restore deleted or edited entries. New accepted research may add a new
record later. Clearing a cache does not alter existing workbooks, manual tables,
model annotations or the independent database structure_cache.json. AI mapping
rows remain review candidates, not automatic step-2 overrides.
