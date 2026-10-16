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
    -> <model>_3_carbon_flow_Enumeration.xlsx
    -> <model>_3_carbon_flow_EMU.xlsx
```

Step 3 includes experimental [flux-weighted single-carbon source contributions](docs/source-contributions.md)
alongside Enumeration. Both reuse the same mappings and preserve prior results;
their assumptions and statistics differ.

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
| S. cerevisiae | [Yeast9 v9.1.1](https://github.com/SysBioChalmers/yeast-GEM/tree/v9.1.1/model) | `data/yeast9/yeast-GEM.mat`, `yeast-GEM.yml` |
| E. coli | [BiGG iML1515](https://bigg.ucsd.edu/models/iML1515) | `data/ecoli/iML1515.mat`, `iML1515.json` |
| Y. lipolytica | [Pinned iYali5 source](https://github.com/bioengtaltech/eciYali5-GEM/tree/3b97d38d92327a099f3b0f6f9cd92e8aeab78c10/model) | `data/yarrowia/iYali5.mat`, `eciYali5-GEM.yml` |

Keep each source model in its own subfolder under `data/`, with its companion
JSON/YAML files in the same folder. Configure the computational model path only
once in step 1. The generated plan records its source path and SHA256; step 2
automatically reads supported native structures and scans that folder for matching
annotations. Companion filenames need not match the computational model.

MAT is an upstream model **data format**, read directly by COBRApy and SciPy.
Model files, large datasets and environments are not bundled. Ordinary iYali5
supplies the Yarrowia pFBA network; its paired eciYali5 YAML supplies chemical
structures only. Do not substitute the enzyme-constrained network.

## 3. Configure and run step 1

Edit the user settings and relevant named model blocks in
[step1_prepare_carbon_tracing_inputs.py](step1_prepare_carbon_tracing_inputs.py):

```bash
python3 step1_prepare_carbon_tracing_inputs.py
```

Step 1 uses the configured `.envs/preparation/bin/python` when it is installed,
even if the launching interpreter also has COBRApy. You can invoke it explicitly.
Paths are resolved from the project directory, independently of the working folder.

The six numbered blocks use standard COBRApy operations directly:
`gene.knock_out()`, `Metabolite`, `Reaction.add_metabolites()`, explicit
`lower_bound` / `upper_bound` assignments, `model.optimize()` and `pfba()`.
Each strain's modifications stay together, even where operations repeat between
strains. Model versions are not enforced; source paths and hashes are recorded.
The non-IUPD growth objective
is set explicitly. Maximum-growth solutions are checked before use; production
solutions are validated by the shared exporter before writing any result.
See the [COBRApy model-building tutorial](https://cobrapy.readthedocs.io/en/latest/building_model.html)
for these operations.

The default blocks prepare non-IUPD/IUPD pairs for Saccharomyces_cerevisiae,
Escherichia_coli and Yarrowia_lipolytica. Each IUPD block independently reloads
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
| Y. lipolytica | 1 | 2 | 2.09 |

Limits are mmol/gDW/h. S. cerevisiae and E. coli inherit native glucose bounds;
iYali5 has glucose closed, so each Yarrowia block explicitly selects glucose
medium with `glucose_uptake = 1.0`. IUPD retains the growth-phase glucose bounds. Native oxygen
settings are inherited. Other carbon uptake is closed in the E. coli and
Yarrowia comparisons. E. coli non-IUPD/IUPD both receive functional SQS activity and
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

Step 2 applies model-scoped overrides from `data/manual_smiles.csv`, preserves
existing input structures, then uses compatible source/companion annotations before
online structure completion. It checks metabolite identity rather than choosing
a file by species name. Conflicting structures and name-selected database results
mark affected reactions for review. Structure sources and candidates are recorded
in the enriched metabolite table. Inputs remain unchanged. Inspect `ReviewRequired=TRUE` and, if necessary,
correct `Final_Mapping` before step 3. Rerunning step 2 can replace those edits.
`AIGenerated` is stored only in `Final_Mapping` and defaults to FALSE.
The three-tool mapping stage does not perform AI review.

## 5. Analyze carbon origins

```bash
python3 step3_carbon_flow.py
```

`ANALYSIS_METHODS = ["Enumeration", "EMU"]` at the top of step 3 runs both
analyses by default. They read step 2 independently and never invoke mappers or
an LLM. EMU runs first, so it does not depend on Enumeration completing.

| Method | Output | Meaning |
|---|---|---|
| Enumeration | `<model>_3_carbon_flow_Enumeration.xlsx` | Possible molecular label states, by propagation round |
| EMU | `<model>_3_carbon_flow_EMU.xlsx` | Flux-weighted single-carbon and mean source contributions at steady state |

AI-reviewed mappings add `_AIGenerated` before `.xlsx` for either method.

EMU here means **size-1 EMU-inspired origin balances**, not full multi-carbon
isotopomer simulation or fitted experimental isotope abundances. It reads original
and optional AI-reviewed mappings separately, requiring the original step-1
input/plan for either version. Missing/ambiguous mappings remain Unknown. The Excel
has only MetaboliteContributions and BoundaryInputs; calculation checks
and provenance remain in `source_contributions_run.json`. Symmetry-related
positions retain their source when all possible precursor positions are fully
known to share it; unresolved mixed-origin ambiguity remains Unknown.
See [method, assumptions and result columns](docs/source-contributions.md).

Enumeration retains its default CO2/bicarbonate label termination. EMU includes
internal CO2 refixation in the fixed-flux source balance. Both outputs record their
method and inorganic-carbon policy; their outputs describe different quantities. Do not interpret them as the same calculation under identical
assumptions. `ROUNDS` and `ALLOW_INORGANIC_CARBON_RECYCLING` control Enumeration
only. To select one method for a run:

```bash
python3 step3_carbon_flow.py --method Enumeration
python3 step3_carbon_flow.py --method EMU
```

Enumeration reads both normal and optional AI-reviewed step-2 mappings for models
in the selected plan, using the same ordered carbon sources for both versions.
It never calls mappers or changes input/mapping workbooks. Each
`<model>_3_carbon_flow_Enumeration.xlsx` contains only `Carbon_Flow_R1` onward and `Summary`.
The default maximum is 10 rounds, with earlier stopping on convergence.

`MAX_EXCEL_SIZE_MB = 50` applies separately to each output. EMU saves a complete
workbook atomically and rejects an oversized replacement. For Enumeration: Generated rows are staged in temporary disk storage rather than retained
as a whole-round DataFrame or an in-memory Excel workbook. `OUTPUT_BATCH_ROWS =
50000` controls size-check batches; checks also occur at layer and round boundaries.
At a checkpoint, the program streams a real candidate XLSX and measures its
compressed size. If it reaches 50,000,000 bytes, enumeration stops and the
unfinished round is discarded. Only complete rounds below the limit atomically
replace the output; the preceding complete result is retained unchanged.

There is no worksheet pagination. A round exceeding 1,048,575 data rows cannot fit
in one sheet and is discarded with `row_limit`. Size stops use `file_size_limit`,
never convergence. If round one cannot be saved, the report records zero completed
rounds and no current output; an existing file is from an earlier run.

Checks occur at batch boundaries, not at every byte: temporary candidates may
exceed the limit before rejection. Rebuilding sorted candidate files costs disk
space and I/O. Label states and their indexes are still needed in memory; this is
not full external-memory propagation, and independent memory safety checks remain.

No aggregate Summary workbook is generated. Enumeration retains its per-file
Summary for manual plotting or fitting. The Enumeration manifest is
`carbon_flow_run.json`; the EMU manifest is
`source_contributions_run.json`. A failed method does not cancel the other method,
but the process returns an error. Inspect both manifests rather than trusting old files.

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

### Optional: propagate using AI-reviewed mapping copies

After external review, place `<model>_2_atom_mapping_AIGenerated.xlsx` copies
in the configured results directory, keeping the original model names and input
hashes. Successful replacements belong in `Final_Mapping.MappedReactionSMILES`,
with `AIGenerated=TRUE`; leave `Mapping_Selection` unchanged as the original tool-selection record. Preserve raw tool predictions and all unchanged reactions.

```bash
python3 step3_carbon_flow.py
```

Both methods automatically process both `<model>_2_atom_mapping.xlsx` and any
`<model>_2_atom_mapping_AIGenerated.xlsx` in the selected plan's output folder.
Each version produces independent Enumeration and EMU outputs, retaining the
`_AIGenerated` suffix for reviewed copies. Either mapping version can be present
alone. EMU still requires the original step-1 input/plan. Excel lock
files, temporary files and unrelated names are ignored. Only models in the plan
are included; model selection applies to both available versions of that model.
Both versions must pass input-hash, direction and structure checks. The suffix
selects a file; it does not certify correctness. Step 3 performs no AI review or
mapper prediction. The separate step-4 entry has been removed.

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

The default commands reproduce the six-model study. A different analysis can
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
`carbon-tracing-ai-review-skill` delegates flagged reactions to
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
