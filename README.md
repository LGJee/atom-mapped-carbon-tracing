# Atom-Mapped Carbon Tracing

Prepare pFBA-active metabolic networks with **COBRApy**, compare three atom-mapping
tools and calculate fixed-flux carbon-source contributions. In
`step1_prepare_carbon_tracing_inputs.py`, each of the four named strain blocks follows
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

Want an AI assistant to install it? Copy the self-contained instruction in
[INSTALL.md](INSTALL.md#ask-an-ai-to-install). Installation does not run the study.
Existing computed workbooks are available in the [results index](results/README.md);
no installation is needed to open them in Excel or LibreOffice.


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
Original model files, full mapper datasets and environments are not bundled.
Selected computed workbooks and their analysis plan are published.

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
step 3 controls EMU numerical settings and resources. Their defaults need no edits.

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
rule-only partials show a reduced unit subreaction in the mapped-SMILES columns.
Their actual positional transfers and molar coefficients are stored in hidden
Mapping_Info records; no Partial_Mapping sheet is generated. These maps are
explicitly Partial and require review. Step 3 uses the original reaction flux
and coefficients, leaving unresolved carbon contributions Unknown.

AI review changes only Final_Mapping in a separate copy. Confirmed partial rules
stay Partial; rejected or scientifically unresolved rules become Unresolved and
are not propagated even though their original audit records remain. A validated
complete AI replacement becomes Complete and supersedes partial records.


## 5. Analyze carbon origins

Run from the project root:

```bash
python3 step3_carbon_flow.py
python3 step3_carbon_flow.py Escherichia_coli_non_IUPD Escherichia_coli_IUPD
```

Step 3 now runs only the existing size-1 EMU-inspired source-balance solver.
It calculates fixed-flux mean carbon-origin fractions, keeping Unknown explicit.
It is not a full isotopomer simulator, isotope-label time course or experimental
measurement. It does not enumerate label combinations or use rounds/ProductionID.

For each model in the selected plan, both the ordinary step-2 file and an
available _AIGenerated copy are analyzed independently. Output names remain
<model>_3_carbon_flow.xlsx and <model>_3_carbon_flow_AIGenerated.xlsx.
Each new workbook contains only EMU_Metabolites and EMU_BoundaryInputs.
Internal validation and structure-scope metadata are retained in the solver but
are not displayed as extra metabolite columns. See [assumptions](docs/source-contributions.md).

At the top of step3_carbon_flow.py, set EMU_ALLOW_INORGANIC_CARBON_RECYCLING
(default True), EMU_TOLERANCE (1e-8), resource limits and the tracing interpreter
only when required. Normal and AI outputs use the same settings. Step 3 invokes
no mappers and no AI. External inputs not named as configured sources receive
separate External categories; unresolved origins are not redistributed.

MAX_CARBON_PROCESSES limits independent model/version workers. CPU/memory-aware
scheduling, live memory checks, a working-disk reserve and atomic workbook saves
remain enabled. MAX_EXCEL_SIZE_MB defaults to 50 decimal MB. A failed or
oversized export does not replace the existing workbook; a resource failure has
a nonzero exit status. No partial EMU workbook is presented as complete.
Allocation checks are safeguards, not a hard guarantee against OS-level OOM.

The current carbon_flow_run.json records the exact mapping hash, selected sources,
mapping version, solver checks, assumptions, emu_solved and emu_completed.
A workbook left by an earlier run is not proof of success in a new run.
Historical combined/Enumeration workbooks are preserved as historical results;
they are not the output contract of the EMU-only version. No new Summary,
Carbon_Flow or carbon_flow_summary.xlsx is generated.

## Optional AI review and custom cases

Ordinary scripts do not require AI. Optional companion skills prepare models,
review structures/correspondences or orchestrate stages:

- [Model setup](https://github.com/LGJee/carbon-tracing-model-setup-skill): existing model/culture configuration and step-1 exports.
- [Mapping review](https://github.com/LGJee/carbon-tracing-review-skill): review missing/flagged mappings and produce _AIGenerated copies.
- [Enzymatic mapping](https://github.com/LGJee/enzymatic-atom-mapping-skill): evidence-supported reaction atom correspondences.
- [Workflow](https://github.com/LGJee/carbon-tracing-workflow-skill): project deployment and stage orchestration.

Keep original predictions and Mapping_Selection unchanged during final mapping
review. Approved edits go into Final_Mapping of a separate _AIGenerated copy.
Keep source hashes and plan identities consistent. Rerun step 3 to analyze both
versions. Older skill instructions that request rounds, Enumeration, Carbon_Flow
or Summary must be adapted to the [current interface](docs/interfaces.md).

For an independent case, copy examples/step1_prepare_custom_model_inputs.py to
results/results_by_AI/<case>/step1_prepare_inputs.py and edit that case's model
blocks. Reuse project environments and data; do not copy the whole installation.
Use --plan results/results_by_AI/<case>/analysis_plan.json with steps 2 and 3.

## Checks, provenance and interpretation

For a bounded, offline workflow check using a deterministic test mapper:

```bash
.envs/localmapper/bin/python -m unittest scripts.test_analysis_plan scripts.test_step3_methods
```

These tests use temporary directories; they do not call real predictors.
The optional real installation example is python3 step2_atom_mapping.py --example
followed by python3 step3_carbon_flow.py --example after all environments/data
are installed. The mapper example can download model weights and is distinct
from the offline check.

See [detailed guide](docs/user-guide.md), [validation](docs/validation.md),
[third-party terms](THIRD_PARTY.md) and [project license](LICENSE).
A single FBA/pFBA solution is one representative flux distribution. More mappings,
broader reachability or fewer Unknown values do not establish biological accuracy.
