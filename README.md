# Atom-Mapped Carbon Tracing

Prepare pFBA-active metabolic networks, compare three atom-mapping tools and
trace possible carbon-source label propagation. Configure the analysis at the
start of **step1_prepare_carbon_tracing_inputs.m**. Steps 2 and 3 read its
automatically generated plan. Normal analyses need no Python edits or handwritten
JSON; advanced execution settings are at the top of the corresponding Python step.

```text
step1_prepare_carbon_tracing_inputs.m
    → <model>_1_pFBA_nonzero_reactions.xlsx + data/analysis_plan.json
step2_atom_mapping.py
    → <model>_2_atom_mapping.xlsx
step3_carbon_flow.py
    → <model>_3_carbon_flow.xlsx
```

## 1. Install MATLAB/COBRA and obtain the models

Install MATLAB and [COBRA Toolbox](https://github.com/opencobra/cobratoolbox).
Follow the [official installation and solver instructions](https://opencobra.github.io/cobratoolbox/stable/installation.html),
initialize COBRA and configure a compatible solver. GLPK is the reference solver.
COBRA is not copied into this project.

Clone this project in Linux/WSL:

```bash
git clone https://github.com/LGJee/atom-mapped-carbon-tracing.git
cd atom-mapped-carbon-tracing
```

Download models into `data/`:

| Species | Source | Local files |
|---|---|---|
| S. cerevisiae | [Yeast9 v9.1.1](https://github.com/SysBioChalmers/yeast-GEM/tree/v9.1.1/model) | `data/yeast9/yeast-GEM.mat`, `yeast-GEM.yml` |
| E. coli | [BiGG iML1515](https://bigg.ucsd.edu/models/iML1515) | `data/ecoli/iML1515.mat`, `iML1515.json` |
| Y. lipolytica | [iYali5 and companion structures](https://github.com/bioengtaltech/eciYali5-GEM/tree/3b97d38d92327a099f3b0f6f9cd92e8aeab78c10/model) | `data/yarrowia/iYali5.mat`, `eciYali5-GEM.yml` |

The default test setup includes all three species. The six named model blocks
run sequentially as written; there are no model run switches. Each WT/IUPD pair fixes both models to 50% of its WT maximum growth rate,
maximizes squalene formation and then minimizes total absolute flux. S. cerevisiae and E. coli
use native glucose bounds for growth optimization. The published iYali5 medium
has glucose closed; `cfg.yarrowiaGlucoseUptake = 1` explicitly selects glucose
and other carbon uptake is closed. WT uptake limits then double after fixing
growth: S. cerevisiae 1 -> 2, E. coli 10 -> 20 and Y. lipolytica 1 -> 2 mmol/gDW/h.
Other medium settings remain
species-specific. Native oxygen bounds and the source biomass objective are
inherited without repeated settings; the WT growth objective is checked before
optimization. IUPD retains the growth-phase glucose bounds and permits isoprenol
uptake up to 2.09 times the configured glucose limit: 2.09, 20.9 and 2.09 mmol/gDW/h,
respectively. Each IUPD block requires its matching WT block to run first. E. coli WT/IUPD both receive a heterologous SQS reaction and
a squalene output, since these are absent from iML1515.
Do not present those simulations as experimental strains or paper reproduction.

Keep MATLAB and Python working on the same project folder. Windows MATLAB may
open the project through its `\\wsl.localhost` directory. If using separate
computers, copy the input workbooks **and `data/analysis_plan.json`**, preserving
the project-relative paths. Windows usernames and WSL distribution names are
not written into the plan.

## 2. Configure and run step 1 in MATLAB

Edit paths and analysis settings at the top of
[step1_prepare_carbon_tracing_inputs.m](step1_prepare_carbon_tracing_inputs.m):
model/output paths, the Yarrowia glucose-medium limit and ordered carbon sources. Culture conditions, knockouts, introduced
reactions and objectives are written as direct COBRA calls in six clearly named
model blocks. Each block loads its original model independently, so IUPD never
inherits an already-fixed WT growth bound.

Each model block uses the short working variables `modelWT` or `modelIUPD`:

```matlab
modelIUPD = deleteModelGenes(modelIUPD,{'YML126C'});
modelIUPD = addMetabolite(modelIUPD,'isoprenol[c]', ...
    'metName','isoprenol','metFormula','C5H10O');
modelIUPD = addReaction(modelIUPD,'Ex_isopentenol_e_t','reactionFormula', ...
    'isoprenol[e] -> isoprenol[c]','lowerBound',0,'upperBound',1000);
```

The reaction and gene IDs above belong to the S. cerevisiae reference block.
Edit the direct COBRA calls in the appropriate block to change its conditions.
Each block saves its model and solution under species-specific names before
the next block reuses the working variables. Workbook export follows all six
blocks; input checks and shared export details are local functions below them.

```matlab
run('step1_prepare_carbon_tracing_inputs.m')
```

Default pairs are `Saccharomyces_cerevisiae_WT`/`Saccharomyces_cerevisiae_IUPD`, `Escherichia_coli_WT`/`Escherichia_coli_IUPD` and
`Yarrowia_lipolytica_WT`/`Yarrowia_lipolytica_IUPD`. Each pair uses species-specific genes and metabolite
IDs, with the same introduced two-step isoprenol-to-IPP utilization chemistry.
Details and limitations are in [model notes](docs/yeast9-example.md).

Each input workbook contains `Nonzero` and `Metabolites` sheets.
The reaction sheet has a short fixed name to respect Excel's 31-character limit.
Step 1 exports identities, formulas, charges and available database identifiers.
Step 2 reads model structures for analyses listed in `MODEL_YAML_FILES` at its top:
Yeast9 annotations for S. cerevisiae and the paired eciYali5 YAML for Yarrowia.
The supported COBRA/RAVEN YAML layouts supply scalar or first-list SMILES.
IDs, names, compartments and available formulas are checked before companion
structures are used. Manual overrides and existing input structures take priority;
PubChem supplies remaining missing structures. Yarrowia pFBA uses ordinary iYali5,
without enzyme constraints. The manual table supplies the introduced isoprenol
structures and explicit `[H+]` for iYali5 metabolites named `H+_p+1`.
Unconfigured analyses never receive another model's
annotations merely because IDs overlap. The generated plan lists exactly the models
from this export, so unrelated old files in `results/` are not analyzed.
It is written only after all inputs succeed; rerun step 1 after changing its settings.
The generated plan contains only model/input records, input hashes and the output
folder. Mapper selection, batch/atom limits, structure
completion and mapper environments are configured in step2_atom_mapping.py.
Propagation rounds, concurrency, memory and interpreter are configured in
step3_carbon_flow.py. Their defaults support normal execution without changes.
After updating an older checkout, rerun step 1 to regenerate the plan.
The plan and mapping diagnostics record each input SHA256. Step 3 rejects a
same-name mapping from a different step-1 export; rerun step 2 before tracing.

## 3. Install the Python environments (once)

In Linux/WSL x86_64 with Conda and Git:

```bash
bash scripts/setup_environments.sh
python3 scripts/download_enzymemap.py
```

Three isolated environments are checked and reused when compatible. Model files,
environment directories and large EnzymeMap data are not bundled. See the
[installation guide](docs/user-guide.md#3-install-python-environments-and-enzymemap-data).

## 4. Predict and select atom mappings

```bash
python3 step2_atom_mapping.py
```

The script reads step 1's plan, runs LocalMapper, RXNMapper and EnzymeMap with
CPU/memory-limited concurrency, then selects final correspondences. Each model
has one `<model>_2_atom_mapping.xlsx` containing the three original predictions,
`Mapping_Selection`, `Final_Mapping` and the enriched metabolite table.

Step 2 applies model-scoped overrides from `data/manual_smiles.csv` before online
structure completion. Name-selected structures mark affected reactions for
review. Inputs remain unchanged. Inspect `ReviewRequired=TRUE` and, if necessary,
correct `Final_Mapping` before step 3. Rerunning step 2 can replace those edits.
`AIGenerated` is stored only in `Final_Mapping` and defaults to FALSE.
The three-tool mapping stage does not perform AI review.

## 5. Propagate carbon labels

```bash
python3 step3_carbon_flow.py
```

The script reads the selected mappings and ordered carbon sources from the plan.
It never calls mappers or changes input/mapping workbooks. Each
`<model>_3_carbon_flow.xlsx` contains only `Carbon_Flow_R1` onward and `Summary`.
The default maximum is 10 rounds, with earlier stopping on convergence.

By default, labels reaching CO2 or bicarbonate are recorded as endpoints but
are not reused in downstream reactions. Carbonic acid and carbonate are treated
the same way in every compartment. Carboxylation reactions still transfer labels
from their other substrates. This excludes inorganic-carbon recycling from the
reachability analysis; it does not assume that all CO2 escapes or that refixation
is biologically negligible. To include recycling or trace an inorganic-carbon
source, set `ALLOW_INORGANIC_CARBON_RECYCLING = True` at the top of
`step3_carbon_flow.py`. Step 4 inherits this setting. Summary still counts labelled
inorganic-carbon endpoints, and its mapped-network denominator is unchanged.

For an optional subset or one-run override:

```bash
python3 step2_atom_mapping.py Escherichia_coli_WT Escherichia_coli_IUPD
python3 step3_carbon_flow.py Escherichia_coli_WT Escherichia_coli_IUPD --rounds 1
```

### Optional: propagate using AI-reviewed mapping copies

After external review, place `<model>_2_atom_mapping_AIGenerated.xlsx` copies
in the configured results directory, keeping the original model names and input
hashes. Successful replacements belong in `Final_Mapping.MappedReactionSMILES`,
with `AIGenerated=TRUE`; keep `Mapping_Selection.SelectedMappedReactionSMILES`
consistent. Preserve raw tool predictions and all unchanged reactions.

```bash
python3 step4_carbon_flow_ai.py
```

Step 4 uses available reviewed copies from the current plan and the same
propagation settings/implementation as step 3. It writes separate
`<model>_3_carbon_flow_AIGenerated.xlsx` files. It performs no AI review or mapper
prediction itself. A filename suffix selects a copy; it does not certify the
mapping's correctness. Missing reviewed copies never fall back to ordinary ones.
Model selection and `--rounds` work as in step 3.

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
Compare WT/IUPD within each species under explicitly stated conditions; differences
in model coverage and medium prevent interpreting cross-species percentages as
direct biological rankings.

## Custom inputs and installation check

For a different model, adapt a named COBRA block and its source list in step 1.
For Excel inputs prepared elsewhere, use
[examples/step1_use_existing_inputs.py](examples/step1_use_existing_inputs.py)
instead of the MATLAB exporter. Set propagation rounds with `ROUNDS` at the top
of step3_carbon_flow.py; the default is 10. See the
[detailed guide](docs/user-guide.md#4-configure-a-custom-network-and-carbon-sources)
for the input schema and mixed-source configuration.

The small installation example does not require MATLAB or an analysis plan:

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
