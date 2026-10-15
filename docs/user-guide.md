# Detailed workflow and custom analyses

This page documents optional customization and technical details. For the default
reference workflow, follow the [quick start](../README.md).

Trace carbon-source atom origins through pFBA-active metabolic networks using
LocalMapper, RXNMapper and EnzymeMap, consensus selection, and multi-round label
propagation. Configure one or more carbon sources for single-source or
mixed-carbon cultures. Each analysis produces three plain Excel workbooks, one per stage.
This release does **not** include LLM-assisted mapping or automatic biochemical corrections.

```text
Metabolic model + MATLAB/COBRA → pFBA Excel
                              ↓
            three independent atom-mapping tools
                              ↓
                  consensus / Final_Mapping
                              ↓
              carbon propagation → Excel results
```

## 1. Install MATLAB and COBRA Toolbox

Install MATLAB and [COBRA Toolbox](https://github.com/opencobra/cobratoolbox),
initialize COBRA and configure a compatible LP solver according to the
[official installation guide](https://opencobra.github.io/cobratoolbox/stable/installation.html).
COBRA is not redistributed here.

Clone this project:

```bash
git clone https://github.com/LGJee/atom-mapped-carbon-tracing.git
cd atom-mapped-carbon-tracing
```

## 2. Prepare input workbooks

Provide one Excel workbook per metabolic network or culture condition, with
nonzero reaction fluxes and metabolite identities. Structures are optional at export. The Python workflow accepts
any model that supplies the following schema; WT/IUPD names and Yeast9 reaction
IDs are not required.

| Sheet | Required columns |
|---|---|
| Reaction sheet (any name) | `ReactionID`, `ReactionName`, `Flux`, `Equation` |
| Metabolite sheet (any name) | `MetaboliteID`, `MetaboliteName` |

No SMILES columns are required. If the model already has structures, optionally
export `SMILES_Model`; the legacy `SMILES_Yeast9` column remains supported.
`SMILES_Manual` takes priority when supplied. Optional `MolecularFormula`,
`Charge`, `InChIKey` and `PubChemCID` aid step-2 structure retrieval.
The supplied MATLAB exporter does not read SMILES annotations. Step 2 reads
model-native structures only for analysis names configured in `MODEL_YAML_FILES`.
The reader supports the supplied COBRA and RAVEN YAML layouts, taking a scalar
SMILES or the first entry in a SMILES list. Yeast9 supplies S. cerevisiae structures;
the companion eciYali5 YAML supplies Yarrowia structures while ordinary iYali5
provides its pFBA network. The manual table also identifies iYali5 protons
(`H+_p+1`) as `[H+]`; no generic lipid/polymer structure is invented.
Names, compartments and available formulas must agree
with each matched metabolite ID. Manual overrides and existing input structures
take precedence; PubChem supplies remaining missing structures. Local YAML/manual
sources still work when `COMPLETE_MISSING_SMILES=False` disables online lookup.
Add a custom analysis name to that dictionary when reusing compatible annotations.
Names must include compartments, such as `D-glucose[c]`, and match the names in
reaction equations. Equations use `=>`, `->` or `<=>`, terms separated by ` + `,
and positive integer stoichiometry, e.g. `A[c] + 2 B[c] => C[c]`.
Flux sign determines the active direction. Fractional reactions and missing
structures are recorded as unmappable.

[step1_prepare_carbon_tracing_inputs.m](../step1_prepare_carbon_tracing_inputs.m)
contains the model, medium, knockout, objective and carbon-source settings.
The default test setup exports three species' WT/IUPD pairs. Both models of each species fix growth at 50% of the WT maximum and maximize squalene
formation before minimizing total absolute flux. The E. coli comparison models
both include introduced SQS activity and a squalene output; medium settings
remain species-specific.
S. cerevisiae and E. coli WT growth optimization preserves native glucose bounds.
iYali5 ships with glucose uptake closed, so step 1 explicitly enables a glucose
medium using `cfg.yarrowiaGlucoseUptake` (default 1) and closes other carbon uptake.
After fixing half of the maximum growth rate, WT glucose uptake limits double:
S. cerevisiae 1 -> 2, E. coli 10 -> 20 and Y. lipolytica 1 -> 2 mmol/gDW/h.
The iYali5 glucose limit permits uptake rather than forcing it.
IUPD retains the growth-phase glucose bounds (limits 1, 10 and 1),
permits isoprenol uptake up to 2.09 times those limits (2.09, 20.9 and 2.09),
and uses the same fixed growth rate as its paired WT, not half its own maximum.
Native oxygen bounds are inherited (lower bound -1000 in all three reference
models). The loaded WT biomass objective is validated and used directly. Set a
different medium bound or objective explicitly when adapting another model.
IUPD uses its paired WT growth reference before selecting the squalene objective.
pFBA supplies one representative solution, not all possible pathways.

The MATLAB exporter also writes `data/analysis_plan.json` with portable paths,
ordered sources, the output folder and input hashes. Do not edit this generated file; change
step 1 and rerun it. Copy it with the Excel inputs when moving between computers.

## 3. Install Python environments and EnzymeMap data

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

Three isolated environments are created in `.envs/`. Existing environments at
those paths are checked and reused, never deleted. Checks include actual
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

## 4. Configure a custom network and carbon sources

For model/input preparation, edit step 1. Paths and ordered sources are at the top;
model-specific medium, knockout, added reactions and objectives appear directly
in the six named COBRA blocks. Each block uses `modelWT` or `modelIUPD`,
then saves a species-specific model and solution for the separate export section.
Validation and export helpers are at the end of the file. Model-specific source lists are
`cfg.carbonSources.<model>`. Set maximum propagation rounds with `ROUNDS` at
the top of step3_carbon_flow.py. Execution settings belong to the Python steps:
step 2 configures tools, batch/atom limits, mapping resources, structure completion
and mapper/data paths; step 3 configures rounds, propagation resources and its interpreter.
The generated plan does not override those execution settings. The atom limit
counts both reaction sides. Optional command-line overrides take precedence for
that run. Older generated plans are rejected with a request to rerun step 1.

To run fewer reference models, omit their model blocks, corresponding export
calls and `cfg.carbonSources` entries together. Keep the matching WT block
when retaining an IUPD block: its WT maximum supplies the shared growth reference.
Each model still loads its original MAT file; IUPD does not inherit WT production
bounds.

For Excel inputs prepared elsewhere, edit `INPUTS` and `OUTPUT_DIR` at the top
of [examples/step1_use_existing_inputs.py](../examples/step1_use_existing_inputs.py),
then run it instead of the MATLAB exporter:

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
Python execution settings do not require a new MATLAB export. If export fails, no usable plan
remains: fix the error and rerun step 1 before proceeding. `--example` bypasses
the plan for the supplied installation test.


## 5. Run carbon propagation

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
rather than silently giving no propagation. Step 4 inherits the same setting.
The default deliberately limits the scope of qualitative reachability; it is
not evidence that CO2 all escapes or that its biological contribution is small.

Each rerun starts from configured sources and replaces old carbon-flow and
Summary sheets. It does not resume a previous unfinished enumeration.

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
set AIGenerated to TRUE only for successful replacements, and synchronize
Mapping_Selection.SelectedMappedReactionSMILES. Append successful formerly
unmapped reactions to Final_Mapping using its existing schema. AI provenance
alone does not clear ReviewRequired or establish biochemical correctness.

```bash
python3 step4_carbon_flow_ai.py
python3 step4_carbon_flow_ai.py MyCulture --rounds 1
```

Place these copies in the plan's output directory (normally results/). Only
available reviewed models from that plan are processed; a requested missing
model or absence of reviewed copies causes an error. Step 4 retains the input
hash/direction/structure checks, uses step 3's source lists, rounds and resource
settings, and writes <model>_3_carbon_flow_AIGenerated.xlsx. Original mapping and
ordinary carbon-flow outputs remain unchanged. It runs no mappers or AI skill;
review and workbook updates happen separately before this optional entry.

## 6. Inspect Excel results

Step 2 writes `<OUTPUT_DIR>/<analysis>_2_atom_mapping.xlsx`; step 3 reads it
and writes `<OUTPUT_DIR>/<analysis>_3_carbon_flow.xlsx`. The default output
directory is `results`. Raw MATLAB inputs use `<model>_1_` filenames.
MATLAB inputs remain unchanged. Close results in Excel before replacement;
if rename is denied, a completed `.tmp.xlsx` remains. No timestamped intermediate
result directories are created.

The step 2 workbook contains copied input tables, three raw tool sheets,
`Mapping_Selection`, `Final_Mapping` and hidden `Mapping_Info` diagnostics.
It is unchanged by step 3. The step 3 workbook contains only `Summary` and
`Carbon_Flow_Rn` sheets (paginated if Excel's row limit is exceeded).
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

Results describe **qualitative carbon-origin reachability under the chosen pFBA
network and mapping assumptions**, not measured 13C abundance, quantitative
carbon allocation, isotope kinetics or one molecule's physical history.
Multi-round combinations can grow rapidly even without rerunning mappers.
Uncertain mappings, incomplete reactions and polymer-as-monomer structures can
limit position-specific interpretation. Coverage concerns the mapped subnetwork.

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
.envs/localmapper/bin/python -m scripts.test_carbon_tracing
.envs/localmapper/bin/python -m scripts.test_mapper_backends
.envs/localmapper/bin/python -m scripts.test_consensus_mapping
.envs/localmapper/bin/python -m scripts.test_parallel_carbon
.envs/localmapper/bin/python -m scripts.test_stages
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
No upstream implementations, weights, COBRA or full EnzymeMap dataset are
redistributed. Follow citation instructions for
[COBRA](https://github.com/opencobra/cobratoolbox),
[Yeast9](https://github.com/SysBioChalmers/yeast-GEM),
[LocalMapper](https://github.com/snu-micc/LocalMapper),
[RXNMapper](https://github.com/rxn4chemistry/rxnmapper) and
[EnzymeMap](https://github.com/hesther/enzymemap).

The reference MAT files store the model as `model` (Yeast9 and iYali5) or
`iML1515` (BiGG E. coli). The six direct MATLAB blocks show these load operations.
When replacing a MAT model, adapt `data.<variable>` as well as the model-specific
exchange, growth and production reaction IDs. Model-building operations are
intentionally explicit; shared functions handle Excel export and provenance.
