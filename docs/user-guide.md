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

Choose [installation with Codex](../INSTALL.md#option-a-install-with-codex) or
[manual installation](../INSTALL.md#option-b-install-manually).
The commands below describe the manual method.

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

Step 1 keeps four numbered strain blocks. Edit the relevant block to change its
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

## 3. Custom cases and mapping

Copy examples/step1_prepare_custom_model_inputs.py into
results/results_by_AI/<case>/step1_prepare_inputs.py. Edit explicit COBRApy
model blocks and run the copy. Keep models and companion annotations together
under data/<full_species_name>/. Generated plans and results stay beside the
custom script, sharing the project environments.

From the project root:

```bash
python3 step2_atom_mapping.py --plan results/results_by_AI/my_case/analysis_plan.json
python3 step3_carbon_flow.py --plan results/results_by_AI/my_case/analysis_plan.json
```

Step 2 retrieves/reuses structures, runs published mapping tools and selects
correspondences. Manual tables take priority. Final_Mapping includes every
reaction with coverage and review flags; partial positional evidence is retained
without pretending unknown regions are mapped. Hidden Mapping_Info stores
provenance and partial records; no Partial_Mapping sheet is produced.
See the [mapping rules](carbon-mapping-rules.md).

Mapper tasks retain CPU/memory-aware scheduling. Identical reaction inputs can
share predictions, and EnzymeMap uses shared database lookup and bounded rule
workers. A shared task scheduler coordinates these parallel jobs.

## 4. Analyze carbon sources

Step 3 runs the size-1 EMU-inspired solver only. It reads ordinary and optional
AI-reviewed mappings from the selected plan, with the same carbon sources for
each version. It never launches AI or atom mappers. Model/version tasks may
run concurrently within CPU/memory limits.

Relevant script settings are EMU_ALLOW_INORGANIC_CARBON_RECYCLING (True),
EMU_TOLERANCE (1e-8), MAX_CARBON_PROCESSES, CARBON_MEMORY_GB,
MIN_AVAILABLE_MEMORY_GB (2), MAX_EXCEL_SIZE_MB (50 decimal MB) and TRACING_PYTHON.
The default EMU setting allows internal inorganic-carbon
recycling; disabling it is an explicit modeling assumption, not a claim that
all biological CO2 escapes.

Each <model>_3_carbon_flow.xlsx contains:

- EMU_Metabolites: metabolite identity, carbon count, inflow/outflow, source-origin
  fractions and Unknown (plus external-source categories when present).
- EMU_BoundaryInputs: actual external inputs and origin categories.

The AI mapping version produces <model>_3_carbon_flow_AIGenerated.xlsx.
The solver retains validation/structure-scope metadata internally but omits
FluxResidual, Status, StructureType, KnownCarbonCount and ContributionScope in
the metabolite sheet. Partial-structure fractions refer only to explicit known
carbons, not to its unspecified group. See [scope](source-contributions.md).

Unknown is not redistributed among supplied sources.

## 5. Completion, safety and troubleshooting

Check the exit code and current carbon_flow_run.json, including mapping hashes,
emu_completed and output_file. emu_solved only confirms that numerical solving
finished. A file from an older run does not establish new success.

Memory checks reserve at least 2 GiB or 10% of Linux/cgroup memory, with estimates
before large allocations. Working-disk checks retain a reserve. Atomic exports
replace outputs only after complete saving and size validation. An oversized
or failed export returns an error and leaves previous output intact. These
safeguards cannot guarantee that OS/native-library allocations never run out
of memory. Free memory or reduce MAX_CARBON_PROCESSES when other applications
consume resources.

Routine output is short and English; --verbose shows diagnostics and tracebacks.
Mapper library messages are saved under the selected results/logs/ directory.
Loading and sparse solving have no reliable ETA. Use output/status files rather
than elapsed time to determine completion.

## 6. Bounded offline checks

Run from the project root:

```bash
.envs/localmapper/bin/python -m unittest scripts.test_analysis_plan scripts.test_step3_methods scripts.test_output_limits scripts.test_source_contributions scripts.test_task_scheduler
.envs/preparation/bin/python -m unittest scripts.test_model_preparation
```

The first set uses deterministic mapper fixtures, temporary workbooks and
hand-solvable networks, not full-model predictions. Tests needing EnzymeMap
dependencies run in its own environment. See [validation](validation.md).

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
only size-1 EMU-inspired balances. EMU reports model-dependent
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
