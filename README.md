# Atom-Mapped Carbon Tracing

Trace carbon-source atom origins through pFBA-active metabolic networks using
LocalMapper, RXNMapper and EnzymeMap, consensus selection, and multi-round label
propagation. Each strain produces one plain Excel workbook. This release does
**not** include LLM-assisted mapping or automatic biochemical corrections.

```text
MATLAB + COBRA + Yeast9 → pFBA Excel
                              ↓
            three independent atom-mapping tools
                              ↓
                  consensus / Final_Mapping
                              ↓
              carbon propagation → Excel results
```

## 1. Install COBRA Toolbox and obtain Yeast9

Install MATLAB and [COBRA Toolbox](https://github.com/opencobra/cobratoolbox),
initialize COBRA and configure a compatible LP solver according to the
[official installation guide](https://opencobra.github.io/cobratoolbox/stable/installation.html).
COBRA is not redistributed here.

Clone this project:

```bash
git clone https://github.com/LGJee/atom-mapped-carbon-tracing.git
cd atom-mapped-carbon-tracing
```

Obtain `model/yeast-GEM.mat` and `model/yeast-GEM.yml` from
[yeast-GEM v9.1.1](https://github.com/SysBioChalmers/yeast-GEM/tree/v9.1.1),
commit `2d594ae1c4a2d550ccef120d96a58c7bbf586255`.
Place both files in this project's `model/`, or configure their folder below.
The files must be from the same version; the YAML reader targets this version's
metabolite ID and SMILES annotation layout.

## 2. Generate input workbooks in MATLAB

Edit the settings at the start of [pFBA_network.m](pFBA_network.m), then run it
after initializing COBRA. Set `cfg.cobraFolder` only if initialization should
be performed from an explicit installation folder. The reference solver is
`glpk`; choose an installed compatible solver when necessary.

The script preserves these **study-specific**, not universal, conditions:

| Setting | Reference value |
|---|---|
| Model / growth objective | Yeast9 9.1.1 / `r_2111` |
| Glucose lower bound during growth optimization | -1 mmol/gDW/h |
| Oxygen lower bound | -1000 mmol/gDW/h |
| IUPD | Delete `YML126C` (ERG13); add isoprenol utilization reactions |
| IUPD isoprenol lower bound | -2.09 mmol/gDW/h |
| Fixed growth | 50% of each strain's initial maximum |
| WT glucose lower bound after fixing growth | -2 mmol/gDW/h |
| Production objective / nonzero tolerance | Squalene `r_1012` / `abs(flux)>1e-8` |

Remaining medium bounds are inherited from the model. Negative exchange flux
denotes uptake; a lower bound permits but does not force uptake.
`optimizeCbModel(model,'max','one')` preserves the optimal production objective,
then minimizes the sum of absolute fluxes. This is one representative pFBA
solution, not every possible pathway or optimal solution.

The script writes two workbooks:

| Workbook | Reaction sheet | Metabolite sheet |
|---|---|---|
| `results/WT_pFBA_nonzero_reactions.xlsx` | `WT_nonzero` | `Metabolites` |
| `results/IUPD_pFBA_nonzero_reactions.xlsx` | `IUPD_nonzero` | `Metabolites` |

Reaction columns: `ReactionID`, `ReactionName`, `Flux`, `Equation`.
Metabolite columns: `MetaboliteID`, `MetaboliteName`, `SMILES_Yeast9`, `SMILES_Manual`.
Names include compartments, such as `D-glucose[c]`; manual SMILES take priority.
The reference supplies manual structures for the added isoprenol metabolites
and the specified FBP(4−) form. No PubChem lookup or Windows Excel COM is used.
Running this MATLAB export again replaces its previous input exports.

Other models may supply the same schema. Python does not require WT/IUPD
reaction IDs. Equations use `=>`, `->` or `<=>`, terms separated by ` + `,
and positive integer stoichiometry, e.g. `A[c] + 2 B[c] => C[c]`.
Fractional reactions and missing structures are recorded as unmappable.

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

Three isolated environments are created in `.envs/`. Existing environments at
those paths are checked and reused, never deleted. Checks include actual
LocalMapper/RXNMapper calls and actual rdchiral rule execution for EnzymeMap's
dependencies. The full EnzymeMap adapter is exercised by the example after
downloading its data. An incompatible environment fails the check.

The LocalMapper/tracing stack pins Python 3.6.13, RDKit 2021.9.4, pandas 1.1.5
and openpyxl 3.0.10 to preserve study serialization. Python 3.6 is end-of-life;
this is a research reproduction environment, not an internet-facing service.
All specifications are in [environments/](environments/).

For a custom prefix, run `bash scripts/setup_environments.sh /path/to/envs`
and set `environment_dir` accordingly in the JSON configuration. Existing
compatible environments elsewhere can be referenced by an optional
`python_overrides` object containing Python paths for `LocalMapper`,
`RXNMapper` and `EnzymeMap`. LocalMapper's interpreter also runs carbon tracing.

The downloader obtains official **EnzymeMap v2** CSV and rules from commit
`c91578856f11312d7941b290f8b773eb41b54973`. The compressed CSV is about 13 MB.
Before unpickling rules, SHA256 must match
`21cdebccdda5ea4308555b6dd8dfb27413950725ace66e7b8cb7ed3dba2c1e1a`.
Only load trusted official rule files. Data and model weights are not bundled.
On a timeout, rerun the downloader or obtain files from the
[official pinned data directory](https://github.com/hesther/enzymemap/tree/c91578856f11312d7941b290f8b773eb41b54973/data)
and configure their paths. Rebuilding the original BRENDA-derived database is
outside this workflow.

## 4. Configure and run atom mapping

Edit [config.json](config.json), or copy it to `config.local.json` and pass
`--config config.local.json`. Relative paths resolve against the configuration
file, not the terminal directory. Configure input/sheet names, ordered carbon
sources, rounds, batch size, maximum reaction atoms, output and resource limits.
Maximum reaction atoms counts **both sides combined**.

```bash
python3 run_carbon_tracing.py --stage map --rounds 1
python3 run_carbon_tracing.py WT --stage map
```

Pure identity-preserving compartment transport bypasses the prediction tools
and size limit. Boundary/exchange reactions require no two-sided mapping.
Other reactions are independently submitted to the three tools. All jobs finish
before selection. Concurrency is constrained by CPU affinity and available
memory; the per-process memory budget is an estimate, not a peak guarantee.
Lower concurrency, batch size or reaction-size limits when necessary.

EnzymeMap first searches matching database records, then uses official
single-step rules, optionally normalizing protonation and restoring the original
input structures. It must preserve supplied participants and multiplicity.
Aggregate reactions or missing structures can remain unavailable.

Predictions are reused only when inputs, structures, direction and adapter
settings match. Run with `--refresh-mappings` after changing model weights,
datasets or rule contents without changing their configuration paths.

### Consensus policy

Compare structurally aligned carbon correspondence graphs, ignoring arbitrary
IDs, SMILES order and equivalent symmetry. Require complete carbon coverage.
Select the largest agreeing group; break ties by EnzymeMap → LocalMapper with
a true template flag → RXNMapper → LocalMapper with a false flag.
Keep one complete original candidate; never splice candidates' atoms.
Transport uses `SelectedMapper=Identity` without tool voting.

`ReviewRequired=FALSE` means at least two tools support the selected carbon
mapping, or the reaction is pure identity transport. No usable mapping,
single-tool support or conflicting single votes require review. Other-heavy-atom
disagreement is retained but does not change this carbon-specific rule.
Native scores are not comparable correctness probabilities: LocalMapper has a
Boolean template flag, RXNMapper attention confidence, and database EnzymeMap
matches relative rule-frequency quality. Agreement does not prove mechanism.
Selected mappings requiring review still participate in propagation.
`AIGenerated` follows `ReviewRequired` and remains `FALSE` in this release.

## 5. Run carbon propagation

```bash
python3 run_carbon_tracing.py --stage trace --rounds 1
```

This verifies the stored inputs and `Final_Mapping`, then propagates without
calling prediction software. For both stages in one command:

```bash
python3 run_carbon_tracing.py --rounds 10
```

Each rerun starts from configured sources and replaces old carbon-flow and
Summary sheets. It does not resume a previous unfinished enumeration.

- Round one introduces sources in configured order: IUPD glucose first, then
  isoprenol with previous states retained. Each reaction is processed at most
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

Glucose initialization expects a six-membered ring form and uses biochemical
C1–C6. Isoprenol is ordered from the hydroxyl end along
the chain, then the terminal alkene carbon and methyl branch. Source ranges
follow list order: glucose then isoprenol gives C1–C6 and C7–C11; isoprenol alone
gives C1–C5. Other sources use carbon order in the supplied single-component
SMILES, explicitly reported rather than claimed as biochemical numbering.
Reaction-local map IDs are bridges, not these carbon-origin labels.

## 6. Inspect Excel results

Outputs are `results/Integrated/<strain>_integrated_carbon_tracing.xlsx`.
MATLAB inputs remain unchanged. Close results in Excel before replacement;
if rename is denied, a completed `.tmp.xlsx` remains. No timestamped intermediate
result directories are created.

| Sheet | Contents |
|---|---|
| `Mapping_Selection` | Tool scores/validation, support, selected correspondence and reason |
| `Final_Mapping` | One usable correspondence per reaction |
| Three tool sheets | Original predictions and native confidence |
| Original reaction sheet / `Metabolites` | Copied input records |
| `Carbon_Flow_Rn` | Concrete combinations and labelled products in round n |
| `Summary` | One row per round, with counts, increments and percentages by source |
| `Mapping_Info` (hidden) | Adapter configuration, status and template diagnostics |

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
python3 run_carbon_tracing.py --config examples/config.json --rounds 1
```

The supplied example uses synthetic nonzero fluxes for acetaldehyde transport
and NADH-dependent reduction. It is not a Yeast9 pFBA result. The reduction is
actually submitted to all three tools. Expected labelled products are cytosolic
acetaldehyde and ethanol; identity transport is not three independent predictions.

```bash
.envs/localmapper/bin/python -m tests.test_carbon_tracing
.envs/localmapper/bin/python -m tests.test_mapper_backends
.envs/localmapper/bin/python -m tests.test_consensus_mapping
.envs/localmapper/bin/python -m tests.test_parallel_carbon
.envs/localmapper/bin/python -m tests.test_stages
.envs/enzymemap/bin/python -m tests.test_enzymemap_rules
```

Scheduling/failure tests use explicit mock adapters; they do not replace the
real-tool example. See [docs/validation.md](docs/validation.md) for observed
release checks and platform limits. Record the repository commit, configuration,
input hashes, tool versions and external data versions with manuscript analyses.

## Attribution

Authored code is under [MIT](LICENSE); dependencies retain their own terms in
[THIRD_PARTY.md](THIRD_PARTY.md). LocalMapper's pinned LICENSE is CC BY-NC-SA 4.0;
this repository's license does not remove its noncommercial restriction.
No upstream implementations, weights, COBRA or full EnzymeMap dataset are
redistributed. Follow citation instructions for
[COBRA](https://github.com/opencobra/cobratoolbox),
[Yeast9](https://github.com/SysBioChalmers/yeast-GEM),
[LocalMapper](https://github.com/snu-micc/LocalMapper),
[RXNMapper](https://github.com/rxn4chemistry/rxnmapper) and
[EnzymeMap](https://github.com/hesther/enzymemap).
