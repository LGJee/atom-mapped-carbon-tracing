# Atom-Mapped Carbon Tracing

Trace carbon-source atom origins through pFBA-active metabolic networks using
LocalMapper, RXNMapper and EnzymeMap, consensus selection, and multi-round label
propagation. Configure one or more carbon sources for single-source or
mixed-carbon cultures. Each analysis produces one plain Excel workbook.
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
nonzero reaction fluxes and metabolite structures. The Python workflow accepts
any model that supplies the following schema; WT/IUPD names and Yeast9 reaction
IDs are not required.

| Sheet | Required columns |
|---|---|
| Reaction sheet (any name) | `ReactionID`, `ReactionName`, `Flux`, `Equation` |
| Metabolite sheet (any name) | `MetaboliteID`, `MetaboliteName`, `SMILES_Yeast9`, `SMILES_Manual` |

Keep these SMILES column names for compatibility even when using another model.
Either structure column may be blank; `SMILES_Manual` takes priority when supplied.
Names must include compartments, such as `D-glucose[c]`, and match the names in
reaction equations. Equations use `=>`, `->` or `<=>`, terms separated by ` + `,
and positive integer stoichiometry, e.g. `A[c] + 2 B[c] => C[c]`.
Flux sign determines the active direction. Fractional reactions and missing
structures are recorded as unmappable.

For MATLAB input generation, [pFBA_network.m](pFBA_network.m) is a Yeast9
reference script. Adapt the model, medium bounds, growth constraint and objective
to your own study before exporting. pFBA supplies one representative optimal
flux distribution, not all possible pathways. See the optional
[Yeast9 WT/IUPD example](docs/yeast9-example.md) for model downloads, exact
reference conditions and output filenames.

Already have compatible input workbooks? Continue with Python setup below.

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

## 4. Configure your network and carbon sources

Edit [config.json](config.json), or copy it to `config.local.json` and pass
`--config config.local.json`. Relative paths resolve against the configuration
file, not the terminal directory. Configure input/sheet names, ordered carbon
sources, rounds, batch size, maximum reaction atoms, output and resource limits.
Maximum reaction atoms counts **both sides combined**.

Replace the study-specific `analyses` object in the supplied configuration
with your own analysis, keeping the other installation and resource settings.
Save the edited full configuration as `config.local.json`:

```json
{
  "analyses": {
    "MyCulture": {
      "input_file": "inputs/my_network.xlsx",
      "reaction_sheet": "Reactions",
      "metabolite_sheet": "Metabolites",
      "carbon_sources": [
        "D-glucose[e]"
      ]
    }
  }
}
```

For a mixed-carbon culture, list every source to trace, for example
`["D-glucose[e]", "ethanol[e]"]`. Each source must have a matching compartment-
qualified name and a valid structure in the input. Use the nonzero reaction
network computed for that culture's medium and constraints.

Source order controls first-round label introduction; it does not represent
physical uptake order or diauxic growth. The medium and uptake constraints
belong to the upstream metabolic model, not this carbon-source list.

Run atom mapping and selection for all configured analyses:

```bash
python3 run_carbon_tracing.py --config config.local.json --stage map
```

To process only WT, add `WT` immediately after `run_carbon_tracing.py` (if WT is
defined in your configuration). Do not add `--rounds` here: this stage performs
mapping and selection only, without carbon propagation.

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

Run atom mapping, selection and carbon propagation together:

```bash
python3 run_carbon_tracing.py --config config.local.json --rounds 10
```

This example runs up to 10 propagation rounds. Later rounds reuse the selected
atom correspondences without calling the prediction tools again.

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
[reference example](docs/yeast9-example.md) for its numbering details.
Reaction-local map IDs are bridges, not carbon-origin labels.

## 6. Inspect Excel results

Outputs are `<output_dir>/<analysis>_integrated_carbon_tracing.xlsx`.
The default output directory is `results/Integrated`.
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
