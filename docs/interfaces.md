# Calling the workflow from automation

The public interface is the existing Python stage scripts plus their generated
files. No service, manual JSON configuration or separate project installation
per analysis is required. Stage scripts do not call an LLM or a skill.

## Prepare an analysis

Keep the study's `step1_prepare_carbon_tracing_inputs.py` unchanged when preparing
a different analysis. Copy `examples/step1_prepare_custom_model_inputs.py` into
`results/results_by_AI/<name>/` and edit its explicit COBRApy blocks. Each block loads a model,
lists source names, sets medium/modifications and objectives, runs pFBA and
registers its model/solution/sources. Use existing model data under `data/<model>/`.
Choosing and downloading an appropriate model is a caller responsibility; the
stage scripts do not select a model from a species name. A same-folder
`source.json`, when present, is recorded by path and SHA256 for provenance.

```bash
python3 results/results_by_AI/my_analysis/step1_prepare_inputs.py
python3 step2_atom_mapping.py --plan results/results_by_AI/my_analysis/analysis_plan.json
python3 step3_carbon_flow.py --plan results/results_by_AI/my_analysis/analysis_plan.json
```

The example discovers the project root and places outputs beside its generated
plan. Step 1 also accepts `--plan` and `--output-dir`. With no arguments, the
study script still generates the four-model reference plan in
`data/analysis_plan.json` and writes to `results/`.

Steps 2 and 3 accept `--plan` with a project-relative or absolute generated plan path,
optional positional model names, and their existing stage-specific overrides.
`--rounds` applies only to propagation. All relative paths stored inside plans
are relative to the project root, not the plan's folder or current directory.
`--example` is an installation check and cannot be combined with `--plan`.
All analyses share the existing project environments and stage execution settings.

## Generated inputs and mappings

Step 1 writes one `<model>_1_pFBA_nonzero_reactions.xlsx` with `Nonzero` and
`Metabolites` sheets per model. It writes the plan only after all exports succeed.
Plan version 3 contains `models` and `output_dir`; each model records
`input_file`, ordered `carbon_sources` and `input_sha256`. When available it also
records `source_model` and `preparation_script` paths/hashes. Preparation metadata
records solver, objectives, exchange bounds and actual fluxes, knocked-out genes
and added reactions. The plan is generated data, not a second configuration file
for users to edit. Failed preparation invalidates only its selected plan.

Current scenarios are static simultaneous supplies. Opening an uptake bound
permits uptake; it does not force it. `preparation.source_uptake` records actual
net imports. A supplied source that is unused by pFBA is retained with a warning;
the nonzero reaction sheet stays unchanged. Source-list order determines label
identifiers and computational introduction order, not biological feeding time.

Step 2 writes `<model>_2_atom_mapping.xlsx`: raw LocalMapper, RXNMapper and
EnzymeMap predictions, `Mapping_Selection`, `Final_Mapping`, enriched metabolites
and provenance. Every reaction is retained in Final_Mapping, with MappingCoverage
Complete/Partial/Unresolved/Boundary (legacy NoCarbon records remain readable). Partial mapper carbon views use
MappedReactionSMILES; rule-derived partials include reduced unit mapped SMILES and hidden Mapping_Info records, whose
SourceSMILES/ProductSMILES and CarbonAtomMap indices are verified against the
stored input and current rule version before propagation. Input workbooks remain unchanged. Missing structures and failed
predictions retain the existing recorded statuses. Confidence and agreement do
not certify biochemical correctness. See the user guide for selection rules.

## AI review handoff

The intended future AI deployment/analysis orchestrator includes this review by
default after step 2. Review is external to the fixed stage scripts:

1. Pass flagged step-2 workbooks to `carbon-tracing-review-skill`; it delegates
   flagged and blank-mapping reactions to `enzymatic-atom-mapping`; explicit
   NoCarbon/Boundary coverage does not require mapping merely because its SMILES is empty.
2. Save `<model>_2_atom_mapping_AIGenerated.xlsx` in the same plan's output folder.
   Update only `Final_Mapping` for successful validated replacements, including
   `MappedReactionSMILES` and `AIGenerated=TRUE`. Preserve raw predictions,
   `Mapping_Selection`, identifiers and input provenance. Failed review retains
   the original mapping and flags, except scientifically rejected/unresolved carbon
   rules: mark their Final_Mapping.MappingCoverage Unresolved and clear mapped
   SMILES so step 3 does not use them. Confirmed partial rules stay Partial without
   setting AIGenerated. Full accepted AI replacements must be Complete. Preserve
   saved partial records in Mapping_Info as audit context. Keep review reports separately.
3. Run step 3 to propagate both available versions independently:

```bash
python3 step3_carbon_flow.py --plan results/results_by_AI/my_analysis/analysis_plan.json
```

Step 3 discovers normal and reviewed copies for each model in the selected
plan. For either method, either version may be present alone and both are processed
when available.
It retains provenance checks for each copy and performs no review. The suffix
selects a file, not a quality grade.
This revision does not create or implement the future orchestrator skill.

## Propagation status

Step 3 writes `<model>_3_carbon_flow.xlsx` for normal mappings and
`<model>_3_carbon_flow_AIGenerated.xlsx`. Each version combines Enumeration and
size-1 EMU-inspired balances when both methods are selected. `Summary` contains
all completed Enumeration rounds; `EXPORT_DETAIL_ROUNDS=3` exports only the
first three rounds of details in one `Carbon_Flow` sheet by default, with
`Carbon_Flow_Round` as its first column. `EMU_Metabolites` and
`EMU_BoundaryInputs` contain the independent source-balance result. EMU-only
runs have no Summary and are omitted from Summary collection. Only the last saved Summary row
contains `<source>_ExclusiveSourceMetabolites`: distinct mapped metabolite IDs whose
recorded states all carry labels from that one source. Separate pure-source states
from different sources disqualify a metabolite. The optional step-3 setting
`EXPORT_EXCLUSIVE_SOURCE_METABOLITES` (default `False`) adds an
`Exclusive_Source_Metabolites` sheet for that round; it does not change propagation.
`carbon_flow_summary.xlsx` copies each saved
Enumeration workbook's Summary to its own sheet in the same output directory,
including AI-reviewed versions. Existing files from earlier runs are included;
use the manifest below to determine which outputs belong to the current run.
Sheet names retain `_AI`, shorten long model names and resolve collisions numerically.

Each batch writes `carbon_flow_run.json`. Version 1 records the finish
timestamp, batch `status` (`finished` or `error`) and one entry per mapping
workbook. Normal and AI entries may share a model name; `mapping_file` and
`output_file` identify the version unambiguously:

| Field | Meaning |
|---|---|
| `model`, `carbon_sources` | Analysis identity and ordered sources |
| `mapping_file`, `mapping_sha256` | Actual selected mapping workbook |
| `requested_rounds` | Requested upper limit |
| `completed_rounds` | Complete rounds saved in this invocation |
| `output_file` | Current saved output, or null when this run saved none |
| `stop_reason` | `converged`, `round_limit`, `file_size_limit`, `row_limit`, `memory_limit`, `error`, `batch_interrupted` or `not_completed` |

Check the process exit status first. Invalid inputs can fail before a new report
exists; do not treat an old report or workbook as a new success. Memory/errors
raise a nonzero exit status and retain saved rounds. Other interrupted workers
may have only their saved output/round count recovered in an error report.
Size and row-limit stops return normally, retaining only complete rounds.
`output_batch_rows` (default 50,000) controls staged row batches; `max_excel_bytes`
records the configured size threshold (default 50,000,000). At batch/layer
checkpoints the actual candidate XLSX is measured. A candidate reaching the
threshold or exceeding the combined `Carbon_Flow` sheet's row limit is discarded, with
`file_size_limit` or `row_limit`; there are no pagination sheets. Only a complete
candidate below the threshold atomically replaces the previous output.
`completed_rounds` and `output_file` describe this invocation's published rounds,
not the rejected candidate. For combined runs, zero completed Enumeration rounds means no new output; the
previous workbook is retained. EMU-only runs can save independently with zero
Enumeration rounds. emu_solved distinguishes a solved balance from emu_completed
(successfully saved tables). Reruns start at round 1, not a checkpoint.
The `emu` manifest field records numerical residual, unresolved positions,
source-balance assumptions, independent recycling setting and mapping version.

Enumeration describes possible carbon-origin label reachability, not measured
isotope abundance, physical time or quantitative carbon allocation. EMU gives
model-dependent mean source fractions, not full isotopomer distributions. Existing
propagation rules, inorganic-carbon settings and Summary denominators are unchanged.

`--output-dir <folder>` optionally changes the propagation destination without
changing the selected plan's mapping locations. Step 3 selects Enumeration, EMU, or both through ANALYSIS_METHODS; method settings stay separate.


## Known-structure mapper inputs

Step 2 retains original equations and signed fluxes. When some saved metabolites
have no SMILES, it may pass only the known-structure participants to the mappers.
At least one retained participant on each side is required. Invalid structures,
unknown metabolite names, noninteger retained coefficients and atom limits still
block that input. No identity or biochemical role is inferred for omitted pools.
Mapping_Info records MappingInputEquation in the active direction and
OmittedComponentsJSON (side, original coefficient, name and reason). Predictions
on reduced inputs have Partial coverage and ReviewRequired=TRUE. Stage 3 uses only
retained carbon correspondences and keeps full original network turnover for EMU;
unresolved positions remain Unknown. A reduced map cannot be relabeled Complete.
AI review is optional and assesses applicability and the omitted roles when used.
