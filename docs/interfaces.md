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
study script still generates the original six-model plan in
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
and provenance. Input workbooks remain unchanged. Missing structures and failed
predictions retain the existing recorded statuses. Confidence and agreement do
not certify biochemical correctness. See the user guide for selection rules.

## AI review handoff

The intended future AI deployment/analysis orchestrator includes this review by
default after step 2. Review is external to the fixed stage scripts:

1. Pass flagged step-2 workbooks to `carbon-tracing-ai-review-skill`; it delegates
   `ReviewRequired=TRUE` reactions to `enzymatic-atom-mapping`.
2. Save `<model>_2_atom_mapping_AIGenerated.xlsx` in the same plan's output folder.
   Update only `Final_Mapping` for successful validated replacements, including
   `MappedReactionSMILES` and `AIGenerated=TRUE`. Preserve raw predictions,
   `Mapping_Selection`, identifiers and input provenance. Failed review retains
   the original mapping and flags. Keep review reports separately.
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

Step 3 writes `<model>_3_carbon_flow_Enumeration.xlsx` for normal mappings and
`<model>_3_carbon_flow_Enumeration_AIGenerated.xlsx`. Enumeration workbooks contain only saved
`Carbon_Flow_R*` and `Summary` tables. No aggregate Summary workbook is generated.

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
threshold or exceeding a single sheet's row limit is discarded, with
`file_size_limit` or `row_limit`; there are no pagination sheets. Only a complete
candidate below the threshold atomically replaces the previous output.
`completed_rounds` and `output_file` describe this invocation's published rounds,
not the rejected candidate. Zero completed rounds means null output even if an
older workbook exists. Reruns start at round 1, not a checkpoint.

These outputs describe possible carbon-origin label reachability, not measured
isotope abundance, physical time or quantitative carbon allocation. Existing
propagation rules, inorganic-carbon settings and Summary denominators are unchanged.

## EMU mean source contributions

Step 3 defaults to `ANALYSIS_METHODS = ["Enumeration", "EMU"]` and accepts
`--method Enumeration`, `--method EMU`, or both names. `--output-dir` changes
only the destination, not the mapping/input locations in the plan. The legacy
`analyze_source_contributions.py` entry delegates to the same EMU runner.

EMU runs first, independently of Enumeration. It requires a generated analysis
plan and the original step-1 workbook. It verifies its hash, mapping provenance
and the copied reaction/Flux table against that input for each original or AI
mapping version. It reads each version's Final_Mapping independently and accepts
AIReviewed/AIGenerated flags. An AI-only model is supported with the same step-1
requirements. A failed method does not cancel the other selected method;
the overall command still returns nonzero for an error.

`<model>_3_carbon_flow_EMU.xlsx` and
`<model>_3_carbon_flow_EMU_AIGenerated.xlsx` contain only MetaboliteContributions
and BoundaryInputs. FluxResidual and Status are omitted from the exported
metabolite table; internal balance/completeness checks remain active. It implements single-carbon origin balances inspired by size-1
EMU, not full multi-carbon isotopomer simulation. It has no Round or ProductionID.
Unknown is retained without renormalizing known contributions.

`source_contributions_run.json` records the current successfully saved EMU
outputs, Original/AI mapping versions, input/mapping hashes, assumptions, balance residuals, internal status
counts, symmetry-recovery counts and error/finish status. EMU output uses the step-3
size limit for atomic complete-workbook saves. The Enumeration batch rollback,
row-limit and round settings do not apply to the EMU numerical solve.

Enumeration Summary includes Method and InorganicCarbonRecycling (default FALSE).
EMU records its method and internal CO2 refixation policy in
`source_contributions_run.json`; it has no Summary sheet. Existing result files
are not renamed or modified automatically. See [assumptions, fields and validation](source-contributions.md).
