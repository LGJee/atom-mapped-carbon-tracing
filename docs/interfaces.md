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
optional positional model names, and their existing stage-specific overrides. All relative paths stored inside plans
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
the nonzero reaction sheet stays unchanged. Source-list order determines output column order, not biological feeding time.

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

## Carbon-source analysis and completion

Step 3 reads Final_Mapping through the current saved-structure/partial-record
validation. It solves fixed-flux size-1 EMU-inspired origin balances only.
The public command is python3 step3_carbon_flow.py [model names] --plan <plan>;
--output-dir can redirect new results while mapping inputs stay in the plan's
output directory. There is no --rounds or analysis-method selection.

The internal entry is run_carbon_flow(jobs, output_dir, max_processes=None,
memory_gb=None). Jobs specify carbon_sources and mapping_file or mapping_files,
with optional EMU numerical/recycling/resource settings. run_analyses and
finish_analysis also no longer take a rounds argument.

Ordinary and AI mappings produce separate <model>_3_carbon_flow.xlsx and
<model>_3_carbon_flow_AIGenerated.xlsx. Each has only EMU_Metabolites and
EMU_BoundaryInputs. There is no Summary or summary collector. Metabolite exports
omit FluxResidual, Status, StructureType, KnownCarbonCount and ContributionScope;
those fields remain internal. Source fractions and Unknown are unchanged.
Partial-structure averages cover explicit known carbon positions only, with
unknown total CarbonCount left blank. Unspecified polymer totals are unavailable.

The batch writes carbon_flow_run.json (version 1, method EMU). Each model/version
records model, mapping_file, mapping_sha256, carbon_sources, output_file,
emu_solved, emu_completed and stop_reason. A successful result has
emu_completed=true, a new output_file and stop_reason=emu_solved. The emu object
records assumptions, tolerance, residual, structure limitations and recycling.
Round counters are no longer part of the active interface.

Check process exit status and the report timestamp/input hashes. A failed read
can happen before a report is written. Resource/export errors return nonzero
and preserve prior workbooks; emu_solved alone does not mean output was saved.
Stop reasons include memory_limit, disk_space_limit, file_size_limit, error and
not_completed. In a failed batch, another worker may have saved a complete file;
its completion is recovered from the changed atomic output. Historical files
are not current completion records. Concurrent runs should use separate output
directories. Existing historical combined results remain untouched.

Older external skills requesting enumeration/rounds/Summary must update their
calls and completion checks to this contract. Normal scripts call no LLM.

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
