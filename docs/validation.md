# Validation record: EMU-only local development version

Validated on 2026-10-24 (Asia/Shanghai), branch codex/major-rework.
The pre-removal local snapshot is a28ec6c4bc908f2939a7a0509f4dcf16e28b41b7.
Historical combined-method validation records remain accessible in that Git
commit; they describe an older interface, not current supported commands.
This records the initial local cleanup checks. See publication-validation.md for the subsequent publication checks.

## Scope

Step 3 now directly runs fixed-flux, single-carbon source balances. Enumeration,
rounds, state catalogs, provenance-row staging and enumeration workers were
removed. The original/AI workbook naming contract remains. New workbooks have
only EMU_Metabolites and EMU_BoundaryInputs; existing historical results are not
rewritten or removed. No Summary collector is introduced.

The AST of solve_source_contributions is identical to the saved version.
Its inorganic-carbon recognition helper was moved out of the deleted enumeration
module without changing the helper logic. No changes were made to step-1/step-2
scientific settings, models, data tables, environments or existing results.

## Checks actually performed

- 138 compatible tracing-environment unit tests passed, including source balances,
  partial structures/rules, original-versus-AI mapping separation, single/multiple
  sources, plan/hash validation, exported columns, serial/parallel equality,
  scheduler admission, initializer support and abrupt worker failure.
- 17 COBRApy preparation tests passed in the preparation environment. Test calls
  now explicitly pass an empty argument list rather than consuming unittest's
  command-line arguments. Production model-preparation code is unchanged.
- Deterministic mapper-backend, consensus/reuse, and SMILES-completion integration
  scripts passed. Fake mappers run in isolated temporary directories. Legacy
  concurrency assumptions were corrected to respect live memory admission.
- Native EnzymeMap single-step rule smoke test passed in its existing environment.
- step2_atom_mapping.py --help and step3_carbon_flow.py --help passed. The old
  --rounds argument is rejected; analysis plans retain their version-3 contract.
- Atomic export tests preserve prior bytes on size, disk or memory failure and
  distinguish emu_solved from emu_completed. Internal metadata stays unchanged
  while the requested five fields remain absent from EMU_Metabolites.
- All active Python source/test/example files compile. Current source was checked
  for removed enumeration imports/functions/options. Git diff checks were run.

Run commands from the project root with the appropriate interpreter. For example:

```bash
.envs/localmapper/bin/python -m unittest scripts.test_source_contributions scripts.test_step3_methods scripts.test_analysis_plan scripts.test_output_limits scripts.test_parallel_carbon scripts.test_task_scheduler
.envs/preparation/bin/python -m unittest scripts.test_model_preparation
.envs/enzymemap/bin/python -m scripts.test_enzymemap_rules
```

Detailed command/exit records are retained locally in .cache/emu_only/checks.json
and accompanying logs; these are temporary validation artifacts, not public
research results. Environment-specific tests cannot all run under one interpreter.

## Boundaries and external compatibility

No full-model recalculation, real three-predictor sweep, new environment install,
server stress test, GitHub push or external skill update was performed. Existing
ordinary/AI files remain historical until explicitly recalculated with this code.
The retained solver is size-1 EMU-inspired, not full multi-carbon isotopomer EMU.
Unknown is not redistributed and full-reaction turnover remains unchanged.

Companion skills were subsequently adapted separately to the EMU-only interface.
Their current handoff compatibility is recorded in publication-validation.md.
