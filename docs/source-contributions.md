# Size-1 EMU-inspired carbon-source contributions (experimental)

Step 3 includes this analysis by default alongside Enumeration. It reuses the step-1 Flux values and step-2 selected atom
mappings. It computes the mean fraction of each carbon position attributable to
each input source. It implements single-carbon origin balances inspired by the
size-1 EMU balances, not a full EMU isotopomer simulator or 13C-MFA fitter.

Run in Linux/WSL from the project root:

```bash
python3 step3_carbon_flow.py --method EMU \
  --plan results/results_by_AI/my_analysis/analysis_plan.json
```

The legacy `analyze_source_contributions.py` command is a thin alias for this
same runner with EMU selected by default. Put positional model names before
`--method`, or omit it to use the script settings. Both methods run by default
from `step3_carbon_flow.py`; EMU runs first and ignores Enumeration rounds.

Optional positional analysis names select models from that plan. `--output-dir`
can isolate trial outputs. Without it, new files are written beside the existing
case outputs. The entry uses the existing localmapper environment and SciPy;
no new dependency or mapper execution is needed. Available `*_2_atom_mapping.xlsx` and `*_2_atom_mapping_AIGenerated.xlsx` files are read separately. Either can be present alone; AIReviewed/AIGenerated flags are accepted as mapping provenance. Origin solving, statistics and exports are deterministic code and make no LLM/mapper calls. Upstream ordinary mappers may themselves use machine learning; that mapping provenance is separate from numerical origin calculation. Earlier input and mapping
files are never changed. Re-running replaces only this analysis's own outputs.

## Calculation

A node is `(metabolite, SMILES fragment index, atom index)`. For each node i and
source s, let p[i,s] be its origin fraction. Selected mappings connect precursor
nodes j to product nodes i. Molecular copies and integer stoichiometry contribute
separate transfers. Negative Flux reverses the reaction before applying its
already direction-validated mapping. Let w[i,j] be the summed carbon transfer
rate, q[i] the total production rate of that position, and b[i,s] the prescribed
external or unknown input rate. Solve the sparse linear system:

```text
q[i] * p[i,s] - sum_j(w[i,j] * p[j,s]) = b[i,s]
```

Fluxes have a common unit, normally mmol/gDW/h. Reaction bounds are not Flux.
Every production and consumption term contributes to metabolite balance checks,
including sinks and reactions lacking mappings. A molecular imbalance larger
than `1e-6 + 1e-7 * max(inflow, outflow)` makes that pool Unknown. A source is
prescribed only on a boundary inflow, never clamped throughout an internal pool.
Other carbon-containing imports receive separate `External:<name>` categories.
Internal CO2 preserves its ancestry and may be refixed where mappings allow it.
The Enumeration-only `ALLOW_INORGANIC_CARBON_RECYCLING` switch does not change
this source-balance assumption. Boundary imports and exports use their actual
model fluxes; no additional gas-exchange rates are invented.
This differs from step 3's default inorganic-carbon termination policy.

Missing mappings send affected product carbons to Unknown. Graph-disconnected
cycles without any origin anchor also become Unknown. Missing/invalid structures
and structure/formula carbon-count mismatches have no computed positions or
fractions; their result cells remain blank. Mapped transfers to these unavailable
positions are omitted from position solving; a transfer from an unavailable
precursor to a usable product position contributes Unknown. Internal diagnostics
record these cases rather than inventing positions or aborting valid downstream
calculations. Mapping validation and stoichiometric checks remain in place.
No inferred biochemical mapping is
added for failed reactions or cofactor recycling.

Carbon positions equivalent under the saved structure's chirality-aware graph
symmetry are detected with RDKit atom ranks. An arbitrary choice between such
positions must not erase a source common to all possible precursor positions.
The solver starts conservatively with Unknown on ambiguous transfers. It restores
a transfer only when every candidate precursor position has the same fully known
single source (within numerical tolerance), then resolves the balances again.
This propagates verified origins through successive symmetric intermediates;
it never uses metabolite names or pathway-specific overrides.

If candidate positions have different origins, incomplete origins, or even the
same mixed-source fraction vector, this conservative rule does not assign an
arbitrary correspondence or assume random scrambling. Those transfers remain
Unknown. Full stereospecific/symmetry modeling and tighter bounds for unresolved
mixed-origin cases remain outside this implementation. Wrong saved mappings can
still give wrong results; this check does not establish biochemical correctness.

Probabilities, row sums and scaled equation residuals are checked to `1e-8`.
Unknown is never redistributed among known sources. Consequently a known fraction
in a partial result is resolved contribution, not a complete point estimate;
additional unknown carbon may also belong to that source. Metabolite means are
arithmetic means over its carbon positions. They are not flux-weighted averages
across different metabolites or fractions of molecules with a particular origin.

## Output

Each original mapping produces `<model>_3_carbon_flow_EMU.xlsx`; each reviewed
copy produces `<model>_3_carbon_flow_EMU_AIGenerated.xlsx`. Both contain:

| Sheet | Contents |
|---|---|
| MetaboliteContributions | Mean fractions, Unknown, carbon count and inflow/outflow |
| BoundaryInputs | Actual carbon imports and assigned origin categories |

The exported MetaboliteContributions table omits FluxResidual and Status. Flux
balance, per-position equations, structure availability and completeness are still
checked internally. The Excel contains no Summary, AtomContributions,
ReactionDiagnostics or Analysis_Info sheets. The existing
`source_contributions_run.json` retains hashes, assumptions, balance residuals,
metabolite-status counts and resolved/unresolved symmetry-transfer counts.

All exported source fractions use 0..1. Per-position fractions are calculated
internally using saved SMILES indices, not biochemical C1/C2 assignments.
The run manifest records the method, internal CO2 refixation policy, size-1 scope
and source-balance assumptions. No aggregate Summary workbook is generated.

`source_contributions_run.json` records only outputs successfully saved in that
invocation, mapping/input hashes, mapping_version (Original/AI) and finish/error
status. Each version uses its own Final_Mapping and retains step-1 hash,
embedded Flux, direction and structure validation. Outputs are atomically
replaced and capped by step 3 `MAX_EXCEL_SIZE_MB` (default 50 MB). A size or
memory error preserves any prior output and returns a nonzero process status.
EMU saves have no round or convergence interpretation. Prior outputs may remain after failure: check the
current invocation's process exit status and manifest, not file existence alone.

Assumptions are metabolic and origin steady state, fully mixed pools, fixed
net-direction fluxes and no isotope-dependent reaction preference. Unknown gross
forward/reverse exchange is not recovered from net Flux. A single pFBA solution
is one model scenario; it does not establish experimental contributions, flux
robustness, molecular source combinations, source essentiality or orthogonality.

## Validation

```bash
.envs/localmapper/bin/python -m unittest scripts.test_source_contributions
```

Hand-solvable cases cover weighted merging and per-position multiple origins,
condensation/cleavage, cycles, negative Flux, repeated molecular copies, missing
mappings, symmetry uncertainty, external imports, internal source-pool mixing,
invalid flux balance, carbon-count mismatches, scientific-notation coefficients,
input preservation, independent original/AI exports, AI-only inputs, and rejection
of mismatched input hashes or altered embedded Flux.

Reference: Antoniewicz MR, Kelleher JK, Stephanopoulos G. (2007).
[Elementary metabolite units (EMU): a novel framework for modeling isotopic distributions](https://doi.org/10.1016/j.ymben.2006.09.001).
