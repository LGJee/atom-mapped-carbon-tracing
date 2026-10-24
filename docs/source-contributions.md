# Size-1 EMU-inspired source balances

Step 3 calculates carbon-source contributions at each carbon position using the
selected fixed flux distribution.
It is not a complete multi-carbon EMU/isotopomer simulator and does not fit fluxes
or experimental isotope measurements. It solves mean source-origin fractions
under metabolic/origin steady state, well-mixed pools, net-direction fluxes and
the selected positional atom mappings.

## Configuration and output

Step 3 runs this solver directly. EMU_ALLOW_INORGANIC_CARBON_RECYCLING
(default True) and EMU_TOLERANCE (1e-8) remain in step3_carbon_flow.py.
Shared model-level scheduling, memory/disk and workbook-size safeguards remain.

Each ordinary or AI-reviewed mapping produces one <model>_3_carbon_flow.xlsx or
<model>_3_carbon_flow_AIGenerated.xlsx. Each contains EMU_Metabolites (carbon
counts, inflow/outflow, source/Unknown fractions) and EMU_BoundaryInputs.
FluxResidual, Status, StructureType, KnownCarbonCount and ContributionScope
remain internal, not exported.
carbon_flow_run.json records solving, export completion, assumptions and numerical
checks.

## Partial mappings and Unknown

Full-reaction turnover is always calculated from the original equations and
fluxes. Saved CarbonRules pairs use abs(Flux) * TransferCoefficient for each
explicit carbon pair, including decimal amounts such as 55.3. Sparse mapper
carbon maps contribute only their shared positive map IDs. Only the remaining
unmapped product-position inflow is Unknown. Over-assigned donor or product
positions fail validation. Reduced mapped views must never replace the original stoichiometry in this solver.

Molecular formulas provide carbon counts and input checks, not positional
correspondences. Missing structures, missing maps, invalid metabolite balances,
and unanchored carbon cycles can produce Unknown. Unknown is not redistributed
or silently removed by renormalizing known fractions. External inputs not among
configured sources get separate External:<metabolite> categories.

The selected validated mapping defines positional representatives for symmetric
structures; graph symmetry alone does not erase a known origin or imply random
scrambling. A mapping prediction can still be chemically wrong. Carrier carbons
and attached groups must not be confused: acetyl-CoA has 21 CoA carbons and two
acetyl carbons, so the whole-molecule mean can be dominated by the carrier.
A terminal biomass drain does not supply downstream atoms, but missing mappings
of cofactor returns within a biomass reaction can contaminate downstream origin
balances. One pFBA solution and these selected maps do not establish experimental
source contributions, source essentiality or all possible flux distributions.

## Reference and tests

Antoniewicz et al. (2007), [Elementary metabolite units](https://doi.org/10.1016/j.ymben.2006.09.001).
The original framework uses atomic transitions; this implementation retains only
the single-carbon origin-balance scope.

```bash
.envs/localmapper/bin/python -m unittest scripts.test_source_contributions scripts.test_step3_methods
```

Output is atomically saved only after a complete solve and export. A rejected
export keeps previous output unchanged; emu_solved can be true while
emu_completed is false. The process returns a nonzero status on such a failure.
