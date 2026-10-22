# Size-1 EMU-inspired source balances

Step 3 restores the fixed-flux single-carbon origin solver from commit fefa1f9.
It is not a complete multi-carbon EMU/isotopomer simulator and does not fit fluxes
or experimental isotope measurements. It solves mean source-origin fractions
under metabolic/origin steady state, well-mixed pools, net-direction fluxes and
the selected positional atom mappings.

## Configuration and output

At the top of step3_carbon_flow.py, ANALYSIS_METHODS selects Enumeration, EMU or
both. EMU_ALLOW_INORGANIC_CARBON_RECYCLING (default True) and EMU_TOLERANCE
(default 1e-8) are separate from Enumeration rounds and recycling settings.
Shared memory/concurrency and workbook-size safeguards apply to each output.

Both methods read the same Final_Mapping version. Each model/version creates one
<model>_3_carbon_flow.xlsx or <model>_3_carbon_flow_AIGenerated.xlsx:

- Summary and Carbon_Flow belong to Enumeration. By default only R1-R3 detail rows are
  exported by default; Summary continues across all completed rounds.
- EMU_Metabolites contains carbon counts, inflow/outflow and source/Unknown
  fractions. FluxResidual and Status remain internal validation fields.
- EMU_BoundaryInputs records actual carbon imports and their origin categories.

There is no EMU Summary, atom-position/diagnostic sheet, or separate method file.
carbon_flow_run.json records EMU completion, assumptions and numerical checks.
EMU-only workbooks are excluded from carbon_flow_summary.xlsx.

## Partial mappings and Unknown

Full-reaction turnover is always calculated from the original equations and
fluxes. Saved CarbonRules pairs use abs(Flux) * TransferCoefficient for each
explicit carbon pair, including decimal amounts such as 55.3. Sparse mapper
carbon maps contribute only their shared positive map IDs. Only the remaining
unmapped product-position inflow is Unknown. Over-assigned donor or product
positions fail validation. Qualitative unit subreactions from Enumeration must
never replace the original stoichiometry in this solver.

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

Combined output is published atomically with each complete Enumeration round.
If the first round is rejected, previous output stays unchanged; emu_solved may
be true while emu_completed is false. Select EMU alone for independent delivery.
