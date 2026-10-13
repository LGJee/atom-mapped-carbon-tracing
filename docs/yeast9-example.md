# Yeast9 WT/IUPD reference example

This optional example reproduces the original study setup. Its gene changes,
medium bounds, growth constraint and squalene objective are not requirements
of the carbon-tracing program. For another organism or mixed-carbon culture,
supply a workbook with the same input schema and configure your own sources.

## Obtain the model

Obtain `model/yeast-GEM.mat` and `model/yeast-GEM.yml` from
[yeast-GEM v9.1.1](https://github.com/SysBioChalmers/yeast-GEM/tree/v9.1.1),
commit `2d594ae1c4a2d550ccef120d96a58c7bbf586255`.
Place both files in this project's `model/`, or configure their folder below.
The files must be from the same version; the YAML reader targets this version's
metabolite ID and SMILES annotation layout.

## Generate reference inputs

Edit the settings at the start of [pFBA_network.m](../pFBA_network.m), then run it
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

## Configure and run the reference analyses

The root `config.json` contains these WT/IUPD example settings. Use them after
running `pFBA_network.m`, or replace its `analyses` object with your own analysis.

```bash
python3 run_carbon_tracing.py --config config.json --stage map
python3 run_carbon_tracing.py --config config.json --stage trace --rounds 10
```

Use `WT` or `IUPD` as a positional argument to run only that reference analysis.

Glucose initialization uses biochemical C1-C6 for the supported cyclic form.
Isoprenol is numbered from the hydroxyl end along the chain, then the terminal
alkene carbon and methyl branch. Glucose followed by isoprenol therefore uses
C1-C6 and C7-C11, respectively; isoprenol alone uses C1-C5.
These labels identify source atom origins, not experiment-derived abundances.

Return to the [main workflow](../README.md).
