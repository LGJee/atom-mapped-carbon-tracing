# Reference non-IUPD/IUPD model conditions

The default workflow implements these reference study conditions. Its gene changes,
medium bounds, growth constraint and squalene objective are not requirements
of the carbon-tracing program. For another organism or mixed-carbon culture,
supply a workbook with the same input schema and configure your own sources.

## Obtain the model

Obtain `data/Saccharomyces_cerevisiae/yeast-GEM.mat` and `data/Saccharomyces_cerevisiae/yeast-GEM.yml` from
[yeast-GEM v9.1.1](https://github.com/SysBioChalmers/yeast-GEM/tree/v9.1.1),
commit `2d594ae1c4a2d550ccef120d96a58c7bbf586255`.
Place both files in this project's `data/Saccharomyces_cerevisiae/`; configure only the MAT path in
step 1 when relocating them. Step 2 discovers matching annotations in that folder.
The files must be from the same version; the YAML reader targets this version's
metabolite ID and SMILES annotation layout.

## Generate reference inputs

Keep the scientific settings unchanged for reference reproduction. Adjust model paths
in the corresponding strain blocks, or the shared solver setting at the top of
[step1_prepare_carbon_tracing_inputs.py](../step1_prepare_carbon_tracing_inputs.py),
then run `python3 step1_prepare_carbon_tracing_inputs.py`. The preparation
environment supplies COBRApy and the reference solver, GLPK.

The script preserves these **study-specific**, not universal, conditions:

| Setting | Reference value |
|---|---|
| Model / growth objective | Yeast9 9.1.1 / `r_2111` |
| Glucose lower bound during growth optimization | -1 mmol/gDW/h |
| Oxygen lower bound | Inherited native value: -1000 mmol/gDW/h |
| IUPD | Delete `YML126C` (ERG13); add isoprenol utilization reactions |
| IUPD isoprenol lower bound | -2.09 mmol/gDW/h |
| Fixed growth for both non-IUPD and IUPD | 50% of non-IUPD's maximum at native glucose bounds |
| non-IUPD glucose lower bound after fixing growth | -2 mmol/gDW/h |
| Production objective / nonzero tolerance | Squalene `r_1012` / `abs(flux)>1e-8` |

Remaining medium bounds are inherited from the model. Negative exchange flux
denotes uptake; a lower bound permits but does not force uptake.
`cobra.flux_analysis.pfba(model)` preserves the optimal production objective,
then minimizes the sum of absolute fluxes. This is one representative pFBA
solution, not every possible pathway or optimal solution.

For the S. cerevisiae pair, the script writes these two workbooks:

| Workbook | Reaction sheet | Metabolite sheet |
|---|---|---|
| `results/Saccharomyces_cerevisiae_non_IUPD_1_pFBA_nonzero_reactions.xlsx` | `Nonzero` | `Metabolites` |
| `results/Saccharomyces_cerevisiae_IUPD_1_pFBA_nonzero_reactions.xlsx` | `Nonzero` | `Metabolites` |

Reaction columns: `ReactionID`, `ReactionName`, `Flux`, `Equation`.
Metabolite columns include `MetaboliteID`, `MetaboliteName`, `MolecularFormula`
and `Charge`, plus available database identifiers. Step 1 does not read SMILES.
Names include compartments, such as `D-glucose[c]`; manual SMILES take priority.
Step 2 reads the pinned Yeast9 YAML by metabolite ID, and applies the verified
isoprenol and FBP(4−) structures from
`data/manual_smiles.csv`. Input export performs no structure retrieval.
Step 2 enables PubChem completion by default; disable `COMPLETE_MISSING_SMILES`
for manual-only reference comparisons.
Running this input export again replaces its previous input exports.

## Configure and run the reference analyses

Step 1 supplies the reference inputs and ordered carbon sources through its
automatically generated plan. Step 3 sets the 10-round limit with `ROUNDS`. Run steps 2 and 3 after the Python
export without editing their default settings. For another study, change step 1's
model blocks and analysis settings; advanced mapper and propagation execution
settings remain at the top of the corresponding Python step.

```bash
python3 step2_atom_mapping.py
python3 step3_carbon_flow.py
```

Use `Saccharomyces_cerevisiae_non_IUPD` or `Saccharomyces_cerevisiae_IUPD` as a positional argument to run only that reference analysis.

Glucose initialization uses biochemical C1-C6 for the supported cyclic form.
Isoprenol is numbered from the hydroxyl end along the chain, then the terminal
alkene carbon and methyl branch. Glucose followed by isoprenol therefore uses
C1-C6 and C7-C11, respectively; isoprenol alone uses C1-C5.
These labels identify source atom origins, not experiment-derived abundances.

Return to the [main workflow](../README.md).

## E. coli comparison models

The reference setup also includes E. coli iML1515. Its IUPD condition knocks out
`b2515` (ispG, reaction `MECDPDH5`) and adds the same two-step isoprenol
phosphorylation chemistry as Yeast9, using cytosolic ATP, ADP and IPP identifiers.
Both comparison models include heterologous SQS activity and a squalene output,
using the stoichiometry of [RHEA:32295](https://www.rhea-db.org/rhea/32295).
These are model-engineering comparisons, not experimentally validated strains.

The non-IUPD condition uses the native glucose limit of 10 mmol/gDW/h to determine
maximum growth, fixes growth at 50% of that value, then permits glucose uptake up
to 20 while maximizing squalene and minimizing total absolute flux. IUPD uses the
same fixed growth reference, keeps the glucose limit at 10, and permits isoprenol
uptake up to 20.9. Run the non-IUPD block before its IUPD block. Native oxygen
bounds are inherited and other carbon boundary uptake is closed.
These settings are specific to this comparison, not universal culture conditions.
