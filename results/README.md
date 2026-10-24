# Published computed results

These are model preparation inputs and computed analysis outputs, not raw fermentation measurements. No study computations were rerun for publication. Four conditions cover Saccharomyces cerevisiae and Escherichia coli, each non-IUPD/IUPD. non-IUPD does not necessarily mean unmodified wild type.

## Configuration and provenance

- [Analysis plan](../data/analysis_plan.json): validated input hashes, source-model hashes, uptake constraints, objectives and solver metadata.
- [SHA-256 manifest](manifest.json): exact identities of all 20 unchanged workbooks.
- [Preparation script](../step1_prepare_carbon_tracing_inputs.py) and [condition notes](../docs/yeast9-example.md).
- [AI review reports](Ai_reports/): saved review rationale, not proof of a new stage-3 run.

## Workbooks

| Download | Stage | Version | Bytes |
|---|---|---|---:|
| [Escherichia_coli_IUPD_1_pFBA_nonzero_reactions.xlsx](Escherichia_coli_IUPD_1_pFBA_nonzero_reactions.xlsx) | 1 | Original | 46588 |
| [Escherichia_coli_IUPD_2_atom_mapping.xlsx](Escherichia_coli_IUPD_2_atom_mapping.xlsx) | 2 | Original | 597303 |
| [Escherichia_coli_IUPD_2_atom_mapping_AIGenerated.xlsx](Escherichia_coli_IUPD_2_atom_mapping_AIGenerated.xlsx) | 2 | AI-reviewed | 586931 |
| [Escherichia_coli_IUPD_3_carbon_flow.xlsx](Escherichia_coli_IUPD_3_carbon_flow.xlsx) | 3 | Original | 33099 |
| [Escherichia_coli_IUPD_3_carbon_flow_AIGenerated.xlsx](Escherichia_coli_IUPD_3_carbon_flow_AIGenerated.xlsx) | 3 | AI-reviewed | 32783 |
| [Escherichia_coli_non_IUPD_1_pFBA_nonzero_reactions.xlsx](Escherichia_coli_non_IUPD_1_pFBA_nonzero_reactions.xlsx) | 1 | Original | 46818 |
| [Escherichia_coli_non_IUPD_2_atom_mapping.xlsx](Escherichia_coli_non_IUPD_2_atom_mapping.xlsx) | 2 | Original | 600295 |
| [Escherichia_coli_non_IUPD_2_atom_mapping_AIGenerated.xlsx](Escherichia_coli_non_IUPD_2_atom_mapping_AIGenerated.xlsx) | 2 | AI-reviewed | 590621 |
| [Escherichia_coli_non_IUPD_3_carbon_flow.xlsx](Escherichia_coli_non_IUPD_3_carbon_flow.xlsx) | 3 | Original | 31781 |
| [Escherichia_coli_non_IUPD_3_carbon_flow_AIGenerated.xlsx](Escherichia_coli_non_IUPD_3_carbon_flow_AIGenerated.xlsx) | 3 | AI-reviewed | 31553 |
| [Saccharomyces_cerevisiae_IUPD_1_pFBA_nonzero_reactions.xlsx](Saccharomyces_cerevisiae_IUPD_1_pFBA_nonzero_reactions.xlsx) | 1 | Original | 47124 |
| [Saccharomyces_cerevisiae_IUPD_2_atom_mapping.xlsx](Saccharomyces_cerevisiae_IUPD_2_atom_mapping.xlsx) | 2 | Original | 537257 |
| [Saccharomyces_cerevisiae_IUPD_2_atom_mapping_AIGenerated.xlsx](Saccharomyces_cerevisiae_IUPD_2_atom_mapping_AIGenerated.xlsx) | 2 | AI-reviewed | 539997 |
| [Saccharomyces_cerevisiae_IUPD_3_carbon_flow.xlsx](Saccharomyces_cerevisiae_IUPD_3_carbon_flow.xlsx) | 3 | Original | 33545 |
| [Saccharomyces_cerevisiae_IUPD_3_carbon_flow_AIGenerated.xlsx](Saccharomyces_cerevisiae_IUPD_3_carbon_flow_AIGenerated.xlsx) | 3 | AI-reviewed | 34382 |
| [Saccharomyces_cerevisiae_non_IUPD_1_pFBA_nonzero_reactions.xlsx](Saccharomyces_cerevisiae_non_IUPD_1_pFBA_nonzero_reactions.xlsx) | 1 | Original | 48072 |
| [Saccharomyces_cerevisiae_non_IUPD_2_atom_mapping.xlsx](Saccharomyces_cerevisiae_non_IUPD_2_atom_mapping.xlsx) | 2 | Original | 555047 |
| [Saccharomyces_cerevisiae_non_IUPD_2_atom_mapping_AIGenerated.xlsx](Saccharomyces_cerevisiae_non_IUPD_2_atom_mapping_AIGenerated.xlsx) | 2 | AI-reviewed | 558831 |
| [Saccharomyces_cerevisiae_non_IUPD_3_carbon_flow.xlsx](Saccharomyces_cerevisiae_non_IUPD_3_carbon_flow.xlsx) | 3 | Original | 31838 |
| [Saccharomyces_cerevisiae_non_IUPD_3_carbon_flow_AIGenerated.xlsx](Saccharomyces_cerevisiae_non_IUPD_3_carbon_flow_AIGenerated.xlsx) | 3 | AI-reviewed | 32392 |

## Reuse and reproduction

Run from the project root after installation. To analyze existing mappings without replacing published outputs:

```bash
python3 step3_carbon_flow.py --plan data/analysis_plan.json --output-dir results/recalculated
```

To recalculate from the models, obtain the exact upstream files in README and compare hashes with the plan; then run steps 1, 2 and 3. These commands replace corresponding outputs: use a separate checkout or a custom case for comparison. Step 2 does not regenerate AI review.

## Limits and historical data

All published stage-3 files have the current two-sheet EMU schema. The original root carbon_flow_run.json was no longer available when assembling this publication; no replacement success report was fabricated. Exact solver version/settings and mapping-to-output lineage cannot be fully certified for those existing exports. Filenames and schema establish the intended conditions/version, not proof of their full run history. Recompute into a separate directory to obtain a fresh verified run record.

This is size-1 EMU-inspired fixed-flux mean carbon-origin analysis, not full isotopomer simulation. Unknown remains explicit, compartments stay distinct, partial-structure fractions cover known carbon only. A single pFBA solution is not all possible pathways. Lower Unknown or larger coverage does not establish higher accuracy. No manuscript DOI or experimental validation is asserted.

Historical enumeration/summary workbooks are excluded. Original models, environment installations, database caches and full AI evidence archives are not redistributed. AI CSV Source fields retain upstream links and local evidence identifiers; cited local archives are not all published. Selected mappings are preserved in workbooks. Fresh AI review or other solver versions may differ. See [third-party notices](../THIRD_PARTY.md).
