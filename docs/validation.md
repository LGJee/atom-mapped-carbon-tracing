# Release validation

Validation performed on 2026-10-13. Results below separate actual software
execution from deterministic tests using mock adapters.

## MATLAB reference export

Executed the published reference script with MATLAB R2025b, the installed
COBRA Toolbox and GLPK, using the unmodified official Yeast9 v9.1.1 model.
The upstream YAML was also retrieved and compared byte-for-byte with the study
copy; it contains 1,787 SMILES annotations. Its SHA256 is
`2ddfaaa73ead41ab243127d1b34fa56db2c86052f6eabcd12bb4eea93947dc17`.

Both new workbook tables matched the original study inputs, including fluxes
within numerical comparison tolerances, equations, metabolite identities and
SMILES overrides:

| Strain | Nonzero reactions | Participating metabolite rows |
|---|---:|---:|
| WT | 493 | 522 |
| IUPD | 486 | 516 |

The original model, code, inputs and result files were not overwritten.

## Scientific behavior before and after refactoring

The original and release implementations read the same full WT/IUPD inputs
and the same previously generated three-tool candidates. No new predictions
were substituted during this comparison. Dataframes were compared, not merely
their sizes. Selected correspondence strings, mapper choice, support groups,
review/AI flags, fluxes, both rounds' state IDs, labels, labelled SMILES and
Summary values matched.

| Strain | Usable mappings | Round 1 rows | Round 2 rows | Mapped-network carbon metabolites |
|---|---:|---:|---:|---:|
| WT | 335 | 252 | 334 | 334 |
| IUPD | 335 | 266 | 342 | 337 |

This verifies preservation of the implementation's scientific behavior. It is
not biochemical validation of all candidate mappings.

## Real three-tool example

Executed the supplied synthetic acetaldehyde example with actual LocalMapper,
RXNMapper and EnzymeMap calls, followed by one round of carbon propagation.
EnzymeMap used the downloaded official v2 CSV and SHA256-verified rule file.
All three produced usable chemical-reaction mappings. The transport was handled
separately by identity correspondence. Consensus selected EnzymeMap for the
reduction; the labelled products were `acetaldehyde[c]` and `ethanol[c]`, each
carrying C1 and C2. `AIGenerated` remained false.

The example is synthetic and does not assert pFBA-derived fluxes or experimental
isotope measurements. The full biological networks were not freshly predicted
again solely for the refactoring comparison.
The separate `tests.test_enzymemap_rules` integration check forced an empty
database lookup and successfully executed the official single-step rule,
restoring the supplied structures with complete carbon correspondence.

## Deterministic regression checks

The five deterministic test programs (separate from the native-rule integration) cover:

- Exact substrate-state pairing, cleavage, multiple product copies and repeated
  carbon-origin labels.
- Structural alignment, symmetric equivalents, complete carbon coverage,
  selection priority and preservation of raw mapper predictions.
- Batch failure and individual retry, isolated errors and validated transport.
- Six concurrent mock prediction processes with a completion barrier.
- Separate strain propagation processes and equivalence to serial execution.
- Prediction cache validation, flux direction and input-overwrite protection.
- Equivalent combined and separate map/trace execution; trace makes no mapper
  calls and removes stale higher-round sheets.
- Rejection of duplicate IDs, zero/nonfinite flux and conflicting structures.
- Source-specific Summary counts and mapped-network denominator rules.

All checks passed. Python files passed syntax, unused-import and unused-variable
checks; the setup shell script passed `bash -n`. Example input/output sheets were
rendered and inspected. Copies retain plain styles, readable column widths,
wrapped equations and no added colour fills, filters or frozen panes.

## Environment and platform scope

The checked execution platform was Linux x86_64 under Ubuntu WSL on Windows.
The three tool environments are isolated. Critical tracing versions are pinned
to the original study stack; exact tool/model/data references are in the
environment files and main configuration.

The first clean RXNMapper installation exposed an incompatibility with
setuptools 84, which removed `pkg_resources`. Both the RXNMapper and EnzymeMap
specifications now pin setuptools 80.9.0. This is required by the selected
legacy tool stack and should not be removed without rechecking imports.
The clean EnzymeMap installation also required an explicit `requests=2.32.5`
dependency, which was implicit in the original workstation. This dependency
is now included in the environment specification.
All three isolated environments were newly created for release validation,
their functional checks passed, and the real three-tool example was rerun
using these environments. Existing workstation environments remained intact.
Rerunning the setup script reused all three environments and passed without
creating replacements. New-environment raw prediction sheets, selections,
final mappings, carbon flow and Summary matched the original-environment
example results.

Native Linux on a separate machine, native Windows Python, macOS, ARM and GPU
execution were not independently tested. Environment specifications support
Linux/WSL x86_64 only. The research tracing environment uses end-of-life Python
3.6 to preserve the checked serialization behavior.

## Code cleanup

Removed machine-specific Windows-to-WSL launch switching, personal interpreter
paths, the legacy nonzero-sheet mapping writer, old split-tool output lookup,
obsolete IUPD split-source sheet cleanup, unused root/output constants and unused
imports. Consolidated atomic workbook replacement and plain column-width
handling. Retained adapters, input/structure validation, failure propagation,
resource limits and the original label-transfer algorithm.

No AI review stage, trained mapper, web service, additional scientific correction
rules or skill package was added.
