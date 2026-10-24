# Installation

Choose **Option A** to let Codex perform the installation, or **Option B** to run
the commands yourself. To view existing Excel results, no installation is needed.

## Option A: Install with Codex

- **Windows:** you may use Codex on Windows and let it operate through a WSL
  Ubuntu terminal. The project and its dependencies run inside WSL Ubuntu.
  Codex does not need to be installed inside Ubuntu.
- **Linux:** use Codex with access to the local Linux terminal. The project runs
  directly on Linux; no WSL is required.

Both paths require Linux x86_64 for the project environments. On Windows, a
working WSL Ubuntu installation is a prerequisite. The project setup script
configures project dependencies; it does **not** install WSL or Ubuntu.

Open Codex with access to local files and a terminal, then paste the following
instruction into a new task. Codex will perform installation and small checks;
starting a scientific analysis is a separate step.

Copy everything in the instruction below:

> Install https://github.com/LGJee/atom-mapped-carbon-tracing on this computer.
> Find and reuse an existing checkout if available; do not overwrite local changes,
> models, plans or results, and do not silently upgrade an existing installation.
> Read this version's README, INSTALL.md and applicable project instructions first.
> Inspect the operating system and architecture. The supplied environments target
> Linux x86_64. On Windows, run project commands inside WSL Ubuntu, not directly
> in PowerShell. Codex may run on Windows with access to the WSL terminal; it does
> not need to be installed inside Ubuntu. If WSL Ubuntu is unavailable, explain
> that prerequisite and help the user set it up before continuing; the project
> setup script does not install WSL. On Linux, use the local terminal directly.
> Do not claim other platforms work.
> Check Git and Conda/Miniforge, install missing prerequisites using official sources
> where authorized, then follow the current setup instructions. Reuse compatible
> environments. Run the bounded installation checks and report exactly what passed,
> what failed, the project location and the next command. Do not start the study's
> four-model analysis, replace existing data, or equate installation with scientific
> analysis. If network access or system permissions block installation, explain the
> actual blocker instead of claiming success.

## Option B: Install manually

Run in a Linux x86_64 terminal (on Windows, inside WSL), not PowerShell.
Git and an initialized Conda installation are required. Install Miniforge from
https://github.com/conda-forge/miniforge if needed. Clone only if no checkout exists:

```bash
git clone https://github.com/LGJee/atom-mapped-carbon-tracing.git
cd atom-mapped-carbon-tracing
bash scripts/setup_environments.sh
python3 scripts/download_enzymemap.py
```

All subsequent commands below run from that project root. The setup script creates
four isolated environments and reuses existing paths after checking them. Its checks
include real small predictor calls; weights may download. EnzymeMap data are obtained
separately by the downloader. Network access, disk space and upstream availability
are required. Existing incompatible environments fail rather than being deleted.
See [the detailed installation guide](docs/user-guide.md#1-install-python-environments-and-enzymemap-data)
for recovery and custom environment paths.

For a bounded **offline** end-to-end check after dependencies are available:

```bash
.envs/localmapper/bin/python -m unittest scripts.test_analysis_plan scripts.test_step3_methods
.envs/preparation/bin/python scripts/check_environment.py preparation
```

The first command creates temporary inputs, uses deterministic mapper fixtures and
checks ordinary/AI EMU exports. It does not test real predictor accuracy or download
models. The second solves a tiny real COBRApy/GLPK pFBA problem.

An optional real-mapper example (may download weights; writes Example files in results):

```bash
python3 step2_atom_mapping.py --example
python3 step3_carbon_flow.py --example
```

Do not confuse this optional predictor run with the offline check. Full study commands
and model downloads are in [README](README.md). Native Windows/macOS/ARM installation
is not supported by these pinned specifications. Python 3.6 in the legacy mapper
stack is end-of-life; this is a reproduction environment, not a network service.
A new network installation and the real three-mapper example were not rerun during
this publication; see [actual verification](docs/publication-validation.md).

## Common problems

- `Conda was not found`: initialize Conda in the Linux shell and reopen it.
- `Tracing interpreter not found`: run setup from the project root; check configured
  interpreter paths if using a custom environment directory.
- Failed package/weight/database download: resolve connectivity and rerun the affected
  installer or downloader. Do not delete unrelated environments or disable rule hashes.
- Input hash mismatch: use the plan corresponding to those inputs; regenerate mapping
  after deliberate step-1 changes. Do not edit hashes to bypass validation.
- Missing source-model file: obtain the official model listed in README. Published
  step-2 workbooks remain self-contained inputs for step 3; step 1/2 source annotation
  discovery may need the original files.
- Resource-limit failure: inspect the current run report, free resources and reduce
  MAX_CARBON_PROCESSES. A previous workbook does not prove the new run succeeded.
