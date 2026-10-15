"""Step 3: propagate carbon labels using the step 2 mapping workbooks."""

import sys
from carbon_flow_tools.cli import main
from step2_atom_mapping import OUTPUT_DIR as MAPPING_DIR

# Propagation settings. Input paths and ordered carbon sources come from step 1.
# Input/output for legacy internal callers; the generated plan supplies real paths.
INPUT_DIR = MAPPING_DIR
OUTPUT_DIR = INPUT_DIR

# Fallback values for callers that do not use an analysis plan.
CARBON_SOURCES = {}  # Loaded from step 1 for every configured model.

# Propagation
ROUNDS = 10  # Maximum rounds; stop earlier when propagation converges.
# Record CO2/bicarbonate endpoints, but exclude their label recycling by default.
# Set True to trace inorganic-carbon sources or include carbon refixation.
ALLOW_INORGANIC_CARBON_RECYCLING = False

# Resources: None selects a CPU/memory-limited capacity.
MAX_CARBON_PROCESSES = None
CARBON_MEMORY_GB = None

# Compatible tracing environment; no mapper is called here.
TRACING_PYTHON = ".envs/localmapper/bin/python"


def settings():
    """Pass only propagation settings; no mapper settings are needed."""
    return dict(
        analysis_plan="data/analysis_plan.json",
        input_dir=INPUT_DIR,
        carbon_sources=CARBON_SOURCES,
        output_dir=OUTPUT_DIR,
        rounds=ROUNDS,
        allow_inorganic_carbon_recycling=ALLOW_INORGANIC_CARBON_RECYCLING,
        max_carbon_processes=MAX_CARBON_PROCESSES,
        carbon_memory_gb=CARBON_MEMORY_GB,
        tracing_python=TRACING_PYTHON,
    )


if __name__ == "__main__":
    sys.exit(main("trace", settings()))
