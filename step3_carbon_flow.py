"""Step 3: fixed-flux, size-1 EMU-inspired carbon-source balances."""

import sys
from carbon_flow_tools.cli import run_cli
from step2_atom_mapping import OUTPUT_DIR as MAPPING_DIR

# Read both <model>_2_atom_mapping.xlsx and its optional _AIGenerated copy.
# Each mapping produces its own carbon-flow workbook.
# Input paths and ordered carbon sources come from step 1.
# Input/output for legacy internal callers; the generated plan supplies real paths.
INPUT_DIR = MAPPING_DIR
OUTPUT_DIR = INPUT_DIR

# Fallback values for callers that do not use an analysis plan.
CARBON_SOURCES = {}  # Loaded from step 1 for every configured model.

# EMU: fixed-flux, single-carbon origin balances, not a full isotopomer simulator.
EMU_ALLOW_INORGANIC_CARBON_RECYCLING = True
EMU_TOLERANCE = 1e-8

# Shared output and execution limits.
MAX_EXCEL_SIZE_MB = 50  # Decimal MB; reject an oversized candidate before replacement.
# Independent model/mapping-version processes.
MAX_CARBON_PROCESSES = None
CARBON_MEMORY_GB = None
# Check this free-memory reserve or 10% of Linux/WSL memory.
MIN_AVAILABLE_MEMORY_GB = 2.0

# Compatible tracing environment; no mapper is called here.
TRACING_PYTHON = ".envs/localmapper/bin/python"


def settings():
    """Pass only propagation settings; no mapper settings are needed."""
    return dict(
        analysis_plan="data/analysis_plan.json",
        input_dir=INPUT_DIR,
        carbon_sources=CARBON_SOURCES,
        output_dir=OUTPUT_DIR,
        emu_allow_inorganic_carbon_recycling=EMU_ALLOW_INORGANIC_CARBON_RECYCLING,
        emu_tolerance=EMU_TOLERANCE,
        max_excel_size_mb=MAX_EXCEL_SIZE_MB,
        max_carbon_processes=MAX_CARBON_PROCESSES,
        carbon_memory_gb=CARBON_MEMORY_GB,
        min_available_memory_gb=MIN_AVAILABLE_MEMORY_GB,
        tracing_python=TRACING_PYTHON,
    )


if __name__ == "__main__":
    sys.exit(run_cli("trace", settings()))
