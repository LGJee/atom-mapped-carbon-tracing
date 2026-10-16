"""Step 3: propagate normal and AI-reviewed step-2 mappings separately."""

import sys
from carbon_flow_tools.cli import main
from step2_atom_mapping import OUTPUT_DIR as MAPPING_DIR

# Read both <model>_2_atom_mapping.xlsx and its optional _AIGenerated copy.
# Each mapping produces its own carbon-flow workbook.
# Input paths and ordered carbon sources come from step 1.
# Input/output for legacy internal callers; the generated plan supplies real paths.
INPUT_DIR = MAPPING_DIR
OUTPUT_DIR = INPUT_DIR

# Fallback values for callers that do not use an analysis plan.
CARBON_SOURCES = {}  # Loaded from step 1 for every configured model.

# Propagation
ROUNDS = 10  # Maximum rounds; stop earlier when propagation converges.
MAX_EXCEL_SIZE_MB = 50  # Decimal MB; discard a candidate round reaching this size.
OUTPUT_BATCH_ROWS = 50000  # Check actual candidate XLSX size after each batch/layer.
# Record CO2/bicarbonate endpoints, but exclude their label recycling by default.
# Set True to trace inorganic-carbon sources or include carbon refixation.
ALLOW_INORGANIC_CARBON_RECYCLING = False

# Resources: None selects a CPU/memory-limited capacity.
MAX_CARBON_PROCESSES = None
CARBON_MEMORY_GB = None
# Stop safely below this free-memory reserve or 10% of Linux/WSL memory.
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
        rounds=ROUNDS,
        max_excel_size_mb=MAX_EXCEL_SIZE_MB,
        output_batch_rows=OUTPUT_BATCH_ROWS,
        allow_inorganic_carbon_recycling=ALLOW_INORGANIC_CARBON_RECYCLING,
        max_carbon_processes=MAX_CARBON_PROCESSES,
        carbon_memory_gb=CARBON_MEMORY_GB,
        min_available_memory_gb=MIN_AVAILABLE_MEMORY_GB,
        tracing_python=TRACING_PYTHON,
    )


if __name__ == "__main__":
    sys.exit(main("trace", settings()))
