"""Step 3: Enumeration and size-1 EMU-inspired source balances in one workbook."""

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

# Methods share the same selected mapping but solve independently.
ANALYSIS_METHODS = ["Enumeration", "EMU"]

# Enumeration
EXPORT_DETAIL_ROUNDS = 3  # Later rounds retain only unique states and Summary, not provenance rows.
ROUNDS = 10  # Maximum rounds; stop earlier when propagation converges.
MAX_ENUMERATION_PROCESSES = None  # Share each model's CPU/memory budget; 1 forces serial candidates.
# Optional list for the last saved round; Summary always includes the counts.
EXPORT_EXCLUSIVE_SOURCE_METABOLITES = False
OUTPUT_BATCH_ROWS = 50000  # Upper row bound; buffers also flush at an estimated 8 MiB.
# Record CO2/bicarbonate endpoints, but exclude their label recycling by default.
# Set True to trace inorganic-carbon sources or include carbon refixation.
ALLOW_INORGANIC_CARBON_RECYCLING = False

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
        rounds=ROUNDS,
        analysis_methods=ANALYSIS_METHODS,
        export_detail_rounds=EXPORT_DETAIL_ROUNDS,
        max_enumeration_processes=MAX_ENUMERATION_PROCESSES,
        emu_allow_inorganic_carbon_recycling=EMU_ALLOW_INORGANIC_CARBON_RECYCLING,
        emu_tolerance=EMU_TOLERANCE,
        export_exclusive_source_metabolites=EXPORT_EXCLUSIVE_SOURCE_METABOLITES,
        max_excel_size_mb=MAX_EXCEL_SIZE_MB,
        output_batch_rows=OUTPUT_BATCH_ROWS,
        allow_inorganic_carbon_recycling=ALLOW_INORGANIC_CARBON_RECYCLING,
        max_carbon_processes=MAX_CARBON_PROCESSES,
        carbon_memory_gb=CARBON_MEMORY_GB,
        min_available_memory_gb=MIN_AVAILABLE_MEMORY_GB,
        tracing_python=TRACING_PYTHON,
    )


if __name__ == "__main__":
    sys.exit(run_cli("trace", settings()))
