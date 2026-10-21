"""Step 2: predict atom mappings and select final correspondences."""

import sys
from carbon_flow_tools.cli import main

# Mapping execution settings. Defaults work with the supplied environments.
# Step 1 supplies the input files and output folder; no paths need repeating.
# Fallback paths for callers that do not use an analysis plan.
INPUT_DIR = "results"  # Step 1 output folder; files and sheets are detected.
OUTPUT_DIR = "results"
# Mapping
MAPPERS = ["LocalMapper", "RXNMapper", "EnzymeMap"]
BATCH_SIZE = 1
MAX_REACTION_ATOMS = 800  # Count both reaction sides together.
REFRESH_MAPPINGS = False

# Structure completion: manual overrides always apply; False disables online lookup.
COMPLETE_MISSING_SMILES = True
MANUAL_SMILES_FILE = "data/manual_smiles.csv"
MANUAL_ATOM_MAPPINGS_FILE = "data/manual_atom_mappings.csv"
STRUCTURE_CACHE_FILE = "data/structure_cache.json"
STRUCTURE_TIMEOUT = 20  # Seconds per PubChem request.
STRUCTURE_MIN_INTERVAL = 0.25  # At most four sequential requests per second.
# Model structures are discovered from step 1 provenance and the source directory.

# Resources
MAX_MAPPER_PROCESSES = None  # None selects a CPU/memory-limited capacity.
MAPPER_MEMORY_GB = 6

# Installation and EnzymeMap data
ENVIRONMENT_DIR = ".envs"
PYTHON_OVERRIDES = {}  # Optional mapper-name/Python-path pairs.
ENZYMEMAP = {
    "dataset": "data/enzymemap/processed_reactions.csv",
    "rules": "data/enzymemap/rules.pkl",
    "rules_sha256": "21cdebccdda5ea4308555b6dd8dfb27413950725ace66e7b8cb7ed3dba2c1e1a",
    "source_commit": "c91578856f11312d7941b290f8b773eb41b54973",
    "normalize_protonation": True,
}


def settings():
    """Use these execution settings with the inputs exported by step 1."""
    return dict(
        analysis_plan="data/analysis_plan.json",
        input_dir=INPUT_DIR,
        manual_mapping_file=MANUAL_ATOM_MAPPINGS_FILE,
        structure_options=dict(
            lookup=COMPLETE_MISSING_SMILES,
            manual_file=MANUAL_SMILES_FILE,
            cache_file=STRUCTURE_CACHE_FILE,
            timeout=STRUCTURE_TIMEOUT,
            min_interval=STRUCTURE_MIN_INTERVAL,
        ),
        output_dir=OUTPUT_DIR,
        mappers=MAPPERS,
        batch_size=BATCH_SIZE,
        max_reaction_atoms=MAX_REACTION_ATOMS,
        max_mapper_processes=MAX_MAPPER_PROCESSES,
        mapper_memory_gb=MAPPER_MEMORY_GB,
        refresh=REFRESH_MAPPINGS,
        environment_dir=ENVIRONMENT_DIR,
        python_overrides=PYTHON_OVERRIDES,
        enzymemap=ENZYMEMAP,
    )


if __name__ == "__main__":
    sys.exit(main("map", settings()))
