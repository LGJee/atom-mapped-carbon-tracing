"""Compatibility alias for step 3's size-1 EMU-inspired source analysis."""
from carbon_flow_tools.cli import main as run_stage
from step3_carbon_flow import settings


def main(argv=None):
    options = settings()
    options["analysis_methods"] = ["EMU"]
    return run_stage("trace", options, argv)


if __name__ == "__main__":
    import sys
    sys.exit(main())
