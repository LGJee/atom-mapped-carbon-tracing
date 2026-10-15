"""Step 4 (optional): propagate carbon labels using AI-reviewed mapping copies."""

import sys

from carbon_flow_tools.cli import main
from step3_carbon_flow import settings as propagation_settings


def settings():
    """Reuse step 3 settings; select only _AIGenerated mapping workbooks."""
    options = propagation_settings()
    options["ai_generated"] = True
    return options


if __name__ == "__main__":
    sys.exit(main("trace", settings()))
