"""Optional step 1: register existing input workbooks without modifying them."""

import hashlib
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from carbon_flow_tools.cli import load_analysis_plan

# Replace the installation example with your own Excel files and ordered sources.
INPUTS = {
    "Example": {
        "input_file": "examples/Example_1_input.xlsx",
        "carbon_sources": ["acetaldehyde[e]"],
    }
}
OUTPUT_DIR = "results"


def main():
    plan = {"version": 3, "models": {}, "output_dir": OUTPUT_DIR}
    for name, entry in INPUTS.items():
        path = (ROOT / entry["input_file"]).resolve()
        plan["models"][name] = dict(
            input_file=path.relative_to(ROOT.resolve()).as_posix(),
            carbon_sources=entry["carbon_sources"],
            input_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
        )
    target = ROOT / "data/analysis_plan.json"
    target.parent.mkdir(parents=True, exist_ok=True)
    temporary = target.with_suffix(".tmp")
    try:
        temporary.write_text(json.dumps(plan, indent=2) + "\n", encoding="utf-8")
        load_analysis_plan(temporary, "map", {})
        temporary.replace(target)
    finally:
        if temporary.exists():
            temporary.unlink()
    print("Analysis plan saved: " + str(target))


if __name__ == "__main__":
    main()
