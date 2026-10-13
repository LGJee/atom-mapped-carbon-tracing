"""Run atom mapping, carbon propagation, or both from one JSON configuration."""

import argparse
import json
import os
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parent
BACKENDS = ("LocalMapper", "RXNMapper", "EnzymeMap")


def load_config(filename):
    """Resolve relative paths against the configuration file, not the terminal."""
    filename = Path(filename).resolve()
    settings = json.loads(filename.read_text(encoding="utf-8"))
    base = filename.parent
    if not isinstance(settings, dict) or not settings.get("analyses"):
        raise ValueError("Configuration must define a nonempty analyses object.")
    environments = base / settings["environment_dir"]
    interpreters = {
        name: environments / folder / "bin/python"
        for name, folder in (
            ("LocalMapper", "localmapper"),
            ("RXNMapper", "rxnmapper"),
            ("EnzymeMap", "enzymemap"),
        )
    }
    for name, interpreter in settings.get("python_overrides", {}).items():
        if name not in BACKENDS:
            raise ValueError("Unknown interpreter override: " + name)
        interpreters[name] = base / interpreter
    configs = {
        name: {
            "module": "carbon_flow_tools.mappers." + module,
            "python": str(interpreters[name].resolve()),
        }
        for name, module in (
            ("LocalMapper", "localmapper_adapter"),
            ("RXNMapper", "rxnmapper_adapter"),
            ("EnzymeMap", "enzymemap_adapter"),
        )
    }
    enzyme = dict(settings.get("enzymemap", {}))
    for key in ("dataset", "rules"):
        if key in enzyme:
            enzyme[key] = str((base / enzyme[key]).resolve())
    configs["EnzymeMap"]["options"] = enzyme
    jobs = {
        name: dict(job, input_file=(base / job["input_file"]).resolve())
        for name, job in settings["analyses"].items()
    }
    return settings, jobs, configs, (base / settings["output_dir"]).resolve()


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "analyses", nargs="*", help="Analysis names; default: all configured analyses"
    )
    parser.add_argument("--config", type=Path, default=ROOT / "config.json")
    parser.add_argument("--stage", choices=("all", "map", "trace"), default="all")
    parser.add_argument(
        "--rounds", type=int, help="Override the configured propagation rounds"
    )
    parser.add_argument(
        "--mapper", choices=BACKENDS, nargs="+", help="Override configured tools"
    )
    parser.add_argument(
        "--refresh-mappings",
        action="store_true",
        help="Do not reuse stored predictions",
    )
    args = parser.parse_args(argv)
    if os.name == "nt":
        parser.error("Run this command inside a Linux or WSL terminal; see the README.")
    settings, jobs, configs, output = load_config(args.config)
    tracing_python = Path(configs["LocalMapper"]["python"])
    if not tracing_python.is_file():
        raise FileNotFoundError(
            "Tracing interpreter not found; run scripts/setup_environments.sh first: "
            + str(tracing_python)
        )
    if Path(sys.executable).resolve() != tracing_python.resolve():
        return subprocess.run(
            [str(tracing_python), str(Path(__file__).resolve())]
            + (sys.argv[1:] if argv is None else list(argv))
        ).returncode
    names = args.analyses or list(jobs)
    if any(name not in jobs for name in names):
        parser.error("Analysis names must be configured: " + ", ".join(jobs))
    rounds = settings["rounds"] if args.rounds is None else args.rounds
    backends = settings["mappers"] if args.mapper is None else args.mapper
    if any(name not in BACKENDS for name in backends):
        parser.error("Unsupported configured mapper")
    from carbon_flow_tools.pipeline import run_analyses

    run_analyses(
        {name: jobs[name] for name in names},
        rounds,
        settings["batch_size"],
        settings["max_reaction_atoms"],
        output,
        backends,
        configs,
        refresh=args.refresh_mappings,
        stage=args.stage,
        **settings.get("resources", {}),
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
