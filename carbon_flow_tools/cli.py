"""Shared command handling for mapping and carbon propagation entries."""

import argparse
import hashlib
import json
import os
import re
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
BACKENDS = ("LocalMapper", "RXNMapper", "EnzymeMap")


def project_path(value):
    return (ROOT / value).resolve()


def build_mapper_configs(settings):
    """Resolve mapper environments and trusted data paths from step 2 settings."""
    environments = project_path(settings["environment_dir"])
    interpreters = {
        name: environments / folder / "bin/python"
        for name, folder in (
            ("LocalMapper", "localmapper"),
            ("RXNMapper", "rxnmapper"),
            ("EnzymeMap", "enzymemap"),
        )
    }
    for name, interpreter in settings["python_overrides"].items():
        if name not in BACKENDS:
            raise ValueError("Unknown interpreter override: " + name)
        interpreters[name] = project_path(interpreter)
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
    enzyme = dict(settings["enzymemap"])
    for key in ("dataset", "rules"):
        enzyme[key] = str(project_path(enzyme[key]))
    configs["EnzymeMap"]["options"] = enzyme
    return configs


def discover_inputs(directory, stage="map", ai_generated=False):
    """Discover numbered input workbooks for either Python stage."""
    number, field, pattern, previous = (
        ("1", "input_file", "*_1_*.xlsx", "step1_prepare_carbon_tracing_inputs.m")
        if stage == "map"
        else ("2", "mapping_file", "*_2_atom_mapping.xlsx", "step2_atom_mapping.py")
    )
    if stage == "trace" and ai_generated:
        pattern = "*_2_atom_mapping_AIGenerated.xlsx"
    models = {}
    for path in sorted(directory.glob(pattern)):
        if (
            not path.is_file()
            or path.name.startswith("~$")
            or path.name.endswith(".tmp.xlsx")
        ):
            continue
        name = re.fullmatch(r"(.+)_" + number + r"_.+\.xlsx", path.name).group(1)
        if name in models:
            raise ValueError(
                "Multiple step {} inputs for model {} in {}".format(
                    number, name, directory
                )
            )
        models[name] = {field: path}
    if not models:
        raise FileNotFoundError(
            "No {} inputs in {}. Run {} first, or set INPUT_DIR to its output folder.".format(
                pattern, directory, previous
            )
        )
    return models


def load_analysis_plan(path, stage, defaults):
    """Read MATLAB's generated plan; never parse or execute MATLAB source."""
    path = project_path(path)
    if not path.is_file():
        raise FileNotFoundError(
            "Run step1_prepare_carbon_tracing_inputs.m first: " + str(path)
        )
    with path.open(encoding="utf-8") as handle:
        plan = json.load(handle)
    if not isinstance(plan, dict) or plan.get("version") != 3:
        raise ValueError(
            "Unsupported analysis plan version; rerun step1_prepare_carbon_tracing_inputs.m"
        )
    if set(plan) != {"version", "models", "output_dir"}:
        raise ValueError(
            "Unexpected analysis plan fields; rerun step1_prepare_carbon_tracing_inputs.m"
        )
    models = plan.get("models")
    if not isinstance(models, dict) or not models:
        raise ValueError("The analysis plan must contain at least one model")
    output_dir = plan["output_dir"]
    if not isinstance(output_dir, str) or not output_dir or "\\" in output_dir:
        raise ValueError("Use a nonempty portable output path in the analysis plan")
    settings = dict(defaults)
    settings["output_dir"] = output_dir
    if stage == "trace" and (
        type(settings["rounds"]) is not int or settings["rounds"] < 1
    ):
        raise ValueError("Propagation rounds must be a positive integer")
    jobs = {}
    for name, entry in models.items():
        if not re.fullmatch(r"[A-Za-z][A-Za-z0-9_]*", name) or not isinstance(
            entry, dict
        ):
            raise ValueError("Invalid analysis name or model settings: " + str(name))
        digest = entry.get("input_sha256", "")
        if not isinstance(digest, str) or not re.fullmatch(r"[0-9a-f]{64}", digest):
            raise ValueError("Run step 1 to record the input SHA256 for " + name)
        input_file = entry.get("input_file")
        sources = entry.get("carbon_sources")
        if not isinstance(input_file, str) or not input_file or "\\" in input_file:
            raise ValueError("Use portable input paths in the analysis plan: " + name)
        if (
            not isinstance(sources, list)
            or not sources
            or any(not isinstance(s, str) or not s.strip() for s in sources)
            or len(set(sources)) != len(sources)
        ):
            raise ValueError(
                "Configure nonempty unique ordered carbon sources for " + name
            )
        if stage == "map":
            jobs[name] = {
                "input_file": project_path(input_file),
                "input_sha256": digest,
            }
        else:
            suffix = "_AIGenerated" if settings.get("ai_generated") else ""
            mapping = project_path(output_dir) / (
                name + "_2_atom_mapping" + suffix + ".xlsx"
            )
            jobs[name] = {
                "mapping_file": mapping,
                "carbon_sources": sources,
                "input_sha256": digest,
            }
    print("Analysis plan: {} ({} models)".format(path, len(jobs)), flush=True)
    return settings, jobs


def main(stage, settings, argv=None):
    ai_generated = settings.get("ai_generated", False)
    parser = argparse.ArgumentParser(
        description="Predict and select atom mappings."
        if stage == "map"
        else "Propagate carbon labels using step 2 workbooks."
    )
    parser.add_argument(
        "analyses",
        nargs="*",
        help="Model names detected from the previous step's output filenames",
    )
    parser.add_argument(
        "--example",
        action="store_true",
        help="Run the supplied synthetic installation example",
    )
    if stage == "trace":
        parser.add_argument(
            "--rounds", type=int, help="Override the configured maximum rounds"
        )
    else:
        parser.add_argument(
            "--mapper", choices=BACKENDS, nargs="+", help="Override this script's tools"
        )
        parser.add_argument(
            "--refresh-mappings",
            action="store_true",
            help="Do not reuse step 2 predictions",
        )
    args = parser.parse_args(argv)
    if os.name == "nt":
        parser.error("Run inside a Linux or WSL terminal; see the README.")
    planned_models = None
    if not args.example and settings.get("analysis_plan"):
        settings, planned_models = load_analysis_plan(
            settings["analysis_plan"], stage, settings
        )
    configs = build_mapper_configs(settings) if stage == "map" else {}
    tracing_python = project_path(
        configs["LocalMapper"]["python"]
        if stage == "map"
        else settings["tracing_python"]
    )
    if not tracing_python.is_file():
        raise FileNotFoundError(
            "Tracing interpreter not found; run scripts/setup_environments.sh: "
            + str(tracing_python)
        )
    if Path(sys.executable).resolve() != tracing_python:
        script = ROOT / (
            "step2_atom_mapping.py"
            if stage == "map"
            else "step4_carbon_flow_ai.py"
            if ai_generated
            else "step3_carbon_flow.py"
        )
        return subprocess.run(
            [str(tracing_python), str(script)]
            + (sys.argv[1:] if argv is None else list(argv))
        ).returncode
    output = project_path(settings["output_dir"])
    if args.example:
        output = ROOT / "scripts/results/example"
        models = (
            {
                "Example": {
                    "input_file": "examples/Example_1_input.xlsx",
                }
            }
            if stage == "map"
            else {
                "Example": {
                    "mapping_file": "scripts/results/example/Example_2_atom_mapping"
                    + ("_AIGenerated" if ai_generated else "")
                    + ".xlsx",
                    "carbon_sources": ["acetaldehyde[e]"],
                }
            }
        )
    else:
        models = (
            planned_models
            if planned_models is not None
            else discover_inputs(
                project_path(settings["input_dir"]), stage, ai_generated
            )
        )
    if ai_generated:
        models = {
            name: job
            for name, job in models.items()
            if project_path(job["mapping_file"]).is_file()
        }
        if not models:
            raise FileNotFoundError(
                "No *_2_atom_mapping_AIGenerated.xlsx files found. Place AI-reviewed "
                "mapping copies in the configured results folder before running step 4."
            )
    names = args.analyses or list(models)
    if (
        not names
        or len(set(names)) != len(names)
        or any(name not in models for name in names)
    ):
        parser.error("Select unique available model names: " + ", ".join(models))
    jobs = {name: dict(models[name]) for name in names}
    if stage == "map":
        import pandas as pd
        from .excel_io import input_sheets

        for name, job in jobs.items():
            job["input_file"] = project_path(job["input_file"])
            if (
                "input_sha256" in job
                and hashlib.sha256(job["input_file"].read_bytes()).hexdigest()
                != job["input_sha256"]
            ):
                raise ValueError("Input changed after step 1; rerun step 1: " + name)
            with pd.ExcelFile(job["input_file"], engine="openpyxl") as excel:
                job["reaction_sheet"], job["metabolite_sheet"] = input_sheets(excel)
            print(
                "{}: input {}; sheets {}, {}".format(
                    name,
                    job["input_file"],
                    job["reaction_sheet"],
                    job["metabolite_sheet"],
                ),
                flush=True,
            )
        backends = args.mapper or settings["mappers"]
        if any(name not in BACKENDS for name in backends):
            parser.error("Unsupported mapper in script settings")
        from .pipeline import run_analyses

        structure_options = dict(settings["structure_options"])
        for field in ("manual_file", "cache_file"):
            structure_options[field] = project_path(structure_options[field])
        structure_options["model_yaml_files"] = {
            name: project_path(path)
            for name, path in structure_options.get("model_yaml_files", {}).items()
        }
        if args.example:
            structure_options["lookup"] = False
        run_analyses(
            jobs,
            1,
            settings["batch_size"],
            settings["max_reaction_atoms"],
            output,
            backends,
            configs,
            stage="map",
            refresh=args.refresh_mappings or settings["refresh"],
            max_mapper_processes=1
            if args.example
            else settings["max_mapper_processes"],
            mapper_memory_gb=settings["mapper_memory_gb"],
            structure_options=structure_options,
        )
    else:
        for name, job in jobs.items():
            if not args.example and planned_models is None:
                if name not in settings["carbon_sources"]:
                    parser.error(
                        "Set CARBON_SOURCES for {} at the top of step3_carbon_flow.py".format(
                            name
                        )
                    )
                job["carbon_sources"] = settings["carbon_sources"][name]
            job["mapping_file"] = project_path(job["mapping_file"])
            job["allow_inorganic_carbon_recycling"] = settings.get(
                "allow_inorganic_carbon_recycling", False
            )
        if planned_models is not None:
            import pandas as pd

            for name, job in jobs.items():
                if not job["mapping_file"].is_file():
                    raise FileNotFoundError(
                        "Run step2_atom_mapping.py first: " + str(job["mapping_file"])
                    )
                with pd.ExcelFile(job["mapping_file"], engine="openpyxl") as excel:
                    info = pd.read_excel(excel, sheet_name="Mapping_Info")
                hashes = (
                    set(info["SourceInputSHA256"].dropna())
                    if "SourceInputSHA256" in info
                    else set()
                )
                if hashes != {job["input_sha256"]}:
                    raise ValueError(
                        "Mapping does not match this step-1 input; rerun step2_atom_mapping.py: "
                        + name
                    )
        rounds = (
            args.rounds
            if args.rounds is not None
            else (1 if args.example else settings["rounds"])
        )
        from .pipeline import run_carbon_flow

        run_carbon_flow(
            jobs,
            rounds,
            output,
            max_processes=1 if args.example else settings["max_carbon_processes"],
            memory_gb=settings["carbon_memory_gb"],
        )
    return 0
