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


def discover_inputs(directory, stage="map"):
    """Discover numbered input workbooks for either Python stage."""
    number, field, pattern, previous = (
        ("1", "input_file", "*_1_*.xlsx", "step1_prepare_carbon_tracing_inputs.py")
        if stage == "map"
        else ("2", "mapping_file", "*_2_atom_mapping*.xlsx", "step2_atom_mapping.py")
    )
    models = {}
    for path in sorted(directory.glob(pattern)):
        if (
            not path.is_file()
            or path.name.startswith("~$")
            or path.name.endswith(".tmp.xlsx")
        ):
            continue
        match = re.fullmatch(
            r"(.+)_2_atom_mapping(?:_AIGenerated)?\.xlsx"
            if stage == "trace" else r"(.+)_1_.+\.xlsx", path.name
        )
        if match is None:
            continue
        name = match.group(1)
        if stage == "trace" and name in models:
            continue
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
    """Read the generated input plan; never execute configuration as code."""
    path = project_path(path)
    if not path.is_file():
        raise FileNotFoundError(
            "Run step1_prepare_carbon_tracing_inputs.py first: " + str(path)
        )
    with path.open(encoding="utf-8") as handle:
        plan = json.load(handle)
    if not isinstance(plan, dict) or plan.get("version") != 3:
        raise ValueError(
            "Unsupported analysis plan version; rerun step1_prepare_carbon_tracing_inputs.py"
        )
    if set(plan) != {"version", "models", "output_dir"}:
        raise ValueError(
            "Unexpected analysis plan fields; rerun step1_prepare_carbon_tracing_inputs.py"
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
            source = entry.get("source_model")
            if source is not None:
                if (
                    not isinstance(source, dict)
                    or not isinstance(source.get("path"), str)
                    or not source["path"]
                    or not isinstance(source.get("sha256"), str)
                    or not re.fullmatch(r"[0-9a-f]{64}", source["sha256"])
                    or (
                        source.get("variable_name") is not None
                        and not isinstance(source["variable_name"], str)
                    )
                ):
                    raise ValueError(
                        "Invalid source model provenance; rerun step 1: " + name
                    )
                jobs[name]["source_model"] = dict(
                    source, path=str(project_path(source["path"]))
                )
        else:
            mapping = project_path(output_dir) / (name + "_2_atom_mapping.xlsx")
            jobs[name] = {
                "mapping_file": mapping,
                "carbon_sources": sources,
                "input_sha256": digest,
            }
    print("Analysis plan: {} ({} models)".format(path, len(jobs)), flush=True)
    return settings, jobs


def main(stage, settings, argv=None):
    parser = argparse.ArgumentParser(
        description="Predict and select atom mappings."
        if stage == "map"
        else "Enumerate possible carbon-source label propagation."
    )
    parser.add_argument(
        "analyses",
        nargs="*",
        help="Model names detected from the previous step's output filenames",
    )
    parser.add_argument(
        "--plan", help="Generated analysis plan for an independent case (project-relative or absolute)",
    )
    parser.add_argument(
        "--example",
        action="store_true",
        help="Run the supplied synthetic installation example",
    )
    if stage == "trace":
        parser.add_argument(
            "--rounds", type=int, help="Override the maximum Enumeration rounds"
        )
        parser.add_argument("--output-dir", help="Optional result directory; inputs stay unchanged")
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
    if args.example and args.plan:
        parser.error("Use either --plan or --example, not both.")
    if os.name == "nt":
        parser.error("Run inside a Linux or WSL terminal; see the README.")
    planned_models = None
    if args.plan:
        settings = dict(settings, analysis_plan=args.plan)
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
            else "step3_carbon_flow.py"
        )
        arguments = sys.argv[1:] if argv is None else list(argv)
        return subprocess.run([str(tracing_python), str(script)] + arguments).returncode
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
                    "mapping_file": "scripts/results/example/Example_2_atom_mapping.xlsx",
                    "carbon_sources": ["acetaldehyde[e]"],
                }
            }
        )
    else:
        models = (
            planned_models
            if planned_models is not None
            else discover_inputs(
                project_path(settings["input_dir"]), stage
            )
        )
    names = args.analyses or list(models)
    if (
        not names
        or len(set(names)) != len(names)
        or any(name not in models for name in names)
    ):
        parser.error("Select unique available model names: " + ", ".join(models))
    if stage == "trace" and args.output_dir:
        output = project_path(args.output_dir)
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
            job.setdefault("input_sha256", hashlib.sha256(job["input_file"].read_bytes()).hexdigest())
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
        structure_options["source_models"] = {
            name: job["source_model"]
            for name, job in jobs.items()
            if "source_model" in job
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
            manual_mapping_file=project_path(settings.get("manual_mapping_file", "data/manual_atom_mappings.csv")),
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
            job["analysis_methods"] = settings.get("analysis_methods", ["Enumeration"])
            job["export_detail_rounds"] = settings.get("export_detail_rounds", 3)
            job["emu_allow_inorganic_carbon_recycling"] = settings.get("emu_allow_inorganic_carbon_recycling", True)
            job["emu_tolerance"] = settings.get("emu_tolerance", 1e-8)
            job["min_available_memory_gb"] = settings.get("min_available_memory_gb", 2)
            job["max_excel_size_mb"] = settings.get("max_excel_size_mb", 50)
            job["output_batch_rows"] = settings.get("output_batch_rows", 50000)
            job["export_exclusive_source_metabolites"] = settings.get("export_exclusive_source_metabolites", False)
            directory = project_path(job["mapping_file"]).parent
            job["mapping_files"] = [
                path for path in (
                    directory / (name + "_2_atom_mapping.xlsx"),
                    directory / (name + "_2_atom_mapping_AIGenerated.xlsx"),
                ) if path.is_file()
            ]
            if not job["mapping_files"]:
                raise FileNotFoundError(
                    "No step-2 mappings for {}; run step2_atom_mapping.py first in {}"
                    .format(name, directory)
                )
            job["mapping_file"] = job["mapping_files"][0]
            job["allow_inorganic_carbon_recycling"] = settings.get(
                "allow_inorganic_carbon_recycling", False
            )
        if planned_models is not None:
            import pandas as pd

            for name, job in jobs.items():
                for mapping_file in job["mapping_files"]:
                    with pd.ExcelFile(mapping_file, engine="openpyxl") as excel:
                        info = pd.read_excel(excel, sheet_name="Mapping_Info")
                    hashes = (
                        set(info["SourceInputSHA256"].dropna())
                        if "SourceInputSHA256" in info
                        else set()
                    )
                    if hashes != {job["input_sha256"]}:
                        raise ValueError(
                            "Mapping does not match this step-1 input; rerun step2_atom_mapping.py: "
                            + str(mapping_file)
                        )
        rounds = (
            args.rounds
            if args.rounds is not None
            else (1 if args.example else settings["rounds"])
        )
        from .pipeline import run_carbon_flow
        run_carbon_flow(jobs, rounds, output,
            max_processes=1 if args.example else settings["max_carbon_processes"],
            memory_gb=settings["carbon_memory_gb"])
    return 0
