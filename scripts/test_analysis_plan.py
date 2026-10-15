"""Check the MATLAB-to-Python plan, six-model isolation and both stages."""

import hashlib
import json
import sys
import tempfile
from pathlib import Path
from unittest.mock import patch

import pandas as pd

from carbon_flow_tools.cli import load_analysis_plan, main
from scripts.test_parallel_carbon import make_input
import step2_atom_mapping as mapping
import step3_carbon_flow as tracing
from examples import step1_use_existing_inputs as external


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run():
    with tempfile.TemporaryDirectory() as temporary:
        folder = Path(temporary)
        (folder / "absent.csv").write_text(
            "Model,MetaboliteID,MetaboliteName,SMILES\n", encoding="utf-8"
        )
        names = (
            "Saccharomyces_cerevisiae_WT",
            "Saccharomyces_cerevisiae_IUPD",
            "Escherichia_coli_WT",
            "Escherichia_coli_IUPD",
            "Yarrowia_lipolytica_WT",
            "Yarrowia_lipolytica_IUPD",
        )
        models = {}
        hashes = {}
        for name in names:
            path = folder / (name + "_1_input.xlsx")
            make_input(path, "B[c]")
            models[name] = {
                "input_file": str(path),
                "input_sha256": digest(path),
                "carbon_sources": ["A[e]"],
            }
            hashes[path] = digest(path)
        # Unrelated files must not enter a plan-driven analysis.
        make_input(folder / "Old_Culture_1_input.xlsx", "B[c]")
        plan = {
            "version": 3,
            "models": models,
            "output_dir": str(folder),
        }
        config = folder / "analysis_plan.json"
        config.write_text(json.dumps(plan), encoding="utf-8")
        ms = mapping.settings()
        ms["analysis_plan"] = str(config)
        ts = tracing.settings()
        ts["analysis_plan"] = str(config)
        # Execution settings are supplied by the Python entry points, not MATLAB.
        ms.update(
            mappers=["LocalMapper"],
            batch_size=7,
            max_reaction_atoms=500,
            max_mapper_processes=None,
            mapper_memory_gb=1,
            python_overrides={"LocalMapper": sys.executable},
            structure_options={
                "lookup": False,
                "manual_file": str(folder / "absent.csv"),
                "cache_file": str(folder / "cache.json"),
                "timeout": 1,
                "min_interval": 0,
            },
        )
        ts.update(
            rounds=2,
            max_carbon_processes=1,
            carbon_memory_gb=None,
            tracing_python=sys.executable,
        )
        resolved, jobs = load_analysis_plan(config, "map", ms)
        assert list(jobs) == list(names) and resolved["max_mapper_processes"] is None
        assert resolved["mappers"] == ["LocalMapper"]
        assert resolved["batch_size"] == 7 and resolved["max_reaction_atoms"] == 500
        assert resolved["structure_options"] == ms["structure_options"]
        assert resolved["python_overrides"] == ms["python_overrides"]
        assert "rounds" not in resolved
        resolved_trace, trace_jobs = load_analysis_plan(config, "trace", ts)
        assert resolved_trace["rounds"] == 2
        with patch.object(tracing, "ROUNDS", 6):
            round_settings, _ = load_analysis_plan(config, "trace", tracing.settings())
        assert round_settings["rounds"] == 6
        assert resolved_trace["max_carbon_processes"] == 1
        assert resolved_trace["tracing_python"] == sys.executable
        assert all(job["mapping_file"].parent == folder for job in trace_jobs.values())
        # Reject the old all-settings plan instead of silently restoring its overrides.
        legacy = dict(plan, version=2, rounds=10)
        config.write_text(json.dumps(legacy), encoding="utf-8")
        try:
            load_analysis_plan(config, "map", ms)
        except ValueError as error:
            assert "rerun step1" in str(error)
        else:
            raise AssertionError("Legacy settings plan accepted")
        config.write_text(json.dumps(plan), encoding="utf-8")
        # No mapper installation is needed for this deterministic stage check.
        with patch(
            "carbon_flow_tools.cli.build_mapper_configs",
            return_value={
                "LocalMapper": {
                    "python": sys.executable,
                    "module": "scripts.test_parallel_carbon",
                    "options": {"test": True},
                }
            },
        ):
            assert main("map", ms, []) == 0
        assert not (folder / "Old_Culture_2_atom_mapping.xlsx").exists()
        mapping_hashes = {
            folder / (name + "_2_atom_mapping.xlsx"): digest(
                folder / (name + "_2_atom_mapping.xlsx")
            )
            for name in names
        }
        with patch(
            "carbon_flow_tools.cli.build_mapper_configs",
            side_effect=AssertionError("Trace invoked mappers"),
        ):
            assert main("trace", ts, []) == 0
        for name in names:
            file = folder / (name + "_3_carbon_flow.xlsx")
            with pd.ExcelFile(file, engine="openpyxl") as excel:
                assert set(excel.sheet_names) == set(
                    [
                        "Carbon_Flow_R1",
                        "Carbon_Flow_R2",
                        "Summary",
                    ]
                )
                flow = pd.read_excel(excel, sheet_name="Carbon_Flow_R1")
                assert set(flow.Production) == {"A[c]", "B[c]"}
            assert file.is_file()
        assert all(digest(path) == value for path, value in hashes.items())
        assert all(digest(path) == value for path, value in mapping_hashes.items())
        with patch("carbon_flow_tools.pipeline.run_carbon_flow") as propagate:
            assert main("trace", ts, ["Escherichia_coli_IUPD", "--rounds", "3"]) == 0
            assert list(propagate.call_args[0][0]) == ["Escherichia_coli_IUPD"]
            assert propagate.call_args[0][1] == 3
        plan["models"]["Saccharomyces_cerevisiae_WT"]["input_sha256"] = "0" * 64
        config.write_text(json.dumps(plan), encoding="utf-8")
        with patch("carbon_flow_tools.pipeline.run_carbon_flow") as propagate:
            try:
                main("trace", ts, ["Saccharomyces_cerevisiae_WT"])
            except ValueError as error:
                assert "rerun step2" in str(error)
            else:
                raise AssertionError("Stale same-name mapping accepted")
            assert not propagate.called
        plan["models"]["Saccharomyces_cerevisiae_WT"]["input_sha256"] = models[
            "Saccharomyces_cerevisiae_WT"
        ]["input_sha256"]
        plan["models"]["Saccharomyces_cerevisiae_WT"]["carbon_sources"] = []
        config.write_text(json.dumps(plan), encoding="utf-8")
        try:
            load_analysis_plan(config, "trace", ts)
        except ValueError as error:
            assert "carbon sources" in str(error)
        else:
            raise AssertionError("Empty sources accepted")
        plan["models"]["Saccharomyces_cerevisiae_WT"]["carbon_sources"] = ["A[e]"]
        for invalid in (0, 1.5, True):
            ts["rounds"] = invalid
            config.write_text(json.dumps(plan), encoding="utf-8")
            try:
                load_analysis_plan(config, "trace", ts)
            except ValueError as error:
                assert "positive integer" in str(error)
            else:
                raise AssertionError("Invalid propagation rounds accepted")
        config.unlink()
        try:
            load_analysis_plan(config, "map", ms)
        except FileNotFoundError as error:
            assert "Run step1" in str(error)
        else:
            raise AssertionError("Missing plan accepted")
        # The optional example registers existing Excel inputs without touching them.
        existing = folder / (names[0] + "_1_input.xlsx")
        example_inputs = {
            "MyCulture": {"input_file": existing.name, "carbon_sources": ["A[e]"]}
        }
        with patch.object(external, "ROOT", folder), patch.object(
            external, "INPUTS", example_inputs
        ), patch("carbon_flow_tools.cli.ROOT", folder):
            external.main()
            example_plan = json.loads((folder / "data/analysis_plan.json").read_text())
            assert set(example_plan) == {"version", "models", "output_dir"}
            assert example_plan["models"]["MyCulture"]["input_sha256"] == digest(
                existing
            )
            unchanged_plan = digest(folder / "data/analysis_plan.json")
            example_inputs["MyCulture"]["carbon_sources"] = []
            try:
                external.main()
            except ValueError as error:
                assert "carbon sources" in str(error)
            else:
                raise AssertionError("Invalid external sources accepted")
            assert digest(folder / "data/analysis_plan.json") == unchanged_plan
            assert digest(existing) == hashes[existing]
    print(
        "PASS: six planned models, both independent stages, stale-file exclusion and immutable inputs"
    )


if __name__ == "__main__":
    run()
