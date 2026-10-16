"""Check the preparation-to-mapping plan, six-model isolation and both stages."""

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
            "Saccharomyces_cerevisiae_non_IUPD",
            "Saccharomyces_cerevisiae_IUPD",
            "Escherichia_coli_non_IUPD",
            "Escherichia_coli_IUPD",
            "Yarrowia_lipolytica_non_IUPD",
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
            model_folder = folder / "data" / name
            model_folder.mkdir(parents=True)
            source = model_folder / "source_model.json"
            source.write_text(
                json.dumps(
                    dict(
                        metabolites=[
                            dict(
                                id="s_0",
                                name="A",
                                compartment="e",
                                annotation=dict(smiles="CCO"),
                            )
                        ]
                    )
                ),
                encoding="utf-8",
            )
            models[name]["source_model"] = dict(path=str(source), sha256=digest(source))
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
        # Execution settings are supplied by the Python entry points, not step 1.
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
        assert all(
            job["source_model"] == models[name]["source_model"]
            for name, job in jobs.items()
        )
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
            file = folder / (name + "_3_carbon_flow_Enumeration.xlsx")
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
        plan["models"]["Saccharomyces_cerevisiae_non_IUPD"]["input_sha256"] = "0" * 64
        config.write_text(json.dumps(plan), encoding="utf-8")
        with patch("carbon_flow_tools.pipeline.run_carbon_flow") as propagate:
            try:
                main("trace", ts, ["Saccharomyces_cerevisiae_non_IUPD"])
            except ValueError as error:
                assert "rerun step2" in str(error)
            else:
                raise AssertionError("Stale same-name mapping accepted")
            assert not propagate.called
        plan["models"]["Saccharomyces_cerevisiae_non_IUPD"]["input_sha256"] = models[
            "Saccharomyces_cerevisiae_non_IUPD"
        ]["input_sha256"]
        plan["models"]["Saccharomyces_cerevisiae_non_IUPD"]["carbon_sources"] = []
        config.write_text(json.dumps(plan), encoding="utf-8")
        try:
            load_analysis_plan(config, "trace", ts)
        except ValueError as error:
            assert "carbon sources" in str(error)
        else:
            raise AssertionError("Empty sources accepted")
        plan["models"]["Saccharomyces_cerevisiae_non_IUPD"]["carbon_sources"] = ["A[e]"]
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


def check_independent_cases():
    from openpyxl import load_workbook

    with tempfile.TemporaryDirectory() as directory:
        root = Path(directory)
        default_plan = root / "data/analysis_plan.json"
        default_plan.parent.mkdir()
        default_plan.write_text('{"version":3,"default":true}')
        protected = default_plan.read_bytes()
        manual = root / "manual.csv"
        manual.write_text("Model,MetaboliteID,MetaboliteName,SMILES\n")
        mapping_settings = mapping.settings()
        mapping_settings.update(
            analysis_plan=str(default_plan), mappers=["LocalMapper"],
            max_mapper_processes=1, max_carbon_processes=1,
            structure_options=dict(lookup=False, manual_file=str(manual),
                                   cache_file=str(root / "cache.json"), timeout=1, min_interval=0),
        )
        trace_settings = tracing.settings()
        trace_settings.update(analysis_plan=str(default_plan), tracing_python=sys.executable,
                              rounds=2, max_carbon_processes=1)
        configs = {"LocalMapper": dict(python=sys.executable,
                                       module="scripts.test_parallel_carbon", options={"test": True})}
        for name in ("CaseOne", "CaseTwo"):
            case = root / "results/results_by_AI" / name
            case.mkdir(parents=True)
            source = case / "input.xlsx"
            make_input(source, "B[c]")
            workbook = load_workbook(source)
            workbook["Metabolites"].append(["unused", "X[e]", "CCO", ""])
            workbook.save(source)
            workbook.close()
            output = case
            plan = case / "analysis_plan.json"
            plan.write_text(json.dumps(dict(version=3, output_dir=str(output), models={
                name: dict(input_file=str(source), input_sha256=digest(source),
                           carbon_sources=["A[e]", "X[e]"])
            })))
            with patch("carbon_flow_tools.cli.ROOT", root), patch(
                    "carbon_flow_tools.cli.build_mapper_configs", return_value=configs):
                assert main("map", mapping_settings, ["--plan", str(plan)]) == 0
                assert main("trace", trace_settings, ["--plan", str(plan)]) == 0
            result = output / (name + "_3_carbon_flow_Enumeration.xlsx")
            summary = pd.read_excel(result, sheet_name="Summary", engine="openpyxl")
            assert summary["X[e]_LabeledCarbonMetabolites"].eq(0).all()
            assert summary.TotalCarbonMetabolites.eq(3).all()
            assert summary["A[e]_LabeledCarbonMetabolites"].eq(2).all()
            report_path = output / "carbon_flow_run.json"
            report = json.loads(report_path.read_text())
            assert report["models"][0]["stop_reason"] == "converged"
            assert report["models"][0]["completed_rounds"] == 2
            assert report["models"][0]["mapping_sha256"] == digest(output / (name + "_2_atom_mapping.xlsx"))
            previous_result = result.read_bytes()
            first_flow = pd.read_excel(result, sheet_name="Carbon_Flow_R1", engine="openpyxl")
            with patch("carbon_flow_tools.cli.ROOT", root):
                assert main("trace", dict(trace_settings, max_excel_size_mb=.001),
                            ["--plan", str(plan)]) == 0
            stopped = json.loads(report_path.read_text())["models"][0]
            assert stopped["stop_reason"] == "file_size_limit"
            assert stopped["completed_rounds"] == 0
            assert stopped["output_file"] is None
            assert result.read_bytes() == previous_result
            assert result.stat().st_size >= 1000
            pd.testing.assert_frame_equal(
                pd.read_excel(result, sheet_name="Carbon_Flow_R1", engine="openpyxl"), first_flow)
            pd.testing.assert_frame_equal(
                pd.read_excel(result, sheet_name="Summary", engine="openpyxl"), summary)
            assert default_plan.read_bytes() == protected
            # The environment relaunch must retain the explicit plan and unified entry.
            with patch("carbon_flow_tools.cli.ROOT", root), patch(
                    "carbon_flow_tools.cli.sys.executable", "/unrelated/python"), patch(
                    "carbon_flow_tools.cli.subprocess.run") as launch:
                launch.return_value.returncode = 0
                assert main("trace", trace_settings, ["--plan", str(plan)]) == 0
                assert Path(launch.call_args[0][0][1]).name == "step3_carbon_flow.py"
                assert launch.call_args[0][0][-2:] == ["--plan", str(plan)]
        for stage, settings in (("map", mapping_settings), ("trace", trace_settings)):
            try:
                main(stage, settings, ["--example", "--plan", str(default_plan)])
            except SystemExit as error:
                assert error.code == 2
            else:
                raise AssertionError("Example silently ignored an explicit case plan")
        assert default_plan.read_bytes() == protected
        assert not list((root / "results").glob("*.xlsx"))
    print("PASS: isolated case plans, explicit --plan, unused source zero counts and AI relaunch")


if __name__ == "__main__":
    run()
    check_independent_cases()
