"""Check mapping-only output and propagation without any mapper calls."""

import hashlib
from pathlib import Path
import sys
import tempfile
from unittest.mock import patch

import pandas as pd
from openpyxl import load_workbook

from carbon_flow_tools.pipeline import run_analyses, read_input
from tests.test_parallel_carbon import make_input


def main():
    with tempfile.TemporaryDirectory() as directory:
        root = Path(directory)
        source = root / "input.xlsx"
        make_input(source, "B[c]")
        before = hashlib.sha256(source.read_bytes()).hexdigest()
        jobs = {
            "Example": dict(
                input_file=source,
                reaction_sheet="Nonzero",
                metabolite_sheet="Metabolites",
                carbon_sources=["A[e]"],
            )
        }
        configs = {
            "ExampleMapper": {
                "module": "tests.test_parallel_carbon",
                "python": sys.executable,
                "options": {"test": True},
            }
        }
        arguments = (jobs, 2, 1, 800)
        options = dict(max_mapper_processes=1, max_carbon_processes=1)
        run_analyses(
            *arguments, root / "combined", ["ExampleMapper"], configs, **options
        )
        run_analyses(
            *arguments,
            root / "split",
            ["ExampleMapper"],
            configs,
            stage="map",
            **options,
        )
        output = root / "split/Example_integrated_carbon_tracing.xlsx"
        with pd.ExcelFile(output, engine="openpyxl") as excel:
            assert "Final_Mapping" in excel.sheet_names
            assert not any(name.startswith("Carbon_Flow") for name in excel.sheet_names)
            for name in ("Mapping_Selection", "Final_Mapping"):
                table = pd.read_excel(excel, sheet_name=name)
                assert (
                    table.columns[table.columns.get_loc("ReviewRequired") + 1]
                    == "AIGenerated"
                )
                assert table.AIGenerated.eq(False).all()
        with patch(
            "carbon_flow_tools.pipeline.run_mapper",
            side_effect=AssertionError("Mapper called during trace"),
        ):
            run_analyses(
                *arguments,
                root / "split",
                ["ExampleMapper"],
                {},
                stage="trace",
                **options,
            )
        with pd.ExcelFile(output, engine="openpyxl") as split, pd.ExcelFile(
            root / "combined/Example_integrated_carbon_tracing.xlsx", engine="openpyxl"
        ) as combined:
            for name in (
                "Mapping_Selection",
                "Final_Mapping",
                "Carbon_Flow_R1",
                "Carbon_Flow_R2",
                "Summary",
            ):
                pd.testing.assert_frame_equal(
                    pd.read_excel(split, sheet_name=name),
                    pd.read_excel(combined, sheet_name=name),
                )
        assert hashlib.sha256(source.read_bytes()).hexdigest() == before
        run_analyses(
            jobs,
            1,
            1,
            800,
            root / "split",
            ["ExampleMapper"],
            {},
            stage="trace",
            **options,
        )
        with pd.ExcelFile(output, engine="openpyxl") as excel:
            assert "Carbon_Flow_R2" not in excel.sheet_names
        for kind in (
            "duplicate_id",
            "zero_flux",
            "nonfinite_flux",
            "conflicting_structure",
        ):
            invalid = root / (kind + ".xlsx")
            workbook = load_workbook(source)
            if kind == "duplicate_id":
                workbook["Nonzero"]["A3"] = "r1"
            elif kind in ("zero_flux", "nonfinite_flux"):
                workbook["Nonzero"]["D2"] = 0 if kind == "zero_flux" else float("nan")
            else:
                workbook["Metabolites"]["B4"] = "A[c]"
                workbook["Metabolites"]["C4"] = "CC"
            workbook.save(invalid)
            workbook.close()
            try:
                read_input(invalid, "Nonzero", "Metabolites")
            except ValueError:
                pass
            else:
                raise AssertionError("Invalid input accepted: " + kind)
    print(
        "PASS: separate stages, no prediction during trace, flags, fewer-round cleanup, input preservation"
    )


if __name__ == "__main__":
    main()
