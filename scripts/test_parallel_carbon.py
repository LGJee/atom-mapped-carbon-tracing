"""Regression checks for preserved scientific behavior."""

import hashlib
import os
import sys
import tempfile
import time
from pathlib import Path
from unittest.mock import patch

import pandas as pd
from openpyxl import Workbook

from carbon_flow_tools.pipeline import run_analyses, finish_analysis
from carbon_flow_tools.config import parallel_capacity


def load_model(options):
    return None


def map_batch(model, smiles):
    assert all(text == "CCO>>CCO" for text in smiles)
    return [
        {"mapped_rxn": "[CH3:1][CH2:2][OH:3]>>[CH3:1][CH2:2][OH:3]", "confidence": 0.8}
        for _ in smiles
    ]


def synchronized_finish(*args):
    strain, job = args[:2]
    events = Path(job["events"])
    (events / (str(os.getpid()) + ".started")).touch()
    deadline = time.monotonic() + 15
    while len(list(events.glob("*.started"))) < 2:
        if time.monotonic() > deadline:
            raise RuntimeError("The two strain propagation workers did not overlap")
        time.sleep(0.01)
    finish_analysis(*args)


def failed_finish(*args):
    if args[0] == "IUPD":
        raise RuntimeError("Test propagation worker failed")
    finish_analysis(*args)


def make_input(path, product):
    workbook = Workbook()
    sheet = workbook.active
    sheet.title = "Nonzero"
    sheet.append(["ReactionID", "ReactionName", "Equation", "Flux"])
    sheet.append(["r1", "transport", "A[e] -> A[c]", 1])
    sheet.append(["r2", "reaction", "A[c] -> " + product, 1])
    sheet = workbook.create_sheet("Metabolites")
    sheet.append(["MetaboliteID", "MetaboliteName", "SMILES_Yeast9", "SMILES_Manual"])
    for i, name in enumerate(("A[e]", "A[c]", product)):
        sheet.append(["s_" + str(i), name, "CCO", ""])
    workbook.save(path)
    workbook.close()


def main():
    with patch(
        "carbon_flow_tools.config.os.sched_getaffinity", return_value=set(range(8))
    ), patch(
        "carbon_flow_tools.config.Path.read_text",
        return_value="MemTotal: 16777216 kB\nMemAvailable: 12582912 kB\n",
    ):
        assert parallel_capacity(
            2, memory_per_process_gb=None, task_kind="carbon propagation"
        ) == (2, 4)
    with patch(
        "carbon_flow_tools.config.os.sched_getaffinity", return_value=set(range(8))
    ), patch(
        "carbon_flow_tools.config.Path.read_text",
        return_value="MemTotal: 67108864 kB\nMemAvailable: 8388608 kB\n",
    ):
        assert (
            parallel_capacity(
                2, memory_per_process_gb=None, task_kind="carbon propagation"
            )[0]
            == 1
        )
    with patch(
        "carbon_flow_tools.config.os.sched_getaffinity", return_value={0}
    ), patch(
        "carbon_flow_tools.config.Path.read_text",
        return_value="MemTotal: 16777216 kB\nMemAvailable: 12582912 kB\n",
    ):
        assert parallel_capacity(
            2, memory_per_process_gb=None, task_kind="carbon propagation"
        ) == (1, 1)
    with patch("carbon_flow_tools.config.Path.read_text", side_effect=OSError):
        assert (
            parallel_capacity(
                2, memory_per_process_gb=None, task_kind="carbon propagation"
            )[0]
            == 1
        )
    with tempfile.TemporaryDirectory() as folder:
        root = Path(folder)
        events = root / "events"
        events.mkdir()
        jobs = {}
        for strain, product in [("non_IUPD", "B[c]"), ("IUPD", "D[c]")]:
            source = root / (strain + ".xlsx")
            make_input(source, product)
            jobs[strain] = dict(
                input_file=source,
                reaction_sheet="Nonzero",
                metabolite_sheet="Metabolites",
                carbon_sources=["A[e]"],
                events=str(events),
            )
        hashes = {
            j["input_file"]: hashlib.sha256(j["input_file"].read_bytes()).hexdigest()
            for j in jobs.values()
        }
        settings = {
            "ExampleMapper": {
                "module": "scripts.test_parallel_carbon",
                "python": sys.executable,
                "options": {"test": True},
            }
        }
        args = (jobs, 2, 1, 800)
        run_analyses(
            *args,
            root / "serial",
            ["ExampleMapper"],
            settings,
            max_carbon_processes=1,
            carbon_memory_gb=0.25,
        )
        with patch(
            "carbon_flow_tools.pipeline.finish_analysis", new=synchronized_finish
        ):
            run_analyses(
                *args,
                root / "parallel",
                ["ExampleMapper"],
                settings,
                max_carbon_processes=2,
                carbon_memory_gb=0.25,
            )
        assert len(list(events.glob("*.started"))) == 2
        assert all(int(path.stem) != os.getpid() for path in events.glob("*.started"))
        for strain, product in [("non_IUPD", "B[c]"), ("IUPD", "D[c]")]:
            name = strain + "_3_carbon_flow.xlsx"
            with pd.ExcelFile(
                root / "serial" / name, engine="openpyxl"
            ) as serial, pd.ExcelFile(
                root / "parallel" / name, engine="openpyxl"
            ) as parallel:
                assert serial.sheet_names == parallel.sheet_names
                for sheet in serial.sheet_names:
                    pd.testing.assert_frame_equal(
                        pd.read_excel(serial, sheet_name=sheet),
                        pd.read_excel(parallel, sheet_name=sheet),
                    )
                flow = pd.read_excel(parallel, sheet_name="Carbon_Flow_R1")
                assert set(flow.Production) == {"A[c]", product}
                assert flow.ProductionID.tolist() == [1, 2]
                assert flow.CarbonLabels.eq("C1,C2").all()
                assert len(pd.read_excel(parallel, sheet_name="Summary")) == 2

        failed_output = root / "failed"
        failed_output.mkdir()
        protected = failed_output / "IUPD_3_carbon_flow.xlsx"
        protected.write_bytes(b"previous output")
        with patch("carbon_flow_tools.pipeline.finish_analysis", new=failed_finish):
            try:
                run_analyses(
                    *args,
                    failed_output,
                    ["ExampleMapper"],
                    settings,
                    refresh=True,
                    max_carbon_processes=2,
                    carbon_memory_gb=0.25,
                )
            except RuntimeError as error:
                assert str(error) == "Test propagation worker failed"
            else:
                raise AssertionError(
                    "Propagation worker error was not returned to parent"
                )
        assert protected.read_bytes() == b"previous output"
        for path, before in hashes.items():
            assert hashlib.sha256(path.read_bytes()).hexdigest() == before
    print(
        "PASS: two simultaneous carbon processes, isolated states/IDs, serial equivalence, "
        "round order, exceptions and input preservation"
    )


if __name__ == "__main__":
    main()
