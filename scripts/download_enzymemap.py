"""Download the pinned official EnzymeMap v2 dataset and verified rule file."""

import argparse
import gzip
import hashlib
from pathlib import Path
import shutil
from urllib.request import urlopen

COMMIT = "c91578856f11312d7941b290f8b773eb41b54973"
RULES_SHA256 = "21cdebccdda5ea4308555b6dd8dfb27413950725ace66e7b8cb7ed3dba2c1e1a"
BASE = "https://raw.githubusercontent.com/hesther/enzymemap/" + COMMIT + "/data/"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path(__file__).resolve().parents[1] / "data/enzymemap",
    )
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    for remote, name in (
        ("rules.pkl", "rules.pkl"),
        ("processed_reactions.csv.gz", "processed_reactions.csv"),
    ):
        destination = args.output_dir / name
        if destination.is_file():
            if (
                name == "rules.pkl"
                and hashlib.sha256(destination.read_bytes()).hexdigest() != RULES_SHA256
            ):
                raise ValueError(
                    "Existing rule file has a different SHA256; it was not overwritten."
                )
            print("Reusing " + str(destination))
            continue
        temporary = destination.with_suffix(destination.suffix + ".tmp")
        print("Downloading " + remote + " from commit " + COMMIT, flush=True)
        try:
            with urlopen(BASE + remote, timeout=60) as response:
                stream = (
                    gzip.GzipFile(fileobj=response)
                    if remote.endswith(".gz")
                    else response
                )
                with temporary.open("wb") as output:
                    shutil.copyfileobj(stream, output)
            if (
                name == "rules.pkl"
                and hashlib.sha256(temporary.read_bytes()).hexdigest() != RULES_SHA256
            ):
                raise ValueError("Downloaded rules failed SHA256 verification.")
            temporary.replace(destination)
        finally:
            if temporary.exists():
                temporary.unlink()
        print("Saved " + str(destination), flush=True)
    print(
        "Official data remain external to this repository; retain the upstream attribution."
    )


if __name__ == "__main__":
    main()
