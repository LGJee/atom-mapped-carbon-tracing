"""Run one adapter in its own environment; no Excel or carbon propagation."""

import importlib
import json
import sys
from pathlib import Path


sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from carbon_flow_tools.console import Progress, detail, verbose, library_output


def checked_batch(adapter, model, batch):
    results = adapter.map_batch(model, batch)
    if not isinstance(results, list) or len(results) != len(batch):
        raise RuntimeError("Adapter result count differs from input count")
    for result in results:
        if not isinstance(result, dict) or not isinstance(
            result.get("mapped_rxn", ""), str
        ):
            raise RuntimeError(
                "Each result must be a dictionary with a string mapped_rxn"
            )
    json.dumps(results)
    return results


def main():
    module, request, response, batch_size, strain, name = sys.argv[1:]
    payload = json.loads(Path(request).read_text(encoding="utf-8"))
    reactions = payload["reactions"]
    size = int(batch_size)
    if size < 1:
        raise ValueError("Batch size must be a positive integer")
    label = "{} / {}".format(strain, name)
    unit = "reaction" if len(reactions) == 1 else "reactions"
    print("{}: starting ({} {}).".format(label, len(reactions), unit), flush=True)
    log_file = payload.get("log_file", str(Path(response).with_suffix(".log")))
    Path(log_file).parent.mkdir(parents=True, exist_ok=True)
    with Path(log_file).open("w", encoding="utf-8") as log:
        log.write("{}: latest mapper run\n".format(label))
    options = payload.get("options", {})
    with library_output(log_file):
        adapter = importlib.import_module(module)
        model = (
            adapter.load_model(dict(options, reactions=reactions, progress_label=label))
            if options
            else adapter.load_model()
        )
    progress = Progress(label + " | mapping" if verbose() else label, len(reactions), unit)
    results, retries = [], 0
    for start in range(0, len(reactions), size):
        batch = reactions[start : start + size]
        detail(
            "[{} {} {}-{}/{}] Mapping atoms...".format(
                strain, name, start + 1, start + len(batch), len(reactions)
            ),
            flush=True,
        )
        try:
            with library_output(log_file):
                mapped = checked_batch(adapter, model, batch)
        except Exception as error:
            retries += 1
            detail("Batch failed; retrying individually: " + str(error), flush=True)
            mapped = []
            for reaction in batch:
                try:
                    with library_output(log_file):
                        mapped.append(checked_batch(adapter, model, [reaction])[0])
                except Exception as exc:
                    mapped.append({"error": str(exc)})
        results.extend(mapped)
        progress.update(len(results), force=len(results) == len(reactions))
    if retries:
        detail("{} | {} batch(es) retried individually; per-reaction outcomes are saved."
              .format(label, retries), flush=True)
    with Path(log_file).open("a", encoding="utf-8") as log:
        for index, result in enumerate(results, 1):
            if result.get("error") or not result.get("mapped_rxn"):
                log.write("Reaction {}: {}\n".format(index, result.get("error", "No prediction")))
    flagged = sum(bool(result.get("error")) or not result.get("mapped_rxn") for result in results)
    if flagged:
        print("{}: {} {} checking; details are in the log.".format(
            label, flagged, "reaction needs" if flagged == 1 else "reactions need"), flush=True)
    Path(response).write_text(json.dumps(results), encoding="utf-8")


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        if verbose():
            raise
        print("Mapper worker failed: " + str(error), file=sys.stderr, flush=True)
        sys.exit(1)
