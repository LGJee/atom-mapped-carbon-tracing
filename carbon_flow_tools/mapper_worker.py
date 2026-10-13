"""Run one adapter in its own environment; no Excel or carbon propagation."""

import importlib
import json
import sys
from pathlib import Path


sys.path.insert(0, str(Path(__file__).resolve().parent.parent))


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
    adapter = importlib.import_module(module)
    options = payload.get("options", {})

    model = (
        adapter.load_model(dict(options, reactions=reactions))
        if options
        else adapter.load_model()
    )
    results = []
    for start in range(0, len(reactions), size):
        batch = reactions[start : start + size]
        print(
            "[{} {} {}-{}/{}] Mapping atoms...".format(
                strain, name, start + 1, start + len(batch), len(reactions)
            ),
            flush=True,
        )
        try:
            mapped = checked_batch(adapter, model, batch)
        except Exception as error:
            print("Batch failed; retrying individually: " + str(error), flush=True)
            mapped = []
            for reaction in batch:
                try:
                    mapped.append(checked_batch(adapter, model, [reaction])[0])
                except Exception as exc:
                    mapped.append({"error": str(exc)})
        results.extend(mapped)
    Path(response).write_text(json.dumps(results), encoding="utf-8")


if __name__ == "__main__":
    main()
