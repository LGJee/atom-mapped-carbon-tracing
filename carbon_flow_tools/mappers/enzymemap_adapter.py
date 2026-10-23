"""Match EnzymeMap records, then use official single-step rules."""

import csv
from functools import lru_cache
from pathlib import Path

from rdkit import Chem, RDLogger

from ..console import Progress, detail


@lru_cache(maxsize=50000)
def structure(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("Invalid SMILES: " + smiles)
    for atom in mol.GetAtoms():
        atom.SetAtomMapNum(0)

    return Chem.MolToSmiles(mol, isomericSmiles=False)


def key(reaction):
    return tuple(
        tuple(sorted(structure(s) for s in side.split(".")))
        for side in reaction.split(">>")
    )


def transfer(input_reaction, mapped_reaction, normalizer=None):
    """Transfer database map IDs while retaining input charge and stereochemistry."""
    output = []
    for given, stored in zip(input_reaction.split(">>"), mapped_reaction.split(">>")):
        available = stored.split(".")
        parts = []
        for smiles in given.split("."):
            base = Chem.MolFromSmiles(smiles)
            query = normalizer(Chem.Mol(base)) if normalizer else base
            if query.GetNumAtoms() != base.GetNumAtoms():
                raise ValueError(
                    "Normalization changed the atom count; cannot restore the input."
                )
            for index, candidate in enumerate(available):
                if structure(Chem.MolToSmiles(query)) != structure(candidate):
                    continue
                mapped = Chem.MolFromSmiles(candidate)
                match = mapped.GetSubstructMatch(query, useChirality=True)
                if len(match) != base.GetNumAtoms():
                    continue
                for atom in base.GetAtoms():
                    atom.SetAtomMapNum(
                        mapped.GetAtomWithIdx(match[atom.GetIdx()]).GetAtomMapNum()
                    )
                parts.append(Chem.MolToSmiles(base))
                available.pop(index)
                break
            else:
                if (
                    normalizer
                    and base.GetNumAtoms() == 1
                    and base.GetAtomWithIdx(0).GetAtomicNum() == 1
                ):
                    parts.append(Chem.MolToSmiles(base))
                    continue
                raise ValueError(
                    "Database molecules do not match input structures or specified stereochemistry."
                )
        if available and not (
            normalizer
            and all(
                Chem.MolFromSmiles(s).GetNumAtoms() == 1
                and Chem.MolFromSmiles(s).GetAtomWithIdx(0).GetAtomicNum() == 1
                for s in available
            )
        ):
            raise ValueError("Database stoichiometry does not match input.")
        output.append(".".join(parts))
    if len(output) != 2:
        raise ValueError("Reaction must contain >>.")
    return ">>".join(output)


def identity(reaction):
    """Map identical structures without claiming a database prediction."""
    left, right = key(reaction)
    if left != right:
        return None
    mapped, number = [], 1
    for smiles in reaction.split(">>")[0].split("."):
        mol = Chem.MolFromSmiles(smiles)
        for atom in mol.GetAtoms():
            atom.SetAtomMapNum(number)
            number += 1
        mapped.append(Chem.MolToSmiles(mol))
    side = ".".join(mapped)
    try:
        return transfer(reaction, side + ">>" + side)
    except ValueError:
        return None


def load_model(options):
    path = Path(options["dataset"])
    if not path.is_file():
        raise FileNotFoundError("EnzymeMap dataset not found: " + str(path))
    requested = {key(r) for r in options["reactions"]}
    matches = {}
    RDLogger.DisableLog("rdApp.warning")
    progress = Progress(options.get("progress_label", "EnzymeMap") + " | scanning database",
                        unit="records")
    with path.open(encoding="utf-8-sig", newline="") as source:
        rows = csv.DictReader(source)
        if not {"mapped", "unmapped"} <= set(rows.fieldnames or []):
            raise ValueError("EnzymeMap CSV requires mapped and unmapped columns.")
        for number, row in enumerate(rows, 2):
            if number % 10000 == 0:
                progress.update(number - 1)
                detail(
                    "[EnzymeMap lookup] scanned {} records; {} reactions matched.".format(
                        number - 1, len(matches)
                    ),
                    flush=True,
                )
            try:
                reaction_key = key(row["unmapped"])
                candidates = [
                    (reaction_key, row["mapped"], "forward"),
                    (
                        reaction_key[::-1],
                        ">>".join(row["mapped"].split(">>")[::-1]),
                        "reverse",
                    ),
                ]
                for candidate, mapped, direction in candidates:
                    if candidate not in requested:
                        continue
                    quality = float(row.get("quality") or 0)
                    old = matches.get(candidate)
                    if old is None or quality > old["quality"]:
                        matches[candidate] = {
                            "mapped": mapped,
                            "quality": quality,
                            "template": "EnzymeMap CSV line {}; {}; rxn_idx={}; rule_id={}; {}".format(
                                number,
                                direction,
                                row.get("rxn_idx", ""),
                                row.get("rule_id", ""),
                                row.get("source", ""),
                            ),
                        }
            except (ValueError, TypeError):
                continue
    detail(
        "[EnzymeMap lookup] {} requested reactions; {} database matches.".format(
            len(requested), len(matches)
        ),
        flush=True,
    )
    predictor = None
    if options.get("rules"):
        from .enzymemap_rules import RuleMapper

        predictor = RuleMapper(options)
    return {"matches": matches, "predictor": predictor,
            "options": {k: v for k, v in options.items()
                        if k not in ("reactions", "progress_label")}}


def map_batch(model, reaction_smiles):
    results = []
    for reaction in reaction_smiles:
        same = identity(reaction)
        if same:
            results.append(
                {
                    "mapped_rxn": same,
                    "template": "Identity transport; not an EnzymeMap prediction",
                }
            )
            continue
        match = model["matches"].get(key(reaction))
        error = (
            "No EnzymeMap record with matching structures, charges and stoichiometry."
        )
        if match is not None:
            try:
                results.append(
                    {
                        "mapped_rxn": transfer(reaction, match["mapped"]),
                        "confidence": match["quality"],
                        "template": match["template"],
                    }
                )
                continue
            except ValueError as exc:
                error = str(exc)
        if model["predictor"] is not None:
            results.append(model["predictor"].map(reaction))
        else:
            results.append({"error": error})
    return results


# A conservative estimate from isolated rule-worker measurements; complex reactions
# can exceed it. Keep the parent plus all children within the outer job's budget.
RULE_WORKER_MEMORY_GB = 0.5
_rule_predictor = None
_rule_log = None
_rule_error = None


def _initialize_rules(options, log_file):
    global _rule_predictor, _rule_log, _rule_error
    from ..console import library_output
    _rule_log = log_file
    _rule_error = None
    with library_output(log_file):
        try:
            from .enzymemap_rules import RuleMapper
            from ..config import MemoryGuard
            MemoryGuard().check("loading EnzymeMap rules", force=True, extra_bytes=512 * 1024 ** 2)
            _rule_predictor = RuleMapper(options)
        except Exception as error:
            # Keep the worker alive so the scheduler receives initialization failures.
            _rule_error = str(error)


def _predict_rule(reaction):
    from ..console import library_output
    from ..config import MemoryGuard, MemoryPressureError
    if _rule_error is not None:
        raise RuntimeError("EnzymeMap rule worker initialization failed: " + _rule_error)
    guard = MemoryGuard()
    guard.check("starting EnzymeMap rule prediction", force=True)
    with library_output(_rule_log):
        try:
            result = _rule_predictor.map(reaction)
            guard.check("finishing EnzymeMap rule prediction", force=True)
            return result
        except (MemoryPressureError, MemoryError):
            raise
        except Exception as error:
            return {"error": str(error)}


def rule_capacity(count, cpu_threads, memory_gb):
    from ..config import available_memory, parallel_capacity, task_memory_budget
    # One parent owns the database matches; only children execute rule searches.
    try:
        available, total = available_memory()
    except (OSError, ValueError, KeyError):
        return 1
    budget = (available - max(2, available * 0.2, total * 0.1)
              if memory_gb is None else task_memory_budget(memory_gb, total, available))
    limit = min(cpu_threads, max(1, int((budget - 1) / RULE_WORKER_MEMORY_GB)))
    return parallel_capacity(count, limit, RULE_WORKER_MEMORY_GB,
                             "EnzymeMap rule prediction")[0]


def map_reactions(model, reactions, cpu_threads, memory_gb, log_file, progress):
    """Resolve database hits once, then parallelize independent rule searches.

    Rule order and protonation attempts inside each reaction are unchanged.
    Results remain in input order regardless of process completion order.
    """
    from ..console import library_output
    from ..task_scheduler import run_tasks
    from ..config import MemoryGuard, MemoryPressureError
    guard = MemoryGuard()
    guard.check("matching EnzymeMap database", force=True)
    results, pending, completed = [], [], [0]
    label = progress.label
    progress.begin(label + " | database matching", len(reactions), estimate=False)

    def record():
        completed[0] += 1
        progress.update(completed[0], force=completed[0] == len(pending))

    lookup_only = {"matches": model["matches"], "predictor": None}
    with library_output(log_file):
        for index, reaction in enumerate(reactions):
            guard.check("matching EnzymeMap database")
            try:
                value = map_batch(lookup_only, [reaction])[0]
            except (MemoryPressureError, MemoryError):
                raise
            except Exception as error:
                # Match the ordinary worker's per-reaction failure isolation.
                results.append({"error": str(error)})
                progress.update(index + 1)
                continue
            results.append(value)
            if not value.get("mapped_rxn") and model["predictor"] is not None:
                pending.append(index)
            progress.update(index + 1)
    progress.update(len(reactions), force=True)
    if not pending:
        return results
    progress.begin(label + " | rule prediction", len(pending), estimate=True)
    workers = rule_capacity(len(pending), cpu_threads, memory_gb)
    print("EnzymeMap: {} unique reactions; {} need rule prediction; {} rule process(es)."
          .format(len(reactions), len(pending), workers), flush=True)
    if workers == 1:
        for index in pending:
            guard.check("starting serial EnzymeMap rule prediction", force=True)
            with library_output(log_file):
                try:
                    results[index] = model["predictor"].map(reactions[index])
                except (MemoryPressureError, MemoryError):
                    raise
                except Exception as error:
                    results[index] = {"error": str(error)}
            record()
    else:
        predicted = run_tasks(_predict_rule, [(reactions[i],) for i in pending],
            workers, RULE_WORKER_MEMORY_GB,
            weights=[len(reactions[i]) for i in pending], cpu_threads=1,
            initializer=_initialize_rules, initargs=(model["options"], log_file),
            on_result=record, on_wait=lambda: progress.update(completed[0]))
        for index, value in zip(pending, predicted):
            results[index] = value
    return results
