"""Normalize LocalMapper output fields without changing its mapping."""


def load_model():
    from localmapper import localmapper

    return localmapper()


def map_batch(model, reaction_smiles):
    results = model.get_atom_map(reaction_smiles, return_dict=True)
    if isinstance(results, dict):
        results = [results]
    return [
        {
            "mapped_rxn": r.get("mapped_rxn", ""),
            "confidence": r.get("confident", ""),
            "template": r.get("template", ""),
        }
        for r in results
    ]
