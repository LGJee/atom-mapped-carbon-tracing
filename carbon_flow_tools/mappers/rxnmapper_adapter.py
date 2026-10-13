"""Pass RXNMapper correspondences to the shared alignment layer."""


def load_model():
    from rxnmapper import RXNMapper

    return RXNMapper()


def map_batch(model, reaction_smiles):
    return model.get_attention_guided_atom_maps(reaction_smiles)
