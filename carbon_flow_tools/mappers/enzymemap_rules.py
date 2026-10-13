"""Apply official EnzymeMap rules without changing supplied reaction participants."""

from pathlib import Path
import hashlib

import pandas as pd
from rdkit import Chem
from enzymemap import helpers_map, helpers_rdkit

from .enzymemap_adapter import transfer, key


class RuleMapper:
    def __init__(self, options):
        path = Path(options["rules"])
        if not path.is_file():
            raise FileNotFoundError(
                "Official EnzymeMap rules.pkl not found: " + str(path)
            )
        checksum = hashlib.sha256(path.read_bytes()).hexdigest()
        if checksum != options.get("rules_sha256"):
            raise ValueError(
                "Rule SHA256 differs from configuration; pickle was not loaded."
            )
        self.rules = pd.read_pickle(path)
        if not {"SMARTS", "Num_R", "Num_P", "Rxn"} <= set(self.rules.columns):
            raise ValueError("EnzymeMap rules are missing required columns.")
        self.normalize = options.get("normalize_protonation", True)
        if type(self.normalize) is not bool:
            raise ValueError("normalize_protonation must be True or False.")
        print(
            "[EnzymeMap rules] loaded {} official rules.".format(len(self.rules)),
            flush=True,
        )

    def map(self, reaction):
        # shortcut: use official single-step rules; aggregate cycles need explicit modelling.
        attempts = [(reaction, False)]
        if self.normalize:
            normalized = ">>".join(
                ".".join(
                    Chem.MolToSmiles(
                        helpers_rdkit.neutralize_atoms(Chem.MolFromSmiles(s))
                    )
                    for s in side.split(".")
                )
                for side in reaction.split(">>")
            )
            if key(normalized) != key(reaction):
                attempts.append((normalized, True))
        errors = []
        for smiles, neutralized in attempts:
            try:
                mapped, templates, ids, _ = helpers_map.map(
                    [smiles], self.rules, single=True
                )
                if not mapped:
                    continue
                tag = "EnzymeMap rule prediction; rule_id={}; protonation={}; SMARTS={}".format(
                    ids[0],
                    "neutralized then restored" if neutralized else "original",
                    templates[0] or "isomerase",
                )
                normalizer = (
                    helpers_rdkit.neutralize_atoms if neutralized else lambda mol: mol
                )
                try:
                    restored = transfer(reaction, mapped[0], normalizer=normalizer)
                    return {"mapped_rxn": restored, "template": tag}
                except ValueError as error:
                    errors.append(
                        {
                            "mapped_rxn": mapped[0],
                            "template": tag,
                            "error": "Rule prediction cannot restore input structures: "
                            + str(error),
                        }
                    )
            except Exception as error:
                errors.append(
                    {"error": "Official rule calculation failed: " + str(error)}
                )
        return (
            errors[0]
            if errors
            else {
                "error": "No database match; official single-step rules did not map the supplied products."
            }
        )
