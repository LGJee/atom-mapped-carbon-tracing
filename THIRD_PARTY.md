# External software, models and data

Install these components separately. Their implementations, weights and complete
databases are not included in this repository.

| Component | Source | Applicable upstream terms |
|---|---|---|
| COBRApy | [Version 0.30.0](https://github.com/opencobra/cobrapy/tree/0.30.0) | LGPL-2.0-or-later OR GPL-2.0-or-later per upstream package metadata; retain upstream notices and citations |
| GLPK / swiglpk | [GLPK](https://www.gnu.org/software/glpk/) and [swiglpk](https://github.com/biosustain/swiglpk) | Solver/binding distributions retain their upstream GPL terms and notices |
| Yeast9 | [v9.1.1](https://github.com/SysBioChalmers/yeast-GEM/tree/v9.1.1) | Retain the model's license and attribution |
| iYali5 / companion eciYali5 structures | [Pinned upstream version](https://github.com/bioengtaltech/eciYali5-GEM/tree/3b97d38d92327a099f3b0f6f9cd92e8aeab78c10) | CC BY 4.0; retain `data/yarrowia/LICENSE.iYali5.md`, pinned hashes and [publication](https://doi.org/10.1007/s00253-026-13791-4). Ordinary iYali5 supplies pFBA; eciYali5 supplies chemical structures only |
| LocalMapper | [Pinned commit](https://github.com/snu-micc/LocalMapper/tree/8ff6b3ff2109d834388d4e08ba69399ed3419cad) | LICENSE is **CC BY-NC-SA 4.0** with noncommercial/share-alike terms. Package metadata still says GPL v3; consult the upstream LICENSE/developer rather than assuming MIT |
| RXNMapper | [Official repository](https://github.com/rxn4chemistry/rxnmapper) | MIT package license; model/data terms remain applicable |
| EnzymeMap | [Pinned commit](https://github.com/hesther/enzymemap/tree/c91578856f11312d7941b290f8b773eb41b54973) | MIT package license; processed biochemical data retain their official provenance and applicable source terms |
| lxml | [Version 4.9.4](https://pypi.org/project/lxml/4.9.4/) | BSD-3-Clause; native XML backend for openpyxl; bundled libraries retain their upstream notices |
| PubChem | [Official PUG REST API](https://pubchem.ncbi.nlm.nih.gov/docs/pug-rest) | Retrieved structures are cached locally with CID provenance; retain database attribution and applicable source terms |
| rdchiral | [Pinned commit](https://github.com/hesther/rdchiral/tree/a96477c4b2e9c4e3ecbfa0c88d7c394e480bd874) | Retain upstream license and attribution |

Adapters call public interfaces without copying mapper implementations.
The environment files identify additional dependencies, whose distributions
supply their own license notices. Our MIT license does not relicense them.

The example is a synthetic test network. Original non-IUPD/IUPD study inputs and large
research outputs are not published in this code release.
