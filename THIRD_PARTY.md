# External software, models and data

Install these components separately. Their implementations, weights and complete
databases are not included in this repository.

| Component | Source | Applicable upstream terms |
|---|---|---|
| COBRApy | [Version 0.30.0](https://github.com/opencobra/cobrapy/tree/0.30.0) | LGPL-2.0-or-later OR GPL-2.0-or-later per upstream package metadata; retain upstream notices and citations |
| GLPK / swiglpk | [GLPK](https://www.gnu.org/software/glpk/) and [swiglpk](https://github.com/biosustain/swiglpk) | Solver/binding distributions retain their upstream GPL terms and notices |
| Yeast9 | [v9.1.1](https://github.com/SysBioChalmers/yeast-GEM/tree/v9.1.1) | Retain the model's license and attribution |
| LocalMapper | [Pinned commit](https://github.com/snu-micc/LocalMapper/tree/8ff6b3ff2109d834388d4e08ba69399ed3419cad) | LICENSE is **CC BY-NC-SA 4.0** with noncommercial/share-alike terms. Package metadata still says GPL v3; consult the upstream LICENSE/developer rather than assuming MIT |
| RXNMapper | [Official repository](https://github.com/rxn4chemistry/rxnmapper) | MIT package license; model/data terms remain applicable |
| EnzymeMap | [Pinned commit](https://github.com/hesther/enzymemap/tree/c91578856f11312d7941b290f8b773eb41b54973) | MIT package license; processed biochemical data retain their official provenance and applicable source terms |
| lxml | [Version 4.9.4](https://pypi.org/project/lxml/4.9.4/) | BSD-3-Clause; native XML backend for openpyxl; bundled libraries retain their upstream notices |
| PubChem | [Official PUG REST API](https://pubchem.ncbi.nlm.nih.gov/docs/pug-rest) | Retrieved structures are cached locally with CID provenance; retain database attribution and applicable source terms |
| rdchiral | [Pinned commit](https://github.com/hesther/rdchiral/tree/a96477c4b2e9c4e3ecbfa0c88d7c394e480bd874) | Retain upstream license and attribution |

Adapters call public interfaces without copying mapper implementations.
The environment files identify additional dependencies, whose distributions
supply their own license notices. Our MIT license does not relicense them.

The example is a synthetic test network. Selected non-IUPD/IUPD computed inputs and results are published in results/.
They are derived analysis artifacts, not experimental fermentation measurements.
Original model files and full evidence/database caches are not redistributed.
See results/README.md for provenance and limitations.


The E. coli reference uses [BiGG iML1515](https://bigg.ucsd.edu/models/iML1515).
Retrieve original files from BiGG and retain their applicable upstream notices.
Model downloads are not relicensed by this repository. Derived workbooks contain
model identifiers, annotations and flux subsets; cite their original models as well
as this workflow. No new blanket license is asserted over third-party annotations.

Official upstream notices are retained in
[Yeast9 LICENSE](docs/upstream-notices/Yeast9-LICENSE.md) and
[BiGG license page](docs/upstream-notices/BiGG-license.html).
BiGG permits educational, research and nonprofit distribution subject to its notices;
commercial use requires upstream permission. These terms are not replaced by MIT.
