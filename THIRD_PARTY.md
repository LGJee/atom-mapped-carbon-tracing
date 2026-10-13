# External software, models and data

Install these components separately. Their implementations, weights and complete
databases are not included in this repository.

| Component | Source | Applicable upstream terms |
|---|---|---|
| COBRA Toolbox | [Official repository](https://github.com/opencobra/cobratoolbox) | Retain upstream license/citation guidance; solver licensing is separate |
| Yeast9 | [v9.1.1](https://github.com/SysBioChalmers/yeast-GEM/tree/v9.1.1) | Retain the model's license and attribution |
| LocalMapper | [Pinned commit](https://github.com/snu-micc/LocalMapper/tree/8ff6b3ff2109d834388d4e08ba69399ed3419cad) | LICENSE is **CC BY-NC-SA 4.0** with noncommercial/share-alike terms. Package metadata still says GPL v3; consult the upstream LICENSE/developer rather than assuming MIT |
| RXNMapper | [Official repository](https://github.com/rxn4chemistry/rxnmapper) | MIT package license; model/data terms remain applicable |
| EnzymeMap | [Pinned commit](https://github.com/hesther/enzymemap/tree/c91578856f11312d7941b290f8b773eb41b54973) | MIT package license; processed biochemical data retain their official provenance and applicable source terms |
| rdchiral | [Pinned commit](https://github.com/hesther/rdchiral/tree/a96477c4b2e9c4e3ecbfa0c88d7c394e480bd874) | Retain upstream license and attribution |

Adapters call public interfaces without copying mapper implementations.
The environment files identify additional dependencies, whose distributions
supply their own license notices. Our MIT license does not relicense them.

The example is a synthetic test network. Original WT/IUPD study inputs and large
research outputs are not published in this code release.
