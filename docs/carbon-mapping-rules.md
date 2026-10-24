# Carbon mapping fallbacks: scientific basis and limits

The code is in `carbon_flow_tools/carbon_mapping_rules.py`. These are conservative
project rules requiring scientific review, not a published universal mapper and
not a guarantee that every modeled reaction has a complete atom map.

## What a biomass reaction represents

A biomass objective combines precursor requirements, macromolecular composition,
and energy/maintenance costs. Its decimal coefficients describe material demand
per biomass unit; they are not a list of individual molecules in a single
elementary chemical event. [BOFdat, Methods](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1006971)
constructs these demands from composition, sequences and abundance data.

Yeast models can pool several species into a protein, RNA, carbohydrate or lipid
pseudometabolite. The yETFL reconstruction explicitly replaces protein/RNA
pseudoreactions with sequence-based synthesis, documenting the difference between
an aggregate pool and a specified macromolecule.
[yETFL, Data collection and Allocation](https://www.nature.com/articles/s41467-021-25158-6).
COBRAme likewise needs gene/product information to construct expression reactions:
[COBRAme, ME-model architecture](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1006302).

A chemical formula or total amino-acid composition alone does not define all
positional atoms of a protein or biomass pool. The present exporter does not
supply explicit polymer sequence-to-atom structures. Therefore these pools do not
receive fabricated SMILES, virtual carbon positions, or maps into arbitrary
small-molecule representatives. PubChem name lookup is skipped for generic
biomass/macromolecule pool names; explicit user structures and exact identifiers
remain subject to normal structure validation. Existing generic-pool annotations
are flagged for review; skeleton fallbacks do not assign their polymer atoms.

## Selection order

Matching user-maintained manual carbon records have highest priority. The order
below applies when no manual record matches; raw tool predictions are always retained.

0. Check explicit model formulas against structure carbon counts. Conflicts block
   whole-reaction mapper selection, even if all tools produced predictions.
   They can indicate a wrong annotation or unresolved repeat-unit/chain-length
   convention; do not automatically rescale stoichiometry or replace structures.
   Independent cofactor pairs using consistent structures may still be recovered.
1. Prefer complete mapper carbon agreement. When tied or conflicting, use
   EnzymeMap > RXNMapper > LocalMapper; native scores are not interchangeable.
2. Before discarding an otherwise aligned prediction, retain its carbon-only
   view: unique carbon map IDs shared across both sides, with noncarbon and
   ambiguous/one-sided IDs cleared. Require exact component identity, charge,
   isotope and stereochemistry. Prefer a complete view; otherwise select one
   partial mapper view using the same priority. Do not merge conflicting tools.
3. If any mapper generated a nonempty prediction, never leave the selected or
   final mapping empty. When no usable correspondence passes validation, retain
   one original nonempty prediction by EnzymeMap > RXNMapper > LocalMapper priority,
   set ReviewRequired=TRUE and MappingCoverage=Unresolved, and explain the failed
   checks in SelectionReason. Unresolved candidates are audit data, not propagation
   inputs. Mapping_Selection and Final_Mapping both retain the selected string.
4. Only when all mappers returned empty predictions, try the rules below. Missing
   information in other participants does not invalidate a uniquely resolved
   small-molecule skeleton pair. No rule may replace a retained mapper candidate.
5. Leave unsupported carbon unresolved and explain why in SelectionReason.

Published mapper status does not certify a particular prediction. A partial
prediction must not be relabeled complete to avoid an empty cell. Consensus and
carbon-only recovery live in `consensus.py`; raw tool SMILES remain untouched.

## Implemented rule scope

| Case | Carbon handled | What stays unresolved |
| --- | --- | --- |
| ATP/ADP/AMP; GTP/GDP/GMP; CTP/CDP/CMP; UTP/UDP/UMP | Corresponding nucleotide carbon skeleton within its verified connectivity family | Different bases, deoxynucleotides not covered by the reference library, modified nucleotides, multiple competing partners |
| NAD/NADH; NADP/NADPH | Carbon skeleton within each separate redox family | NAD-to-NADP conversions, degradation or unmatched structures |
| Free CoA and thioester-bound CoA | Common uniquely positioned CoA-core carbons after recognizing the full core | Acyl carbons, uncertain symmetric positions, multiple possible CoA partners |
| Identical specified structures | Unique carbon-position retention | Ambiguous symmetry or competing donor/product pairing |
| ATP use in biomass/maintenance or aminoacyl-tRNA charging | ATP-to-ADP/AMP skeleton, if uniquely identifiable with equal coefficients | Unknown protein/tRNA atoms and the rest of the biomass pool |
| RNA/DNA/protein/lipid pool formation | Any independently established cofactor correspondence | Positional polymer mapping without explicit sequence/structure and assembly rules |
| Exchange, sink, demand | Existing boundary handling; retain the reaction record | No invented internal donor or product for an open boundary |

References supporting these reaction classes include
[IUBMB alanyl-tRNA synthetase](https://iubmb.qmul.ac.uk/enzyme/EC6/1/1/1.html),
[IUBMB translation GTPase](https://iubmb.qmul.ac.uk/enzyme/EC3/6/5/3.html), and
[Rhea nucleoside diphosphate kinase](https://www.rhea-db.org/rhea/27688).
These establish biochemical context, not the validity of every application of a
rule. Conserved molecular groups require atom-level evidence beyond stoichiometry:
[Haraldsdottir and Fleming, 2016](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004999).
The CoA carrier must be distinguished from attached carbon, a problem also
highlighted in [NICEpath, section 2.2](https://pmc.ncbi.nlm.nih.gov/articles/PMC8545321/).
Our implementation does not reproduce NICEpath's pathway scoring algorithm.

Rules require valid supplied structures, compatible specified carbon
stereochemistry, equal positive finite amounts in the active direction, and a
unique donor/product pairing. Internal graph normalization never changes stored
structures. Matching is bounded; ambiguity or reaching the match cap prevents a
rule assignment. CoA retains only position pairs common to every valid alignment;
its two equivalent pantothenate methyl positions are not arbitrarily assigned.
No maximum-common-substructure guess or redistribution of unknown carbon is used.

For the reference biomass reaction, ATP-to-ADP keeps 10 carbons with coefficient
55.3, and NADPH-to-NADP keeps 21 with coefficient 3. The unresolved biomass pool
is not declared atom-balanced or fully mapped. Fractional amounts are saved
exactly; qualitative unit transfers avoid rounding or creating 55 copies.

## Saved contract and review

Final_Mapping includes every reaction and MappingCoverage is Complete, Partial,
Unresolved or Boundary. Legacy NoCarbon records remain readable. Mapper partial views use sparse MappedReactionSMILES.
Rule partials show a reduced unit mapped subreaction in both selected/final
mapped-SMILES fields. Hidden Mapping_Info records (RecordType=PartialMapping) retain exact source/product structures,
zero-based atom-index pairs, direction, original coefficient and rule version.
The loader recomputes those records and checks the reduced string; it never
compares a reduced map as if it described the full original equation.
No full-reaction SMILES is claimed for a reaction with unknown participants.

Rules and partial mapper views have ReviewRequired=TRUE and AIGenerated=FALSE.
Step 3 uses only established positional transfers in its source balances, with
original fluxes and coefficients. It revalidates saved structures and rules,
never reruns a mapper, and retains unresolved carbon inflows as Unknown.

The existing reference workbooks illustrate why this matters: THF has formula
C19H21N7O6 but model annotation C1CCOC1 contains four carbons. Consequently r_0727
requires upstream THF structure review rather than blind mapper selection.
[ChEBI:57453](https://www.ebi.ac.uk/chebi/CHEBI%3A57453) independently confirms the
C19 formula and charge -2 for the tetrahydrofolate species; it is not the
four-carbon structure currently supplied by the model annotation.
Glycogen/mannan and dolichyl species also have repeat-unit/chain-length count
conflicts requiring explicit model interpretation. These are input consistency
issues, not opportunities to manufacture new carbon correspondences. Step 3
rejects complete final mappings that still use such inconsistent structures.

The AI-review skill must assess biochemical applicability, not just graph match.
It leaves original Mapping_Selection, raw predictions and saved partial correspondence records intact.
In an AI copy, confirmed partial rules remain Partial; rejected/scientifically
unresolved rules become Unresolved and are disabled. A validated complete AI
replacement becomes Complete and takes precedence. Neither a rule match nor AI
review proves experimental isotope abundance or quantitative carbon allocation.

## Carbon-free and terminal cases

carbon_presence uses valid structures and explicit elemental formulas to return
True/False/unknown. Co is cobalt, not a carbon token; R/X placeholders cannot
prove carbon absence. Missing SMILES is never equivalent to zero carbons.
Automatic NoCarbon classification is not enabled. Reactions without a usable
mapper result remain reviewable; carbon-free formulas alone do not bypass mapping.
Open boundaries are Boundary. Legacy NoCarbon records are still validated against
the saved inputs when read. Source and metabolite sheets stay visible.

Unknown generic product pools with no internal consumer in the original nonzero
network are described as terminal. A missing mapped downstream edge alone never
establishes terminality. This diagnostic does not fabricate pool atom positions
or change fluxes or original equations.

## Withdrawn additions

Automatic NoCarbon classification, FAD/FMN-family skeleton recovery and the
MEcPP-to-HMBPP (IspG) template are not used. Existing output files are not rewritten.
Historical validation outputs created with these additions are retained only for
comparison. Workbooks containing their withdrawn partial rules must be reselected
from the saved mapper predictions before using the current loader; removed rules
are not silently accepted. Earlier nucleotide, NAD(P) and CoA rules remain active.


## Reduced inputs to published mappers

Before prediction, missing-SMILES participants may be left out of the mapper-only
input view. Both retained sides must be nonempty; original equations and fluxes
are not edited. Successful carbon correspondences are provisional Partial maps,
not proof that excluded participants have no carbon role. The omission record is
kept in Mapping_Info and SelectionReason. This general input projection does not
restore the withdrawn IspG-specific template. Optional AI review examines the
whole reaction and can retain or disable the reduced mapping. Ordinary analysis
does not require prior manual or AI approval of the projection.


## User-maintained manual data has priority

The two user-maintained CSVs have only three columns:

- data/manual_smiles.csv: Species, MetaboliteName, SMILES.
- data/manual_atom_mappings.csv: Species, ReactionID, MappedReactionSMILES.

Species is the full binomial shared by all conditions, for example Escherichia_coli.
Step 2 resolves it from source-model provenance or the analysis-name prefix.
Manual records have priority over software predictions and automatic fallbacks;
raw mapper sheets remain unchanged. Atom-mapping rows follow the direction of the
original written equation. A negative selected flux automatically reverses the
mapped reaction before use. Duplicate species/reaction keys and invalid records
raise an error rather than silently selecting a row.

Matched SMILES are written directly into Mapping_Selection.SelectedMappedReactionSMILES
and Final_Mapping.MappedReactionSMILES. Complete validated mappings are Complete;
carbon-only or reduced subreactions remain Partial with ReviewRequired=TRUE.
The program derives the required positional/stoichiometric bookkeeping and saves
it in hidden Mapping_Info. Users do not enter atom indices or transfer coefficients.
Step 3 can reload these records without the external CSV, preserving the supplied
mapped text. Original equations and fluxes are unchanged.

Reduced maps must uniquely identify the participating original molecular structures
(including stereochemistry/charge). Single-fragment reduced components are supported;
full maps use the normal component validator. Indistinguishable named participants
or inconsistent molecular-copy ratios require explicit resolution, not guessing.
Mapped positions cannot exceed the original directional stoichiometric budgets.
Manual provenance is a supplied scientific assertion, not proof provided by software.

The user-designated MECDPDH5/IspG five-carbon record was migrated without changing
its carbon pairs. Its biochemical reference is
[IUBMB EC 1.17.7.3](https://iubmb.qmul.ac.uk/enzyme/EC1/17/7/3.html).
The previous detailed CSV and its evidence remain in the local migration backup.

AI-assisted records are separate four-column indexes:

- data/ai_smiles.csv: Species, MetaboliteName, SMILES, Source.
- data/ai_atom_mappings.csv: Species, ReactionID, MappedReactionSMILES, Source.

Source includes evidence URLs and links to immutable review archives. Only accepted
AI replacements enter the mapping index; failed proposals remain in detailed review
reports. Repeated identical records merge their sources; different mapping candidates
are retained for inspection. AI mapping indexes are not automatic step-2 overrides.
The review skill revalidates exact-context archive candidates before reuse. A user may
promote a checked entry into the corresponding manual table deliberately.
AI structures can fill missing input structures after formula/charge checks and remain
StructureReviewRequired=TRUE; they never overwrite manual or existing input structures.
AI-assisted structure discovery alone does not set atom-mapping AIGenerated=TRUE.


Manual SMILES use exactly three columns: Species, MetaboliteName, SMILES. Species
is the full underscore-separated binomial used by result prefixes, without the
non_IUPD/IUPD condition suffix. For example Escherichia_coli shares one isoprenol
record across conditions and compartments. Terminal compartment suffixes are
ignored for name matching; names are otherwise exact. Local model metabolite IDs
are not manual-table keys. Conflicting structures for one species/name, or an
explicit model formula/charge incompatible with an override, raise an error unless
the protonation-compatible warning policy below applies.
The original detailed workbook IDs and review archives remain available for audit.


## Partial structures and unspecified repeats

`structure_kind` in `carbon_mapping_rules.py` is the common policy used by manual
validation, source annotations, mapper preparation, selection and step 3. It reports
`complete`, `partial`, `polymer`, `invalid`, or `missing`; original structures are
retained. Here `complete` means a syntactically explicit graph, not experimentally
verified identity or fully specified stereochemistry.

- `*` and `[*]` are legal dummy atoms (atomic number zero). They are not known
  carbon, hydrogen, methyl, or evidence that two variable groups are identical.
  For example, [ChEBI phosphatidylcholine](https://www.ebi.ac.uk/chebi/CHEBI:64482)
  has two R groups and an explicit glycerophosphocholine/acyl-carbonyl skeleton.
- A formula such as [KEGG C00110](https://www.kegg.jp/entry/C00110),
  `C20H37O4P(C5H8)n`, describes variable repeats. Carbon count depends on n; no
  positional graph is inferred, and n is never set to one. A supplied repeat-unit
  SMILES does not resolve an unspecified whole polymer. Such a participant is
  omitted only from the mapper input view, with its reason saved in Mapping_Info;
  the original equation, stoichiometry and flux remain unchanged. Other, independently
  mapped participants can still contribute partial carbon transfers.
- SMILES aromatic `n` and `[nH]` are actual nitrogen atoms, not polymer repeats.
  Ring numbers, charge, stereochemistry, isotopes and dot-separated components are
  parsed normally by RDKit; none is stripped to make a structure pass validation.
  See the [RDKit Book](https://www.rdkit.org/docs/RDKit_Book.html).
- Genuine syntax/valence errors in applicable manual rows still fail with the file,
  line, species, compound and original string. Runtime manual/AI structure validation
  checks the current species and compounds present in its input; standalone
  `load_manual(path)` still validates the entire table. Conflicting applicable rows
  and complete-structure formula/charge conflicts remain errors except for the documented
  protonation-compatible warning below. Partial fragments
  cannot be compared with whole-compound formula/charge as if they were complete;
  a known fragment carbon count exceeding a fixed formula count is still inconsistent.

LocalMapper 0.1.1 and RXNMapper 0.4.3 returned predictions for both an ACP dehydration
and a two-variable-group phosphatidylcholine hydrolysis probe. The pinned EnzymeMap
single-step rules returned an ACP prediction but no PC prediction. These small tests
establish parser/execution support for those inputs, not universal chemical accuracy.
The EnzymeMap probe used an empty isolated lookup CSV and the installed official rules;
coverage of the full dataset was not tested. Versioned probe outputs are local validation
artifacts. Published mapper preference and raw predictions are unchanged.

A wildcard prediction can never be `Complete`, even if the tool numbers every dummy.
The existing carbon-only selection keeps only unique explicit C-to-C IDs, clears dummy
and other noncarbon IDs, sets `Partial` and `ReviewRequired=TRUE`, and records the scope
in SelectionReason. Two R groups are not merged or assigned carbons. An unresolved
nonempty software candidate remains visible for review and is excluded from propagation.

Automatic skeleton rules do not infer ACP/lipid transformations from generic graph
similarity. The narrowly supported exception is a pure, same-identity compartment
transport with equal amounts and the same partial graph: only explicit carbon pairs
invariant over every exact alignment are retained. Dummy and symmetry-ambiguous
positions stay unmapped. Other partial transfers require a usable software prediction
or independently reviewed manual correspondence. They use the existing hidden
Mapping_Info records; no Partial_Mapping worksheet is generated.

Step 3 revalidates sparse maps and balances only explicit carbon positions.
A partially described metabolite does not have resolved origins for its unspecified
region. Variable-repeat
polymers have no automatically constructed positions. A configured carbon source must
have a complete explicit structure; generic carbon sources require resolution first.
Fixed-flux source-contribution exports omit StructureType, KnownCarbonCount and
ContributionScope; the solver retains this structure-scope metadata internally.
Fractions for a partial graph are means over its explicitly known carbon positions
only; CarbonCount is blank because the whole count is unknown. Polymer totals/fractions are unavailable, not repeat-unit
fractions or zero contribution. Unknown regions are not redistributed among sources.

Unsupported cases include unspecified polymer lengths, variable-group identity or
transfer without evidence, generic protein/biomass positional structures, ambiguous
carbon alignments, and parser-invalid chemical notation. An explicit chain/sequence,
exact small-molecule structure, or reviewed positional correspondence is needed before
those regions can enter atom-level tracing.


## Protonation-compatible manual structures

A complete manual/AI-table structure may differ from the model formula/charge only
by a protonation-compatible amount: all non-hydrogen element counts must match,
both charges must be finite integers, and delta H must equal delta formal charge.
Such a mismatch warns and continues without editing the supplied SMILES, formula,
charge, stoichiometry or flux. StructureStatus records both forms and
StructureReviewRequired=TRUE makes affected selected/final mappings reviewable.
Formula/charge compatibility does not establish compound identity or a preferred
pH microstate. Carbon/other-element conflicts, missing comparison metadata, and
hydrogen changes without matching charge changes do not qualify for this exception.
Other manual-table checks remain in force. No pH adjustment or Open Babel dependency
is introduced; model annotations/database candidates retain their existing checks.
