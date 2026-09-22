# Directed atom-level centre reconstruction: mapping stage

Start from all 1,341 frozen directed MAIN transformations. Retain the original
component lists, stereochemistry, charge, multiplicity and reaction keys.
Concatenate each side only for atom-mapping input; do not rebalance reactions,
insert inferred oxygen donors or remove carbon coproducts. Missing auxiliaries
mean a newly incorporated product atom may have no main-substrate counterpart.

Use official RXNMapper 0.4.3 with its bundled ALBERT model and documented defaults,
in a new isolated project environment (Torch 2.1.2 module, Transformers 4.38.2,
RDKit 2026.3.5). Preserve package download hashes and installed model/code hashes.
Run offline CPU inference on Slurm. Do not send project reactions to a web demo.
The model's chemical pretraining exposure is not certified independent of the
evaluation reactions. Mapping is an unsupervised representation operation,
not a new catalytic label or a ground-truth mechanism.

For every candidate, retain successful output or an explicit failure. Do not
truncate reactions beyond the model's token limit. Validate each mapped side
against the input's canonical, map-free, stereochemical component multiset.
Require parseability, unique nonzero map identifiers within each side, and
matching element/isotope identity for shared maps. Distinguish unmapped atoms,
one-sided map identifiers, added/removed atoms and changed shared-atom bonds.
Zero mapping confidence is retained, not interpreted as impossible catalysis;
no performance-selected confidence threshold is applied.

The next feature stage will identify changed bond connectivity/order and atom
attributes, record centres on both sides, and encode centre environments for
interaction with ordered protein sites. Stereochemical-only changes require
explicit treatment: a raw chiral-tag difference depends on atom order and is
not by itself a valid inversion detector. Centres without a reliable mapping
must remain unavailable. Symmetry-related alternative maps and localization
stability must be assessed before claiming biologically precise centres.

Official sources: [RXNMapper code and usage](https://github.com/rxn4chemistry/rxnmapper),
[PyPI release metadata](https://pypi.org/project/rxnmapper/0.4.3/), and
[Schwaller et al., Science Advances (2021)](https://advances.sciencemag.org/content/7/15/eabe4166).
The code repository declares an MIT license. This record does not authorize
redistribution of underlying biological source data.
