# Source attribution and transformations

## Normalized human development labels and frozen model

Holmer M, de Bruyn Kops C, Stork C, Kirchmair J. CYPstrate: A Set of Machine
Learning Models for the Accurate Classification of Cytochrome P450 Enzyme
Substrates and Non-Substrates. Molecules 2021, 26, 4678.
https://doi.org/10.3390/molecules26154678

Article/supplement source: https://www.mdpi.com/1420-3049/26/15/4678
Open-access licence identified by the publisher: CC BY 4.0,
https://creativecommons.org/licenses/by/4.0/

Included derivatives: `data/development_labels.json` and the bundled human
kNN model. Project processing normalized structures and scaffold groups,
excluded unresolved structures and conflicting labels, and consolidated labels
into 14,955 isoform-compound pairs. Original train/test designations were treated
as development-exposed, not as a new independent test. The bundle is this
project's similarity-weighted kNN representation, NOT the authors' CYPstrate
random-forest/SVM implementation or model weights. Source and bundle checksums
are retained. Changes and additional constraints are explained in the paper.

## External classifications and saved predictions

Ni et al., Scientific Data (2025), DOI:
https://doi.org/10.1038/s41597-025-05753-8
Data: Curated CYP450 Interaction Dataset: Covering the Majority of Phase I Drug
Metabolism, Figshare record 26630515, version 4:
https://doi.org/10.6084/m9.figshare.26630515.v4

The frozen source acquisition manifest identifies CC BY 4.0 and all 26 source
file identifiers and hashes. Included derivatives are normalized source labels,
saved project predictions, and the 3,035-label evaluation table. These are
source-reported classifications, not assay-harmonized catalytic truth. Original
structure/source/scaffold filters are separate from the later parent-identity
and steroid-topology sensitivity flags. No source labels were revised by those
sensitivities. A flag is not proof of a labelling error.

## General Figure 2 evaluation fixture

The fixture contains opaque candidate identifiers, labels, applicability masks,
scores, protein-group membership and fixed ordering keys. It contains no raw
database export, sequences or structures. It is a numerical saved-score
reproduction subset, not an independently reusable biological evidence atlas.
Its inclusion does not establish permission to redistribute upstream databases.
Redistribution of the upstream biological atlas requires a separate source-rights
decision; that atlas is not in this package.

No licence on this repository overrides third-party source conditions.

## V20 reproduction derivatives

The development-fold table preserves assignments already present in the
normalized CYPstrate labels. Fingerprints are regenerated from those structures
and the normalized Figshare structures; provenance records compare them with
the original arrays. These human-data derivatives retain the source-specific
CC BY 4.0 attribution above. New evaluator/packaging code is project-owned MIT
code. No additional third-party raw dataset or model weight is added by V20.
