# Attribution and scope of the open biological subset

`Source_data/open_biological_core/` contains selected UniProtKB/Swiss-Prot
sequence and catalytic-activity records and Rhea directional reaction structures.
UniProt and Rhea make their database content available under **CC BY 4.0**:

- UniProt Consortium: https://www.uniprot.org/help/license
  (official dataset licence listing: https://sparql.uniprot.org/.well-known/void).
- Rhea: https://www.rhea-db.org/help/license-disclaimer.
- Licence: https://creativecommons.org/licenses/by/4.0/.

Source identifiers, source-file hashes and original reaction evidence comments
are retained. Files are selected and reorganized from the study's frozen
extraction, not represented as unchanged provider bulk exports. Reaction
participants were parsed, canonicalized and separated into main participants
and explicitly recorded auxiliaries. Sequence identities are SHA-256 hashes of
the exact amino-acid strings. The supplied protein provenance records trace
every included sequence to UniProt. Publication identifiers are citation
metadata, not redistributed article text.

The subset comprises **375 sequences, 639 main reactions and 1,232 edges**.
Swiss-Prot and UniProt API are two representations of the **same lineage**;
2,490 admitted assertions do not mean 2,490 independent experiments. The
subset belongs to the published development core. It is not an unseen test
set and does not replace the 600-sequence, 2,304-edge mixed-source analysis.

`rebuild_open_core.py` verifies exact sequences, reaction-level experimental
PubMed links, declared Rhea identifiers, chemical parsing, participant
conservation and normalized edge identities. It then rebuilds a reusable
three-file biological core for the evidence-lookup software. Admission
decisions remain frozen: the original ambiguity checks against the complete
mixed-source assertion collection are not rerun by this command.

These source terms are separate from the MIT licence on project-owned code.
The general metric files contain study-derived measurements and hashed
identities; they are not P450Rdb, SABIO-RK or CLEAN database/model exports.
