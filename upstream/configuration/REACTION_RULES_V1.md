# P02: auxiliary participants and source-linked main transformations

This is an explicit amendment after viewing the initial full-component baseline.
It is a semantic definition change motivated by source documentation, not an
unseen test or a favorable split. The original full reactions and all previous
outputs remain intact. No performance result selects a role or membership rule.

1. Preserve complete source reaction, direction, structures and component counts.
   Create a separate `MAIN` identity; do not overwrite `RXN` identities.
2. Separate only H2O, H+, O2 and H2O2 by exact structure, plus predeclared redox
   carrier pairs present on opposite sides with equal multiplicity. Never remove
   all carbon-free species or take the largest molecule as the main substrate.
   Carbon coproducts (formate, formaldehyde, CO2, etc.) remain in the transformation.
3. Helper identities come from original, formula-consistent P450Rdb component
   records. Multiple normalized identities for one helper role make that role
   unavailable unless an independently documented chemical reference adjudicates
   it. FMNH2/FMN are anchored to CHEBI:57618/58210 as displayed by Rhea 46308
   (consulted 2026-09-09). Conflicting source name/structure records and their
   parent assertions are quarantined, not silently corrected. Reference selection
   is not majority voting. Charge/tautomer equivalence is used **only to recognize helpers**;
   main-molecule charge, stereochemistry, isotopes and stoichiometry are retained.
4. A residual known but unpaired helper, empty side, identity transformation, or
   a side without a carbon-containing compound prevents primary admission.
   Multi-substrate/multi-product reactions remain multicomponent. Single-pair
   chemistry is a separate view, not a forced simplification.
5. Rhea master IDs are undirected. Use only explicitly linked directional IDs;
   reject conflicting or multiple main directions. Protein/macromolecule
   participants other than the documented hemoprotein-reductase, adrenodoxin and
   2Fe-2S-ferredoxin reactive parts require additional adjudication.
6. Retain the initial P450Rdb membership gate. UniProt additions require a valid
   exact sequence with PF00067, unambiguous accession/taxonomy, a reaction-level
   ECO:0000269 plus a linked PubMed reference, an explicit directional link and
   complete parseable main chemistry. Variant/engineered text in the catalytic
   assertion is quarantined. The absence of such words does not prove wild type.
7. This extension is a **source-linked reference-sequence development cohort**.
   It does not certify the exact assayed construct, physiological use or CYP-domain
   attribution within multifunctional proteins. Keep those qualification flags
   explicit. A strict construct-verified truth set remains a separate unmet gate.
8. UniProt API and Swiss-Prot are one source lineage. Deduplicate the same sequence
   and main transformation; retain all source assertions and publication IDs.
   Evidence count means distinct recorded literature identifiers, not databases.
9. BRENDA/SABIO named chemistry is initially an auxiliary matching/adjudication
   layer, not automatically added training positives. Keep substrate, inhibitor,
   kinetic and multi-protein contexts separate; do not infer structure from EC.
10. Full-component baseline scores and new MAIN-label scores use different
    estimands. They cannot be presented as direct performance improvement.

Authoritative source definitions:

- https://www.rhea-db.org/help/reaction-participant — protein carrier reactive
  parts are not free small-molecule substrates; RHEA-COMP:11964/11965 denote
  reduced/oxidized NADPH--hemoprotein reductase with FMNH2/FMN reactive parts.
- https://www.rhea-db.org/help/reaction-side-direction — left/right sides do not
  by themselves define substrate/product direction; the master ID is undirected.
- https://www.rhea-db.org/help/chebi-in-rhea — differing microspecies conventions
  explain some cross-source structure differences; no automatic main-label merge.
- https://www.ebi.ac.uk/chebi/CHEBI:57618 and
  https://www.rhea-db.org/rhea/46308 — reduced FMNH2 reactive part has formula
  C17H21N4O9P and charge -2; the oxidized FMN reactive part has formula
  C17H18N4O9P and charge -3. This adjudicates a mislabeled FMNH2 source structure.
