# Source-quality registry

Implementation assumptions should primarily derive from primary studies, peer-reviewed reviews, and peer-reviewed conceptual analyses. Popular, legal, media, ideological, or unverified sources can suggest claims to test but cannot establish scientific requirements.

The machine-readable registry is [`source_quality_registry.json`](../source_quality_registry.json). The labels reflect the source descriptions supplied with the v2 specification; they are not a substitute for a full bibliographic review.

| ID | Short description | Quality class | Implementation authority |
| --- | --- | --- | --- |
| S01 | Current DarkDNA-Observer repository | `primary_computational` | Direct evidence of current software behavior. |
| S02 | Havstad and Palazzo difference-maker analysis | `peer_reviewed_perspective` | Conceptual constraint. |
| S03 | Fagundes et al. junk/spam DNA | `peer_reviewed_perspective` | Conceptual constraint. |
| S04 | Linquist causal-role myopia | `peer_reviewed_perspective` | Conceptual constraint for stopping rules. |
| S05 | Graur rubbish/junk/deleterious distinction | `peer_reviewed_perspective` | Evolutionary classification constraint. |
| S06 | Moran book review / junk-filled-genome null | `peer_reviewed_review` | Benchmark and null inspiration, not primary evidence. |
| S07 | Adey et al. / Puffin-D benchmark | `primary_computational` | Benchmark inspiration. |
| S08 | GBE dark-versus-junk commentary | `peer_reviewed_perspective` | Interpretive constraint. |
| S09 | Walter ncRNA essay | `peer_reviewed_perspective` | Evidence-vector motivation. |
| S10 | Mattick and Amaral chapter | `book_or_chapter` | Hypothesis source requiring independent tests. |
| S11 | 2026 ncRNA chapter | `book_or_chapter` | Taxonomy/background. |
| S12 | Plant transcriptional-junk review | `peer_reviewed_review` | Screening-layer motivation. |
| S13 | TE DNA/RNA/protein review | `peer_reviewed_review` | TE evidence decomposition. |
| S14 | MER11 phylogenetic subfamily study | `primary_experimental` | Copy/subfamily benchmark inspiration. |
| S15 | Dynamic alternative DNA structures | `peer_reviewed_review` | Conformation evidence-level constraint. |
| S16 | Flipons framework | `peer_reviewed_perspective` | Testable framework, not proof. |
| S17 | Information/entropy analysis | `primary_computational` | Statistical-structure methods only. |
| S18 | Agoni mutation-probability hypothesis | `preprint` | Speculative model only. |
| S19 | Budinsky anti-junk/design document | `ideological_or_advocacy_source` | No implementation assumptions; claims may be tested under severe nulls. |
| S20 | LegalServiceIndia article | `popular_article` | No scientific authority. |
| S21 | New Scientist article | `news_report` | Must resolve to primary study before use. |
| S22a | IEEE 11537839 | `unverified` | No use pending manual full-text review. |
| S22b | ScienceDirect S1084952123001726 | `unverified` | No use pending manual full-text review. |

