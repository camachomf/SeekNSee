# SeekNSee

**Traceable, rule-driven selection and prioritization of peptide candidates derived from in silico enzymatic digestion.**

SeekNSee is an R/Shiny application for exploratory peptide selection from protein sequences. The software combines enzyme-aware in silico digestion with configurable rule evaluation, candidate scoring, traceable blocking decisions, reference-proteome specificity assessment, and structured Excel export.

## Current version

- SeekNSee: `0.2.5`
- Algorithm: `2026-09-enzyme-aware-missed-cleavages-p1prime-1`
- Repository status: active development and repository preparation

## Scientific scope

SeekNSee supports enzyme-aware digestion, up to two missed cleavages, candidate filtering, rule-level traceability, peptide ranking, and export of run metadata and evaluation results.

The current scientific core contains twelve configurable rules:

- Specific against reference proteome
- Unique accession group
- Length
- Dibasic motif
- Glycosylation motif
- Chemical risk
- D/E near tryptic cleavage
- Terminal fragment
- Missed cleavages
- Proline at P1-prime
- Isobaric I/L
- Heuristic charge

## Platform status

- Windows 10/11 x64: built, scientifically validated, manually accepted, and frozen
- Ubuntu 24.04 LTS x86_64: built, scientifically validated, manually accepted, and frozen
- macOS Apple Silicon: planned experimental build; native automated and manual validation not yet completed

Validated release archives will be attached to the GitHub release rather than committed to the source history.

## Development and authorship status

SeekNSee is being developed through a collaborative scientific and software-development effort.

Contributor identification, software authorship, author order, contribution roles, affiliations, and citation metadata are currently under review. A `CITATION.cff` file will be added after this information has been reviewed and approved by the project contributors.

## Repository status and collaboration

This repository is currently private while the source tree, documentation, licensing, third-party notices, and release materials are reviewed. Contributions should be proposed through focused branches and pull requests after the initial repository audit is completed.

## License

A public software license has not yet been finalized. Until a license is added, the source code remains subject to the default protections of copyright law. Do not redistribute or reuse the source outside the authorized collaboration without permission.

## Documentation

Documentation is being prepared in the `docs/` directory. Installation instructions, rule definitions, algorithm documentation, validation evidence, packaging procedures, and troubleshooting guidance will be added before the public release.

## Reference proteome

The reviewed human reference proteome used for validation is not committed to this source repository. Distribution, provenance, versioning, and licensing information will be documented separately before public release.

## Contact

Project-maintainer and institutional contact information are under review and will be added before the repository is made public.
