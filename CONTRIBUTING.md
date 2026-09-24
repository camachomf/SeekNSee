# Contributing to SeekNSee

SeekNSee is currently in a controlled repository-preparation phase.

## Before contributing

1. Coordinate the proposed change in an issue or with the repository maintainer.
2. Use a focused branch with a descriptive name.
3. Do not modify the validated scientific algorithm and packaging logic in the same change.
4. Do not commit reference proteomes, installed R libraries, bundled runtimes, build outputs, exported workbooks, credentials, personal data, or private governance records.
5. Document the purpose, scope, and validation of every scientific or packaging change.

## Branch naming

Examples:

- `docs/rule-documentation`
- `test/scientific-regression`
- `fix/windows-launcher`
- `ci/macos-apple-silicon`

## Pull requests

Pull requests should describe:

- the problem or objective;
- files and components changed;
- scientific impact;
- validation performed;
- known limitations;
- whether release artifacts must be rebuilt.

## Scientific integrity

Changes to digestion, candidate generation, rule definitions, scoring, blocking, specificity logic, reference matching, export structure, or default configuration require explicit regression testing and review.

## Authorship and contribution roles

Git activity documents technical participation but does not independently determine software authorship or author order. Formal contributor roles and citation metadata will be maintained through a separate review process.
