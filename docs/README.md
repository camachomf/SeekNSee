# SeekNSee documentation

This directory documents the scientific behavior, implementation, testing, and distribution of **SeekNSee v0.2.5**.

## Scientific methods

- [Algorithm](ALGORITHM.md): digestion, candidate generation, specificity, scoring, and evaluation sequence.
- [Rules](RULES.md): twelve configurable selection criteria, applicability, scoring, and blocking.
- [Architecture](ARCHITECTURE.md): code organization and component responsibilities.
- [Function contracts](FUNCTION_CONTRACTS.md): inputs, outputs, assumptions, and limitations of core functions.

## Reproducibility and validation

- [Validation](VALIDATION.md): available tests, limitations, and evidence status.
- [Reproducibility](REPRODUCIBILITY.md): baseline provenance, checksums, and scientific regression artifacts.
- [Environment and dependency management](ENVIRONMENT.md): source installation and creation of a project lockfile.

## Installation and distribution

- [Ubuntu installation](INSTALL_UBUNTU.md): running the source application on Linux.
- [Windows installation](INSTALL_WINDOWS.md): running the source application on Windows.
- [Packaging](PACKAGING.md): distinction between source code, prepared binaries, and release verification.
- [Troubleshooting](TROUBLESHOOTING.md): common issues and diagnostic commands.

**Platform claims:** The repository README describes previously built and manually accepted Windows and Ubuntu bundles. The underlying build scripts, artifacts, manifests, and validation evidence are not yet included in the repository; availability of release downloads must be checked when a release is prepared. macOS Apple Silicon remains experimental.
