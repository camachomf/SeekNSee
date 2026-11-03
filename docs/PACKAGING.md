# Packaging and release evidence

This page distinguishes **source availability**, **reported validated bundles**, and **reproducible source builds** for SeekNSee v0.2.5.

## Status at repository preparation

| Target | Repository claim | Available in current source tree |
| --- | --- | --- |
| Windows 10/11 x64 | Built, scientifically validated, manually accepted, frozen | `packaging/windows_x64/README.md`; build scripts and release binary not committed |
| Ubuntu 24.04 LTS x86_64 | Built, scientifically validated, manually accepted, frozen | `packaging/ubuntu_x86_64/README.md`; build scripts and release binary not committed |
| macOS Apple Silicon | Experimental/planned | `packaging/macos_apple_silicon/README.md`; no validated public archive |

The first two entries are **project-maintainer statements from the repository README**, not an independent reproduction of the packaging process during this preparation effort. Do not silently downgrade those historical claims; also do not present them as freshly reproduced until the evidence is reviewed.

## Minimum evidence package per platform

For each release archive, record:

1. Archive filename, size, SHA-256, architecture, OS version, R version, and bundled package versions.
2. Build script(s) with local paths replaced by configurable parameters and any third-party notices.
3. Installation and launch steps performed on a clean machine or isolated environment.
4. Scientific regression and edge-case test results on the target platform.
5. Manual acceptance evidence (including FASTA input, rule controls, output inspection, and Excel download).
6. License and citation metadata approved by the project maintainers.

## Known release risks

- A `renv.lock` documents R package dependencies for **source-based** installs; it is not proof that bundled platform-specific R executables were built from the same dependencies.
- Excel worksheet row capacity (1,048,576 rows) is not explicitly guarded in the current export path; large-run exports require boundary testing and a documented policy.
- Long reference-identifier fields may be abbreviated when written to Excel. The scientific evaluation uses all identifiers, but displayed lists can omit entries; omission counts are reported.
- The presence of a CI workflow is not evidence of a successful remote run until logs are available.

## Public release rule

Do not publish a binary archive without matching checksums, platform acceptance records, third-party license review, and maintainer approval. Do not upload proprietary validation FASTA files unless their distribution is authorized.
