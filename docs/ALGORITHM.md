# SeekNSee scientific algorithm (v0.2.5)

> **Documentation draft — requires maintainer review and execution-based verification.**
> This description is derived from the previously provided SeekNSee `app.R` v0.2.5 (algorithm identifier `2026-09-enzyme-aware-missed-cleavages-p1prime-1`). Confirm that the repository copy matches this source before replacing the placeholder. This document describes implemented behavior, not experimentally established peptide detectability.

## Purpose and scope

SeekNSee is a deterministic, configurable **peptide prioritization and decision-support** workflow for targeted proteomics. It enumerates in silico digestion candidates from target protein sequences, evaluates twelve sequence- and context-based rules against an optional reference proteome, and produces a rule-level audit trail. Its heuristic score is **not a calibrated probability of MS detection, proteotypicity, or assay success**.

## Inputs and normalization

- A reference-proteome FASTA; targets may be chosen from this proteome or supplied in an external FASTA.
- FASTA text is read as UTF-8, any byte-order mark and whitespace are removed, sequences are uppercased, and trailing `*` stop characters are removed. Internal `*`, missing headers, empty sequences, and unsupported residues are rejected.
- Accepted sequence characters are `ACDEFGHIKLMNPQRSTVWYXBZJUO` (plus trailing `*` before normalization). Ambiguity is not resolved during matching.
- The first whitespace-delimited token after `>` is the original identifier. UniProt-style `sp|ACCESSION|...` and `tr|ACCESSION|...` headers allow base-accession and isoform parsing; other identifiers are retained without equivalent biological normalization. Duplicate original IDs receive a distinct internal suffix (`__duplicate_...`). A `GN=` field, when present, is recorded.
- Per-peptide coordinates are **one-based and inclusive**.

## 1. Enzyme-aware digestion

Nine configurations are exposed: `Trypsin`, `Trypsin/P`, `LysC`, `ArgC`, `AspN`, `Chymotrypsin`, `Elastase`, `GluC`, and `ProAlanase`.

The current residue-based cleavage definitions are:

| Configuration | Implemented cleavage boundary |
| --- | --- |
| Trypsin | After K or R, excluding K/R followed by P |
| Trypsin/P | After K or R, including K/R followed by P |
| LysC | After K |
| ArgC | After R |
| AspN | Before D |
| Chymotrypsin | After F, W, or Y, excluding sites followed by P |
| Elastase | After A, V, L, I, S, or G |
| GluC | After E |
| ProAlanase | After P or A |

The definitions above are the **software's current deterministic cleavage rules**; they are not a comprehensive empirical model of protease activity or cleavage-site efficiency.

The algorithm determines cleavage boundaries and enumerates concatenations of adjacent digestion fragments with **0, 1, or 2 missed cleavages**, as configured. Every generated candidate stores sequence, target-relative start/end positions, length and number of missed cleavages. The `digest_protein()` function itself does not apply the candidate length filter; `make_target_candidates()` retains **4–25 aa** for evaluation. The separate `length` prioritization rule passes **8–25 aa**. Consequently, 4–7 aa candidates may enter the rule engine and then be blocked or penalized according to their rule state.

## 2. Reference matching and target identity

For each peptide sequence, `find_reference_matches()` checks exact, case-normalized **substring presence** across reference protein sequences (`grepl(..., fixed = TRUE)`) and records matching entries, total substring occurrences, gene annotations, and normalized base accession groups. Repeated peptide lookups can be cached within a run.

This is **sequence-occurrence specificity, not enzyme-constrained peptide uniqueness**: a substring match counts even if the reference protein would not produce that peptide under the chosen enzyme. I/L-equivalent masses and biological sequence variants are not collapsed.

The `unique_proteome` rule permits:

- For a target selected from the reference proteome: matches to that selected internal target ID only.
- For an external target under `strict_absence`: no matching reference entries.
- For an external target under `allow_exact_target_match`: reference entries whose **entire protein sequence** equals the supplied target sequence; other matches are disallowed.

The `unique_protein` rule uses the number of **normalized accession groups** among matched entries; it is not guaranteed to represent gene-level uniqueness.

## 3. Rule evaluation and scoring

Each candidate is evaluated against the twelve documented rules in [`RULES.md`](RULES.md). Rule configuration can be `ignore`, `apply` or `block`:

- `ignore`: failed checks have no penalty and never block.
- `apply`: a failed applicable check reduces the heuristic score by 1, except `miscleavage`, whose penalty equals the candidate's missed-cleavage count (1 or 2).
- `block`: a failed applicable check makes the candidate **blocked**, without by itself subtracting points from the score.
- A nonapplicable rule produces `not_applicable` and has neither a penalty nor a blocking effect.

The implemented score is:

`score = 12 + sum(rule_score)`

where each `rule_score` is 0 or a nonpositive penalty. **Eligibility (`retained` or `blocked`) is a separate decision from score.** A blocked candidate can therefore retain a high score. Scores should be interpreted only with their rule configuration, enzyme and input/reference context; they are not comparable as calibrated probabilities across different runs.

The default configuration is `block` for `unique_proteome` and `length`, and `apply` for the other ten rules.

## 4. Candidate trace and reporting

Each evaluated candidate instance is identified by the composite `candidate_key` (target ID, coordinates, missed-cleavage count, peptide sequence). An instance is not the same as a globally distinct peptide sequence: the same peptide can occur at multiple positions or in different targets.

`make_target_candidates()` creates **one rule-trace row per candidate and rule** (12 rows per evaluated instance), including applicability, pass/fail/nonapplicable, selected state, score contribution, and whether a rule blocked the candidate.

`build_seeknsee_workbook()` assembles seven sheets:

1. `Run_summary`
2. `Rule_configuration`
3. `Target_proteins`
4. `Ranked_peptides` (retained, score-sorted)
5. `Rule_trace`
6. `Blocked_peptides`
7. `Rule_definitions`

The workbook records the software and algorithm versions, R version, FASTA *filenames*, settings, and candidate counts. **Input file hashes, dependency lockfile identifiers, and Git commit SHA are not recorded by the examined version.**

Long semicolon-delimited identifier lists may be abbreviated at a 30,000-character safe threshold, with omission counts and summary audit metadata. **The examined code does not pre-check Excel's 1,048,576-rows-per-worksheet limit**; this is a release-relevant limitation, especially for `Rule_trace`. An Excel workbook alone is not a lossless archive if identifier lists have been abbreviated.

## 5. Reproducibility boundaries

The score is heuristic and explainable, not experimentally calibrated. The glycosylation motif, chemical-risk flags and `charge_hat = 1 + nK + nR + nH` are **sequence heuristics**, not measurements of modification occupancy, chemical instability, ionization, or observed precursor charge. The acidic-cleavage and P1-prime rules are enzyme- and position-specific as documented in `RULES.md`.

Claims about measured DIA/DDA performance, proteome-scale counts, runtime, and the number of automated assertions require **separately versioned input data, test scripts and benchmark artifacts**. This document by itself does not verify those reported results.
