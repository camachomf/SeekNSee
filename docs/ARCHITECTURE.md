# SeekNSee — Code architecture and function inventory

**Status:** Working draft for maintainer review. **Application version:** 0.2.5. **Algorithm identifier:** `2026-09-enzyme-aware-missed-cleavages-p1prime-1`. **Source:** Static review of the `app/app.R` function inventory and parsed R expression supplied on 3 November 2026, together with the successful 23-check baseline run reported by the collaborator. No modifications to the application are proposed in this document.

## 1. Purpose and release constraints

SeekNSee is a local R/Shiny application for rule-based, traceable prioritization of digestion-derived peptide candidates. The immediate publication-readiness objective is to make the existing algorithm **understandable, testable, installable, reproducible, and citable** without changing its scientific decisions. The code is not a validated predictor of experimental peptide detectability.

**Non-negotiable constraints**

- Preserve the original `app/app.R` and record its SHA-256 before making code changes. The collaborator reported `3cd007d023e6e0ddd59908a74accb156b709c8693d3dd9e7567694db0b907c27`.
- Preserve the nine enzyme configurations, up to two missed cleavages, generation of candidate instances 4–25 aa, twelve configurable rules, and existing retention/scoring semantics unless a separately reviewed scientific change is approved.
- Preserve the seven-sheet workbook schema, or explicitly version and document any output-schema change.
- Distinguish syntax parsing, passing unit checks, deterministic regression, workload benchmarking, Shiny testing, and experimental validation; none substitutes for all the others.
- Document the original method, its assumptions and limitations before changing implementation boundaries.

## 2. Current architecture (`app/app.R`)

The repository already separates application, documentation, packaging, and tests at the directory level. **Within `app/app.R`, however, data handling, scientific logic, report generation, test helpers, Shiny UI, and server orchestration remain in one source file.** The currently observed direct R packages are `shiny`, `dplyr`, `stringr`, and `writexl`.

| Conceptual layer | Responsibilities today | Primary functions/objects | Documentation and regression focus |
|---|---|---|---|
| Data and identifiers | Parse reference/target FASTA, normalize identifiers, represent empty inputs | `empty_proteome`, `extract_accession`, `parse_identifier`, `read_fasta` | FASTA grammar, duplicate headers, accession normalization, acceptable residues, invalid inputs |
| Scientific core — digestion | Locate enzymatic cleavage positions and enumerate fragments | `find_regex_positions`, `cleavage_positions`, `digest_protein`, `has_internal_cleavage_site` | Nine enzymes, terminal boundaries, P1/P1′ exceptions, missed cleavages, coordinates |
| Scientific core — specificity and contextual checks | Locate sequence occurrences in the reference and assess cleavage context | `count_fixed_occurrences`, `find_reference_matches`, `exact_reference_ids`, `valid_tryptic_cut`, `acidic_window_at_cut`, `assess_acidic_near_tryptic_cleavage`, `assess_proline_p1prime` | Per-entry versus per-occurrence match counts, normalized accession groups, exact-target policy, enzyme applicability |
| Scientific core — rules and candidates | Calculate rule outcomes, heuristic score, eligibility, per-rule trace and candidate instances | `evaluate_peptide`, `default_rules`, `validate_rule_configuration`, `make_target_candidates` | Apply/ignore/block, 12 trace records per candidate, 4–25 aa generation vs 8–25 aa length rule, `candidate_key` uniqueness |
| Presentation | Configure inputs, show candidates and errors, control run lifecycle | `build_rules`, `ui`, `server`, nested `recompute`, nested `close_progress` | UI–core decoupling, reactive state transitions, session lifecycle, error reporting |
| Reporting and export | Build seven-sheet workbook, adapt long identifiers and write xlsx | `make_export_filename`, `empty_candidate_table`, `empty_rule_trace_table`, `build_seeknsee_workbook`, `truncate_delimited_identifier_list`, `adapt_sheet_identifier_lists_for_excel`, `adapt_workbook_for_excel`, `empty_excel_preflight_issues`, `validate_excel_workbook`, `export_seeknsee_workbook` | Seven-sheet schema, Excel character and row limits, omissions audit, deterministic contents excluding timestamps |
| Verification and validation helpers | Headless analysis, summaries, accounting checks and convenience test helpers | `run_validation_case`, `summarize_validation_result`, `validate_internal_identities`, `evaluate_test_peptide`, `compare_digestion` | Executable fixtures, expected outputs, edge cases, empirical validation boundaries |

**Ownership note:** The table above is a conceptual responsibility map, *not* evidence that these functions already live in separate modules.

## 3. Detailed function inventory

Line numbers refer to the collaborator's `grep` output for the unmodified `app/app.R` and should be regenerated after any edit.

### 3.1 Data and identifiers

| Function | Approx. line | Responsibility | Critical checks |
|---|---:|---|---|
| `empty_proteome` | 48 | Return a typed empty reference-proteome table | Consistent column names/types |
| `extract_accession` | 57 | Obtain a filename-safe accession token | UniProt vs generic headers, punctuation, empty input |
| `parse_identifier` | 72 | Parse accession, base accession and isoform suffix | UniProt `sp|`/`tr|`, isoforms, custom accessions |
| `read_fasta` | 94 | Read/validate sequences, headers, gene metadata and unique instance IDs | BOM, blank lines, stop characters, invalid residues, duplicated records, large files |

### 3.2 Digestion and reference matching

| Function | Approx. line | Responsibility | Critical checks |
|---|---:|---|---|
| `find_regex_positions` | 311 | Return motif positions | No-match/one-match/multiple-match behavior |
| `cleavage_positions` | 316 | Locate enzyme-specific cut boundaries | Each enzyme, N/C termini, Trypsin and Trypsin/P difference |
| `digest_protein` | 348 | Enumerate fragments with 0–2 missed cleavages and coordinates | Three-fragment joining, contiguous fragments, input validation |
| `has_internal_cleavage_site` | 404 | Detect an internal enzyme cleavage motif | Verify live call sites before deciding whether it is removable |
| `count_fixed_occurrences` | 447 | Count occurrences of a fixed peptide string | Repeated motif behavior, exact versus overlapping matches |
| `find_reference_matches` | 452 | Scan the reference sequences and optionally cache matches | Multiple proteins, occurrences vs entries, cache isolation/invalidation, scale |
| `exact_reference_ids` | 592 | Find entries with full-length sequence identity to target | Identical isoforms and external-target policy |

### 3.3 Cleavage context, rules and candidate engine

| Function | Approx. line | Responsibility | Critical checks |
|---|---:|---|---|
| `valid_tryptic_cut` | 495 | Decide whether a K/R cleavage boundary is valid under Trypsin | K/R followed by P and edge boundaries |
| `acidic_window_at_cut` | 505 | Identify D/E near a valid cleavage boundary | P3–P3′ context, sequence edges, exclusion of scissile residue |
| `assess_acidic_near_tryptic_cleavage` | 522 | Assess whether a candidate borders acidic-risk sites | Trypsin-only applicability, both ends, no valid site |
| `assess_proline_p1prime` | 556 | Evaluate internal K/R–P boundary for Trypsin/P | P1′ orientation, N-terminus, non-Trypsin/P |
| `evaluate_peptide` | 599 | Evaluate 12 rules and return score, block status and trace | Rule states, nonapplicable rules, `score = 12 + sum(trace)`, blocked score independent of status |
| `build_rules` | 724 | Translate Shiny controls into a rule configuration | Missing/invalid controls; avoid importing Shiny into the core |
| `make_target_candidates` | 731 | Generate candidates, evaluate them, build candidate and trace tables | 4–25-aa enumeration versus 8–25-aa length check, one row per candidate instance, 12 rule rows |
| `default_rules` | 903 | Provide default `ignore/apply/block` states | Defaults exactly match interface, documentation and tests |
| `validate_rule_configuration` | 923 | Check names and state values for all 12 rules | Unknown rules, missing names and invalid states |

### 3.4 Workbook and Excel support

| Function | Approx. line | Responsibility | Critical checks |
|---|---:|---|---|
| `make_export_filename` | 85 | Generate an accession- and time-based xlsx filename | Valid names, collision considerations |
| `empty_candidate_table` | 963 | Define typed schema for empty candidate output | Schema matches nonempty output |
| `empty_rule_trace_table` | 996 | Define typed schema for empty rule trace | Schema matches nonempty trace |
| `build_seeknsee_workbook` | 1016 | Assemble seven named workbook sheets | Schema, sorting, counts and score consistency |
| `truncate_delimited_identifier_list` | 1103 | Abbreviate long semicolon-delimited ID lists for Excel | Character budget, count and audit of omitted IDs |
| `adapt_sheet_identifier_lists_for_excel` | 1183 | Adapt selected workbook columns and collect audit records | Only permitted columns changed; accurate audit |
| `adapt_workbook_for_excel` | 1272 | Adapt all sheets and validate cell lengths | Preserve full scientific evaluation; warn of display truncation |
| `empty_excel_preflight_issues` | 1367 | Provide typed empty error table | Column contracts |
| `validate_excel_workbook` | 1381 | Detect cells exceeding the configured character limit | Boundary 32,767 characters, very large sheets |
| `export_seeknsee_workbook` | 1441 | Create output directory and write workbook with `writexl` | File integrity, writexl failures, sheet size limits |

### 3.5 Validation helpers and interactive application

| Function/object | Approx. line | Responsibility | Critical checks |
|---|---:|---|---|
| `run_validation_case` | 1452 | Run headless reference/target workflow with explicit parameters | Input selection, deterministic candidate and trace outputs |
| `summarize_validation_result` | 1576 | Report descriptive output counts, scores and rule failures | Empty outputs, counting semantics |
| `validate_internal_identities` | 1669 | Check candidate, retained/blocked, trace and score identities | Multiple candidate instances and zero-output cases |
| `evaluate_test_peptide` | 1718 | Inspect one candidate's complete 12-rule evaluation | Individual rules and state combinations |
| `compare_digestion` | 1765 | Compare a digestion result with expected peptides | Every enzyme and missed-cleavage setting where supported |
| `ui` | ~1790–2303 | Render Shiny page, about/rules sections and controls | Controls match documented algorithm and defaults |
| `server` | 2304 | Coordinate uploads, validation, reactive state, visualization and downloads | Session lifecycle, stale state, file upload errors, download parity |
| Nested `close_progress` | 2330 | Close upload progress safely | Error/early exit cleanup |
| Nested `recompute` | 2523 | Recompute candidate results from reactive inputs | Same result as headless runner for identical inputs |

Two nested `contextual_value` helpers occur in the export section (reported around lines 1240 and 1411); they are implementation details, not public functions.

## 4. Scientific contracts to document explicitly

1. **Peptide instances:** reference-based matching is a literal substring search across reference sequences (`grepl(..., fixed=TRUE)`), not a full digest-aware identification model. The implementation distinguishes matching entries from within-entry occurrence counts.
2. **Digestion:** the selected enzyme defines cleavage boundaries; contiguous cleavage fragments can be combined to create 0, 1 or 2 missed-cleavage candidates. The candidate instance records one-based inclusive start and end coordinates.
3. **Length:** `make_target_candidates` retains generated instances of length **4–25 aa** for rule evaluation; the configurable `length` rule passes only **8–25 aa**. The distinction is intentional and must not be collapsed in a refactor.
4. **Score and eligibility:** all 12 rules are evaluated. A failed `apply` rule normally contributes −1 (or −1/−2 for one/two missed cleavages); a failed `block` rule sets `status = blocked` but does not directly decrease the numerical score. `ignore` still records the check outcome.
5. **Applicability:** the D/E tryptic-boundary assessment is Trypsin-specific; dibasic check applies to Trypsin and Trypsin/P; the P1′ proline check is Trypsin/P- and context-specific.
6. **Charge:** `z_hat = 1 + count(K) + count(R) + count(H)` is a heuristic and should not be described as measured charge or a physical charge-state predictor.
7. **Reference specificity:** the interpretation depends on whether the target comes from the reference proteome or an external FASTA with strict absence versus full-sequence identity allowance.
8. **Output:** the workbook contains `Run_summary`, `Rule_configuration`, `Target_proteins`, `Ranked_peptides`, `Rule_trace`, `Blocked_peptides`, and `Rule_definitions`.

## 5. Release blockers and targeted verification

These are **risks requiring verification**, not claims that a defect has been experimentally demonstrated.

| Priority | Finding | Why it matters | Verification/action before release |
|---|---|---|---|
| P0 | The 23 passing starter checks do not reproduce the manuscript's reported 198 assertions | Published validation claim must map to accessible code/evidence | Locate original test harness or label starter checks accurately; archive reproducible verification results |
| P0 | Excel adaptation checks the **cell-character limit**, but the supplied implementation does not visibly check the **worksheet-row limit** | The benchmark reports millions of rule-trace rows; Excel has a maximum of 1,048,576 rows per worksheet | Add boundary tests; establish explicit behavior for over-limit exports (fail with explanation, split or use an alternative format) in a separately reviewed change |
| P0 | `README.md`, install instructions, packaging folders and `tests/packaging` are not yet verified by clean-machine runs | A public release must be installable without maintainer help | Test installation on at least one clean supported environment and record steps; document other platforms as verified or unverified |
| P1 | `find_reference_matches` rescans the full reference vector for each uncached peptide | Potentially expensive for a large reference and many candidates | Establish memory/time benchmark before optimization; compare output hashes after any change |
| P1 | `make_target_candidates` grows per-candidate and per-rule lists inside nested loops | High memory allocation overhead for large datasets | Benchmark memory and output before considering batching/preallocation |
| P1 | Shiny server uses shared `active_sessions` and calls `stopApp()` when the last session ends | Could be suitable for an intentionally local single-user launcher, but unsafe for multi-user/server hosting | Document deployment model; test local launcher and hosted Shiny separately before altering lifecycle semantics |
| P1 | `renderUI` search uses `grepl(q, ...)` without `fixed=TRUE` | User input is currently interpreted as a regular expression and may cause regex errors | Document search semantics or test literal search and modify only with explicit UI behavior acceptance |
| P1 | Reproducibility metadata include versions/file names/time but no displayed FASTA digests or dependency lock | Same filenames do not guarantee identical input bytes | Consider input SHA-256 recording and `renv.lock`; test restoring the environment |
| P2 | Nested helpers / `has_internal_cleavage_site` may be candidates for consolidation | Avoid dead code, but deletion could affect unseen callers | Search real references and regression tests before removal |

**Important:** Excel identifiers being truncated in *displayed lists* is already explicitly audited by current helper functions. This is distinct from row-count overflow and should not be described as silent loss of the scientific score or evaluation.

## 6. Proposed modular layout (after reference snapshots)

Keep the existing repository's top-level layout and `app/app.R` entry point; extract stable, Shiny-independent helpers only after tests exist:

```text
SeekNSee/
├── app/
│   ├── app.R                  # Shiny entry point and light wiring
│   ├── R/
│   │   ├── data.R            # FASTA, identifiers, typed empty tables
│   │   ├── digestion.R       # Enzyme boundaries and missed cleavages
│   │   ├── specificity.R     # Reference matches and context checks
│   │   ├── rules.R           # Definitions, defaults, rule evaluation
│   │   ├── candidates.R      # Candidate and trace generation
│   │   ├── reporting.R       # Workbook construction and safe export
│   │   └── validation.R      # Headless runner and accounting checks
│   └── www/
├── docs/
├── packaging/
├── tests/
└── .github/
```

Do **not** move `build_rules` into the scientific core: it directly reads Shiny inputs and belongs to presentation/orchestration. Keep an explicit deterministic source order or package namespace; document it rather than relying on implicit attach behavior.

## 7. Regression program before extraction

- **Reference case A:** synthetic FASTA with two distinct reference entries and an external target; verify exact candidate table and per-rule trace with stable sorting.
- **Reference case B:** Trypsin/Trypsin/P discriminating K–P and R–P boundaries, including 0/1/2 missed-cleavage enumeration and coordinate mapping.
- **Reference case C:** all rule states (`ignore`, `apply`, `block`) and non-applicability, including simultaneous block reasons, terminal peptides and 4–7 aa candidates.
- **Reference case D:** identical target/reference sequence, duplicated normalized accessions, multiple sequence occurrences, and cache reuse/invalidation.
- **Reference case E:** zero candidates, multiple targets, Excel cell over-limit and row-count boundary behavior.
- **UI smoke:** load FASTA, select target, run analysis, update rules, export workbook, and end/reopen a session.

For deterministic comparisons, exclude or normalize `Run_summary` fields such as run timestamp and elapsed time, and normalize file-output paths. The **candidate**, **rule trace**, **status**, **score** and **rule definitions** should be compared directly and strictly. Save fixture FASTA, expected outputs, their digests and test commands in version control.

## 8. Immediate next step

1. Review this draft with the maintainer, especially intentional semantics, deployment model and publication claims.
2. Preserve `app/app.R` unchanged and extend the baseline test suite to check full, ordered outputs for small synthetic fixtures.
3. Add R documentation comments (`#'`) to the highest-level scientific functions, focusing on parameters, returns, invariants and why the method behaves as it does; only then consider extraction to `app/R/`.
4. After each isolated change, run syntax check, all baseline checks and deterministic output comparisons; stage logical changes, but do not commit or push until jointly reviewed.

---

*This file is a review draft and must not be interpreted as test certification, completed refactoring or maintainer approval.*
