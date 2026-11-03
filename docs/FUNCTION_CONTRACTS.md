# SeekNSee — Scientific function contracts (v0.2.5)

**Status:** Documentation draft, for review before adding inline comments. **Source:** The unchanged `app/app.R` v0.2.5 function definitions inspected in the collaboration. This file does not modify the application or validate its empirical predictive performance.

## Data ingestion — `read_fasta(path, progress_callback = NULL)`

**Purpose.** Read FASTA/FAA data into a normalized protein-entry table.

**Inputs.** `path`: readable FASTA file. `progress_callback`: optional function called with named arguments `completed` and `total` during parsing.

**Returns.** Data frame with `id`, `original_id`, `accession`, `base_id`, `isoform`, `gene`, `name`, `sequence`, `sequence_length`.

**Contract and interpretation.** Removes a UTF-8 BOM, blank lines and whitespace from sequence text; uppercases sequences; accepts protein alphabet symbols plus `X B Z J U O` and terminal `*`; rejects internal stop characters, invalid residues, missing headers and empty records. Extracts UniProt-style accession/isoform when possible and disambiguates duplicate `id` values using a suffix. Do not equate normalized identifier uniqueness with biological sequence uniqueness.

**Regression obligations.** Valid/invalid FASTA; UniProt and custom headers; duplicate identifiers; isoforms; terminal `*`; single entry; multiple entries; progress callback. Verify exact output columns and order.

## Digestion — `digest_protein(sequence, enzyme, missed_cleavages = 0L)`

**Purpose.** Enumerate enzyme-bounded fragments and contiguous joins spanning 0–2 missed cleavage events.

**Inputs.** A protein sequence, one supported enzyme name, and a scalar integer missed-cleavage maximum from 0 to 2.

**Returns.** Data frame with `peptide`, `start`, `end`, `length`, `missed_cleavages`.

**Contract and interpretation.** One-based inclusive coordinates. Builds boundaries from `cleavage_positions()` then joins neighboring fragments up to the requested count. **This function does not itself apply the 4–25 aa candidate window**; the window is applied in `make_target_candidates()`.

**Regression obligations.** Nine enzymes, N/C terminal boundaries, Trypsin vs Trypsin/P K/R–P behavior, AspN cuts preceding D, no sites, 0/1/2 missed cleavages, exact order and coordinates.

## Reference matches — `find_reference_matches(peptide, reference_proteome, match_cache = NULL)`

**Purpose.** Count and identify exact peptide-substring occurrences in a reference protein table.

**Inputs.** Candidate sequence; reference-proteome data frame; optional cache environment keyed by peptide sequence.

**Returns.** List with `entry_count`, `occurrence_count`, `ids`, `genes`, `bases`.

**Contract and interpretation.** `entry_count` counts reference entries containing the substring; `occurrence_count` counts occurrences, including multiple occurrences within a protein. This is **substring matching**, not a separate enzyme-aware proteotypicity prediction. The cache is valid only for the specific reference proteome from which it was populated; a new reference requires a fresh cache.

**Regression obligations.** No hits, one hit, repeated substring, shared peptides in different reference entries, empty reference, stable order of IDs, cache equivalence and invalidation.

## Rule evaluation — `evaluate_peptide(...)`

**Purpose.** Evaluate all 12 selection rules for one candidate and return both numeric priority and eligibility.

**Principal inputs.** Candidate peptide; named rule-state list (`ignore`, `apply`, `block`); reference proteome; target ID/source/sequence; enzyme; one-based peptide boundaries; missed cleavage count; reference-specificity policy; optional cache and precomputed allowed IDs.

**Returns.** Named list: `valid`, `status`, `score`, `trace`, `pass`, `applicable`, `blocked_by`, `charge_hat`, `missed_cleavages`, detailed acidic/P1′ assessments, reference matches and allowed/disallowed IDs.

**Contract and interpretation.** All 12 rules are evaluated even if one blocks the peptide. Base score is 12; failed `apply` rules normally subtract 1, except missed cleavage subtracts its count (0–2). Failed `block` rules set eligibility but do **not directly** lower the numeric score. `ignore` records evidence without penalties. An inapplicable rule does not penalize or block. `charge_hat = 1 + nK + nR + nH` is only a heuristic. Reference-proteome specificity changes with reference-vs-external target source and external policy.

**Regression obligations.** All three rule states; failed and non-applicable rules; simultaneous blockers; length thresholds; exact target allowance; Trypsin vs Trypsin/P applicability; score/trace identity; no changes to numeric types.

## Candidate generation — `make_target_candidates(...)`

**Purpose.** Generate target-specific candidate instances, annotate them, and preserve the full rule-level audit trail.

**Principal inputs.** Target table, source label, enzyme, rule-state list, reference table and specificity policy, maximum missed cleavages, optional shared matching cache.

**Returns.** List of **complete** data frames `candidates`, `rule_trace`, `target_summary`.

**Contract and interpretation.** Retains 4–25 aa candidates from digestion, including 4–7 aa instances that subsequently fail the *length* rule (8–25 aa). Candidate keys incorporate target ID, start/end coordinates, missed cleavage count and sequence. There must be 12 rule-trace rows for each candidate. `target_summary` counts candidates and retained/blocked statuses. Candidate instances are not the same as unique peptide sequences.

**Regression obligations.** Empty and multi-target inputs; candidate count vs retention; 12× trace cardinality; exact full candidate/trace schemas; scores, statuses, positions, identifiers and row order.

## Review and scope control

These contracts are **descriptions of the observed implementation**, not new proposed scientific behavior. Once the `v0.2.5_full_outputs.rds` snapshot has been captured and verified, migrate the approved text into `#'` comments above each respective function in `app/app.R`, preserving all executable statements. Re-run the 23-check suite and complete comparison after every documented edit. If any semantics in these descriptions are disputed, correct documentation or create an explicitly reviewed behavioral-change issue rather than silently modifying code.
