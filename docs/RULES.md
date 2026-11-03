# SeekNSee rules (v0.2.5)

> **Documentation draft — derived from the previously shared `app.R` v0.2.5.** Confirm against the Git checkout and obtain maintainer review before publication. Rules are deterministic selection heuristics, not quantitative MS detectability predictors.

## Rule-state semantics

| State | Failed, applicable criterion | Score impact | Candidate eligibility |
| --- | --- | --- | --- |
| `ignore` | Recorded but disregarded for decisions | 0 | Unchanged |
| `apply` | Recorded and penalized | Usually −1; `miscleavage`: −1 or −2 | Unchanged |
| `block` | Recorded as an exclusion | 0 from this rule | `blocked` |

Nonapplicable criteria do not penalize or block under any state. The score starts at **12**, the number of rules, and eligibility is determined independently from the score. Default `block` states: `unique_proteome`, `length`. All other rules default to `apply`.

## The twelve rules

| Rule key | Pass condition | Applicability / important interpretation |
| --- | --- | --- |
| `unique_proteome` | No disallowed reference entries contain the peptide as an exact substring | Allowed matches depend on reference target vs external target and external `strict_absence` vs `allow_exact_target_match` policy. No enzyme-specific reference digest |
| `unique_protein` | Matches belong to at most one normalized reference base-accession group | Based on parsed FASTA identifiers, **not** guaranteed biological gene uniqueness |
| `length` | Peptide has **8–25 residues** | Evaluated candidates are already restricted to **4–25 residues** |
| `dibasic` | No `KK`, `KR`, `RK` or `RR` motif | Applicable to `Trypsin` and `Trypsin/P` only |
| `glyco` | No `N-X-S/T` motif where X ≠ P (`N[^P][ST]`) | Sequence-level motif flag; does not measure glycosylation |
| `chemical` | No M, W, C, `NG`, `QG`, `DP` or N-terminal Q/E | Heuristic flags for potentially problematic sequence features, not observed chemistry |
| `acidic` | No D/E within the implemented P3–P3′ vicinity of a valid K/R cleavage boundary delimiting the candidate | Applicable only to `Trypsin`; K/R–P is not a valid tryptic boundary. The implementation inspects positions cut−2, cut−1, cut+1, cut+2, cut+3 (the cleavage residue itself is K/R) |
| `terminal` | Peptide is neither N-terminal nor C-terminal in its target protein | Fails when candidate start is 1 or end equals protein length |
| `miscleavage` | Candidate spans **zero** missed cleavage boundaries | Applicable to every candidate. Under `apply`, penalty equals the number of missed cleavages (1 or 2) |
| `proline_p1prime` | P1′ residue is not P at an internal K/R boundary | Applicable only to `Trypsin/P` candidates starting after an internal K/R; otherwise not applicable |
| `isobaric` | Peptide contains neither I nor L | Flags I/L ambiguity; not a spectral or isotopic model |
| `charge` | `charge_hat <= 2`, where `charge_hat = 1 + nK + nR + nH` | A charge heuristic, **not** an experimentally observed charge state |

## Typical interpretation pitfalls

- **A high score does not mean the candidate passed exclusion rules.** `block` does not decrease score; always inspect `candidate_status` and `blocked_by`.
- **A nonapplicable rule is not a passed rule.** Its `pass_fail` value is `not_applicable`.
- **Candidate length has two thresholds.** Enumeration filter: 4–25 aa; prioritization rule: 8–25 aa.
- **Specificity is based on sequence occurrence.** It does not model whether a shared substring is produced by the selected protease in the reference entry.
- **The score does not represent a probability** and should not be compared across enzyme, rule-state or reference-policy configurations without retaining the full settings and trace.

## Example

For a synthetic 8-aa peptide `AAAAAAAK` evaluated against an empty reference with default rules and `Trypsin`, the peptide is retained, but the `terminal` rule fails (the target sequence and peptide coincide). Its resulting score is **11**, not 12. A 7-aa terminal peptide may be **blocked by `length`** while having a high score; this demonstrates the separation of eligibility and scoring.

## Scientific review requirements

The rule definitions must be reconciled with the Methods, supplementary rule table, GUI labels, and exported `Rule_definitions` sheet. Any change to cleavage behavior, reference-matching semantics, rule thresholds, penalty weights or default states is a **scientific change**, not a documentation-only edit, and requires updated regression expectations and an algorithm-version decision.
