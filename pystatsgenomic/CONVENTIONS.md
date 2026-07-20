# PyStatsGenomic Conventions (adopt-and-extend)

This document governs the **pystatsgenomic public API**. It does not restate the
law from scratch — it **adopts the pystatistics constitution as binding** and adds
a small set of genomic-specific amendments (G-series) for this library's domain
(statistical and population genetics: HWE, linkage disequilibrium, association,
population structure).

## 0. Adoption

The binding base is **`pystatistics/CONVENTIONS.md`**. Everything in it applies
to pystatsgenomic verbatim: the naming law S0–S6, the selector taxonomy, the
backend & precision convention, the result-object conventions
(`…Solution` wrapping `Result[…Params]`, `core.result.SolutionReprMixin`), the
exception conventions, and amendments **A1–A15**.

**A15 matters here more than anywhere else in the ecosystem,** because genomics
*uses* an unusually large amount of general statistics (regression per locus,
multiple-testing correction, PCA) while contributing genetic framing. The pull to
re-implement upstream machinery "because GWAS needs it" is the single most likely
way this package goes wrong.

When this document and the base disagree, **this document wins for
pystatsgenomic**; where this document is silent, the base governs. Like the base,
this document is **self-amending**: a new ambiguity is resolved once, here, as a
numbered amendment (G-series), not re-litigated per occurrence.

Reuse, don't fork: pystatsgenomic imports `pystatistics.core.exceptions`,
`pystatistics.core.result`, and `pystatistics.core.compute.backend` rather than
defining its own parallels. One hierarchy across the ecosystem.

---

## pystatsgenomic amendments (G-series)

### G1 — Every public return is a `…Solution`

`HWEResult` (`hwe/_common.py`) becomes **`HWESolution`** and `AlleleFrequencies`
becomes **`AlleleFrequenciesSolution`** — both wrapping frozen `…Params` payloads
in `core.result.Result` with the uniform metadata accessors, `summary()`, and
`_repr_html_` via `SolutionReprMixin`. `AlleleFrequencies` is a *top-level public
return* (`allele_frequencies()` returns it directly), so it takes the envelope; a
value object that only ever appears nested inside another Solution would not.

Note `HWEResult` carries the `…Result` suffix the base renames to `…Solution` —
fix it now, before there are users.

### G2 — Exceptions come from `pystatistics.core.exceptions`

No bare `raise ValueError` / `raise TypeError` for validation. The current `hwe/`
modules raise three `ValueError` and one `TypeError`; all become
`ValidationError`. Fit-quality / non-convergence failures →`ConvergenceError`.

### G3 — Dependency peg

`pyproject.toml` currently pegs **`pystatistics>=3.0.1`** — badly stale (the
current base is 5.1, and 5.0 hard-renamed the exception taxonomy and result
objects this document depends on). It declares **`pystatistics>=5.1`**.

### G4 — What delegates, and what is genuinely genomic (A15 applied)

Delegates upstream — **this list is the whole point of the amendment**:

- **Multiple-testing correction (Benjamini–Hochberg / FDR / Bonferroni)** —
  `hypothesis.p_adjust`. This is the classic trap: FDR is general statistics and
  does **not** become genomic because GWAS is its most famous consumer. Never
  implement it here.
- **Single-locus association** — `regression` (linear or logistic per SNP). A GWAS
  scan orchestrates many upstream fits; it does not contain its own.
- **PCA for population stratification** — `multivariate`. The genomic value is the
  ancestry *interpretation* of the components, not the decomposition.
- **Chi-square goodness-of-fit** — `hypothesis.chisq_test` where it fits; the
  expected genotype counts are computed here, the test is not re-derived.

Stays local — the genetic contribution:

- The **Hardy–Weinberg exact test** (Wigginton, Cui & Abecasis 2005) — a
  domain-specific combinatorial test with no upstream equivalent.
- Allele and genotype frequency computation, linkage disequilibrium (D, D′, r²),
  Fst, genomic inflation (λ_GC), diversity and neutrality statistics, TDT.

Not a violation: `scipy.stats.chi2.sf` for a p-value quantile (A15).

**Cochran–Armitage trend** is a promotion candidate, not local work: it is a
general categorical trend test whose dominant *use* is per-SNP. Promote upstream,
wrap here.

**DYR (don't reimplement) also binds sideways:** sequence bioinformatics (GC
content, codon usage, motifs) is Biopython's job, and differential expression has
mature tooling. Neither is claimed here merely because it is genomic.

### G5 — Genotype-count vocabulary exemption (documented S1/A1 carve-out)

The genotype-count parameters **`n_AA`, `n_Aa`, `n_aa`** keep their genetics
spelling as a **documented exemption** to S1/A1, on exactly the A14 reasoning:
they *are* the canonical vocabulary a geneticist already knows, the capitalization
is semantic (major/minor allele, not shouting), and any "corrected" spelling would
be less clear, not more. The same applies to established locus/allele symbols
where they are the reference vocabulary.

This carve-out is **closed**: it covers genotype counts and established allele
notation only. A new abbreviated public identifier that is not established genetics
vocabulary is an S1/A1 violation to fix, not an exemption to assume.

### G6 — First-pass scope discipline

The first pass implements **Tier 1 of `ROADMAP.md` only**. Items in the roadmap's
*"Ambiguous — discuss before building"* section stay parked and are not decided by
implementation.
