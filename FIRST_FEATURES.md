# First Features — pystatsgenomic

This package is currently a **reserved `0.0.1` skeleton**. This file specifies the
first 1–2 features a future session should implement to turn it into a real
`0.1.0` (the "harden later" step of SGC-Bio roadmap item B-3).

## Ground rules (read first)

- **Stay in this layer.** Only implement things that are genuinely *genomic /
  computational-biology*. General statistics belong in `pystatistics`. The
  classic trap: **multiple-testing correction (Benjamini–Hochberg / FDR) is
  general stats and belongs in `pystatistics`, NOT here** — don't grab it just
  because genomics uses it. Apply the "which-of-the-4-layers" test.
- **Build on `pystatistics`** for general tests/CIs; don't reimplement them.
- **Follow `pystatsbio`'s conventions** — hatchling packaging, the Coding Bible
  (`CLAUDE.md` there: fail loud, one job per module, tests first, deterministic),
  typed, ruff/mypy clean.
- **Ship with tests** (normal / edge / failure) and bump to `0.1.0` once landed.

## Feature 1 (primary): Hardy–Weinberg equilibrium

Unambiguously genomic, self-contained, and a natural first function.

- Suggested API: `hardy_weinberg(n_AA, n_Aa, n_aa, *, method="exact") ->
  HWEResult`
- Compute:
  - allele frequencies p (major) and q (minor)
  - expected genotype counts under HWE (`p²`, `2pq`, `q²` × N)
  - **chi-square test** (1 df) of observed vs expected, and/or
  - **exact test** (Wigginton, Cui & Abecasis 2005) — preferred for small samples
    / rare alleles; make `method` select `"exact"` vs `"chisq"`.
- Result object fields: `p, q, observed (AA/Aa/aa), expected (AA/Aa/aa),
  statistic, p_value, method, summary`.
- Failure behavior (fail loud): reject negative/non-integer counts and an
  all-zero sample with explicit errors.
- Reference: Wigginton JE, Cui J, Abecasis GR, "A Note on Exact Tests of
  Hardy–Weinberg Equilibrium," Am J Hum Genet 2005.
- Tests: a known in-equilibrium sample → high p-value; strong disequilibrium →
  low p-value; the exact-test reference example; the failure cases.

## Feature 2 (optional secondary)

The allele/genotype-frequency computation above is reusable on its own — expose
it as a small public helper (`allele_frequencies(n_AA, n_Aa, n_aa)`) if useful.
Defer anything heavier (GWAS scans, LD) — out of scope for `0.1.0`.
