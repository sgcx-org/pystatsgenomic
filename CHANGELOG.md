# Changelog

## 0.1.0

### Added
- `hwe` subpackage for Hardy-Weinberg equilibrium analysis of a biallelic locus:
  - `hwe.hardy_weinberg(...)`: equilibrium test via the Pearson chi-square test
    (1 df) or the Wigginton, Cui & Abecasis (2005) exact test (the default,
    preferred for small samples / rare alleles).
  - `hwe.allele_frequencies(...)`: major/minor allele-frequency helper.
  Fails loud on negative/non-integer counts and an empty sample.

### Changed
- Added `pystatistics` as a runtime dependency: the chi-square goodness-of-fit
  statistic is computed by `pystatistics.hypothesis.chisq_test` rather than
  re-implemented. The HWE-specific 1-df p-value is computed locally (a generic
  goodness-of-fit test would use 2 df).

First feature release, promoting the package from a reserved skeleton.
