# PyStatsGenomic

**Genomics and computational-biology statistical computing for Python.**

> **Status: `0.1.0` — first module shipped.** Hardy–Weinberg equilibrium is
> available now; more statistical-genetics methods are on the roadmap below.

PyStatsGenomic is part of the open-core PyStatistics family:

| Package | Layer |
|---|---|
| [`pystatistics`](https://github.com/sgcx-org/pystatistics) | Fundamental, general statistics |
| [`pystatsbio`](https://github.com/sgcx-org/pystatsbio) | Biotech / pharma statistics |
| [`pystatsclinical`](https://github.com/sgcx-org/pystatsclinical) | Clinical-trial / clinical-research statistics |
| **`pystatsgenomic`** | Genomics / computational-biology statistics |
| [`pystatsfinance`](https://github.com/sgcx-org/pystatsfinance) | Financial / quantitative statistics |
| [`pystatsinsurance`](https://github.com/sgcx-org/pystatsinsurance) | Actuarial / insurance statistics |

Like its siblings, it builds on `pystatistics` for the general statistical layer
and adds methods specific to genomics.

## What's available now (0.1.0)

**Hardy–Weinberg equilibrium** for a biallelic locus, via the exact test
(Wigginton, Cui & Abecasis 2005, the default) or the Pearson chi-square test:

```python
from pystatsgenomic import hwe

result = hwe.hardy_weinberg(n_AA=298, n_Aa=489, n_aa=213)
print(result.summary())

# Allele frequencies on their own:
freqs = hwe.allele_frequencies(n_AA=298, n_Aa=489, n_aa=213)
```

`hardy_weinberg` reports allele frequencies, expected genotype counts, the
chi-square statistic, and the p-value (1 df — HWE estimates an allele
frequency). It fails loud on negative/non-integer counts and an empty sample.

## Roadmap (candidates, not commitments)

- Linkage disequilibrium (D, D′, r²) and genotyping QC metrics.
- Single-locus association and GWAS-scan orchestration (with genomic control).
- Population-genetics summaries: Fst, nucleotide diversity, Tajima's D.

## Installation

```bash
pip install pystatsgenomic
```

## License

MIT © Hai-Shuo. Part of the [SGCX](https://sgcx.org) open-core ecosystem.
