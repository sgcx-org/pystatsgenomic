# PyStatsGenomic

**Genomics and computational-biology statistical computing for Python.**

> **Status: early / reserved (`0.0.1`).** This package reserves the
> `pystatsgenomic` name within the **PyStatistics open-core ecosystem** and will
> grow into a full genomics statistics library. The public API is forthcoming.

PyStatsGenomic is part of the open-core PyStatistics family:

| Package | Layer |
|---|---|
| [`pystatistics`](https://github.com/sgcx-org/pystatistics) | Fundamental, general statistics |
| [`pystatsbio`](https://github.com/sgcx-org/pystatsbio) | Biotech / pharma statistics |
| **`pystatsgenomic`** | Genomics / computational-biology statistics |

Like its siblings, it builds on `pystatistics` for the general statistical layer
and adds methods specific to genomics.

## Planned scope (candidates, not commitments)

- Hardy–Weinberg equilibrium test for a SNP, with allele/genotype frequencies.
- Genomics-aware helpers for high-dimensional differential-abundance summaries.

## Installation

```bash
pip install pystatsgenomic
```

## License

MIT © Hai-Shuo. Part of the [SGCX](https://sgcx.org) open-core ecosystem.
