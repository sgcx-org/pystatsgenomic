# Roadmap — pystatsgenomic

A working **checklist** of features that plausibly belong at the *genomic /
statistical-genetics / population-genetics* layer. Tick items off as they land,
in roughly the way `pystatistics/docs/ROADMAP.md` tracks its modules.

**This list is not a commitment.** It captures *what kind of thing belongs here*
and in roughly what order we'd build it. Re-order, drop, or promote/demote
freely. Tiers are priority bands, not deadlines.

## Layer rules (the gate every item must pass)

- **Stay in this layer.** Only genuinely *genomic / statistical-genetics*
  methods. **The classic trap: multiple-testing correction (Benjamini–Hochberg /
  FDR) is general stats and belongs in `pystatistics`, NOT here** — don't grab it
  just because GWAS uses it.
- **Build on `pystatistics`** for general tests/CIs/regression — don't
  reimplement them. The genomic contribution is the *genetic framing*
  (genotypes, alleles, loci, populations), not the underlying test.
- **`pystatsgenomic` vs `pystatsbio`:** statistical/population genetics (HWE, LD,
  GWAS, Fst, diversity) is genomic; wet-lab / sequence-bioinformatics
  (assembly, alignment, motif/composition) leans `pystatsbio` — and where a
  public library already does it (Biopython), **DYR**: don't wrap it.
- Conventions: `pystatsbio` Coding Bible, typed, ruff/mypy clean, tests, version
  bump via the release flow.

## Tier 1 — self-contained, unmistakably genomic

- [ ] **Hardy–Weinberg equilibrium** — χ² and exact (Wigginton 2005) tests.
  *v0.1.0 target — full spec in [FIRST_FEATURES.md](FIRST_FEATURES.md).*
- [ ] **Allele / genotype frequencies** — small public helper. *See
  FIRST_FEATURES.md.*
- [ ] **Linkage disequilibrium** — D, D′, r² between locus pairs.
- [ ] **Genotyping QC metrics** — per-SNP / per-sample call rate, MAF, HWE-based
  filtering thresholds.

## Tier 2 — association & population structure

- [ ] **Single-locus association** — allelic / genotypic / Cochran–Armitage
  trend tests, framed per-SNP (effect size + p-value). *(Trend test home is
  ambiguous — see below.)*
- [ ] **Genomic inflation factor (λ_GC)** and basic genomic-control adjustment.
- [ ] **GWAS scan orchestration** — run association across many SNPs, collect
  per-SNP results, apply λ_GC; hand multiple-testing correction off to
  `pystatistics`.
- [ ] **Population differentiation (Fst)** — Weir & Cockerham / Nei.

## Tier 3 — population genetics & relatedness

- [ ] **Diversity / neutrality stats** — nucleotide diversity π, Watterson's θ,
  Tajima's D.
- [ ] **Relatedness / kinship** — IBS/IBD estimates, genomic relationship matrix
  (GRM).
- [ ] **Transmission disequilibrium test (TDT)** — family-based association.
- [ ] **Effective number of independent tests** (LD-aware) and/or PC-based
  population-stratification adjustment. *(PCA wrapper is ambiguous — see below.)*

## Ambiguous — discuss before building

- **Cochran–Armitage trend test** — a general categorical trend test, but its
  dominant use is per-SNP genotype association. General (`pystatistics`) with a
  genomic wrapper, or native here?
- **Differential expression (RNA-seq, DESeq2/edgeR-style NB GLM + shrinkage)** —
  heavy, and mature public tools already exist (**DYR**). If we do it at all, is
  it genomic or `pystatsbio`? Likely out of scope for now.
- **Sequence composition (GC content, codon usage, motif counts)** — this is
  bioinformatics, and Biopython covers it (**DYR**). Probably *not* here, maybe
  `pystatsbio`, maybe nowhere.
- **PCA for population stratification** — PCA itself lives in `pystatistics`; the
  genomic value is thin (ancestry interpretation of PCs). Worth a wrapper, or
  just call `pystatistics` directly?
- **Heritability / polygenic risk scores** — genomic, but heavy and
  model-dependent; park until the basics land and we know the use case.
