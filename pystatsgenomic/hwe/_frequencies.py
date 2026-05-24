"""Allele and genotype frequencies for a biallelic locus."""

from __future__ import annotations

from pystatsgenomic.hwe._common import AlleleFrequencies


def _validate_counts(n_AA: int, n_Aa: int, n_aa: int) -> int:
    """Validate genotype counts and return the total number of individuals.

    Raises
    ------
    TypeError
        If any count is not an int.
    ValueError
        If any count is negative or the sample is all-zero.
    """
    for name, value in (("n_AA", n_AA), ("n_Aa", n_Aa), ("n_aa", n_aa)):
        if not isinstance(value, int) or isinstance(value, bool):
            raise TypeError(f"{name} must be an int, got {type(value).__name__}")
        if value < 0:
            raise ValueError(f"{name} must be non-negative, got {value}")
    n = n_AA + n_Aa + n_aa
    if n == 0:
        raise ValueError("genotype counts are all zero; no individuals to analyze")
    return n


def allele_frequencies(n_AA: int, n_Aa: int, n_aa: int) -> AlleleFrequencies:
    """Compute major/minor allele frequencies for a biallelic locus.

    Parameters
    ----------
    n_AA, n_Aa, n_aa : int
        Counts of the three genotypes (homozygous A, heterozygous, homozygous a).

    Returns
    -------
    AlleleFrequencies
        ``p`` is the larger (major) allele frequency, ``q`` the minor; p + q = 1.

    Raises
    ------
    TypeError
        If any count is not an int.
    ValueError
        If any count is negative or the sample is empty.
    """
    n = _validate_counts(n_AA, n_Aa, n_aa)
    freq_a_allele = (2 * n_AA + n_Aa) / (2 * n)
    freq_b_allele = (2 * n_aa + n_Aa) / (2 * n)
    p = max(freq_a_allele, freq_b_allele)
    q = min(freq_a_allele, freq_b_allele)
    return AlleleFrequencies(p=p, q=q, n=n)
