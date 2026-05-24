"""Result types for Hardy-Weinberg equilibrium analysis."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class AlleleFrequencies:
    """Allele frequencies for a biallelic locus.

    Attributes
    ----------
    p : float
        Frequency of the major (more common) allele.
    q : float
        Frequency of the minor (less common) allele. p + q == 1.
    n : int
        Number of individuals genotyped.
    """

    p: float
    q: float
    n: int

    def summary(self) -> str:
        """Human-readable one-line summary."""
        return f"Allele frequencies (n={self.n}): p={self.p:.4f}, q={self.q:.4f}"


@dataclass(frozen=True)
class HWEResult:
    """Result of a Hardy-Weinberg equilibrium test for a biallelic locus.

    Attributes
    ----------
    p : float
        Frequency of the major allele.
    q : float
        Frequency of the minor allele.
    observed : tuple of int
        Observed genotype counts (AA, Aa, aa) as supplied.
    expected : tuple of float
        Expected genotype counts under HWE (AA, Aa, aa).
    statistic : float
        Pearson chi-square goodness-of-fit statistic (1 df), reported for both
        methods as a descriptive measure of departure from HWE.
    p_value : float
        p-value. Source depends on ``method``: the exact test (Wigginton et al.
        2005) for 'exact', the chi-square test (1 df) for 'chisq'.
    method : str
        'exact' or 'chisq'.
    n : int
        Number of individuals genotyped.
    """

    p: float
    q: float
    observed: tuple[int, int, int]
    expected: tuple[float, float, float]
    statistic: float
    p_value: float
    method: str
    n: int

    def summary(self) -> str:
        """Human-readable summary."""
        return "\n".join(
            [
                f"Hardy-Weinberg Equilibrium ({self.method})",
                "=" * 50,
                f"n             : {self.n}",
                f"Allele freq   : p={self.p:.4f}, q={self.q:.4f}",
                f"Observed      : AA={self.observed[0]}, "
                f"Aa={self.observed[1]}, aa={self.observed[2]}",
                f"Expected      : AA={self.expected[0]:.2f}, "
                f"Aa={self.expected[1]:.2f}, aa={self.expected[2]:.2f}",
                f"Chi-square    : {self.statistic:.4f} (1 df)",
                f"p-value       : {self.p_value:.4g}",
            ]
        )
