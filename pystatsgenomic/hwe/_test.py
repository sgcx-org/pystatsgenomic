"""Hardy-Weinberg equilibrium test for a biallelic locus.

Provides both the classical Pearson chi-square goodness-of-fit test (1 df) and
the exact test of Wigginton, Cui & Abecasis (2005), which is preferred for
small samples or rare alleles.

Reference
---------
Wigginton JE, Cui J, Abecasis GR. "A Note on Exact Tests of Hardy-Weinberg
Equilibrium." Am J Hum Genet 2005; 76:887-893.
"""

from __future__ import annotations

from typing import Literal

from pystatistics.hypothesis import chisq_test
from scipy import stats

from pystatsgenomic.hwe._common import HWEResult
from pystatsgenomic.hwe._frequencies import _validate_counts


def _chi_square(
    n_AA: int, n_Aa: int, n_aa: int, n: int,
) -> tuple[float, tuple[float, float, float]]:
    """Pearson chi-square statistic and expected genotype counts under HWE.

    The genotype proportions expected under HWE (p^2, 2pq, q^2) are
    domain-specific; the goodness-of-fit statistic itself is delegated to
    ``pystatistics.hypothesis.chisq_test`` rather than re-implemented here. The
    p-value is *not* taken from that call: a generic 3-category GOF uses 2 df,
    whereas HWE uses 1 df (an allele frequency is estimated), so the caller
    computes the p-value separately (see ``hardy_weinberg``).
    """
    freq_a = (2 * n_AA + n_Aa) / (2 * n)
    freq_b = 1.0 - freq_a
    props = (freq_a * freq_a, 2.0 * freq_a * freq_b, freq_b * freq_b)

    if freq_a == 0.0 or freq_b == 0.0:
        # Monomorphic locus: an expected proportion is zero and so is the
        # matching observed count, giving a statistic of 0. chisq_test would
        # divide by a zero expected count, so handle this directly.
        return 0.0, (props[0] * n, props[1] * n, props[2] * n)

    sol = chisq_test([n_AA, n_Aa, n_aa], p=list(props))
    expected = (float(sol.expected[0]), float(sol.expected[1]), float(sol.expected[2]))
    return float(sol.statistic), expected


def _exact_p_value(n_AA: int, n_Aa: int, n_aa: int, n: int) -> float:
    """Exact HWE p-value via the Wigginton, Cui & Abecasis (2005) recurrence.

    The conditional distribution of heterozygote counts given the allele count
    is computed by an exact recurrence; the p-value is the summed probability
    of all heterozygote configurations no more likely than the observed one.
    """
    obs_hets = n_Aa
    rare = 2 * min(n_AA, n_aa) + n_Aa  # number of rare-allele copies
    if rare == 0:
        # Monomorphic locus: trivially in equilibrium.
        return 1.0

    probs = [0.0] * (rare + 1)

    # Start at the most probable heterozygote count with matching parity.
    mid = (rare * (2 * n - rare)) // (2 * n)
    if mid % 2 != rare % 2:
        mid += 1
    probs[mid] = 1.0
    total = 1.0

    # Recurse downward in heterozygote count (steps of 2).
    curr_hets = mid
    curr_hom_r = (rare - mid) // 2
    curr_hom_c = n - curr_hets - curr_hom_r
    while curr_hets > 1:
        probs[curr_hets - 2] = (
            probs[curr_hets]
            * curr_hets
            * (curr_hets - 1.0)
            / (4.0 * (curr_hom_r + 1.0) * (curr_hom_c + 1.0))
        )
        total += probs[curr_hets - 2]
        curr_hom_r += 1
        curr_hom_c += 1
        curr_hets -= 2

    # Recurse upward in heterozygote count.
    curr_hets = mid
    curr_hom_r = (rare - mid) // 2
    curr_hom_c = n - curr_hets - curr_hom_r
    while curr_hets <= rare - 2:
        probs[curr_hets + 2] = (
            probs[curr_hets]
            * 4.0
            * curr_hom_r
            * curr_hom_c
            / ((curr_hets + 2.0) * (curr_hets + 1.0))
        )
        total += probs[curr_hets + 2]
        curr_hom_r -= 1
        curr_hom_c -= 1
        curr_hets += 2

    p_obs = probs[obs_hets]
    tail = sum(pr for pr in probs if pr <= p_obs)
    return min(tail / total, 1.0)


def hardy_weinberg(
    n_AA: int,
    n_Aa: int,
    n_aa: int,
    *,
    method: Literal["exact", "chisq"] = "exact",
) -> HWEResult:
    """Test a biallelic locus for Hardy-Weinberg equilibrium.

    Parameters
    ----------
    n_AA, n_Aa, n_aa : int
        Counts of the three genotypes (homozygous A, heterozygous, homozygous a).
    method : {'exact', 'chisq'}
        'exact' (default) uses the Wigginton et al. (2005) exact test, preferred
        for small samples or rare alleles. 'chisq' uses the Pearson chi-square
        test (1 df). The chi-square statistic is reported either way.

    Returns
    -------
    HWEResult

    Raises
    ------
    TypeError
        If any count is not an int.
    ValueError
        If any count is negative, the sample is empty, or method is unknown.
    """
    if method not in ("exact", "chisq"):
        raise ValueError(f"method must be 'exact' or 'chisq', got {method!r}")
    n = _validate_counts(n_AA, n_Aa, n_aa)

    freq_a = (2 * n_AA + n_Aa) / (2 * n)
    freq_b = 1.0 - freq_a
    p = max(freq_a, freq_b)
    q = min(freq_a, freq_b)

    statistic, expected = _chi_square(n_AA, n_Aa, n_aa, n)

    if method == "exact":
        p_value = _exact_p_value(n_AA, n_Aa, n_aa, n)
    else:
        p_value = float(stats.chi2.sf(statistic, df=1))

    return HWEResult(
        p=p,
        q=q,
        observed=(n_AA, n_Aa, n_aa),
        expected=expected,
        statistic=statistic,
        p_value=p_value,
        method=method,
        n=n,
    )
