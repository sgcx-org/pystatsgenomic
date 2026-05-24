"""Tests for hwe.hardy_weinberg and hwe.allele_frequencies."""

import pytest

from pystatsgenomic import hwe


# ---------------------------------------------------------------------------
# allele_frequencies
# ---------------------------------------------------------------------------

def test_allele_frequencies_basic():
    # freq(A) = (2*36 + 48) / 200 = 0.6
    f = hwe.allele_frequencies(36, 48, 16)
    assert f.p == pytest.approx(0.6)
    assert f.q == pytest.approx(0.4)
    assert f.n == 100
    assert f.p + f.q == pytest.approx(1.0)


def test_allele_frequencies_monomorphic():
    f = hwe.allele_frequencies(100, 0, 0)
    assert f.p == pytest.approx(1.0)
    assert f.q == pytest.approx(0.0)


def test_allele_frequencies_all_zero_raises():
    with pytest.raises(ValueError, match="all zero"):
        hwe.allele_frequencies(0, 0, 0)


def test_allele_frequencies_negative_raises():
    with pytest.raises(ValueError, match="non-negative"):
        hwe.allele_frequencies(-1, 0, 0)


# ---------------------------------------------------------------------------
# hardy_weinberg — normal cases
# ---------------------------------------------------------------------------

def test_in_equilibrium_high_p_value():
    # Exact HWE counts for p=0.6: chi-square is 0, p-values are 1.0.
    r = hwe.hardy_weinberg(36, 48, 16, method="chisq")
    assert r.statistic == pytest.approx(0.0, abs=1e-9)
    assert r.p_value == pytest.approx(1.0)
    assert r.expected == pytest.approx((36.0, 48.0, 16.0))

    r_exact = hwe.hardy_weinberg(36, 48, 16, method="exact")
    assert r_exact.p_value == pytest.approx(1.0)


def test_strong_disequilibrium_low_p_value():
    # No heterozygotes at all -> strong departure from HWE.
    r = hwe.hardy_weinberg(50, 0, 50, method="chisq")
    assert r.statistic == pytest.approx(100.0)
    assert r.p_value < 1e-10

    r_exact = hwe.hardy_weinberg(50, 0, 50, method="exact")
    assert r_exact.p_value < 1e-3


def test_exact_small_sample_reference():
    # (AA=1, Aa=0, aa=1): given 2 A and 2 a alleles, P(het=0) = 1/3,
    # P(het=2) = 2/3; the exact p-value for the observed het=0 is 1/3.
    r = hwe.hardy_weinberg(1, 0, 1, method="exact")
    assert r.p_value == pytest.approx(1.0 / 3.0, abs=1e-9)


def test_chisq_statistic_small_sample():
    # (1, 0, 1): expected (0.5, 1.0, 0.5) -> chi2 = 0.5 + 1.0 + 0.5 = 2.0
    r = hwe.hardy_weinberg(1, 0, 1, method="chisq")
    assert r.statistic == pytest.approx(2.0)
    assert r.p_value == pytest.approx(0.15729920705, abs=1e-6)


# ---------------------------------------------------------------------------
# Edge cases
# ---------------------------------------------------------------------------

def test_monomorphic_locus():
    r = hwe.hardy_weinberg(100, 0, 0, method="exact")
    assert r.p == pytest.approx(1.0)
    assert r.q == pytest.approx(0.0)
    assert r.statistic == pytest.approx(0.0)
    assert r.p_value == pytest.approx(1.0)


def test_default_method_is_exact():
    r = hwe.hardy_weinberg(36, 48, 16)
    assert r.method == "exact"


def test_summary_is_string():
    r = hwe.hardy_weinberg(36, 48, 16)
    assert isinstance(r.summary(), str)
    assert "Hardy-Weinberg" in r.summary()


# ---------------------------------------------------------------------------
# Failure cases
# ---------------------------------------------------------------------------

def test_negative_count_raises():
    with pytest.raises(ValueError, match="non-negative"):
        hwe.hardy_weinberg(-5, 48, 16)


def test_all_zero_raises():
    with pytest.raises(ValueError, match="all zero"):
        hwe.hardy_weinberg(0, 0, 0)


def test_non_int_raises():
    with pytest.raises(TypeError):
        hwe.hardy_weinberg(36.0, 48, 16)


def test_bad_method_raises():
    with pytest.raises(ValueError, match="method"):
        hwe.hardy_weinberg(36, 48, 16, method="bogus")
