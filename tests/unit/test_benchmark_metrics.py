"""Unit tests for validation/benchmarks/metrics.py."""

import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from validation.benchmarks.metrics import (  # noqa: E402
    aafe,
    fold_error,
    within_n_fold,
)


class TestFoldError:
    def test_equal_values(self):
        assert fold_error(10.0, 10.0) == pytest.approx(1.0)

    def test_2_fold_over(self):
        assert fold_error(20.0, 10.0) == pytest.approx(2.0)

    def test_2_fold_under(self):
        assert fold_error(5.0, 10.0) == pytest.approx(2.0)

    def test_zero_predicted_raises(self):
        with pytest.raises(ValueError):
            fold_error(0.0, 10.0)

    def test_zero_observed_raises(self):
        with pytest.raises(ValueError):
            fold_error(10.0, 0.0)

    def test_negative_raises(self):
        with pytest.raises(ValueError):
            fold_error(-1.0, 10.0)


class TestAAFE:
    def test_identity(self):
        predicted = [1.0, 2.0, 3.0]
        observed = [1.0, 2.0, 3.0]
        assert aafe(predicted, observed) == pytest.approx(1.0)

    def test_uniform_2_fold(self):
        predicted = [2.0, 4.0, 8.0]
        observed = [1.0, 2.0, 4.0]
        assert aafe(predicted, observed) == pytest.approx(2.0)

    def test_mixed(self):
        predicted = [2.0, 1.0]  # 2-fold over, 2-fold under
        observed = [1.0, 2.0]
        # geometric mean of (2, 2) = 2
        assert aafe(predicted, observed) == pytest.approx(2.0)

    def test_length_mismatch_raises(self):
        with pytest.raises(ValueError, match="length"):
            aafe([1.0, 2.0], [1.0])

    def test_empty_raises(self):
        with pytest.raises(ValueError, match="empty"):
            aafe([], [])


class TestWithinNFold:
    def test_all_within(self):
        predicted = [1.1, 1.5, 0.8]
        observed = [1.0, 1.0, 1.0]
        assert within_n_fold(predicted, observed, n=2.0) == pytest.approx(1.0)

    def test_half_within(self):
        predicted = [1.1, 3.5]   # 1.1-fold, 3.5-fold
        observed = [1.0, 1.0]
        assert within_n_fold(predicted, observed, n=2.0) == pytest.approx(0.5)

    def test_custom_n(self):
        predicted = [2.5]
        observed = [1.0]
        assert within_n_fold(predicted, observed, n=2.0) == pytest.approx(0.0)
        assert within_n_fold(predicted, observed, n=3.0) == pytest.approx(1.0)
