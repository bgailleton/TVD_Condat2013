import numpy as np
"""Smoke tests for the TVDCondat2013 extension."""

import os
import sys
import numpy as np

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
import TVDCondat2013 as tvd


def test_tvd_preserves_dtype_and_values():
    """Lambda zero should yield an identical array with the same dtype."""
    x32 = np.array([0, 1, 2, 3], dtype=np.float32)
    for func in (tvd.tvd_2013, tvd.tvd_2017):
        y32 = func(x32, np.float32(0.0))
        assert y32.dtype == np.float32
        np.testing.assert_array_equal(y32, x32)

    x64 = np.array([0.0, 1.0, 2.0, 3.0], dtype=np.float64)
    for func in (tvd.tvd_2013, tvd.tvd_2017):
        y64 = func(x64, 0.0)
        assert y64.dtype == np.float64
        np.testing.assert_array_equal(y64, x64)
def _total_variation(arr: np.ndarray) -> float:
    """Compute the total variation of a 1-D array."""
    return np.sum(np.abs(np.diff(arr)))


def test_tvd_reduces_variation():
    """Denoising should not increase total variation."""
    x32 = np.array([0, 2, 1, 3, 0], dtype=np.float32)
    for func in (tvd.tvd_2013, tvd.tvd_2017):
        y32 = func(x32, np.float32(0.5))
        assert _total_variation(y32) <= _total_variation(x32) + 1e-6

    x64 = np.array([0, 2, 1, 3, 0], dtype=np.float64)
    for func in (tvd.tvd_2013, tvd.tvd_2017):
        y64 = func(x64, 0.5)
        assert _total_variation(y64) <= _total_variation(x64) + 1e-12

