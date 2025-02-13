import numpy as np
import pytest

from bioemu_benchmarks.eval.md_emulation.plot import _get_axis_ranges

NUM_DATA = 5
BUFFER = 0.0


@pytest.mark.parametrize("test_dim", [1, 2, 5])
def test_get_axis_ranges(test_dim):
    test_data = [np.random.randn(NUM_DATA, test_dim) for _ in range(test_dim)]

    axis_ranges = _get_axis_ranges(test_data, buffer=BUFFER)

    assert len(axis_ranges) == test_dim

    for k in range(len(axis_ranges)):
        target_min = min([np.min(x[:, k]) for x in test_data])
        target_max = max([np.max(x[:, k]) for x in test_data])
        target_buffer = (target_max - target_min) * BUFFER
        target_min -= target_buffer
        target_max += target_buffer

        np.testing.assert_allclose(target_min, axis_ranges[k][0])
        np.testing.assert_allclose(target_max, axis_ranges[k][1])
