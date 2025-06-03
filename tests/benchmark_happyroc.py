"""
Microbenchmarks for ROC computation in Haplo.happyroc
"""

import pytest
from Haplo import compute_roc_points


@pytest.mark.parametrize("n", [100, 1000, 5000])
def test_compute_roc_points_benchmark(benchmark, n):
    # Prepare synthetic counts
    tp = list(range(1, n + 1))
    fp = list(range(n, 0, -1))
    fn = [n // 2] * n
    # Benchmark the compute_roc_points function
    benchmark(compute_roc_points, tp, fp, fn)
