import pytest

# skip entire module if scipy not installed
pytest.importorskip(
    "scipy", reason="scipy is required for confidence interval calculations"
)
from Tools.ci import jeffreysCI


@pytest.mark.parametrize(
    "x,n,alpha,expected",
    [
        (0, 0, 0.05, (0.0, 0.0, 1.0)),  # zero trials
        (0, 10, 0.05, (0.0, 0.0, pytest.approx(1 - (0.05 / 2) ** (1 / 10)))),
        (10, 10, 0.05, (1.0, pytest.approx((0.05 / 2) ** (1 / 10)), 1.0)),
        (
            5,
            10,
            0.05,
            pytest.approx((0.5,)),  # symmetric case p=0.5
        ),
    ],
)
def test_jeffreysCI_basic(x, n, alpha, expected):
    # Test known cases for Jeffreys CI
    res = jeffreysCI(x, n, alpha)
    # Compare proportion
    assert res[0] == pytest.approx(expected[0], rel=1e-6)
    # Check bounds are within [0,1]
    assert 0.0 <= res[1] <= 1.0
    assert 0.0 <= res[2] <= 1.0


def test_jeffreysCI_cache():
    # Ensure repeated calls hit the cache without error
    v1 = jeffreysCI(3, 15, 0.1)
    v2 = jeffreysCI(3, 15, 0.1)
    assert v1 == v2
