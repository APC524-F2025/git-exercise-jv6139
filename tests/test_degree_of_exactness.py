from decimal import Decimal

import pytest

from gauss import compute_gauss_quadrature


@pytest.mark.parametrize("npoints", [3, 5, 8])
def test_polynomial_exactness(npoints, precision):
    tol = Decimal(10) ** (-precision // 2)
    x, w = compute_gauss_quadrature(npoints, precision)

    for k in range(2 * npoints):
        integ = sum(wi * xi**k for wi, xi in zip(w, x, strict=True)) - (
            1 + (-1) ** k
        ) / Decimal(k + 1)
        assert abs(integ) < tol, (
            f"Should be exact to degree {2 * npoints - 1}, but achieves error {integ} for x**{k}"
        )

    k = 2 * npoints
    integ = sum(wi * xi**k for wi, xi in zip(w, x, strict=True))
    assert abs(integ) > tol, (
        f"Should not be exact to degree {2 * npoints}, but achieves error {integ} for x**{k}"
    )
