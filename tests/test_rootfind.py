from decimal import Decimal

from gauss.legendre.legroot import poly_rootfind


def _assert_roots_close(x, xtrue, prec):
    tol = Decimal(10) ** (-prec)
    for i, roots in enumerate(zip(x, xtrue, strict=True)):
        xi, xtruei = roots
        assert abs(xi - xtruei) < tol, (
            f"knot {i} not equal:\ncomputed: {xi}\n  actual: {xtruei}"
        )


def test_easy_polynomials(precision):
    # x^2 - 1/2
    _assert_roots_close(
        poly_rootfind(
            [-Decimal("0.25"), Decimal(0), Decimal(1)],
            -1,
            1,
            Decimal(10) ** (-precision),
        ),
        [-Decimal("0.5"), Decimal("0.5")],
        precision,
    )

    # x(x^2 - 1/2)
    _assert_roots_close(
        poly_rootfind(
            [Decimal(0), -Decimal("0.25"), Decimal(0), Decimal(1)],
            -1,
            1,
            Decimal(10) ** (-precision),
        ),
        [-Decimal("0.5"), 0, Decimal("0.5")],
        precision,
    )
