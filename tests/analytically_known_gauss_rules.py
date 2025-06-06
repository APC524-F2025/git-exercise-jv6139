from decimal import Decimal

from gauss import compute_gauss_quadrature


def _assert_quadratures_close(x, w, xtrue, wtrue, prec):
    tol = Decimal(10) ** (-prec)
    npoints = len(xtrue)
    assert len(x) == len(w), (
        f"number of knots ({len(x)}) should equal number of weights ({len(x)})."
    )
    assert len(x) == len(xtrue), f"number of knots ({len(x)}) should be {npoints}"
    for i in range(npoints):
        assert abs(x[i] - xtrue[i]) < tol, (
            f"knot {i} not equal:\ncomputed: {x[i]}\n  actual: {xtrue[i]}"
        )
        assert abs(w[i] - wtrue[i]) < tol, (
            f"weight {i} not equal:\ncomputed: {w[i]}\n  actual: {wtrue[i]}"
        )


# https://en.wikipedia.org/wiki/Gauss%E2%80%93Legendre_quadrature


def test_1point(precision):
    xtrue = [0]
    wtrue = [2]

    npoints = len(xtrue)
    _assert_quadratures_close(
        *compute_gauss_quadrature(npoints, precision), xtrue, wtrue, precision
    )


def test_2point(precision):
    one_over_root_three = 1 / Decimal(3).sqrt()
    xtrue = [-one_over_root_three, one_over_root_three]
    wtrue = [1, 1]

    npoints = len(xtrue)
    _assert_quadratures_close(
        *compute_gauss_quadrature(npoints, precision), xtrue, wtrue, precision
    )


def test_3point(precision):
    root_three_fifths = (Decimal(3) / Decimal(5)).sqrt()
    xtrue = [-root_three_fifths, 0, root_three_fifths]
    five_nineths = Decimal(5) / Decimal(9)
    wtrue = [five_nineths, Decimal(8) / Decimal(9), five_nineths]

    npoints = len(xtrue)
    _assert_quadratures_close(
        *compute_gauss_quadrature(npoints, precision), xtrue, wtrue, precision
    )


def test_4point(precision):
    tmp0 = (Decimal(6) / Decimal(5)).sqrt() * Decimal(2) / Decimal(7)
    tmp1 = (Decimal(3) / Decimal(7) - tmp0).sqrt()
    tmp2 = (Decimal(3) / Decimal(7) + tmp0).sqrt()
    xtrue = [-tmp2, -tmp1, tmp1, tmp2]
    tmp3 = Decimal(30).sqrt() / Decimal(36)
    one_half = 1 / Decimal(2)
    wtrue = [one_half - tmp3, one_half + tmp3, one_half + tmp3, one_half - tmp3]

    npoints = len(xtrue)
    _assert_quadratures_close(
        *compute_gauss_quadrature(npoints, precision), xtrue, wtrue, precision
    )


def test_5point(precision):
    tmp0 = (Decimal(10) / Decimal(7)).sqrt() * Decimal(2)
    tmp1 = (5 - tmp0).sqrt() / Decimal(3)
    tmp2 = (5 + tmp0).sqrt() / Decimal(3)
    xtrue = [-tmp2, -tmp1, 0, tmp1, tmp2]
    tmp3 = Decimal(13) * Decimal(70).sqrt() / Decimal(900)
    tmp4 = Decimal(322) / Decimal(900)
    wtrue = [
        tmp4 - tmp3,
        tmp4 + tmp3,
        Decimal(128) / Decimal(225),
        tmp4 + tmp3,
        tmp4 - tmp3,
    ]

    npoints = len(xtrue)
    _assert_quadratures_close(
        *compute_gauss_quadrature(npoints, precision), xtrue, wtrue, precision
    )
