import pytest
from fractions import Fraction
from decimal import Decimal, getcontext
from itertools import zip_longest

from gauss.legendre.bonnet import legendre_polynomial

@pytest.fixture(params=[5, 10, 100, 1000])
def precision(request):
    getcontext().prec = request.param + 1  # guard digit
    return request.param

def _to_decimal(x):
    # Convert ints/floats/Decimals/Fractions to Decimal
    if isinstance(x, Fraction):
        # Decimal can't take Fraction directly; convert exactly
        return Decimal(x.numerator) / Decimal(x.denominator)
    return Decimal(x)

def _assert_poly_close(p, ptrue, prec):
    """Compare coefficient lists to tolerance 10^(-prec)."""
    tol = Decimal(10) ** (-prec)
    for k, (a, atrue) in enumerate(zip_longest(p, ptrue, fillvalue=0)):
        a = _to_decimal(a)
        atrue = _to_decimal(atrue)
        assert abs(a - atrue) <= tol, (
            f"degree {k} coefficient not within tol:\n"
            f"computed: {a} x**{k}\n  actual: {atrue} x**{k}\n  tol: {tol}"
        )

def test_leg0(precision):
    _assert_poly_close(legendre_polynomial(0), [1], precision)

def test_leg1(precision):
    _assert_poly_close(legendre_polynomial(1), [0, 1], precision)

def test_leg2(precision):
    _assert_poly_close(
        legendre_polynomial(2),
        [-Decimal("0.5"), 0, Decimal("1.5")],
        precision,
    )

def test_leg3(precision):
    _assert_poly_close(
        legendre_polynomial(3),
        [0, -Decimal("1.5"), 0, Decimal("2.5")],
        precision,
    )

def test_leg4(precision):
    eighth = Decimal(1) / Decimal(8)
    _assert_poly_close(
        legendre_polynomial(4),
        [Decimal(3) * eighth, 0, -Decimal(30) * eighth, 0, Decimal(35) * eighth],
        precision,
    )
