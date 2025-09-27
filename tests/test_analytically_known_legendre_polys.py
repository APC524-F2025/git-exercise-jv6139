import pytest
from fractions import Fraction
from decimal import Decimal, getcontext
from itertools import zip_longest

from gauss.legendre.bonnet import legendre_polynomial


# Frank’s fixture + helpers
@pytest.fixture(params=[5, 10, 100, 1000])
def precision(request):
    prec = request.param
    getcontext().prec = prec + 1
    return prec


def poly_fraction_to_decimal(p: list[Fraction]) -> list[Decimal]:
    return [a.numerator / Decimal(a.denominator) for a in p]


def _assert_poly_close(p, ptrue, prec):
    tol = Decimal(10) ** (-prec)
    p = poly_fraction_to_decimal(p)
    for k, (a, atrue) in enumerate(zip_longest(p, ptrue, fillvalue=0)):
        assert abs(a - atrue) <= tol, (
            f"degree {k} coefficient not within tol:\n"
            f"computed: {a} x**{k}\n  actual: {atrue} x**{k}\n  tol: {tol}"
        )


# Iris’s equaltiy helper 
def _assert_poly_equals(p, ptrue):
    for k, coefs in enumerate(zip_longest(p, ptrue, fillvalue=0)):
        a, atrue = coefs
        assert a == atrue, (
            f"degree {k} coefficient not equal:\n"
            f"computed: {a} x**{k}\n  actual: {atrue} x**{k}"
        )
