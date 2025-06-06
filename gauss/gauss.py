from collections.abc import Callable
from decimal import Decimal
from decimal import getcontext as decimalcontext

from .legendre.bonnet import legendre_polynomial
from .legendre.legroot import poly_rootfind, polyderiv, polyeval


def compute_gauss_quadrature(n: int, prec: int) -> tuple[list[Decimal], list[Decimal]]:
    """Computes the knots x[i] and weights w[i] for the gaussian quadrature rule
       int(f) = sum (  w[i] * f(x[i])  )

    Args:
        n (int): the number of quadrature points of the quadrature rule
        prec (int): the precision (number of decimal places) to set to

    Returns:
        tuple[list[Decimal], list[Decimal]]: the lists x and w of knots and
                weights, respectively.
    """
    # add 2 decimal places just for safety
    decimalcontext().prec = prec + 2

    Pn = legendre_polynomial(n)
    # used for weights
    Pnprime = polyderiv(Pn)

    # add 1 decimal place just for safety
    knots = poly_rootfind(Pn, -1, 1, Decimal("0.1") ** (prec + 1))

    weights = [Decimal(2) / ((1 - x**2) * polyeval(Pnprime, x) ** 2) for x in knots]
    return knots, weights


def compute_integral(
    f: Callable[[Decimal | float], Decimal | float], n: int, prec: int
):
    """Computes the n-node quadrature of the function f, on the interval [-1,1] with the precision `prec`.

    Args:
        f (Callable[[Decimal  |  float],Decimal  |  float]): The function to integrate
        n (int): The number of nodes for the Gaussian Quadrature
        prec (int): the decimal precision of the quadrature rule
    """
    x, w = compute_gauss_quadrature(n, prec)

    sol = Decimal(0)
    for wi, xi in zip(w, x, strict=True):
        sol += wi * Decimal(f(xi))

    return sol
