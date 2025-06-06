from decimal import Decimal
from fractions import Fraction
from typing import TypeVar

# int * coeff works for both Decimal and Fraction
T = TypeVar("T", Decimal, Fraction)


def polyderiv(p: list[T]) -> list[T]:
    """Computes the derivative of the given polynomial."""
    return [k * a for k, a in enumerate(p)][1:]


def polyeval(p: list[Decimal], x: Decimal | float) -> Decimal:
    """Computes p(x) for Decimal x and Decimal coefficients of p."""
    if isinstance(x, float):
        x = Decimal(x)
    if len(p) == 1:
        return p[0]
    result = p[0] + x * p[1]
    xk = x
    for a in p[2:]:
        xk *= x
        result += xk * a
    return result


def poly_rootfind(
    p: list[Decimal], lb: Decimal | float, ub: Decimal | float, tol: Decimal | float
) -> list[Decimal]:
    """Numerically computes the roots of the given polynomial `p` using
    a sequence of bisection methods. Assumes all of the roots are real,
    and lie in the range `(lb,ub)`.

    The true values (at the current decimal accuracy) will be within a
    distance `tol` of the given values.

    Args:
        p (list[Decimal]): the list of coefficients of p in increasing degree.
            Assumes p[-1] != 0.
        lb (Decimal): Lower bound of roots to find
        ub (Decimal): Upper bound of roots to find
        tol (Decimal): the tolerance for the roots

    Returns:
        list[Decimal]: the list of roots of p, in increasing value
    """
    degree = len(p) - 1
    if degree == 2:
        vert = -p[1] / (2 * p[2])
        off = (p[1] ** 2 - 4 * p[2] * p[0]).sqrt() / (2 * p[2])
        return [vert - off, vert + off]

    if degree == 1:
        return [p[0] / p[1]]

    deriv = polyderiv(p)
    # since roots are all real and unique, derivative roots lie between each root.
    root_bounds = [lb, *poly_rootfind(deriv, lb, ub, tol), ub]

    roots = []
    for i in range(degree):
        # bisection method: assume one root falls between a and b
        # therefore (p(a)*p(b) < 0)
        a = root_bounds[i]
        b = root_bounds[i + 1]

        pa = polyeval(p, a)

        while True:
            c = (a + b) / 2
            if c - a < tol:
                break
            pc = polyeval(p, c)
            if pa * pc < 0:
                # root is between a and c
                b = c
            else:
                # root is between c and b
                a = c

        roots.append((a + b) / 2)

    return roots
