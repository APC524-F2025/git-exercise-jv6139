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
    p: list[Decimal],
    lb: Decimal | float,
    ub: Decimal | float,
    tol: Decimal | float,
    pprime: list[Decimal] | None = None,
    root_guesses: list[Decimal | None | float] | None = None,
    max_iters: int = 20
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
        pprime (list[Decimal] | None): the polynomial p'. If not given, then it is computed.
        root_guesses (list[Decimal|None] | None) A list of guesses for the roots.
        max_iters (int): specifies the maximum number of Newton iterations.

    Returns:
        list[Decimal]: the list of roots of p, in increasing value
    """
    degree = len(p) - 1
    if pprime is None:
        pprime = polyderiv(p)
    if root_guesses is None:
        root_guesses = [None for _ in range(degree)]

    # initial guess for first root can be given
    if root_guesses[0] is None:
        root_guesses[0] = lb

    if degree == 2:
        vert = -p[1] / (2 * p[2])
        off = (p[1] ** 2 - 4 * p[2] * p[0]).sqrt() / (2 * p[2])
        return [vert - off, vert + off]

    if degree == 1:
        return [p[0] / p[1]]

    roots = []
    # newton's method + deflation. It is known that Newton step is an underestimate
    # from outside the convex hull of roots
    for i in range(degree):
        if root_guesses[i] is None:
            x = roots[-1]
            # using LH rule, but we need second derivative; just use an approximation for now.
            eval_eps = (Decimal(ub) - Decimal(lb)) / Decimal(10) ** 6
            ppx = polyeval(pprime, x)
            pppx = (polyeval(pprime, x + eval_eps) - polyeval(pprime, x - eval_eps)) / (
                2 * eval_eps
            )
            step = -ppx / (pppx / 2 - sum(ppx / (x - r) for r in roots[:-1]))
        else:
            x = Decimal(root_guesses[0])
            step = -polyeval(p, x) / polyeval(pprime, x)
        
        num_iters = 0
        while num_iters < max_iters and abs(step) > tol:
            x += step
            px = polyeval(p,x)
            step = -px/(polyeval(pprime,x) - sum(px/(x - r) for r in roots))
            num_iters += 1

        roots.append(x + step)


    return roots
