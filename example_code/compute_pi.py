from decimal import Decimal

from gauss import compute_composite_integral

PRECISION = 80

for n in [5, 10, 20, 40]:
    pi_approx = (
        compute_composite_integral(
            lambda x: (1 - Decimal(x) ** 2).sqrt(), 20, n, PRECISION
        )
        * 2
    )
    print(pi_approx)
