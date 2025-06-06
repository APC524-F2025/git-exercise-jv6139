from decimal import Decimal

from gauss import compute_integral

PRECISION = 40

for n in [5, 10, 20, 40, 80]:
    pi_approx = (
        compute_integral(
            lambda x: (1 - ((Decimal(x) + 1) / 2) ** 2).sqrt(), n, PRECISION
        )
        * 2
    )
    print(pi_approx)
