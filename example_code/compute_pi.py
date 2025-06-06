from gauss import compute_gauss_quadrature

x, w = compute_gauss_quadrature(20, 100)


# approximate pi by integrating the semicircle sqrt(1 - x^2)
half_pi_approx = 0
for wi, xi in zip(w, x, strict=True):
    half_pi_approx += wi * (1 - xi**2).sqrt()

print(half_pi_approx * 2)