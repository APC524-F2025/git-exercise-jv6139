from gauss.legendre.bonnet import legendre_polynomial, next_bonnet

MAX_DEGREE = 50


def polyprod(p, q):
    result = [0] * (len(p) + len(q))
    for i, pi in enumerate(p):
        for j, qj in enumerate(q):
            result[i + j] += pi * qj
    return result


def polydot(p, q):
    pq = polyprod(p, q)

    # integrate pq from -1 to 1: power rule
    result = 0
    for i, ai in enumerate(pq):
        if i % 2 == 1:
            # odd function: integral 0
            continue
        result += 2 * ai / (i + 1)
    return result


def test_legendre_orthogonality():
    polys = [legendre_polynomial(0), legendre_polynomial(1)]
    assert polydot(*polys) == 0
    for n in range(2, MAX_DEGREE):
        Pn = next_bonnet(n - 1, polys[-1], polys[-2])
        for m, Pm in enumerate(polys):
            dot = polydot(Pn, Pm)
            assert dot == 0, (
                f"Dot product between P[{m}] and P[{n}] should be zero."
                f" Got {dot} instead."
            )
        polys.append(Pn)
