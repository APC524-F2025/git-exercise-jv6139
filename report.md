# CSXXX Project Report: Gaussian Integration

Carla Character,
Frank Fictional,
Iris Imaginary

<sup><sub>(Spring XXXX, Professor Pat Pretend)</sub></sup>

Gaussian quadrature is a technique in numerical integration that computes an integral with the highest degree of exactness for a given number of sample points (typically called "knots").
It is known [TODO cite] that an $n$-point quadrature rule

$$\int_{-1}^1 f(x)~dx\approx \sum_{k=1}^n w_k f(x_k)$$

cannot be exact for all polynomials up to degree $2n$, by examining the polynomial $f(x)=\prod_{k=1}^n (x-x_k)^2$, which has a nonzero integral, but a quadrature result of 0. However, by the use of orthogonal polynomials, Gaussian quadrature is exact for polynomials up to degree $2n-1$, which is therefore optimal. In this project, we create a code that can produce these $w_k$, $x_k$ up to an arbitrary degree.

## Legendre Polynomials

(Iris)

## Finding the Roots of $P_n$

Since the roots of $P_n$ are all real and unique, it is easy to numerically compute them. We employ Newton's method to retrieve each root in order.

### Exterior Newton's Method

Suppose we have a degree $n > 0$ polynomial $p$, with distinct real roots $\xi_1 < \dots < \xi_n$. We notice that outside of $(\xi_1,\xi_n)$, $p$ is convex or concave.

- This is because $p'$ and $p''$ have all real roots inside $(\xi_1,\xi_n)$ (since they must be between the roots of $p$ and $p'$, respectively), so $p'$ and $p''$ cannot change sign. Without loss of generality, assume $\lim_{x\to -\infty} p(x) = -\infty$. Then, if $x < \xi_1$, we know that $p'(x) > 0$. If $p''(x) > 0$, then $\lim_{x\to -\infty} p''(x) > 0$ by nature of polynomials (it cannot converge to zero). Hence, at some point $p'$ must change sign to the left of $x$, which contradicts the assertion that $p'$ has all of its roots to the right of $\xi_1$.

This means that if we start our Newton iteration at $x < \xi_1$, then each step will stay below $\xi_1$, so we will always retrieve the smallest root.

### Deflation

Once we found a root $\xi$, we can repeat the process for $q(x) = \frac{p(x)}{(x-\xi)}$. To avoid difficulties of stacking numerical errors on $\xi$, we notice by the quotient rule:
$$q'(x)= \frac{(x-\xi)p'(x) - p(x)}{(x-\xi)^2}$$

$$\frac{q(x)}{q'(x)} = \frac{q(x)(x-\xi)^2}{(x-\xi)p'(x) - p(x)} = \frac{p(x)}{p'(x) - \frac{p(x)}{x-\xi}}$$

For induction, if $p_k(x) = p(x) / \prod_{j=1}^k (x - \xi_j)$ and
$$\frac{p_k(x)}{p_k'(x)} = \frac{p(x)}{p'(x) - \sum_{j=1}^k \frac{p(x)}{x-\xi_j}}$$
then defining $p_{k+1}(x) = p_{k}(x) / (x - \xi_{k+1})$ yields
$$\frac{p_{k+1}(x)}{p_{k+1}'(x)} = \frac{p_k(x)}{p_k'(x) - \frac{p_k(x)}{x-\xi_{k+1}}} = \frac{1}{\frac{p'(x) - \sum_{j=1}^k \frac{p(x)}{x-\xi_j}}{p(x)} - \frac{1}{x-\xi_{k+1}}}= \frac{p(x)}{p'(x) - \sum_{j=1}^{k+1} \frac{p(x)}{x-\xi_j}}$$

Thus, our algorithm is as follows:

- Start at $k=0$, setting $p_0(x) = p(x)$, and choosing some $x < \xi_1$.
- Assume we know $\xi_1,\dots,\xi_k < x$. (At $k=0$, we know none of the roots.)
- Using the formulas we found above, compute the Newton iteration (from $x$) to find $\xi_{k+1}$, which is the smallest root of $p_k$.
- Increment $k$. Set $x = \xi_{k+1}$, and go to step 2.

When $x = \xi_k$, the formula for $p_k(x)/p_k'(x)$ has a zero in the denominator. Hence, on the first step of each Newton iteration, when we start at a root,
we use a different formula derived from L'Hospital's rule:

$$\lim_{x\to \xi_k} \frac{p(x)}{p'(x)- \sum_{j=1}^k \frac{p(x)}{x-\xi_j}} = \lim_{x\to \xi_k} \frac{p(x)(x-\xi_k)}{p'(x)(x - \xi_k) - p(x) - \sum_{j=1}^{k-1} \frac{p(x)(x - \xi_k)}{x-\xi_j}}=\lim_{x\to \xi_k} \frac{p(x) + p'(x)(x-\xi_k)}{p'(x) + p''(x)(x - \xi_k) - p'(x) - \sum_{j=1}^{k-1} \frac{p(x) + p'(x)(x-\xi_k)}{x-\xi_j}} = \frac{2p'(\xi_k)}{p''(\xi_k) - \sum_{j=1}^{k-1} \frac{2p'(\xi_k)}{\xi_k-\xi_j}}$$

## Gaussian Quadrature

(Frank)
