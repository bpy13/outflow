"""Cash-Karp Runge-Kutta step with the paper code's error measure and step-size control.

The paper code (``CRline.py``, function ``RK5``) takes the embedded fourth-order solution as
the result and the difference between the fifth- and fourth-order solutions as the error
estimate (Press et al., Numerical Recipes, §16.2). Steps are accepted when the largest relative
error is at most the tolerance, and the step size is updated as

    h <- 0.9 h (tol / delta)^0.2   after an accepted step,
    h <- 0.9 h (tol / delta)^0.25  after a rejected step.

States are tuples of floats; the derivative returns a tuple of the same length.
"""

from collections.abc import Callable, Sequence
from operator import mul

Derivative = Callable[[float, tuple[float, ...]], Sequence[float]]

# Cash-Karp tableau (nodes A, coefficients B, fifth-order weights C5, fourth-order weights C4).
_A = (0.0, 1.0 / 5.0, 3.0 / 10.0, 3.0 / 5.0, 1.0, 7.0 / 8.0)
_B = (
    (),
    (1.0 / 5.0,),
    (3.0 / 40.0, 9.0 / 40.0),
    (3.0 / 10.0, -9.0 / 10.0, 6.0 / 5.0),
    (-11.0 / 54.0, 5.0 / 2.0, -70.0 / 27.0, 35.0 / 27.0),
    (1631.0 / 55296.0, 175.0 / 512.0, 575.0 / 13824.0, 44275.0 / 110592.0, 253.0 / 4096.0),
)
_C5 = (37.0 / 378.0, 0.0, 250.0 / 621.0, 125.0 / 594.0, 0.0, 512.0 / 1771.0)
_C4 = (2825.0 / 27648.0, 0.0, 18575.0 / 48384.0, 13525.0 / 55296.0, 277.0 / 14336.0, 1.0 / 4.0)
_ERROR_WEIGHTS = tuple(c5 - c4 for c5, c4 in zip(_C5, _C4, strict=True))  # fifth - fourth order

# Growth factor used when the error estimate is exactly zero (the formula would divide by
# zero). The paper code never met this case in its runs.
_GROWTH_WHEN_ERROR_IS_ZERO = 5.0


def cash_karp_step(
    derivative: Derivative,
    h: float,
    r: float,
    y: Sequence[float],
    first_slope: Sequence[float],
) -> tuple[tuple[float, ...], tuple[float, ...]]:
    """One Cash-Karp step of size h from (r, y). Returns (new state, error estimate).

    ``first_slope`` is dy/dr at (r, y). It is passed in rather than computed so that the
    caller can supply a special value at r = 0, where the paper code uses the L'Hopital limit
    dv/dr = q / (3 rho_0) for the first stage only.

    As in the paper code, the new state is the embedded fourth-order solution and the error is
    (fifth order - fourth order).
    """
    # stages[i] holds the stage increments k1, k2, ... of component i. The weighted sums use
    # sum(map(mul, weights, increments)), which runs its loop in C: the solver spends most of
    # its time in this function. The products and their order are those of the plain form
    # sum(c * k for c, k in zip(weights, increments)), so the results are bit-for-bit the same.
    stages = [[h * s] for s in first_slope]
    for stage in range(1, 6):
        state = [yi + sum(map(mul, _B[stage], ki)) for yi, ki in zip(y, stages, strict=True)]
        slope = derivative(r + _A[stage] * h, state)
        for ki, s in zip(stages, slope, strict=True):
            ki.append(h * s)
    new_y = tuple(yi + sum(map(mul, _C4, ki)) for yi, ki in zip(y, stages, strict=True))
    error = tuple(sum(map(mul, _ERROR_WEIGHTS, ki)) for ki in stages)
    return new_y, error


def relative_error(values: Sequence[float], errors: Sequence[float]) -> float:
    """Largest |error / value| over the components, skipping components whose value is zero.

    The paper code computes max(|dy / y|); a component that is identically zero (for example
    P_C when there is no cosmic-ray injection) gives 0/0. Skipping it explicitly makes the
    result independent of where that component sits (docs/decisions.md).
    """
    largest = 0.0
    for value, error in zip(values, errors, strict=True):
        if value != 0.0:
            ratio = abs(error / value)
            if ratio > largest:
                largest = ratio
    return largest


def next_step_size(h: float, delta: float, tolerance: float, accepted: bool) -> float:
    """The paper code's step update: 0.9 h (tol/delta)^0.2 if accepted, ^0.25 if rejected."""
    if delta == 0.0:
        return h * _GROWTH_WHEN_ERROR_IS_ZERO
    exponent = 0.2 if accepted else 0.25
    return 0.9 * h * (tolerance / delta) ** exponent
