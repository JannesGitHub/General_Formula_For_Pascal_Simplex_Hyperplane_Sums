"""Exact non-recursive formulas for Pascal-simplex hyperplane sums.

Indices start at zero. For dimension d >= 1 and offset m >= 0,
f(0) = ... = f(m) = 1 and f(n) = f(n-1) + (d-1)*f(n-m-1).
See NON_RECURSIVE_FORMULA.md for proofs and the eigenvalue formula.
Only the Python standard library is required.
"""

from math import comb
from operator import index


def _integer_at_least(value, minimum, name):
    if isinstance(value, bool):
        raise TypeError(f"{name} must be an integer, not bool")
    value = index(value)
    if value < minimum:
        raise ValueError(f"{name} must be at least {minimum}")
    return value


def calculate_term(dimension, m, n):
    """Return f(n) exactly, without calculating any earlier sequence terms.

    f(n) = sum(comb(n - m*j, j) * (dimension - 1)**j,
               j = 0, ..., floor(n/(m+1))).

    The sum has floor(n/(m+1)) + 1 terms, with arbitrary-precision
    integer arithmetic. dimension=1 gives 1; m=0 gives dimension**n.
    """
    dimension = _integer_at_least(dimension, 1, "dimension")
    m = _integer_at_least(m, 0, "m")
    n = _integer_at_least(n, 0, "n")
    if dimension == 1:
        return 1
    if m == 0:
        return dimension**n
    a = dimension - 1
    return sum(comb(n - m*j, j) * a**j for j in range(n // (m + 1) + 1))


def calculateSequenceClosedForm(dimension, m, k):
    """Return the first k terms, matching the notebook's argument order."""
    dimension = _integer_at_least(dimension, 1, "dimension")
    m = _integer_at_least(m, 0, "m")
    k = _integer_at_least(k, 0, "k")
    return [calculate_term(dimension, m, n) for n in range(k)]


def companion_matrix(dimension, m):
    """Return the integer matrix C with f(n) = (C**n)[0, 0].

    For m >= 1 the first row is [1, 0, ..., 0, dimension-1]
    and the subdiagonal is one. For m=0 both coefficients occupy
    the same entry, so C = [[dimension]].
    """
    dimension = _integer_at_least(dimension, 1, "dimension")
    m = _integer_at_least(m, 0, "m")
    size = m + 1
    matrix = [[0] * size for _ in range(size)]
    matrix[0][0] = 1
    matrix[0][-1] += dimension - 1
    for i in range(1, size):
        matrix[i][i - 1] = 1
    return tuple(tuple(row) for row in matrix)


def calculate_term_matrix(dimension, m, n):
    """Evaluate (C**n)[0,0] exactly using iterative binary powering.

    This avoids generating preceding f values. It uses O((m+1)**3 log n)
    integer arithmetic operations for n >= 2, not counting integer bit costs.
    """
    power = companion_matrix(dimension, m)
    n = _integer_at_least(n, 0, "n")
    size = len(power)
    state = [1] + [0] * (size - 1)
    exponent = n
    while exponent:
        if exponent & 1:
            state = [sum(a*b for a, b in zip(row, state)) for row in power]
        exponent //= 2
        if exponent:
            columns = tuple(zip(*power))
            power = tuple(
                tuple(sum(a*b for a, b in zip(row, column)) for column in columns)
                for row in power
            )
    return state[0]
