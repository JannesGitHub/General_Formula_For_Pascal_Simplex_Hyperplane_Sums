"""Independent checks against geometry, the original notebook, and matrices.

Run: python -m unittest -v
NumPy, when installed, enables additional floating-point spectral checks.
"""

import ast
import itertools
import json
import math
from pathlib import Path
import unittest

from pascal_closed_form import (
    calculate_term,
    calculate_term_matrix,
    calculateSequenceClosedForm,
    companion_matrix,
)

try:
    import numpy as np
except ImportError:
    np = None


NOTEBOOK = Path(__file__).with_name(
    "Generalizing_Formula_for_Pascal's_Simplex_Hyperplane_Sums.ipynb"
)


def load_notebook_functions():
    """Load only inspected function definitions; do not run notebook examples."""
    notebook = json.loads(NOTEBOOK.read_text(encoding="utf-8"))
    wanted = {
        "createPascalSimplex", "calculate_sum_nested_list",
        "calculateHypersumsInPascalSimplex", "calculateSequence",
    }
    namespace = {"factorial": math.factorial}
    for cell in notebook["cells"]:
        if cell["cell_type"] != "code":
            continue
        parsed = ast.parse("".join(cell["source"]))
        for node in parsed.body:
            if isinstance(node, ast.FunctionDef) and node.name in wanted:
                module = ast.Module(body=[node], type_ignores=[])
                exec(compile(module, str(NOTEBOOK), "exec"), namespace)
    if not wanted <= namespace.keys():
        raise AssertionError("An original notebook function is missing")
    return namespace


def direct_multinomial_hyperplane(dimension, m, n):
    """Sum individual simplex entries on k1+(m+1)*(k2+...+kd)=n."""
    total = 0
    for coordinates in itertools.product(range(n + 1), repeat=dimension):
        if coordinates[0] + (m + 1) * sum(coordinates[1:]) != n:
            continue
        numerator = math.factorial(sum(coordinates))
        denominator = math.prod(math.factorial(k) for k in coordinates)
        total += numerator // denominator
    return total


class ExactFormulaTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.original = load_notebook_functions()

    def test_matches_notebook_recurrence(self):
        for d in range(1, 13):
            for m in range(13):
                with self.subTest(d=d, m=m):
                    self.assertEqual(
                        calculateSequenceClosedForm(d, m, 121),
                        self.original["calculateSequence"](d, m, 121),
                    )

    def test_matches_original_nested_simplex_geometry(self):
        for d in range(1, 7):
            simplex = self.original["createPascalSimplex"](12, d)
            for m in range(1, 9):
                with self.subTest(d=d, m=m):
                    self.assertEqual(
                        calculateSequenceClosedForm(d, m, 12),
                        self.original["calculateHypersumsInPascalSimplex"](simplex, m),
                    )

    def test_matches_independent_multinomial_enumeration(self):
        for d in range(1, 5):
            for m in range(5):
                for n in range(11):
                    with self.subTest(d=d, m=m, n=n):
                        self.assertEqual(
                            calculate_term(d, m, n),
                            direct_multinomial_hyperplane(d, m, n),
                        )

    def test_exact_matrix_powers(self):
        for d in range(1, 8):
            for m in range(8):
                for n in [*range(15), 31, 64, 100]:
                    with self.subTest(d=d, m=m, n=n):
                        self.assertEqual(calculate_term(d, m, n), calculate_term_matrix(d, m, n))

    def test_known_sequences_and_indexing(self):
        self.assertEqual(calculateSequenceClosedForm(2, 1, 10), [1, 1, 2, 3, 5, 8, 13, 21, 34, 55])
        self.assertEqual(calculateSequenceClosedForm(2, 2, 10), [1, 1, 1, 2, 3, 4, 6, 9, 13, 19])
        self.assertEqual(calculateSequenceClosedForm(2, 4, 15), [1, 1, 1, 1, 1, 2, 3, 4, 5, 6, 8, 11, 15, 20, 26])
        for n in range(101):
            self.assertEqual(calculate_term(3, 1, n), (2**(n + 1) + (-1)**n) // 3)
        for d in range(1, 8):
            for m in range(8):
                self.assertEqual(calculate_term(d, m, m), 1)
                self.assertEqual(calculate_term(d, m, m + 1), d)

    def test_degenerate_cases_and_large_integer(self):
        self.assertEqual(calculateSequenceClosedForm(2, 1, 0), [])
        self.assertEqual(calculate_term(1, 50, 10000), 1)
        self.assertEqual(calculate_term(7, 0, 100), 7**100)
        self.assertEqual(companion_matrix(7, 0), ((7,),))
        self.assertEqual(calculate_term(2, 1, 1000), calculate_term_matrix(2, 1, 1000))
        self.assertGreater(calculate_term(2, 1, 1000), 10**200)

    def test_invalid_inputs(self):
        for function in (calculate_term, calculate_term_matrix, calculateSequenceClosedForm):
            for args in ((0, 1, 2), (2, -1, 2), (2, 1, -1)):
                with self.subTest(function=function.__name__, args=args):
                    with self.assertRaises(ValueError):
                        function(*args)
            for args in ((2.5, 1, 2), (2, 1.5, 2), (2, 1, 2.5), (True, 1, 2)):
                with self.subTest(function=function.__name__, args=args):
                    with self.assertRaises(TypeError):
                        function(*args)


@unittest.skipIf(np is None, "NumPy not installed; exact tests still run")
class SpectralChecks(unittest.TestCase):
    def test_binet_formula(self):
        for d in range(2, 8):
            for m in range(1, 9):
                coefficients = [1, -1] + [0] * (m - 1) + [-(d - 1)]
                roots = np.roots(coefficients)
                for n in (0, 1, m, m + 1, 15, 30):
                    with self.subTest(d=d, m=m, n=n):
                        spectral = sum(root**(n + 1) / ((m + 1)*root - m) for root in roots)
                        exact = calculate_term(d, m, n)
                        self.assertLess(abs(spectral - exact), 1e-9 * max(1, exact))

    def test_analytic_singular_values(self):
        for d in range(1, 9):
            for m in range(1, 11):
                with self.subTest(d=d, m=m):
                    a = d - 1
                    delta = math.sqrt(a**4 + 4)
                    expected = sorted([
                        math.sqrt((a*a + 2 + delta) / 2),
                        math.sqrt((a*a + 2 - delta) / 2),
                        *([1.0] * (m - 1)),
                    ], reverse=True)
                    actual = np.linalg.svd(np.array(companion_matrix(d, m), dtype=float), compute_uv=False)
                    np.testing.assert_allclose(actual, expected, atol=1e-12, rtol=1e-12)

    def test_singular_values_cannot_be_powered_like_eigenvalues(self):
        matrix = np.array(companion_matrix(2, 2), dtype=float)
        u, singular, vh = np.linalg.svd(matrix)
        correct = matrix @ matrix
        incorrect = u @ np.diag(singular**2) @ vh
        self.assertGreater(np.linalg.norm(correct - incorrect), 0.1)


if __name__ == "__main__":
    unittest.main()
