# General_Formula_For_Pascal_Simplex_Hyperplane_Sums
This repository contains the implementation of a general formula for calculating hyperplane sums in a Pascal Simplex. 

## Exact non-recursive formulas

With zero-based indices, the implemented sequence has $f_0=\cdots=f_m=1$
and $f_n=f_{n-1}+(d-1)f_{n-m-1}$ for $n\geq m+1$. Its exact formula is

$$
f_n(d,m)=\sum_{j=0}^{\lfloor n/(m+1)\rfloor}\binom{n-mj}{j}(d-1)^j.
$$

For $d\geq2$, $m\geq1$, let $\lambda_1,\ldots,\lambda_{m+1}$ be all roots
of $t^{m+1}-t^m-(d-1)=0$. They are distinct and give the Binet formula

$$
f_n(d,m)=\sum_{r=1}^{m+1}\frac{\lambda_r^{n+1}}{(m+1)\lambda_r-m}.
$$

See [the complete derivation](NON_RECURSIVE_FORMULA.md) for the geometric
proof, companion matrix, explicit eigenvalue coefficients, the distinction
between eigenvalues and singular values, and general positive integer
hyperplane directions. For $d=1$, the sequence is identically $1$.

```python
from pascal_closed_form import calculate_term, calculateSequenceClosedForm

calculate_term(3, 1, 10)               # 683, computed without earlier terms
calculateSequenceClosedForm(2, 1, 8)   # [1, 1, 2, 3, 5, 8, 13, 21]
```

The implementation uses only the Python standard library. Run
`python -m unittest -v` for exact checks against the notebook and independent
multinomial sums. NumPy, if available, enables additional numerical checks
of the eigenvalue and singular-value formulas. Run the notebook from the
repository directory so its new examples can import `pascal_closed_form.py`.

## Background

Pascal’s Triangle is one of the most fascinating mathematical structures, appearing in numerous areas of mathematics, from combinatorics to number theory. It is built by arranging binomial coefficients in a triangular pattern, revealing surprising relationships between numbers. One such connection, perhaps lesser known, links Pascal’s Triangle to the Fibonacci sequence.

The Fibonacci sequence is a simple yet powerful sequence of numbers, where each term is the sum of the two preceding ones. It appears in nature, art, and even financial models.
Remarkably, hidden within Pascal’s Triangle, Fibonacci numbers can be discovered along specific diagonal paths, as you can see in the picture below.

![](https://aperiodical.com/wp-content/uploads/2021/12/image-3.png)

As you can see, the sequence of diagonal sums equals the Fibonacci sequence, which is defined as:

$$
F_n =
\begin{cases}
    0, & \text{if } n = 0 \\
    1, & \text{if } n = 1 \\
    F_{n-1} + F_{n-2}, & \text{if } n \geq 2
\end{cases}
$$

In this notebook, we will explore this connection in greater depth and extend it beyond two dimensions. By generalizing Pascal’s Triangle into higher-dimensional structures, known as Pascal simplices, we will uncover a formula for summing elements along specific hyperplanes. This exploration will not only deepen our understanding of Pascal’s Triangle and the Fibonacci sequence but also reveal elegant mathematical patterns that emerge in higher dimensions.

