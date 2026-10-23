# Non-recursive formulas for Pascal-simplex hyperplane sums

This note proves a finite binomial formula and an explicit eigenvalue formula
for the hyperplane sums computed by this repository. It also explains why
singular values cannot generally replace eigenvalues in the matrix argument.

The starting point is the implementation in
`Generalizing_Formula_for_Pascal's_Simplex_Hyperplane_Sums.ipynb` at commit
`316b8d1d98192fe4bd43853fe5ab7a633e09fc97`.

## 1. Precise statement and indexing

Let $d\geq1$ be an integer dimension and $m\geq1$ the integer offset used by
the original hyperplane-sum function. Use **zero-based indexing**, exactly as
the Python lists do:

$$
f_0=f_1=\cdots=f_m=1,\qquad
f_n=f_{n-1}+(d-1)f_{n-m-1}\quad(n\geq m+1).
$$

In particular, $f_{m+1}=d$, not $1$ unless $d=1$. With one-based indexing
$g_N=f_{N-1}$, the initial conditions instead read $g_1=\cdots=g_{m+1}=1$.
The notebook's original docstring had the two recurrence coefficients swapped;
its actual implementation already used the recurrence above.

**The finite, exact, non-recursive formula is**

$$
\boxed{
f_n(d,m)=\sum_{j=0}^{\lfloor n/(m+1)\rfloor}
\binom{n-mj}{j}(d-1)^j.
}
$$

It needs no earlier sequence values, no matrix diagonalization, and no
floating-point arithmetic. The term with $j=0$ is $1$, including when $d=1$.
In one-based notation, replace every $n$ in the formula by $N-1$.

For $d\geq2$, let $\lambda_1,\ldots,\lambda_{m+1}$ be **all** complex roots of

$$
p(t)=t^{m+1}-t^m-(d-1).
$$

They are distinct, and the corresponding **Binet formula** is

$$
\boxed{
f_n(d,m)=\sum_{r=1}^{m+1}\frac{\lambda_r^{n+m}}{p'(\lambda_r)}
=\sum_{r=1}^{m+1}\frac{\lambda_r^{n+1}}{(m+1)\lambda_r-m}.
}
$$

For fixed $d,m$, the roots and coefficients are constants independent of $n$.
The number of summands is exactly $m+1$, regardless of $n$. Conjugate roots
give conjugate summands, so the result is real; the first formula proves it
is an integer. For $d=1$, use $f_n=1$ directly instead of dividing by
derivatives at possibly repeated zero roots.

## 2. Proof directly from the simplex

An entry in a layer of the $d$-dimensional Pascal simplex is the multinomial
coefficient

$$
\binom{k_1+\cdots+k_d}{k_1,\ldots,k_d}
=\frac{(k_1+\cdots+k_d)!}{k_1!\cdots k_d!},\qquad k_i\geq0.
$$

The notebook stores the first inner index in reverse order: position $j$ in
layer $\ell$ fixes $k_1=\ell-j$ and leaves $k_2+\cdots+k_d=j$.
Summing all entries of that inner list gives, by the multinomial theorem,

$$
\begin{aligned}
S_{\ell,j}
&=\sum_{k_2+\cdots+k_d=j}
\frac{\ell!}{(\ell-j)!k_2!\cdots k_d!}\\
&=\binom{\ell}{j}
\sum_{k_2+\cdots+k_d=j}\frac{j!}{k_2!\cdots k_d!}\\
&=\binom{\ell}{j}(d-1)^j.
\end{aligned}
$$

The hyperplane function selects `simplex[n - m*j][j]` and sums its contents.
Set $\ell=n-mj$. The validity condition $j\leq\ell$ is equivalent to
$j\leq\lfloor n/(m+1)\rfloor$, giving the boxed finite sum.
For $d=2$ the inner list is already a scalar, and the same identity holds.
The $d=1$ function returns $1$ directly.

This also identifies precisely the family of hyperplanes handled by the
code, in multinomial coordinates:

$$
\boxed{k_1+(m+1)(k_2+\cdots+k_d)=n.}
$$

Thus the formula is proved for the actual geometric sums, rather than merely
being a fit to the recurrence. These are a particular family of hyperplanes;
arbitrary normal vectors require the generalization in Section 8.

## 3. Generating function and a second proof

Write $a=d-1$ and $F(z)=\sum_{n\geq0}f_nz^n$. Extend the sequence by $f_n=0$
for negative $n$. The recurrence now holds for all $n\geq1$, with $f_0=1$.
Consequently,

$$
\boxed{F(z)=\frac{1}{1-z-a z^{m+1}}.}
$$

This identity can equally be obtained from the multinomial generating
function $1/(1-x_1-\cdots-x_d)$ by substituting $x_1=z$ and
$x_2=\cdots=x_d=z^{m+1}$.

Expand as a formal power series, so no analytic convergence assumptions are
needed:

$$
F(z)=\sum_{\ell\geq0}(z+a z^{m+1})^\ell
=\sum_{\ell\geq0}\sum_{j=0}^{\ell}
\binom{\ell}{j}a^j z^{\ell+mj}.
$$

The coefficient of $z^n$ has $\ell=n-mj$, recovering the finite sum.
This is also a combinatorial interpretation: $f_n$ counts ordered tilings of
length $n$ by unit tiles and length-$(m+1)$ tiles in $a$ colors.
With $j$ long tiles there are $n-mj$ tiles total, giving
$\binom{n-mj}{j}a^j$ choices.

## 4. Matrix equation and explicit diagonalization

For $m\geq1$, define the $(m+1)\times(m+1)$ companion matrix

$$
C=\begin{pmatrix}
1&0&\cdots&0&a\\
1&0&\cdots&0&0\\
0&1&\cdots&0&0\\
\vdots&&\ddots&&\vdots\\
0&0&\cdots&1&0
\end{pmatrix}.
$$

Equivalently, $C_{1,1}=1$, $C_{1,m+1}=a$, $C_{i,i-1}=1$ for
$2\leq i\leq m+1$, and all other entries are zero. Thus for $m=1$ this
means $C=\left(\begin{smallmatrix}1&a\\1&0\end{smallmatrix}\right)$.

With $v_n=(f_n,f_{n-1},\ldots,f_{n-m})^T$ and negative indices set to zero,

$$
v_0=e_1,\qquad v_{n+1}=Cv_n,\qquad
\boxed{f_n=e_1^T C^n e_1.}
$$

If one prefers to use only the original nonnegative-index initial terms,
then $v_m=(1,\ldots,1)^T$ and $v_n=C^{n-m}v_m$ for $n\geq m$.

The characteristic polynomial is $\det(tI-C)=p(t)=t^{m+1}-t^m-a$.
For a root $\lambda_r$, an eigenvector is

$$
w_r=(\lambda_r^m,\lambda_r^{m-1},\ldots,1)^T.
$$

For $a>0$ all roots are simple: a repeated root must also annihilate

$$
p'(t)=t^{m-1}((m+1)t-m).
$$

The only possible candidates are $0$ (when $m>1$) and $m/(m+1)$, but

$$
p(0)=-a\ne0,\qquad
p\!\left(\frac{m}{m+1}\right)
=-\frac{m^m}{(m+1)^{m+1}}-a<0.
$$

Therefore the Vandermonde matrix $V=(w_1\ \cdots\ w_{m+1})$ is invertible,
$C=V\operatorname{diag}(\lambda_r)V^{-1}$, and
$C^n=V\operatorname{diag}(\lambda_r^n)V^{-1}$.

There is no need to leave the coefficients as an unspecified inverse matrix.
Lagrange interpolation of $t^q$ for $0\leq q\leq m$ gives

$$
\sum_{r=1}^{m+1}\frac{\lambda_r^q}{p'(\lambda_r)}
=\begin{cases}0,&q<m,\\1,&q=m.\end{cases}
$$

To see this, compare the coefficient of $t^m$ on both sides of
$t^q=\sum_r\lambda_r^q p(t)/((t-\lambda_r)p'(\lambda_r))$.
Consequently,

$$
e_1=\sum_r\frac{w_r}{p'(\lambda_r)},\qquad
e_1^TC^ne_1=\sum_r\frac{\lambda_r^{n+m}}{p'(\lambda_r)}.
$$

Canceling $\lambda_r^{m-1}$ in the derivative proves the stated Binet formula.
Equivalently, partial fractions of $F$ give

$$
F(z)=\sum_r\frac{\lambda_r}{(m+1)\lambda_r-m}
\frac{1}{1-\lambda_r z}.
$$

## 5. Why singular values do not give the same power formula

An SVD writes $C=U\Sigma V^T$ with orthogonal $U,V$. In general,

$$
C^2=U\Sigma(V^TU)\Sigma V^T\ne U\Sigma^2V^T.
$$

The middle factor cannot be discarded. More generally, for $n\geq1$,
$C^n=U\Sigma(V^TU\Sigma)^{n-1}V^T$. This retains a matrix power and
depends on the singular vectors as well as the singular values.
Eigenvalue diagonalization instead has $V^{-1}V=I$ between every factor,
which is exactly what makes the scalar powers $\lambda_r^n$ work.

For this particular matrix the singular values can even be found explicitly.
Up to reordering coordinates,

$$
C^TC\sim
\begin{pmatrix}2&a\\a&a^2\end{pmatrix}\oplus I_{m-1}.
$$

Hence

$$
\boxed{\sigma_\pm=
\sqrt{\frac{a^2+2\pm\sqrt{a^4+4}}{2}}},
$$

with $m-1$ additional singular values equal to $1$.
The two nontrivial values are independent of $m$. For example, for $d=2$,
the largest singular value is always the golden ratio, whereas the dominant
eigenvalue is about $1.618034$ when $m=1$ and $1.465571$ when $m=2$.
Thus replacing the eigenvalues by singular values would predict the wrong
growth rate. The dimension and the full SVD still contain information about
$m$; the claim here is specifically that powering the singular values does
not solve the recurrence.

## 6. Special cases and asymptotics

For $m=1$, put

$$
\lambda_\pm=\frac{1\pm\sqrt{4d-3}}2.
$$

Then the general expression reduces to

$$
\boxed{f_n(d,1)=
\frac{\lambda_+^{n+1}-\lambda_-^{n+1}}{\sqrt{4d-3}}.}
$$

Examples, with every row beginning at $n=0$:

| Parameters | Formula or polynomial | Initial terms |
| --- | --- | --- |
| $d=1$, any $m$ | $f_n=1$ | $1,1,1,1,1,\ldots$ |
| $d=2,m=1$ | $f_n=F_{n+1}$ | $1,1,2,3,5,8,\ldots$ |
| $d=3,m=1$ | $f_n=(2^{n+1}+(-1)^n)/3$ | $1,1,3,5,11,21,\ldots$ |
| $d=2,m=2$ | roots of $t^3-t^2-1$ | $1,1,1,2,3,4,6,9,\ldots$ |
| $d=2,m=4$ | roots of $t^5-t^4-1$ | $1,1,1,1,1,2,3,4,5,6,8,\ldots$ |

For $d\geq2$, there is a unique positive root $\rho>1$ of
$\rho^m(\rho-1)=d-1$. It strictly dominates the moduli of the other roots.
Indeed, any root satisfies $|\lambda|^{m+1}\leq|\lambda|^m+a$.
For $|\lambda|\geq1$, the function $t^m(t-1)$ is strictly increasing, giving
$|\lambda|\leq\rho$. Equality in the triangle inequality forces
$\lambda^m$ and then $\lambda^{m+1}$ to be positive real, so $\lambda=\rho$.

Writing $\beta=\max_{\lambda_r\ne\rho}|\lambda_r|<\rho$ gives

$$
f_n=\frac{\rho^{n+1}}{(m+1)\rho-m}+O(\beta^n),\qquad
\frac{f_{n+1}}{f_n}\longrightarrow\rho.
$$

The error need not tend to zero: some nondominant roots can have modulus
at least $1$. Therefore rounding only the dominant term is not a general
exact formula.

## 7. Executable implementation and verification

`pascal_closed_form.py` provides `calculate_term(dimension, m, n)` using
the exact finite sum, `calculateSequenceClosedForm(dimension, m, k)`,
the companion matrix, and exact binary matrix powering.
No external dependency is required for these functions.

```python
from pascal_closed_form import calculate_term, calculateSequenceClosedForm

calculate_term(3, 1, 10)                  # 683
calculateSequenceClosedForm(2, 2, 10)     # [1, 1, 1, 2, 3, 4, 6, 9, 13, 19]
```

The formulas and new implementation also support the extension $m=0$,
where $F(z)=1/(1-dz)$, $f_n=d^n$, and the matrix is the one-entry matrix
$[d]$. The original geometry routine assumes $m\geq1$ and must not be
called with $m=0$; its stopping condition would not terminate for $d\geq2$.
The tests use independent multinomial enumeration to check that extension.

Run the checks from the repository root:

```sh
python -m unittest -v
```

The checks cover:

- 18,876 exact comparisons with the notebook's implemented recurrence:
  $1\leq d\leq12$, $0\leq m\leq12$, $0\leq n\leq120$.
- 576 exact comparisons with the original nested-simplex geometry:
  $1\leq d\leq6$, $1\leq m\leq8$, $0\leq n\leq11$.
- 220 exact comparisons with independently enumerated multinomial entries:
  $1\leq d\leq4$, $0\leq m\leq4$, $0\leq n\leq10$.
- 1,008 exact comparisons with binary matrix powers, known sequences,
  boundary indices, invalid inputs, and a value at $n=1000$.
- When NumPy is available: 288 numerical Binet checks across degrees 2--9,
  80 checks of the derived singular values, and a counterexample to powering
  the singular values as if they were eigenvalues.

Numerical root evaluation is illustrative, not a certified way to return
large exact integers: cancellation and finite-precision root errors can
grow with $n$. Use the integer finite sum or exact matrix powering for exact
answers. The proofs establish the identities for all stated parameters;
the tests additionally check their implementation and agreement with the
repository's indexing and geometry.

## 8. Extension to other hyperplane directions

For positive integer weights $w_1,\ldots,w_d$, define

$$
H_n=\sum_{w_1k_1+\cdots+w_dk_d=n}
\frac{(k_1+\cdots+k_d)!}{k_1!\cdots k_d!}.
$$

The same multinomial substitution gives

$$
\sum_{n\geq0}H_nz^n=\frac{1}{1-\sum_{i=1}^d z^{w_i}}.
$$

Thus arbitrary **positive integer normal vectors** also give rational
generating functions and constant-coefficient recurrences. This positivity
condition keeps every slice finite; it does not cover arbitrary mixed-sign
or zero weights without additional restrictions.

Let $L=\max_i w_i$ and $p_w(t)=t^L-\sum_i t^{L-w_i}$. If its roots are
simple, the same partial-fraction argument gives the explicit formula

$$
H_n=\sum_{p_w(\lambda)=0}\frac{\lambda^{n+L-1}}{p_w'(\lambda)}.
$$

With repeated roots the partial fractions instead yield polynomial factors
in $n$ multiplying $\lambda^n$; simplicity must be checked for these more
general weights. The repository's family is $w_1=1$ and
$w_2=\cdots=w_d=m+1$, for which simplicity was proved above when $d\geq2$.

## Method references

The derivations above apply standard generating-function and linear-algebra
methods to the repository's specific index convention. They are not a claim
that the general method is new.

- Robert Sedgewick and Philippe Flajolet, [Generating Functions, An Introduction
  to the Analysis of Algorithms](https://aofa.cs.princeton.edu/30gf/): rational
  generating functions, linear recurrences, and partial fractions.
- Gilbert Strang, [Singular Value Decomposition, MIT OpenCourseWare](https://ocw.mit.edu/courses/18-06sc-linear-algebra-fall-2011/pages/positive-definite-matrices-and-applications/singular-value-decomposition/):
  the decomposition $U\Sigma V^T$ and the eigenvalues of $C^TC$.
