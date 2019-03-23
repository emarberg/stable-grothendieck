# stable-grothendieck

Python code for computations involving stable Grothendieck polynomials, their shifted versions, and related symmetric polynomials.

## Setup
1. Install Python 3
1. Install pytest: `pip3 install pytest`

## Tests
1. Run the tests at the command line: `pytest`
1. If there's time to spare, run the slow tests with `pytest --runslow` 

## Usage

Use functions in the `python3` interpreter.

```python3
from utils import *

# number of variables
n = 6

# define a partition
mu = (3, 2, 1)

# compute Schur function in variables x_1, x_2, ..., x_n
schur(mu, n)
s(mu, n)

# compute Schur decomposition
f = schur(mu, n) * schur(mu, n)
SymmetricPolynomial.schur_expansion(f, n)

# compute Schur P-function in variables x_1, x_2, ..., x_n
schur_P(mu, n)
P(mu, n)

# compute Schur Q-function in variables x_1, x_2, ..., x_n
schur_Q(mu, n)
Q(mu, n)

# compute Schur S-function in variables x_1, x_2, ..., x_n
schur_S(mu, n)
S(mu, n)

# compute stable Grothendieck polynomial in variables x_1, x_2, ..., x_n
grothendieck(mu, n)
G(mu, n)

# compute K-theoretic Schur P-function in variables x_1, x_2, ..., x_n
grothendieck_P(mu, n)
GP(mu, n)

# compute K-theoretic Schur Q-function in variables x_1, x_2, ..., x_n
grothendieck_Q(mu, n)
GQ(mu, n)

# compute K-theoretic Schur S-function in variables x_1, x_2, ..., x_n
grothendieck_S(mu, n)
GS(mu, n)

# compute dual stable Grothendieck polynomial in variables x_1, x_2, ..., x_n
g(mu, n)
dual_grothendieck(mu, n)

# compute shifted dual stable Grothendieck polynomials in variables x_1, x_2, ..., x_n
gp(mu, n)
dual_grothendieck_P(mu, n)

gq(mu, n)
dual_grothendieck_Q(mu, n)

# add, subtract, and multiply polynomials
f = P(mu, n)
g = Q(mu, n)
8 * f == g

alpha, beta = (3,), (4,)
f = GP(alpha, n) - GP(beta, n)
g = GQ(alpha, n)
f == g

# generate tableaux and reverse plane partitions
from tableaux import Tableau
mu = (4, 2, 1)

a = Tableau.semistandard(mu, 3)
b = Tableau.semistandard_setvalued(mu, 3)

c = Tableau.semistandard_shifted(mu, 3, diagonal_primes=True)
d = Tableau.semistandard_shifted_setvalued(mu, 4, diagonal_primes=False)

e = Tableau.semistandard_rpp(mu, 3)
f = Tableau.semistandard_shifted_rpp(mu, 3, diagonal_nonprimes=False)
```