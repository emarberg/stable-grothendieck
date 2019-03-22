# stable-grothendieck

Python code for computations involving stable Grothendieck polynomials, their shifted versions, and related symmetric polynomials.

## Setup
1. Install Python 3
1. Install pytest: `pip3 install pytest`

## Tests
Run the tests at the command line: `pytest`

## Usage

1. `python3`
```
from utils import *

# number of variables
n = 6

# define a partition
mu = (3, 2, 1)

# compute Schur function
schur(mu, n)
s(mu, n)

# compute Schur decomposition
f = schur(mu, n) * schur(mu, n)
SymmetricPolynomial.schur_expansion(f, n)

# compute Schur P-function
schur_P(mu, n)
P(mu, n)

# compute Schur Q-function
schur_Q(mu, n)
Q(mu, n)

# compute Schur S-function
schur_S(mu, n)
S(mu, n)

# compute stable Grothendieck polynomials
grothendieck(mu, n)
G(mu, n)

# compute K-theoretic Schur P-function
grothendieck_P(mu, n)
GP(mu, n)

# compute K-theoretic Schur Q-function
grothendieck_Q(mu, n)
GQ(mu, n)

# compute K-theoretic Schur S-function
grothendieck_S(mu, n)
GS(mu, n)

# compute dual stable Grothendieck polynomial
g(mu, n)
dual_grothendieck(mu, n)

# compute shifted dual stable Grothendieck polynomials
gp(mu, n)
dual_grothendieck_P(mu, n)

gq(mu, n)
dual_grothendieck_Q(mu, n)

```