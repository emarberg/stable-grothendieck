from symmetric import SymmetricPolynomial
from permutations import Permutation


def schur(mu, n):
    return SymmetricPolynomial.schur(mu, n)


def schur_Q(mu, n): # noqa
    return SymmetricPolynomial.schur_q(mu, n)


def schur_P(mu, n): # noqa
    return SymmetricPolynomial.schur_p(mu, n)


def schur_S(mu, n): # noqa
    return SymmetricPolynomial.schur_s(mu, n)


def grothendieck(mu, n, degree_bound=None): # noqa
    if type(mu) == Permutation:
        return mu.stable_grothendieck(n, degree_bound)
    else:
        return SymmetricPolynomial.stable_grothendieck(mu, n, degree_bound)


def grothendieck_Q(mu, n, degree_bound=None): # noqa
    return SymmetricPolynomial.stable_grothendieck_q(mu, n, degree_bound)


def grothendieck_P(mu, n, degree_bound=None): # noqa
    if type(mu) == Permutation:
        return mu.symplectic_stable_grothendieck(n, degree_bound)
    else:
        return SymmetricPolynomial.stable_grothendieck_p(mu, n, degree_bound)

def grothendieck_S(mu, n, degree_bound=None): # noqa
    if type(mu) == Permutation:
        return mu.signed_involution_stable_grothendieck(n, degree_bound)
    else:
        return SymmetricPolynomial.stable_grothendieck_s(mu, n, degree_bound)


def dual_grothendieck(mu, n): # noqa
    return SymmetricPolynomial.dual_stable_grothendieck(mu, n)


def dual_grothendieck_Q(mu, n): # noqa
    return SymmetricPolynomial.dual_stable_grothendieck_q(mu, n)


def dual_grothendieck_P(mu, n): # noqa
    return SymmetricPolynomial.dual_stable_grothendieck_p(mu, n)


s = schur
P = schur_P
Q = schur_Q
S = schur_S

G = grothendieck
GP = grothendieck_P
GQ = grothendieck_Q
GS = grothendieck_S

gp = dual_grothendieck_P
gq = dual_grothendieck_Q


