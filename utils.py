from symmetric import SymmetricPolynomial
from permutations import Permutation


def schur(num_variables, mu):
    return SymmetricPolynomial.schur(num_variables, mu)


def schur_Q(num_variables, mu): # noqa
    return SymmetricPolynomial.schur_q(num_variables, mu)


def schur_P(num_variables, mu): # noqa
    return SymmetricPolynomial.schur_p(num_variables, mu)


def schur_S(num_variables, mu): # noqa
    return SymmetricPolynomial.schur_s(num_variables, mu)


def grothendieck(num_variables, mu, degree_bound=None): # noqa
    if type(mu) == Permutation:
        return mu.stable_grothendieck(num_variables, degree_bound)
    else:
        return SymmetricPolynomial.stable_grothendieck(num_variables, mu, degree_bound)


def grothendieck_Q(num_variables, mu, degree_bound=None): # noqa
    return SymmetricPolynomial.stable_grothendieck_q(num_variables, mu, degree_bound)


def grothendieck_P(num_variables, mu, degree_bound=None): # noqa
    if type(mu) == Permutation:
        return mu.symplectic_stable_grothendieck(num_variables, degree_bound)
    else:
        return SymmetricPolynomial.stable_grothendieck_p(num_variables, mu, degree_bound)

def grothendieck_S(num_variables, mu, degree_bound=None): # noqa
    if type(mu) == Permutation:
        return mu.signed_involution_stable_grothendieck(num_variables, degree_bound)
    else:
        return SymmetricPolynomial.stable_grothendieck_s(num_variables, mu, degree_bound)


def dual_grothendieck(num_variables, mu): # noqa
    return SymmetricPolynomial.dual_stable_grothendieck(num_variables, mu)


def dual_grothendieck_Q(num_variables, mu): # noqa
    return SymmetricPolynomial.dual_stable_grothendieck_q(num_variables, mu)


def dual_grothendieck_P(num_variables, mu): # noqa
    return SymmetricPolynomial.dual_stable_grothendieck_p(num_variables, mu)


s = schur
P = schur_P
Q = schur_Q
S = schur_S

G = grothendieck
GP = grothendieck_P
GQ = grothendieck_Q
GS = grothendieck_S

g = dual_grothendieck
gp = dual_grothendieck_P
gq = dual_grothendieck_Q
