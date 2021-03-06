from symmetric import SymmetricPolynomial
from permutations import Permutation
from polynomials import beta # noqa


def schur(num_variables, mu, nu=()):
    return SymmetricPolynomial.schur(num_variables, mu, nu)


def schur_Q(num_variables, mu, nu=()): # noqa
    return SymmetricPolynomial.schur_q(num_variables, mu, nu)


def schur_P(num_variables, mu, nu=()): # noqa
    return SymmetricPolynomial.schur_p(num_variables, mu, nu)


def schur_S(num_variables, mu, nu=()): # noqa
    return SymmetricPolynomial.schur_s(num_variables, mu, nu)


def grothendieck(num_variables, mu, nu=(), degree_bound=None): # noqa
    if type(mu) == Permutation:
        assert nu == ()
        return mu.stable_grothendieck(num_variables, degree_bound=degree_bound)
    else:
        return SymmetricPolynomial.stable_grothendieck(num_variables, mu, nu, degree_bound=degree_bound)


def grothendieck_Q(num_variables, mu, nu=(), degree_bound=None): # noqa
    return SymmetricPolynomial.stable_grothendieck_q(num_variables, mu, nu, degree_bound=degree_bound)


def grothendieck_P(num_variables, mu, nu=(), degree_bound=None): # noqa
    if type(mu) == Permutation:
        assert nu == ()
        return mu.symplectic_stable_grothendieck(num_variables, degree_bound)
    else:
        return SymmetricPolynomial.stable_grothendieck_p(num_variables, mu, nu, degree_bound)

def grothendieck_S(num_variables, mu, nu=(), degree_bound=None): # noqa
    if type(mu) == Permutation:
        assert nu == ()
        return mu.signed_involution_stable_grothendieck(num_variables, degree_bound)
    else:
        return SymmetricPolynomial.stable_grothendieck_s(num_variables, mu, nu, degree_bound)


def grothendieck_doublebar(num_variables, mu, nu=(), degree_bound=None): # noqa
    return SymmetricPolynomial.stable_grothendieck_doublebar(num_variables, mu, nu, degree_bound=degree_bound)


def grothendieck_Q_doublebar(num_variables, mu, nu=(), degree_bound=None): # noqa
    return SymmetricPolynomial.stable_grothendieck_q_doublebar(num_variables, mu, nu, degree_bound=degree_bound)


def grothendieck_P_doublebar(num_variables, mu, nu=(), degree_bound=None): # noqa
    return SymmetricPolynomial.stable_grothendieck_p_doublebar(num_variables, mu, nu, degree_bound=degree_bound)


def dual_grothendieck(num_variables, mu, nu=()): # noqa
    return SymmetricPolynomial.dual_stable_grothendieck(num_variables, mu, nu)


def dual_grothendieck_Q(num_variables, mu, nu=()): # noqa
    return SymmetricPolynomial.dual_stable_grothendieck_q(num_variables, mu, nu)


def dual_grothendieck_P(num_variables, mu, nu=()): # noqa
    return SymmetricPolynomial.dual_stable_grothendieck_p(num_variables, mu, nu)


s = schur
P = schur_P
Q = schur_Q
S = schur_S

G = grothendieck
GP = grothendieck_P
GQ = grothendieck_Q
GS = grothendieck_S

G_doublebar = grothendieck_doublebar
GP_doublebar = grothendieck_P_doublebar
GQ_doublebar = grothendieck_Q_doublebar

g = dual_grothendieck
gp = dual_grothendieck_P
gq = dual_grothendieck_Q

schur_expansion = SymmetricPolynomial.schur_expansion

G_expansion = SymmetricPolynomial.grothendieck_expansion
g_expansion = SymmetricPolynomial.dual_grothendieck_expansion

gp_expansion = SymmetricPolynomial.gp_expansion
gq_expansion = SymmetricPolynomial.gq_expansion

GP_expansion = SymmetricPolynomial.GP_expansion
GQ_expansion = SymmetricPolynomial.GQ_expansion

P_expansion = SymmetricPolynomial.P_expansion
Q_expansion = SymmetricPolynomial.Q_expansion
