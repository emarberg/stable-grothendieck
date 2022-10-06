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


def mp_dual_grothendieck(num_variables, mu, nu=()): # noqa
    return SymmetricPolynomial.mp_dual_stable_grothendieck(num_variables, mu, nu)


def mp_dual_grothendieck_Q(num_variables, mu, nu=()): # noqa
    return SymmetricPolynomial.mp_dual_stable_grothendieck_q(num_variables, mu, nu)


def mp_dual_grothendieck_P(num_variables, mu, nu=()): # noqa
    return SymmetricPolynomial.mp_dual_stable_grothendieck_p(num_variables, mu, nu)


def shifted_ribbon(alpha):
    nu = []
    mu = []
    for i in range(1, len(alpha)):
        nu.append(sum(alpha[i:]))
        mu.append(nu[-1] + alpha[i - 1])
    mu.append(alpha[-1])
    return tuple(mu), tuple(nu)


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

j = SymmetricPolynomial._slow_transposed_dual_stable_grothendieck

g = dual_grothendieck
gp = dual_grothendieck_P
gq = dual_grothendieck_Q

mp_g = mp_dual_grothendieck
mp_gp = mp_dual_grothendieck_P
mp_gq = mp_dual_grothendieck_Q

jp = SymmetricPolynomial._slow_transposed_dual_stable_grothendieck_p
jq = SymmetricPolynomial._slow_transposed_dual_stable_grothendieck_q

schur_expansion = SymmetricPolynomial.schur_expansion

G_expansion = SymmetricPolynomial.grothendieck_expansion

g_expansion = SymmetricPolynomial.dual_grothendieck_expansion
j_expansion = SymmetricPolynomial.transposed_dual_grothendieck_expansion
gp_expansion = SymmetricPolynomial.gp_expansion
gq_expansion = SymmetricPolynomial.gq_expansion

mp_g_expansion = SymmetricPolynomial.mp_dual_grothendieck_expansion
mp_gp_expansion = SymmetricPolynomial.mp_gp_expansion
mp_gq_expansion = SymmetricPolynomial.mp_gq_expansion


gp_free_expansion = SymmetricPolynomial.gp_free_expansion

jp_expansion = SymmetricPolynomial.jp_expansion
jq_expansion = SymmetricPolynomial.jq_expansion

GP_expansion = SymmetricPolynomial.GP_expansion
GQ_expansion = SymmetricPolynomial.GQ_expansion

P_expansion = SymmetricPolynomial.P_expansion
Q_expansion = SymmetricPolynomial.Q_expansion
