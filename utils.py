from symmetric import SymmetricPolynomial
from permutations import Permutation
from polynomials import beta # noqa


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


def shifted_ribbon(alpha):
    nu = []
    mu = []
    for i in range(1, len(alpha)):
        nu.append(sum(alpha[i:]))
        mu.append(nu[-1] + alpha[i - 1])
    mu.append(alpha[-1])
    return tuple(mu), tuple(nu)


s = SymmetricPolynomial.schur
P = SymmetricPolynomial.schur_p
Q = SymmetricPolynomial.schur_q
S = SymmetricPolynomial.schur_s

G = grothendieck
GP = grothendieck_P
GQ = SymmetricPolynomial.stable_grothendieck_q
GS = grothendieck_S

mn_G = SymmetricPolynomial.mn_stable_grothendieck
mn_GP = SymmetricPolynomial.mn_stable_grothendieck_p
mn_GQ = SymmetricPolynomial.mn_stable_grothendieck_q

G_doublebar = SymmetricPolynomial.stable_grothendieck_doublebar
GP_doublebar = SymmetricPolynomial.stable_grothendieck_p_doublebar
GQ_doublebar = SymmetricPolynomial.stable_grothendieck_q_doublebar
GS_doublebar = SymmetricPolynomial.stable_grothendieck_s_doublebar

g = SymmetricPolynomial.dual_stable_grothendieck
gp = SymmetricPolynomial.dual_stable_grothendieck_p
gq = SymmetricPolynomial.dual_stable_grothendieck_q
gs = SymmetricPolynomial.dual_stable_grothendieck_s

mp_g = SymmetricPolynomial.mp_dual_stable_grothendieck
mp_gp = SymmetricPolynomial.mp_dual_stable_grothendieck_p
mp_gq = SymmetricPolynomial.mp_dual_stable_grothendieck_q

j = SymmetricPolynomial.slow_transposed_dual_stable_grothendieck
jp = SymmetricPolynomial.slow_transposed_dual_stable_grothendieck_p
jq = SymmetricPolynomial.slow_transposed_dual_stable_grothendieck_q
js = SymmetricPolynomial.slow_transposed_dual_stable_grothendieck_s

schur_expansion = SymmetricPolynomial.schur_expansion

G_expansion = SymmetricPolynomial.grothendieck_expansion

g_expansion = SymmetricPolynomial.dual_grothendieck_expansion
j_expansion = SymmetricPolynomial.j_expansion
gp_expansion = SymmetricPolynomial.gp_expansion
gq_expansion = SymmetricPolynomial.gq_expansion
gs_expansion = SymmetricPolynomial.gs_expansion

mp_g_expansion = SymmetricPolynomial.mp_dual_grothendieck_expansion
mp_gp_expansion = SymmetricPolynomial.mp_gp_expansion
mp_gq_expansion = SymmetricPolynomial.mp_gq_expansion

mn_G_expansion = SymmetricPolynomial.mn_grothendieck_expansion
mn_GP_expansion = SymmetricPolynomial.mn_GP_expansion
mn_GQ_expansion = SymmetricPolynomial.mn_GQ_expansion

gp_free_expansion = SymmetricPolynomial.gp_free_expansion

jp_expansion = SymmetricPolynomial.jp_expansion
jq_expansion = SymmetricPolynomial.jq_expansion
js_expansion = SymmetricPolynomial.js_expansion

GP_expansion = SymmetricPolynomial.GP_expansion
GQ_expansion = SymmetricPolynomial.GQ_expansion
GS_expansion = SymmetricPolynomial.GS_expansion

P_expansion = SymmetricPolynomial.P_expansion
Q_expansion = SymmetricPolynomial.Q_expansion
S_expansion = SymmetricPolynomial.S_expansion
