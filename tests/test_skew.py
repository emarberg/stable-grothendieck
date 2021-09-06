from utils import GP, GQ, SymmetricPolynomial
from partitions import Partition
from vectors import Vector


gp_skew = {}
gq_skew = {}


def substitute(vec, b):
    return Vector({key: value.substitute(0, b) for (key, value) in vec.items()})


def gp_expand(f, n, mu, nu):
    return SymmetricPolynomial.GP_expansion(f(n, mu, nu))


def search(n, k):
    partitions = list(Partition.all(k, strict=True))
    for mu in partitions:
        for nu in partitions:
            if Partition.contains(mu, nu) and len(nu) > 0 and nu[0] != mu[0] and not any(nu[i] > mu[i + 1] for i in range(min(len(nu), len(mu) - 1))):
                e = len([(i, j) for (i, j) in Partition.skew(mu,  nu) if i == j])
                print('\ntarget:', gp_expand(GQ, n, mu, nu), '     ', mu, nu)
                print()
                seen = set()
                s = []
                for a in Partition.successors(mu, strict=True):
                    for b in Partition.successors(nu, strict=True):
                        if Partition.contains(a, b) and Partition.skew_key(a, b, True) not in seen:
                            s += [(str(gp_expand(GP, n, a, b)), a, b)]
                            seen.add(Partition.skew_key(a, b, True))
                m = max([len(x) for x, _, _ in s])
                for x, a, b in s:
                    print('       ', x, ' ' * (m - len(x)), '       ', a, b)
                input('')


def expand_search(n, k):
    if n not in gp_skew:
        gp_skew[n] = {}
    if n not in gq_skew:
        gq_skew[n] = {}

    partitions = list(Partition.all(k, strict=True))
    for mu in partitions:
        for nu in partitions:
            if Partition.contains(mu, nu) and not (nu and nu[0] == mu[0]):
                key = (nu, mu)
                gp_skew[n][key] = substitute(SymmetricPolynomial.GP_expansion(GP(n, mu, nu)), -2)
                gq_skew[n][key] = substitute(SymmetricPolynomial.GP_expansion(GQ(n, mu, nu)), -2)





