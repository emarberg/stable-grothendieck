from utils import GP, GQ, SymmetricPolynomial
from partitions import Partition
from vectors import Vector


gp_skew = {}
gq_skew = {}


def substitute(vec, b=None):
    if b is None:
        return vec
    return Vector({key: value.substitute(0, b).constant_term() for (key, value) in vec.items()})


def gp_expand(f, n, mu, nu):
    return substitute(SymmetricPolynomial.GP_expansion(f(n, mu, nu)), -2)


def search(n, k):
    partitions = list(Partition.all(k, strict=True))
    for mu in partitions:
        for nu in partitions:
            if Partition.contains(mu, nu) and len(nu) > 0 and nu[0] != mu[0] and not any(nu[i] > mu[i + 1] for i in range(min(len(nu), len(mu) - 1))):
                target = gp_expand(GQ, n, mu, nu)
                print('\ntarget:', target, '     ', mu, nu)
                print()

                seen = set()
                s = []
                for aa in Partition.successors(mu, strict=True):
                    for a in Partition.successors(aa, strict=True):
                        for bb in Partition.successors(nu, strict=True):
                            for b in Partition.successors(bb, strict=True):
                                if Partition.contains(a, b) and Partition.skew_key(a, b, True) not in seen:
                                    s += [(gp_expand(GP, n, a, b), a, b)]
                                    seen.add(Partition.skew_key(a, b, True))

                display = [(str(x), a, b) for (x, a, b) in s]
                width = 48
                display = [(x if len(x) < width else x[:width] + '...', a, b) for (x, a, b) in display]
                m = max([len(x) for x, _, _ in display])
                for i, (x, a, b) in enumerate(display):
                    print('       ', x, ' ' * (m - len(x)), '       ', a, b, ':', i + 1)

                assert target.is_expressable(*[x[0] for x in s])


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
                gp_skew[n][key] = substitute(SymmetricPolynomial.GP_expansion(GP(n, mu, nu)), None)
                gq_skew[n][key] = substitute(SymmetricPolynomial.GQ_expansion(GQ(n, mu, nu)), None)





