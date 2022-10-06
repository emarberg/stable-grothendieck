from partitions import Partition
from utils import *


def partition_iterator(rows):
    n = 0
    while True:
        for mu in Partition.all(n, max_row=rows):
            yield mu, ()
        n += 1


def strict_partition_iterator(rows):
    n = 0
    while True:
        for mu in Partition.all(n, max_row=rows, strict=True):
            yield mu, ()
        n += 1


def skew_iterator(rows):
    n = 0
    while True:
        for mu in Partition.all(n, max_row=rows):
            for nu in Partition.subpartitions(mu):
                yield mu, nu
        n += 1


def strict_skew_iterator(rows):
    n = 0
    while True:
        for mu in Partition.all(n, max_row=rows, strict=True):
            for nu in Partition.subpartitions(mu, strict=True):
                yield mu, nu
        n += 1


# s, P, Q, g, j, G, J, gp, gq, jp, gq, GP, GQ, gp*gp, gq*gq, GP*GP, GQ*GQ, GS?, skew
data = {
    's': (partition_iterator, s, schur_expansion),
    'P': (strict_partition_iterator, P, P_expansion),
    'Q': (strict_partition_iterator, Q, Q_expansion),
    'g': (partition_iterator, g, g_expansion),
    'G': (partition_iterator, G, G_expansion),
    'gp': (strict_partition_iterator, gp, gp_expansion),
    'gq': (strict_partition_iterator, gq, gq_expansion),
    'GP': (strict_partition_iterator, GP, GP_expansion),
    'GQ': (strict_partition_iterator, GQ, GQ_expansion),
    'mp_g': (partition_iterator, mp_g, mp_g_expansion),
    'mp_gp': (strict_partition_iterator, mp_gp, mp_gp_expansion),
    'mp_gq': (strict_partition_iterator, mp_gq, mp_gq_expansion),
    'skew_G': (skew_iterator, G, None),
    'skew_GP': (strict_skew_iterator, GP, None),
    'skew_GQ': (strict_skew_iterator, GQ, None),
    'ss_skew_G': (skew_iterator, G_doublebar, None),
    'ss_skew_GP': (strict_skew_iterator, GP_doublebar, None),
    'ss_skew_GQ': (strict_skew_iterator, GQ_doublebar, None),
}


def decompose(n, iterator, function, decomp):
    mu, nu = next(iterator)
    f = function(n, mu, nu)
    try:
        expansion = decomp(f)
        return expansion.is_nonnegative()
    except:
        return False


def update(results, trials_left):
    print()
    for (x, y) in results:
        if results[x, y] and not any(x != z != y and results[x, z] and results[z, y] for z in data):
            print('(', x, ')', '-->', '(', y, ')')
    print()
    # for (x, y) in results:
    #     if not results[x, y]:
    #         print('(', x, ')', '-/->', '(', y, ')')
    print()
    print('trials left:', trials_left)
    print()


def test_positivity(n, trials=1000):
    iterators = {name: val[0](n) for name, val in data.items()}
    results = {(x, y): True for x in data for y in data if x != y}
    for i in range(trials):
        for x in data:
            it = iterators[x]
            _, fn, _ = data[x]
            for y in data:
                if x == y:
                    continue
                _, _, dec = data[y]
                results[x, y] = results[x, y] and decompose(n, it, fn, dec)
        update(results, trials - i)
                    


