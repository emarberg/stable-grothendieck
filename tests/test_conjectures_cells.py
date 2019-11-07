from permutations import Permutation
from insertion import InsertionAlgorithm
from tableaux import Tableau
from partitions import Partition


def rsk(pi):
    return InsertionAlgorithm.hecke(pi.oneline)


def des_m(pi):
    return pi.right_descent_set


def des_n(pi):
    ans = set()
    for i in range(1, pi.rank):
        s = Permutation.s_i(i)
        if len(s * pi * s) >= len(pi):
            ans.add(i)
    return ans


def dual_equivalence(tab, i):
    word = tab.row_reading_word()

    subword = tuple(a for a in word if a in [i -1, i, i + 1])
    if subword in [(i, i + 1, i - 1), (i - 1, i + 1, i)]:
        word = tuple(a + (1 if a == i - 1 else -1 if a == i else 0) for a in word)

    elif subword in [(i, i - 1, i + 1), (i + 1, i - 1, i)]:
        word = tuple(a + (1 if a == i else -1 if a == i + 1 else 0) for a in word)

    else:
        assert subword in [(i - 1, i, i + 1), (i + 1, i, i - 1)]

    return Tableau.from_row_reading_word(word)


def representative_m(mu):
    a = 0
    w = Permutation()
    for b in Partition.transpose(mu):
        assert b % 2 == 0
        for i in range(b // 2):
            w *= Permutation.transposition(a + i + 1, a + b - i)
        a += b
    return w

def representative_n(mu):
    if mu == (2, 2, 2, 2):
        return Permutation(5, 6, 7, 8, 1, 2, 3, 4)
    if mu == (3, 3, 2, 2):
        return Permutation(5, 6, 9, 10, 1, 2, 8, 7, 3, 4)
    if mu == (2, 2, 2, 2, 1, 1):
        return Permutation(2, 1, 7, 8, 9, 10, 3, 4, 5, 6)
    if mu == (4, 4, 2, 2):
        return Permutation(5, 6, 11, 12, 1, 2, 9, 10, 7, 8, 3, 4)
    if mu == (3, 3, 3, 3):
        return Permutation(9, 10, 11, 12, 6, 5, 8, 7, 1, 2, 3, 4)
    if mu == (3, 3, 2, 2, 1, 1):
        return Permutation(2, 1, 7, 8, 11, 12, 3, 4, 10, 9, 5, 6)
    if mu == (2, 2, 2, 2, 2, 2):
        return Permutation(7, 8, 9, 10, 11, 12, 1, 2, 3, 4, 5, 6)
    if mu == (2, 2, 2, 2, 1, 1, 1, 1):
        return Permutation(2, 1, 4, 3, 9, 10, 11, 12, 5, 6, 7, 8)
    if mu == (5, 5, 2, 2):
        return Permutation(5, 6, 13, 14, 1, 2, 11, 12, 10, 9, 7, 8, 3, 4)
        q = sorted([w for w in Permutation.fpf_involutions(12) if des_n(w) == des_m(representative_m(mu))])
        x = q[0]
        print(q)
        return x


    a = 0
    w = Permutation()
    assert len(mu) % 2 == 0
    mu = [mu[i] for i in range(0, len(mu), 2)]
    for b in mu:
        for i in range(b):
            w *= Permutation.transposition(a + i + 1, a + 2 * b - i)
            if i % 2 != 0:
                s = Permutation.s_i(a + i)
                w = s * w * s
        a += 2 * b
    return w.star()


def bidirected_edges_m(w):
    assert w.is_fpf_involution()
    n = w.rank

    for i in range(2, n):
        s = Permutation.s_i(i - 1)
        t = Permutation.s_i(i)

        y = s * w * s
        if len(t * w * t) <= len(w) < len(y) < len(t * y * t):
            yield y, i

        if len(t * w * t) > len(w) > len(y) >= len(t * y * t):
            yield y, i

        y = t * w * t
        if len(s * w * s) <= len(w) < len(y) < len(s * y * s):
            yield y, i

        if len(s * w * s) > len(w) > len(y) >= len(s * y * s):
            yield y, i


def molecule_m(mu):
    w = representative_m(mu)
    ans = set()
    level = {w}
    while level:
        ans |= level
        level = {z for y in level for z, _ in bidirected_edges_m(y) if z not in ans}
    return ans


def get_molecules_m(n, verbose=False):
    ans = {}
    nn = double_fac(n - 1)
    for i, w in enumerate(Permutation.fpf_involutions(n)):
        mu = rsk(w)[0].shape()
        if mu not in ans:
            ans[mu] = set()
        ans[mu].add(w)

        if verbose:
            a = nn - i
            if a % 100 == 0:
                print(a)
    return ans


def bidirected_edges_n(w):
    assert w.is_fpf_involution()
    n = w.rank

    for i in range(2, n):
        s = Permutation.s_i(i - 1)
        t = Permutation.s_i(i)

        y = s * w * s
        if len(t * w * t) < len(w) < len(y) <= len(t * y * t):
            yield y, i

        if len(t * w * t) >= len(w) > len(y) > len(t * y * t):
            yield y, i

        y = t * w * t
        if len(s * w * s) < len(w) < len(y) <= len(s * y * s):
            yield y, i

        if len(s * w * s) >= len(w) > len(y) > len(s * y * s):
            yield y, i


def molecule_n(mu):
    w = representative_n(mu)
    ans = set()
    level = {w}
    while level:
        ans |= level
        level = {z for y in level for z, _ in bidirected_edges_n(y) if z not in ans}
    return ans


def get_molecules_n(n, verbose=False):
    ans = {}
    for mu in Partition.generate(n, even_parts=True):
        mu = Partition.transpose(mu)
        molecule = molecule_n(mu)
        ans[mu] = molecule
    return ans


def double_fac(n):
    if n <= 2:
        return 1
    return double_fac(n - 2) * n


def construct_molecular_correspondence(mu):
    mapping = {}
    w = representative_m(mu)
    mapping[w] = representative_n(mu)
    level = {w}
    a = 0
    while level:
        nextlevel = set()
        for w in level:
            edges_m = {i: y for y, i in bidirected_edges_m(w)}
            edges_n = {i: y for y, i in bidirected_edges_n(mapping[w])}
            if set(edges_m) != set(edges_n):
                print()
                print('shape:', mu)
                print('level:', a)
                print()
                print(w, set(edges_m), '?=', set(edges_n), mapping[w])
                print()
                assert set(edges_m) == set(edges_n)
            for i in set(edges_m) & set(edges_n):
                y = edges_m[i]
                z = edges_n[i]
                if y not in mapping:
                    nextlevel.add(y)
                    mapping[y] = z
                else:
                    assert mapping[y] == z
        level = nextlevel
        a += 1
    return mapping


def test_molecular_correspondence(k=10):
    for n in range(2, k + 1, 2):
        for mu in Partition.generate(n, even_parts=True):
            construct_molecular_correspondence(Partition.transpose(mu))
        # a = get_molecules_m(n)
        # b = get_molecules_n(n)
        # print(a)
        # print(b)
        # print(n)
        # print()
        # assert set(a) == set(b)
        # assert all(len(a[mu]) == len(b[mu]) for mu in a)


def test_get_molecules_m(n=10):
    ans = {}
    fpf = set(Permutation.fpf_involutions(n))
    while fpf:
        w = fpf.pop()
        mu = rsk(w)[0].shape()
        molecule = molecule_m(mu)
        fpf -= molecule
        assert mu not in ans
        ans[mu] = molecule
    bns = get_molecules_m(n)
    assert all(ans[mu] == bns[mu] for mu in bns)
    assert all(representative_m(mu) in bns[mu] for mu in bns)


def test_get_molecules_n(n=10):
    ans = {}
    fpf = set(Permutation.fpf_involutions(n))
    bns = get_molecules_n(n)
    for mu in bns:
        molecule = molecule_n(mu)
        fpf -= molecule
        ans[mu] = molecule
    assert len(fpf) == 0
    assert all(ans[mu] == bns[mu] for mu in bns)


def test_bidirected_edges_m(n=10):
    for w in Permutation.fpf_involutions(n):
        p, _ = rsk(w)
        for y, i in set(bidirected_edges_m(w)):
            print(w, '<--', i, '-->', y)
            q, _ = rsk(y)
            print(p)
            print(q)
            print()
            assert q == dual_equivalence(p, i)

