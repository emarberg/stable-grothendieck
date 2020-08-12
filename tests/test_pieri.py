from tableaux import Tableau
from partitions import Partition
from vectors import Vector
import pytest
import utils
import time
from collections import defaultdict


def pairs_lhs(n):
    lhs = defaultdict(list)
    for nu in Partition.all(n, strict=True):
        for lam in Partition.decompose_shifted_shape_by_rims(nu):
            for mu in Partition.decompose_shifted_shape_by_vertical_strips(lam):
                if len(mu) == len(lam):
                    lhs[(nu, mu)].append(lam)
    return lhs


def pairs_rhs(n):
    rhs = defaultdict(list)
    for nu in Partition.all(n, strict=True):
        for lam in Partition.decompose_shifted_shape_by_vertical_strips(nu):
            if len(lam) == len(nu):
                for mu in Partition.decompose_shifted_shape_by_rims(lam):
                    rhs[(nu, mu)].append(lam)
    return rhs


def seconds(t):
    return str(int(1000 * t) / 1000.0)


def pairs_io(n=5):
    assert 0 <= n < 256

    def tobytes(nu, lam, mu):
        return bytes(nu + (0,) + lam + (0,) + mu + (0,))

    def frombytes(b, i):
        nu = []
        while i < len(b) and b[i] != 0:
            nu.append(b[i])
            i += 1
        i += 1
        lam = []
        while i < len(b) and b[i] != 0:
            lam.append(b[i])
            i += 1
        i += 1
        mu = []
        while i < len(b) and b[i] != 0:
            mu.append(b[i])
            i += 1
        i += 1
        return tuple(nu), tuple(lam), tuple(mu), i

    def read(filename):
        print()
        print('Trying to read from `%s`' % filename)
        t0 = time.time()
        with open(filename, 'rb') as file:
            b = bytearray(file.read())
        print('* Opened file in %s seconds' % seconds(time.time() - t0))
        ans = defaultdict(list)
        i = 0
        while i < len(b):
            nu, lam, mu, i = frombytes(b, i)
            ans[(nu, mu)].append(lam)
        print('* Succeeded in %s seconds' % seconds(time.time() - t0))
        print()
        return ans

    def write(file, ans):
        t0 = time.time()
        b = bytearray()
        for nu, mu in ans:
            for lam in ans[(nu, mu)]:
                b += tobytes(nu, lam, mu)
        print('* Writing to file `%s`' % file)
        with open(file, 'wb') as f:
            f.write(b)
        print('* Succeeded in %s seconds' % seconds(time.time() - t0))
        print()

    directory = '/Users/emarberg/examples/gpq/'
    lfile = directory + 'rims_then_strips_%s.b' % n
    rfile = directory + 'strips_then_rims_%s.b' % n

    try:
        lhs = read(lfile)
    except FileNotFoundError:
        print('* Failed, computing instead')
        lhs = pairs_lhs(n)
        write(lfile, lhs)

    try:
        rhs = read(rfile)
    except FileNotFoundError:
        print('* Failed, computing instead')
        rhs = pairs_rhs(n)
        write(rfile, rhs)

    return lhs, rhs


def test_rims():
    assert set(Partition.rims((), 0)) == {()}
    assert set(Partition.rims((), 1)) == {(), ((1, 1),)}
    assert set(Partition.rims((), 2)) == {(), ((1, 1),), ((1, 1), (1, 2))}
    assert set(Partition.rims((1,), 1)) == {(), ((1, 2),), ((1, 2), (2, 2))}
    assert set(Partition.rims((3, 2))) == {(), ((1, 4),), ((1, 4), (2, 4)), ((3, 3),), ((1, 4), (3, 3)), ((1, 4), (2, 4), (3, 3)),  ((1, 4), (2, 4), (3, 3), (3, 4))} # noqa


def test_KOG():  # noqa
    mu = (6, 4, 1)
    assert set(Tableau.KOG(0, mu)) == {Tableau()}

    assert set(Tableau.KOG(1, ())) == {Tableau({(1, 1): 1})}

    nu = (7, 6, 3, 1)
    assert set(Tableau.KOG_by_shape(5, mu)[nu]) == {
        Tableau({
            (1, 7): 1,
            (2, 6): 1,
            (2, 7): 5,
            (3, 4): 2,
            (3, 5): 4,
            (4, 4): 3,
        }),
        Tableau({
            (1, 7): 1,
            (2, 6): 2,
            (2, 7): 5,
            (3, 4): 2,
            (3, 5): 4,
            (4, 4): 3,
        }),
        Tableau({
            (1, 7): 1,
            (2, 6): 2,
            (2, 7): 5,
            (3, 4): 3,
            (3, 5): 5,
            (4, 4): 4,
        }),
        Tableau({
            (1, 7): 1,
            (2, 6): 2,
            (2, 7): 5,
            (3, 4): 3,
            (3, 5): 4,
            (4, 4): 4,
        }),
        Tableau({
            (1, 7): 1,
            (2, 6): 4,
            (2, 7): 5,
            (3, 4): 1,
            (3, 5): 3,
            (4, 4): 2,
        }),
        Tableau({
            (1, 7): 1,
            (2, 6): 4,
            (2, 7): 5,
            (3, 4): 2,
            (3, 5): 4,
            (4, 4): 3,
        }),
        Tableau({
            (1, 7): 1,
            (2, 6): 4,
            (2, 7): 5,
            (3, 4): 2,
            (3, 5): 3,
            (4, 4): 3,
        })
    }


def test_KLG():  # noqa
    mu = (6, 4, 1)
    assert set(Tableau.KLG(0, mu)) == {Tableau()}

    assert set(Tableau.KLG(1, ())) == {Tableau({(1, 1): 1})}

    nu = (7, 6, 3, 1)
    assert set(Tableau.KLG_by_shape(5, mu)[nu]) == {
        Tableau({
            (1, 7): -1,
            (2, 6): -1,
            (2, 7): 5,
            (3, 4): -2,
            (3, 5): 4,
            (4, 4): 3,
        }),
        Tableau({
            (1, 7): -1,
            (2, 6): -2,
            (2, 7): 5,
            (3, 4): -2,
            (3, 5): 4,
            (4, 4): 3,
        }),
        Tableau({
            (1, 7): -1,
            (2, 6): -2,
            (2, 7): 5,
            (3, 4): -3,
            (3, 5): 5,
            (4, 4): 4,
        }),
        Tableau({
            (1, 7): -1,
            (2, 6): -2,
            (2, 7): 5,
            (3, 4): -3,
            (3, 5): 4,
            (4, 4): 4,
        }),
        Tableau({
            (1, 7): -1,
            (2, 6): 4,
            (2, 7): 5,
            (3, 4): -1,
            (3, 5): 3,
            (4, 4): 2,
        }),
        Tableau({
            (1, 7): -1,
            (2, 6): 4,
            (2, 7): 5,
            (3, 4): -2,
            (3, 5): 4,
            (4, 4): 3,
        }),
        Tableau({
            (1, 7): -1,
            (2, 6): 4,
            (2, 7): 5,
            (3, 4): -2,
            (3, 5): 3,
            (4, 4): 3,
        }),
        Tableau({
            (1, 7): -1,
            (2, 6): 4,
            (2, 7): 5,
            (3, 4): -2,
            (3, 5): 3,
            (4, 4): 2,
        }),
        Tableau({
            (1, 7): -1,
            (2, 6): -2,
            (2, 7): 5,
            (3, 4): -3,
            (3, 5): 4,
            (4, 4): 3,
        })
    }


def test_fast_KOG(): # noqa
    for mu in Partition.all(7, strict=True):
        for p in range(7):
            counts = Tableau.KOG_counts_by_shape(p, mu)
            tabs = Tableau.KOG_by_shape(p, mu)
            print(mu, p)
            print(counts)
            print(tabs)
            print()
            assert set(counts) == set(tabs)
            for nu in counts:
                assert counts[nu] == len(tabs[nu])


def test_fast_KLG(): # noqa
    for mu in Partition.all(7, strict=True):
        for p in range(7):
            counts = Tableau.KLG_counts_by_shape(p, mu)
            tabs = Tableau.KLG_by_shape(p, mu)
            print(mu, p)
            print(counts)
            print(tabs)
            print()
            assert set(counts) == set(tabs)
            for nu in counts:
                print('*', nu)
                assert counts[nu] == len(tabs[nu])


def test_GP_pieri(): # noqa
    for mu in Partition.all(6, strict=True):
        for p in [0, 1, 2]:
            ans = Vector()
            for nu, tabs in Tableau.KOG_by_shape(p, mu).items():
                ans += Vector({nu: utils.beta**(sum(nu) - sum(mu) - p) * len(tabs)})
            f = utils.GP(len(mu) + 1, mu) * utils.GP(len(mu) + 1, (p,))
            assert ans == utils.GP_expansion(f)


def test_GQ_pieri(): # noqa
    for mu in Partition.all(5, strict=True):
        for p in [0, 1, 2]:
            ans = Vector()
            for nu, tabs in Tableau.KLG_by_shape(p, mu).items():
                ans += Vector({nu: utils.beta**(sum(nu) - sum(mu) - p) * len(tabs)})
            f = utils.GQ(len(mu) + 1, mu) * utils.GQ(len(mu) + 1, (p,))
            assert ans == utils.GQ_expansion(f)


@pytest.mark.slow
def test_GP_pieri_slow(): # noqa
    for mu in Partition.all(10, strict=True):
        for p in [1, 2, 3]:
            ans = Vector()
            for nu, tabs in Tableau.KOG_by_shape(p, mu).items():
                ans += Vector({nu: utils.beta**(sum(nu) - sum(mu) - p) * len(tabs)})
            f = utils.GP(len(mu) + 1, mu) * utils.GP(len(mu) + 1, (p,))
            assert ans == utils.GP_expansion(f)


@pytest.mark.slow
def test_GQ_pieri_slow(): # noqa
    for mu in Partition.all(10, strict=True):
        for p in [1, 2, 3]:
            ans = Vector()
            for nu, tabs in Tableau.KLG_by_shape(p, mu).items():
                ans += Vector({nu: utils.beta**(sum(nu) - sum(mu) - p) * len(tabs)})
            f = utils.GQ(len(mu) + 1, mu) * utils.GQ(len(mu) + 1, (p,))
            assert ans == utils.GQ_expansion(f)
