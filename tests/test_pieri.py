from tableaux import Tableau
from partitions import Partition
from vectors import Vector
import pytest
import utils


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
