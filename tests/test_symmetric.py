from symmetric import SymmetricMonomial, SymmetricPolynomial
from tableaux import Tableau, Partition
import itertools
import pytest


def test_destandardize():
    mu = (2, 1)
    n = 3
    assert set(SymmetricMonomial._destandardize(mu, n)) == {
        (1, 2, 0),
        (1, 0, 2),
        (0, 1, 2),
        (2, 1, 0),
        (2, 0, 1),
        (0, 2, 1),
    }
    assert len(list(SymmetricMonomial._destandardize(mu, n))) == 6

    mu = (2, 2)
    n = 3
    assert set(SymmetricMonomial._destandardize(mu, n)) == {
        (2, 2, 0),
        (2, 0, 2),
        (0, 2, 2),
    }
    assert len(list(SymmetricMonomial._destandardize(mu, n))) == 3


def test_basic():
    x = SymmetricMonomial(1, (1,))
    assert type(x) in [SymmetricMonomial, SymmetricPolynomial]

    f = x * x
    print(x)
    print(f)
    print(type(f), f, SymmetricMonomial(1, (2,)), f - SymmetricMonomial(1, (2,)))
    assert type(f) in [SymmetricMonomial, SymmetricPolynomial]
    assert f == SymmetricMonomial(1, (2,))

    x = SymmetricMonomial(2, (1,))
    print(x)
    print(x * x)
    assert x * x == SymmetricMonomial(2, (2,)) + 2 * SymmetricMonomial(2, (1, 1))

    f = x + 1
    assert type(f) in [SymmetricMonomial, SymmetricPolynomial]

    f = x * 1
    assert type(f) in [SymmetricMonomial, SymmetricPolynomial]
    assert f == x

    f = x * 0
    assert type(f) in [SymmetricMonomial, SymmetricPolynomial]

    f = x * -1
    assert type(f) in [SymmetricMonomial, SymmetricPolynomial]

    f = x * -1 + x
    assert type(f) in [SymmetricMonomial, SymmetricPolynomial]
    assert f.is_zero()

    x = SymmetricMonomial(2, (1,))
    y = SymmetricMonomial(2, (1, 1))
    assert type(y) in [SymmetricMonomial, SymmetricPolynomial]

    f = x * y
    assert type(f) in [SymmetricMonomial, SymmetricPolynomial]
    assert f == SymmetricMonomial(2, (2, 1))


def test_slow_symmetric_functions():
    for mu in [(), (1,), (1, 1), (2,), (2, 1)]:
        for n in range(4):
            f = SymmetricPolynomial.slow_schur(mu, n)
            g = SymmetricPolynomial.slow_stable_grothendieck(mu, n)
            h = SymmetricPolynomial.slow_schur_s(mu, n)
            k = SymmetricPolynomial.slow_stable_grothendieck_s(mu, n)
            assert g.lowest_degree_terms() == f
            assert k.lowest_degree_terms() == h


def test_slow_dual_stable_grothendieck():
    mu = (1,)
    g = SymmetricPolynomial.slow_dual_stable_grothendieck(mu, 3)
    print(g)
    assert g == SymmetricMonomial(3, (1,))

    mu = (2,)
    g = SymmetricPolynomial.slow_dual_stable_grothendieck(mu, 3)
    print(g)
    assert g == SymmetricMonomial(3, (1, 1)) + SymmetricMonomial(3, (2,))

    mu = (2, 1)
    for t in Tableau.semistandard_rpp(mu, 3):
        print(t)
    g = SymmetricPolynomial.slow_dual_stable_grothendieck(mu, 3)
    print(g)
    assert g == SymmetricMonomial(3, (1, 1)) + SymmetricMonomial(3, (2,)) + 2 * SymmetricMonomial(3, (1, 1, 1)) + SymmetricMonomial(3, (2, 1))


def test_grothendieck_cauchy():
    n = 6
    v = 2
    s = SymmetricPolynomial()
    for k in range(n + 1):
        for mu in Partition.generate(k):
            s += SymmetricPolynomial.slow_dual_stable_grothendieck(mu, v) * SymmetricPolynomial.stable_grothendieck(mu, v)
            s = s.truncate(n)
            print(s)
            print()


@pytest.mark.slow
def test_symmetric_functions():
    nn = 6
    for mu in itertools.chain(*[Partition.generate(n + 1) for n in range(nn)]):
        for n in range(nn):
            print(mu, n)
            print()

            f = SymmetricPolynomial.schur(mu, n)
            g = SymmetricPolynomial.stable_grothendieck(mu, n)

            fs = SymmetricPolynomial.slow_schur(mu, n)
            gs = SymmetricPolynomial.slow_stable_grothendieck(mu, n)

            print(f)
            print(fs)
            print()
            print(g)
            print(gs)
            print()
            print()
            print()

            assert f == fs
            assert g == gs

            h = SymmetricPolynomial.schur_s(mu, n)
            k = SymmetricPolynomial.stable_grothendieck_s(mu, n)
            assert g.lowest_degree_terms() == f
            assert k.lowest_degree_terms() == h


@pytest.mark.slow
def test_strict_symmetric_functions():
    nn = 6
    for mu in itertools.chain(*[Partition.generate(n + 1, strict=True) for n in range(nn)]):
        for n in range(nn):
            print(mu, n)
            print()

            # Schur-P and GP

            f = SymmetricPolynomial.schur_p(mu, n)
            g = SymmetricPolynomial.stable_grothendieck_p(mu, n)

            fs = SymmetricPolynomial.slow_schur_p(mu, n)
            gs = SymmetricPolynomial.slow_stable_grothendieck_p(mu, n)

            print(f)
            print(fs)
            print()
            print(g)
            print(gs)
            print()
            print()
            print()

            assert f == fs
            assert g == gs

            # Schur-Q and GQ

            f = SymmetricPolynomial.schur_q(mu, n)
            g = SymmetricPolynomial.stable_grothendieck_q(mu, n)

            fs = SymmetricPolynomial.slow_schur_q(mu, n)
            gs = SymmetricPolynomial.slow_stable_grothendieck_q(mu, n)

            print(f)
            print(fs)
            print()
            print(g)
            print(gs)
            print()
            print()
            print()

            assert f == fs
            assert g == gs
