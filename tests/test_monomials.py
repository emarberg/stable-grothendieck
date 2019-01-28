from monomials import Monomial, Polynomial
from tableaux import Tableau, Partition


def test_destandardize():
    mu = (2, 1)
    n = 3
    assert set(Monomial._destandardize(mu, n)) == {
        (1, 2, 0),
        (1, 0, 2),
        (0, 1, 2),
        (2, 1, 0),
        (2, 0, 1),
        (0, 2, 1),
    }
    assert len(list(Monomial._destandardize(mu, n))) == 6

    mu = (2, 2)
    n = 3
    assert set(Monomial._destandardize(mu, n)) == {
        (2, 2, 0),
        (2, 0, 2),
        (0, 2, 2),
    }
    assert len(list(Monomial._destandardize(mu, n))) == 3


def test_basic():
    x = Monomial(1, (1,))
    assert type(x) in [Monomial, Polynomial]

    f = x * x
    print(x)
    print(f)
    print(type(f), f, Monomial(1, (2,)), f - Monomial(1, (2,)))
    assert type(f) in [Monomial, Polynomial]
    assert f == Monomial(1, (2,))

    x = Monomial(2, (1,))
    print(x)
    print(x * x)
    assert x * x == Monomial(2, (2,)) + 2 * Monomial(2, (1, 1))

    f = x + 1
    assert type(f) in [Monomial, Polynomial]

    f = x * 1
    assert type(f) in [Monomial, Polynomial]
    assert f == x

    f = x * 0
    assert type(f) in [Monomial, Polynomial]

    f = x * -1
    assert type(f) in [Monomial, Polynomial]

    f = x * -1 + x
    assert type(f) in [Monomial, Polynomial]
    assert f.is_zero()

    x = Monomial(2, (1,))
    y = Monomial(2, (1, 1))
    assert type(y) in [Monomial, Polynomial]

    f = x * y
    assert type(f) in [Monomial, Polynomial]
    assert f == Monomial(2, (2, 1))


def test_slow_symmetric_functions():
    for mu in [(), (1,), (1, 1), (2,), (2, 1)]:
        for n in range(4):
            f = Polynomial.slow_schur(mu, n)
            g = Polynomial.slow_stable_grothendieck(mu, n)
            h = Polynomial.slow_schur_s(mu, n)
            k = Polynomial.slow_stable_grothendieck_s(mu, n)
            assert g.lowest_degree_terms() == f
            assert k.lowest_degree_terms() == h


def test_symmetric_functions():
    for mu in Partition.generate(6):
        for n in range(6):
            f = Polynomial.schur(mu, n)
            g = Polynomial.stable_grothendieck(mu, n)

            fs = Polynomial.slow_schur(mu, n)
            gs = Polynomial.slow_stable_grothendieck(mu, n)

            print(mu, n)
            print(f)
            print(fs)
            print()

            assert f == fs
            assert g == gs

            h = Polynomial.schur_s(mu, n)
            k = Polynomial.stable_grothendieck_s(mu, n)
            assert g.lowest_degree_terms() == f
            assert k.lowest_degree_terms() == h
