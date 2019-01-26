from monomials import Monomial, Polynomial
from tableaux import Tableau


def test_basic():
    x = Monomial('x')
    assert type(x) in [Monomial, Polynomial]

    f = x * x
    assert type(f) in [Monomial, Polynomial]

    f = x + 1
    assert type(f) in [Monomial, Polynomial]

    f = x * 1
    assert type(f) in [Monomial, Polynomial]

    f = x * 0
    assert type(f) in [Monomial, Polynomial]

    f = x * -1
    assert type(f) in [Monomial, Polynomial]

    f = x * -1 + x
    assert type(f) in [Monomial, Polynomial]

    y = Monomial('y')
    assert type(y) in [Monomial, Polynomial]

    f = x * y
    assert type(f) in [Monomial, Polynomial]

    f = x / y
    assert type(f) in [Monomial, Polynomial]

    f = x / y + y / x
    assert type(f) in [Monomial, Polynomial]

    f = (x / y + y / x) * (x / y + y / x) * (x / y + y / x)
    assert type(f) in [Monomial, Polynomial]

    f = (x / y + y / x)**3
    assert type(f) in [Monomial, Polynomial]


def test_from_tableau():
    x1 = Monomial('x1')
    x2 = Monomial('x2')

    t = Tableau({(1, 1): 1, (1, 2): (1, 2)})
    assert Monomial.from_tableau(t) == x1 ** 2 * x2

    t = Tableau({(1, 1): (1, 2), (1, 2): 2})
    assert Monomial.from_tableau(t) == x1 * x2 ** 2


def test_symmetrize():
    f = Polynomial.base(Monomial('x1'))
    assert f.symmetrize(3) == Monomial('x1') + Monomial('x2') + Monomial('x3')


def test_slow_symmetric_functions():
    for mu in [(), (1,), (1, 1), (2,), (2, 1)]:
        for n in range(4):
            f = Polynomial.slow_schur(mu, n)
            g = Polynomial.slow_stable_grothendieck(mu, n)
            h = Polynomial.slow_schur_s(mu, n)
            k = Polynomial.slow_stable_grothendieck_s(mu, n)
            assert g.lowest_degree_terms() == f
            assert k.lowest_degree_terms() == h
            assert f.is_symmetric(n)
            assert g.is_symmetric(n)
            assert h.is_symmetric(n)
            assert k.is_symmetric(n)


def test_symmetric_functions():
    for mu in Tableau.generate_partitions(6):
        for n in range(6):
            f = Polynomial.schur(mu, n)
            g = Polynomial.stable_grothendieck(mu, n)

            fs = Polynomial.slow_schur(mu, n)
            gs = Polynomial.slow_stable_grothendieck(mu, n)
            print(mu, n)
            print(f)
            print(fs)
            print('*', f.symmetrize(n))
            print()
            h = Polynomial.schur_s(mu, n)
            k = Polynomial.stable_grothendieck_s(mu, n)
            assert g.lowest_degree_terms() == f
            assert k.lowest_degree_terms() == h
            assert f.symmetrize(n) == fs
            assert g.symmetrize(n) == gs
