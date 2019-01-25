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


def test_symmetric_functions():
    for mu in [(), (1,), (2,), (2, 1), (3, 2, 2, 1), (3, 3, 2, 1, 1)]:
        for n in range(7):
            f = Polynomial.schur(mu, n)
            g = Polynomial.stable_grothendieck(mu, n)
            print(mu, n)
            print(f)
            print(g)
            print()
            assert g.lowest_degree_terms() == f
            # assert g.is_symmetric(n)
            # assert f.is_symmetric(n)
