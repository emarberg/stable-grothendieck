from symmetric import SymmetricMonomial, SymmetricPolynomial
from tableaux import Tableau, Partition
from polynomials import Polynomial
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


def test_polynomials():
    m = SymmetricMonomial(5, ())
    assert m.polynomial() == Polynomial.one()

    m = SymmetricMonomial(5, (1,))
    x = Polynomial.x
    assert m.polynomial() == x(1) + x(2) + x(3) + x(4) + x(5)

    m = SymmetricMonomial(3, (2, 1,))
    x = Polynomial.y
    assert m.polynomial('y') == \
        x(1) * x(2)**2 + x(1)**2 * x(2) + \
        x(1) * x(3)**2 + x(1)**2 * x(3) + \
        x(3) * x(2)**2 + x(3)**2 * x(2)


def test_slow_symmetric_functions():
    for mu in [(), (1,), (1, 1), (2,), (2, 1)]:
        for n in range(4):
            f = SymmetricPolynomial._slow_schur(mu, n)
            g = SymmetricPolynomial._slow_stable_grothendieck(mu, n)
            h = SymmetricPolynomial._slow_schur_s(mu, n)
            k = SymmetricPolynomial._slow_stable_grothendieck_s(mu, n)
            assert g.lowest_degree_terms() == f
            assert k.lowest_degree_terms() == h


def test_slow_dual_stable_grothendieck():
    mu = (1,)
    g = SymmetricPolynomial._slow_dual_stable_grothendieck(mu, 3)
    print(g)
    assert g == SymmetricMonomial(3, (1,))

    mu = (2,)
    g = SymmetricPolynomial._slow_dual_stable_grothendieck(mu, 3)
    print(g)
    assert g == SymmetricMonomial(3, (1, 1)) + SymmetricMonomial(3, (2,))

    mu = (2, 1)
    for t in Tableau.semistandard_rpp(mu, 3):
        print(t)
    g = SymmetricPolynomial._slow_dual_stable_grothendieck(mu, 3)
    print(g)
    assert g == SymmetricMonomial(3, (1, 1)) + SymmetricMonomial(3, (2,)) + 2 * SymmetricMonomial(3, (1, 1, 1)) + SymmetricMonomial(3, (2, 1))


def test_slow_dual_stable_grothendieck_pq():
    mu = (1,)
    g = SymmetricPolynomial._slow_dual_stable_grothendieck_p(mu, 3)
    h = SymmetricPolynomial._slow_dual_stable_grothendieck_q(mu, 3)
    print(g)
    print(h)
    print()
    assert g == SymmetricMonomial(3, (1,))
    assert h == 2 * SymmetricMonomial(3, (1,))

    mu = (2,)
    for t in Tableau.semistandard_marked_rpp(mu, 3, False):
        print(t)
    g = SymmetricPolynomial._slow_dual_stable_grothendieck_p(mu, 3)
    print(g)
    assert g == \
        SymmetricPolynomial._slow_dual_stable_grothendieck((1, 1), 3) + \
        SymmetricPolynomial._slow_dual_stable_grothendieck((2,), 3)

    mu = (2, 1)
    for t in Tableau.semistandard_marked_rpp(mu, 3):
        print(t)
    g = SymmetricPolynomial._slow_dual_stable_grothendieck_p(mu, 3)
    h = SymmetricPolynomial._slow_dual_stable_grothendieck_q(mu, 3)
    print(g)
    print(h)
    print()
    assert g == SymmetricMonomial(3, (1, 1)) + SymmetricMonomial(3, (2,)) + 2 * SymmetricMonomial(3, (1, 1, 1)) + SymmetricMonomial(3, (2, 1))

    mu = (3, 2)
    g = SymmetricPolynomial._slow_dual_stable_grothendieck_p(mu, 5)
    h = SymmetricPolynomial._slow_dual_stable_grothendieck_q(mu, 5)

    mu = (3, 2, 1)
    assert SymmetricPolynomial._slow_dual_stable_grothendieck_p(mu, 6) == \
        SymmetricPolynomial._slow_dual_stable_grothendieck(mu, 6)


def test_decompose():
    schur = SymmetricPolynomial.schur
    schur_P = SymmetricPolynomial.schur_p  # noqa
    dual_grothendieck = SymmetricPolynomial.dual_stable_grothendieck
    grothendieck = SymmetricPolynomial.stable_grothendieck
    GP = SymmetricPolynomial.stable_grothendieck_p  # noqa
    GQ = SymmetricPolynomial.stable_grothendieck_q  # noqa

    n = 6
    mu = (3, 2, 1)
    nu = (4, 2, 1)

    f = schur(mu, n) * schur(mu, n)
    exp = SymmetricPolynomial.schur_expansion(f)
    assert f == sum([schur(a, n) * coeff for a, coeff in exp.items()])

    f = schur_P(nu, n)
    exp = SymmetricPolynomial.schur_expansion(f)
    assert f == sum([schur(a, n) * coeff for a, coeff in exp.items()])

    f = GP(nu, n)
    exp = SymmetricPolynomial.grothendieck_expansion(f)
    assert f == sum([grothendieck(a, n) * coeff for a, coeff in exp.items()])

    f = GQ(nu, n)
    exp = SymmetricPolynomial.grothendieck_expansion(f)
    assert f == sum([grothendieck(a, n) * coeff for a, coeff in exp.items()])

    f = SymmetricPolynomial.dual_stable_grothendieck_p(nu, n)
    exp = SymmetricPolynomial.dual_grothendieck_expansion(f)
    assert f == sum([dual_grothendieck(a, n) * coeff for a, coeff in exp.items()])

    f = SymmetricPolynomial.dual_stable_grothendieck_q(nu, n)
    exp = SymmetricPolynomial.dual_grothendieck_expansion(f)
    assert f == sum([dual_grothendieck(a, n) * coeff for a, coeff in exp.items()])


@pytest.mark.slow
def test_symmetric_functions():
    nn = 6
    for mu in itertools.chain(*[Partition.generate(n + 1) for n in range(nn)]):
        for n in range(nn):
            print(mu, n)
            print()

            f = SymmetricPolynomial.schur(mu, n)
            g = SymmetricPolynomial.stable_grothendieck(mu, n)
            h = SymmetricPolynomial.dual_stable_grothendieck(mu, n)

            fs = SymmetricPolynomial._slow_schur(mu, n)
            gs = SymmetricPolynomial._slow_stable_grothendieck(mu, n)
            hs = SymmetricPolynomial._slow_dual_stable_grothendieck(mu, n)

            print(f)
            print(fs)
            print()
            print(g)
            print(gs)
            print()
            print(h)
            print(hs)
            print()
            print()

            assert f == fs
            assert g == gs
            assert h == hs

            h = SymmetricPolynomial.schur_s(mu, n)
            k = SymmetricPolynomial.stable_grothendieck_s(mu, n)
            assert g.lowest_degree_terms() == f
            assert k.lowest_degree_terms() == h


@pytest.mark.slow
def test_strict_symmetric_functions():
    nn = 5
    for mu in itertools.chain(*[Partition.generate(n + 1, strict=True) for n in range(nn)]):
        for n in range(nn):
            print(mu, n)
            print()

            # Schur-P and GP

            f = SymmetricPolynomial.schur_p(mu, n)
            g = SymmetricPolynomial.stable_grothendieck_p(mu, n)
            h = SymmetricPolynomial.dual_stable_grothendieck_p(mu, n)

            fs = SymmetricPolynomial._slow_schur_p(mu, n)
            gs = SymmetricPolynomial._slow_stable_grothendieck_p(mu, n)
            hs = SymmetricPolynomial._slow_dual_stable_grothendieck_p(mu, n)

            print(f)
            print(fs)
            print()
            print(g)
            print(gs)
            print()
            print(h)
            print(hs)
            print()
            print()
            print()

            assert f == fs
            assert g == gs
            assert h == hs

            # Schur-Q and GQ

            f = SymmetricPolynomial.schur_q(mu, n)
            g = SymmetricPolynomial.stable_grothendieck_q(mu, n)
            h = SymmetricPolynomial.dual_stable_grothendieck_q(mu, n)

            fs = SymmetricPolynomial._slow_schur_q(mu, n)
            gs = SymmetricPolynomial._slow_stable_grothendieck_q(mu, n)
            hs = SymmetricPolynomial._slow_dual_stable_grothendieck_q(mu, n)

            print(f)
            print(fs)
            print()
            print(g)
            print(gs)
            print()
            print(h)
            print(hs)
            print()
            print()
            print()

            assert f == fs
            assert g == gs
            assert h == hs
