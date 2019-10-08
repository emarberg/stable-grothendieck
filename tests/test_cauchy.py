from polynomials import Polynomial
from symmetric import SymmetricPolynomial
from tableaux import Partition
import pytest


def _test(n_max, v_max):
    for n in range(n_max + 1):
        for v in range(1, v_max + 1):
            for lam1 in Partition.all(n):
                for lam2 in Partition.all(n):
                    print()
                    print()
                    print('* n =', n, ', mu =', lam1, ', nu =', lam2)
                    print()

                    print('Computing LHS . . .')
                    print()

                    s = Polynomial()
                    for mu in Partition.all(max(n, 1)):
                        a = SymmetricPolynomial.dual_stable_grothendieck(v, mu, lam2).truncate(n)
                        b = SymmetricPolynomial.stable_grothendieck_doublebar(v, mu, lam1).truncate(n)
                        s += (a.polynomial('y') * b.polynomial('x')).truncate(n)
                        print('  ', mu, ':', s)
                        print()
                    print('LHS =', s)
                    print()
                    print()

                    print('Computing RHS . . .')
                    print()

                    t = Polynomial()
                    for kappa in Partition.subpartitions(lam2):
                        a = SymmetricPolynomial.dual_stable_grothendieck(v, lam1, kappa).truncate(n)
                        b = SymmetricPolynomial.stable_grothendieck_doublebar(v, lam2, kappa).truncate(n)
                        t += (a.polynomial('y') * b.polynomial('x')).truncate(n)
                        print('  ', kappa, ':', t)
                        print()

                    x = Polynomial.x
                    y = Polynomial.y
                    for i in range(1, v + 1):
                        for j in range(1, v + 1):
                            a = x(i) * y(j)
                            term = Polynomial.one()
                            for e in range(1, n + 1):
                                term += a**e
                            t = (t * term).truncate(n)
                    print('RHS =', t)
                    print()
                    print()
                    print('diff =', s - t)
                    print()
                    assert s == t


def test_grothendieck_cauchy():
    _test(6, 3)


@pytest.mark.slow
def test_grothendieck_cauchy_slow():
    _test(9, 6)


def kernel(n, v):
    t = Polynomial.one()
    x = Polynomial.x
    y = Polynomial.y
    for i in range(1, v + 1):
        for j in range(1, v + 1):
            a = x(i) * y(j)
            base = Polynomial.one() - x(i) + x(i) * y(j)
            term = Polynomial.one()
            for e in range(1, n + 1):
                term += a**e
            base *= term
            term = Polynomial.one()
            for e in range(1, n + 1):
                term += x(i)**e
            base *= term
            t = (t * base).truncate(n)
    return t


def _test_shifted_q(n_max, v_max):
    for n in range(n_max + 1):
        for v in range(1, v_max + 1):
            a, b = [], []
            for k in range(max(n, 1)):
                for mu in Partition.generate(k, strict=True):
                    a += [SymmetricPolynomial.dual_stable_grothendieck_p(v, mu).truncate(n)]
                    b += [SymmetricPolynomial.stable_grothendieck_q(v, mu).truncate(n)]

            s = Polynomial()
            for i in range(len(a)):
                s += (a[i].polynomial('y') * b[i].polynomial('x')).truncate(n)

            print(n, v)
            print(s)
            print()

            t = kernel(n, v)

            print(t)
            print()
            print('*', s - t)
            print()
            print()
            assert s == t


def test_grothendieck_cauchy_q():
    _test_shifted_q(5, 3)


@pytest.mark.slow
def test_grothendieck_cauchy_q_slow():
    _test_shifted_q(8, 5)


def _test_shifted_p(n_max, v_max):
    for n in range(n_max + 1):
        for v in range(1, v_max + 1):
            a, b = [], []
            for k in range(max(n, 1)):
                for mu in Partition.generate(k, strict=True):
                    a += [SymmetricPolynomial.dual_stable_grothendieck_q(v, mu).truncate(n)]
                    b += [SymmetricPolynomial.stable_grothendieck_p(v, mu).truncate(n)]

            s = Polynomial()
            for i in range(len(a)):
                s += (a[i].polynomial('y') * b[i].polynomial('x')).truncate(n)

            print(n, v)
            print(s)
            print()

            t = kernel(n, v)

            print(t)
            print()
            print('*', s - t)
            print()
            print()
            assert s == t


def test_grothendieck_cauchy_p():
    _test_shifted_q(5, 3)


@pytest.mark.slow
def test_grothendieck_cauchy_p_slow():
    _test_shifted_q(8, 5)
