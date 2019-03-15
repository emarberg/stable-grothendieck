from polynomials import Polynomial
from symmetric import SymmetricPolynomial
from tableaux import Partition
import pytest


def _test(n_max, v_max):
    for n in range(n_max + 1):
        for v in range(1, v_max + 1):
            a, b = [], []
            for k in range(max(n, 1)):
                for mu in Partition.generate(k):
                    if mu and mu[0] + k > n:
                        continue
                    a += [SymmetricPolynomial.dual_stable_grothendieck(mu, v).truncate(n)]
                    b += [SymmetricPolynomial.stable_grothendieck(mu, v).truncate(n)]

            s = Polynomial()
            for i in range(len(a)):
                s += (a[i].polynomial('x') * b[i].polynomial('y')).truncate(n)
                print(s)
                print()

            t = Polynomial.one()
            x = Polynomial.x
            y = Polynomial.y
            for i in range(1, v + 1):
                for j in range(1, v + 1):
                    a = x(i) * y(j)
                    term = Polynomial.one()
                    for e in range(1, n + 1):
                        term += a**e
                    t = (t * term).truncate(n)
            print(t)
            print()
            print(s - t)
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
                    a += [SymmetricPolynomial.dual_stable_grothendieck_p(mu, v).truncate(n)]
                    b += [SymmetricPolynomial.stable_grothendieck_q(mu, v).truncate(n)]

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
    _test_shifted_q(6, 3)


@pytest.mark.slow
def test_grothendieck_cauchy_q_slow():
    _test_shifted_q(8, 5)


def _test_shifted_p(n_max, v_max):
    for n in range(n_max + 1):
        for v in range(1, v_max + 1):
            a, b = [], []
            for k in range(max(n, 1)):
                for mu in Partition.generate(k, strict=True):
                    a += [SymmetricPolynomial.dual_stable_grothendieck_q(mu, v).truncate(n)]
                    b += [SymmetricPolynomial.stable_grothendieck_p(mu, v).truncate(n)]

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
    _test_shifted_q(6, 3)


@pytest.mark.slow
def test_grothendieck_cauchy_p_slow():
    _test_shifted_q(8, 5)
