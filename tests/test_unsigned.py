from signed import SignedPermutation
from symmetric import SymmetricPolynomial
import pytest


def _test(r, upper):
    w = SignedPermutation.longest_unsigned(r)
    mu = tuple(i for i in range(r - 1, 0, -2))
    for n in range(upper):
        f = SymmetricPolynomial.stable_grothendieck_q(mu, n, n)
        g, rows = w.marked_stable_grothendieck(degree_bound=n)
        print(f)
        print()
        print(g)
        print()
        print(f - g)
        print()
        print()
        print()
        assert f == g


def test_one():
    _test(1, 8)


def test_two():
    _test(2, 8)


@pytest.mark.slow
def test_three():
    _test(3, 6)


@pytest.mark.slow
def test_four():
    _test(4, 6)


@pytest.mark.slow
def test_five():
    _test(5, 6)


@pytest.mark.slow
def test_six():
    _test(6, 6)
