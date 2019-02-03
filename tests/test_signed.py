from signed import SignedPermutation
from symmetric import Polynomial
import pytest


def _test(r, upper):
    w = SignedPermutation.longest_element(r)
    mu = tuple(i for i in range(r, 0, -1))
    for n in range(upper):
        f = Polynomial.stable_grothendieck_s(mu, n)
        f = Polynomial({
            mon: f[mon] for mon in f if mon.degree() <= n
        })
        g = w.signed_involution_stable_grothendieck(degree_bound=n)
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
    _test(1, 12)


def test_two():
    _test(2, 8)


@pytest.mark.slow
def test_three():
    _test(3, 10)


@pytest.mark.slow
def test_four():
    _test(4, 10)


@pytest.mark.slow
def test_five():
    _test(5, 8)


@pytest.mark.slow
def test_six():
    _test(6, 6)
