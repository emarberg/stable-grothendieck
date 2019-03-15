from permutations import SignedPermutation
from symmetric import SymmetricPolynomial
import pytest


def test_generate():
    assert set(SignedPermutation.fpf_involutions(2)) == {SignedPermutation(2, 1)}
    assert set(SignedPermutation.fpf_involutions(3)) == set()
    assert set(SignedPermutation.fpf_involutions(4)) == {
        SignedPermutation(2, 1, 4, 3),
        SignedPermutation(3, 4, 1, 2),
        SignedPermutation(4, 3, 2, 1),
    }


def test_fpf_involution_shape():
    w = SignedPermutation(2, 1)
    assert w.fpf_involution_shape() == ()

    w = SignedPermutation(3, 4, 1, 2)
    assert w.fpf_involution_shape() == (1,)

    w = SignedPermutation(4, 3, 2, 1)
    assert w.fpf_involution_shape() == (2,)


def _test(it, upper):
    for w in it:
        nu = w.fpf_involution_shape()
        if nu:
            for n in range(upper):
                f = SymmetricPolynomial.stable_grothendieck_p(nu, n, n)
                g = w.symplectic_stable_grothendieck(degree_bound=n)
                print(w, nu)
                print()
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
    _test(SignedPermutation.fpf_involutions(1), 12)


def test_two():
    _test(SignedPermutation.fpf_involutions(2), 8)


def test_three():
    _test(SignedPermutation.fpf_involutions(3), 10)


@pytest.mark.slow
def test_four():
    _test(SignedPermutation.fpf_involutions(4), 8)


@pytest.mark.slow
def test_six():
    _test(SignedPermutation.fpf_involutions(6), 6)
