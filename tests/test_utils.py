from utils import G, GP, GS
from permutations import Permutation


def test_grothendieck():
    w = Permutation(4, 3, 2, 1)
    f = G((3, 2, 1), 3)
    print(f)
    g = G(w, 3)
    print(g)
    assert f == g


def test_grothendieck_p():
    w = Permutation(4, 3, 2, 1)
    f = GP((2,), 3)
    print(f)
    g = GP(w, 3)
    print(g)
    assert f == g


def test_grothendieck_s():
    w = Permutation(-1)
    f = GS((1,), 3)
    print(f)
    g = GS(w, 3)
    print(g)
    assert f == g

