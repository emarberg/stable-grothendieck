from permutations import Permutation
from partitions import Partition


def test_grassmannian():
    w = Permutation.get_grassmannian(4, 4, 3, 2, 2, 2, 1)
    assert w.shape() == (4, 4, 3, 2, 2, 2, 1)

    shapes = set(Partition.subpartitions([4, 3, 2, 1]))
    gr = list(Permutation.grassmannians(5))
    assert len(gr) == len(shapes)
    assert {w.shape() for w in gr} == shapes


def test_inv_grassmannian():
    w = Permutation.get_inv_grassmannian(4, 3, 1)
    t = Permutation.transposition
    assert w == t(1, 5) * t(2, 6) * t(4, 7)
    assert w.involution_shape() == (4, 3, 1)

    shapes = set(Partition.subpartitions([5, 3, 1], True))
    gr = list(Permutation.inv_grassmannians(6))
    assert len(gr) == len(shapes)
    assert {w.involution_shape() for w in gr} == shapes


def test_fpf_grassmannian():
    w = Permutation.get_fpf_grassmannian(4, 3, 1)
    t = Permutation.transposition

    assert w == t(1, 6) * t(2, 7) * t(4, 8) * t(3, 5)
    assert w.fpf_involution_shape() == (4, 3, 1)

    shapes = set(Partition.subpartitions([6, 4, 2], True))
    gr = list(Permutation.fpf_grassmannians(8))
    assert len(gr) == len(shapes)
    assert {w.fpf_involution_shape() for w in gr} == shapes
