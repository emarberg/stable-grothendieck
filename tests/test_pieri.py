from tableaux import Tableau
from partitions import Partition


def test_rims():
    assert set(Partition.rims((), 0)) == {()}
    assert set(Partition.rims((), 1)) == {(), ((1, 1),)}
    assert set(Partition.rims((), 2)) == {(), ((1, 1),), ((1, 1), (1, 2))}
    assert set(Partition.rims((1,), 1)) == {(), ((1, 2),), ((2, 2), (1, 2))}
    assert set(Partition.rims((3, 2))) == {(), ((1, 4),), ((2, 4), (1, 4)), ((3, 3),), ((3, 3), (1, 4)), ((3, 3), (2, 4), (1, 4)), ((3, 3), (3, 4), (2, 4), (1, 4))} # noqa
