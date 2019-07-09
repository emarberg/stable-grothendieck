from tableaux import Tableau, Partition, MarkedReversePlanePartition


def test_transpose():
    mu = ()
    assert Partition.transpose(mu) == mu

    mu = (1,)
    assert Partition.transpose(mu) == mu

    mu = (2,)
    assert Partition.transpose(mu) == (1, 1)

    mu = (5, 3, 3, 2, 2, 1, 1, 1, 1)
    assert Partition.transpose(mu) == (9, 5, 3, 1, 1)


def test_standardize():
    t = Tableau({
        (1, 1): 1, (1, 2): (-3, 2), (1, 3): (4, 5), (1, 4): (-7, -6, 7), (1, 5): (9, 10),
        (2, 2): (4, 4), (2, 3): (6, -7), (2, 4): -8,
        (3, 3): 8
    })
    assert t.standardize() == Tableau({
        (1, 1): 1, (1, 2): (-3, 2), (1, 3): (6, 7), (1, 4): (-10, -8, 12), (1, 5): (15, 16),
        (2, 2): (4, 5), (2, 3): (9, -11), (2, 4): -13,
        (3, 3): 14
    })


def test_generate_partitions():
    assert set(Partition.generate(0)) == {()}
    assert set(Partition.generate(1)) == {(1,)}
    assert set(Partition.generate(2)) == {(1, 1), (2,)}
    assert set(Partition.generate(3)) == {(1, 1, 1), (2, 1), (3,)}
    assert set(Partition.generate(4)) == {(1, 1, 1, 1), (2, 1, 1), (2, 2), (3, 1), (4,)}


def test_horizontal_strips():
    mu = tuple()
    assert list(Tableau._horizontal_strips(mu)) == [
        (tuple(), set(), []),
    ]

    mu = (1,)
    assert list(Tableau._horizontal_strips(mu)) == [
        (tuple(), {(1, 1)}, []),
        ((1,), set(), [(1, 1)]),
    ]

    mu = (2,)
    assert list(Tableau._horizontal_strips(mu)) == [
        (tuple(), {(1, 1), (1, 2)}, []),
        ((1,), {(1, 2)}, [(1, 1)]),
        ((2,), set(), [(1, 2)]),
    ]

    mu = (1, 1)
    assert list(Tableau._horizontal_strips(mu)) == [
        ((1,), {(2, 1)}, []),
        ((1, 1), set(), [(2, 1)]),
    ]

    mu = (2, 1)
    assert list(Tableau._horizontal_strips(mu)) == [
        ((1,), {(1, 2), (2, 1)}, []),
        ((2,), {(2, 1)}, [(1, 2)]),
        ((1, 1), {(1, 2)}, [(2, 1)]),
        ((2, 1), set(), [(1, 2), (2, 1)]),
    ]


def test_shifted_horizontal_strips():
    mu = tuple()
    assert list(Tableau._shifted_horizontal_strips(mu)) == [
        (tuple(), set(), []),
    ]

    mu = (1,)
    assert list(Tableau._shifted_horizontal_strips(mu)) == [
        (tuple(), {(1, 1)}, []),
        ((1,), set(), [(1, 1)]),
    ]

    mu = (2,)
    assert list(Tableau._shifted_horizontal_strips(mu)) == [
        (tuple(), {(1, 1), (1, 2)}, []),
        ((1,), {(1, 2)}, [(1, 1)]),
        ((2,), set(), [(1, 2)]),
    ]

    mu = (2, 1)
    assert list(Tableau._shifted_horizontal_strips(mu)) == [
        ((2,), {(2, 2)}, []),
        ((2, 1), set(), [(2, 2)]),
    ]

    mu = (3, 1)
    assert list(Tableau._shifted_horizontal_strips(mu)) == [
        ((2,), {(2, 2), (1, 3)}, []),
        ((3,), {(2, 2)}, [(1, 3)]),
        ((2, 1), {(1, 3)}, [(2, 2)]),
        ((3, 1), set(), [(1, 3), (2, 2)]),
    ]


def test_shifted_vertical_strips():
    mu = tuple()
    assert list(Tableau._shifted_vertical_strips(mu)) == [
        (tuple(), set(), []),
    ]

    mu = (1,)
    assert list(Tableau._shifted_vertical_strips(mu)) == [
        (tuple(), {(1, 1)}, []),
        ((1,), set(), [(1, 1)]),
    ]

    mu = (2,)
    assert list(Tableau._shifted_vertical_strips(mu)) == [
        ((1,), {(1, 2)}, []),
        ((2,), set(), [(1, 2)]),
    ]

    mu = (2, 1)
    assert list(Tableau._shifted_vertical_strips(mu)) == [
        ((1,), {(1, 2), (2, 2)}, []),
        ((2,), {(2, 2)}, [(1, 2)]),
        ((2, 1), set(), [(2, 2)]),
    ]

    mu = (3, 1)
    print(list(Tableau._shifted_vertical_strips(mu)))
    assert list(Tableau._shifted_vertical_strips(mu)) == [
        ((2,), {(2, 2), (1, 3)}, []),
        ((3,), {(2, 2)}, [(1, 3)]),
        ((2, 1), {(1, 3)}, [(2, 2)]),
        ((3, 1), set(), [(1, 3), (2, 2)]),
    ]



def test_shifted_rpp_horizontal_strips():
    mu = ()
    assert Tableau._shifted_rpp_horizontal_strips(mu) == [
        ((), set())
    ]

    mu = (1,)
    assert Tableau._shifted_rpp_horizontal_strips(mu) == [
        ((), {(1, 1)})
    ]

    mu = (2, 1)
    assert {nu for nu, _ in Tableau._shifted_rpp_horizontal_strips(mu)} == {
        (2,),
        (1,),
        (),
    }

    mu = (3, 1)
    assert {nu for nu, _ in Tableau._shifted_rpp_horizontal_strips(mu)} == {
        (3,),
        (2,),
        (1,),
        (),
        (2, 1),
    }


def test_rpp_horizontal_strips():
    mu = ()
    assert Tableau._rpp_horizontal_strips(mu) == [
        ((), set())
    ]

    mu = (1,)
    assert Tableau._rpp_horizontal_strips(mu) == [
        ((), {(1, 1)}),
        ((1,), set()),
    ]

    mu = (1, 1)
    assert {nu for nu, _ in Tableau._rpp_horizontal_strips(mu)} == {
        (),
        (1, 1),
        (1,),
    }

    mu = (2, 1)
    assert {nu for nu, _ in Tableau._rpp_horizontal_strips(mu)} == {
        (2,),
        (1, 1),
        (1,),
        (),
        (2, 1),
    }

    mu = (3, 1)
    assert {nu for nu, _ in Tableau._rpp_horizontal_strips(mu)} == {
        (3,),
        (2,),
        (1,),
        (),
        (2, 1),
        (1, 1),
        (3, 1),
    }

    mu = (2, 1)
    assert {nu for nu, _ in Tableau._rpp_vertical_strips(mu)} == {
        (2,),
        (1, 1),
        (1,),
        (),
        (2, 1),
    }


def test_shifted_rpp_horizontal_strips():
    mu = ()
    assert Tableau._shifted_rpp_horizontal_strips(mu) == [
        ((), set())
    ]

    mu = (1,)
    assert Tableau._shifted_rpp_horizontal_strips(mu) == [
        ((), {(1, 1)}),
        ((1,), set())
    ]

    mu = (2, 1)
    assert {nu for nu, _ in Tableau._shifted_rpp_horizontal_strips(mu)} == {
        (2, 1),
        (2,),
        (1,),
        (),
    }

    mu = (3, 1)
    assert {nu for nu, _ in Tableau._shifted_rpp_horizontal_strips(mu)} == {
        (3, 1),
        (3,),
        (2,),
        (1,),
        (),
        (2, 1),
    }


def test_shifted_rpp_verticle_strips():
    mu = ()
    assert Tableau._shifted_rpp_vertical_strips(mu) == [
        ((), set())
    ]

    mu = (1,)
    assert Tableau._shifted_rpp_vertical_strips(mu) == [
        ((), {(1, 1)}),
        ((1,), set())
    ]

    mu = (2, 1)
    assert {nu for nu, _ in Tableau._shifted_rpp_vertical_strips(mu)} == {
        (2, 1),
        (2,),
        (1,),
        (),
    }

    mu = (3, 1)
    assert {nu for nu, _ in Tableau._shifted_rpp_vertical_strips(mu)} == {
        (3, 1),
        (3,),
        (2,),
        (1,),
        (),
        (2, 1),
    }


def test_semistandard_marked_rpp():
    mu = (1,)
    print(Tableau.semistandard_marked_rpp(mu, 1, True))
    print()
    assert Tableau.semistandard_marked_rpp(mu, 1, True) == {
        MarkedReversePlanePartition({(1, 1): 1}),
        MarkedReversePlanePartition({(1, 1): -1}),
    }
    print(Tableau.semistandard_marked_rpp(mu, 1, False))
    print()
    assert Tableau.semistandard_marked_rpp(mu, 1, False) == {
        MarkedReversePlanePartition({(1, 1): -1}),
    }


def test_semistandard():
    mu = ()
    assert Tableau.semistandard(mu, 0) == {Tableau()}
    assert Tableau.semistandard(mu, 1) == {Tableau()}
    assert Tableau.semistandard(mu, 2) == {Tableau()}

    mu = (1,)
    assert Tableau.semistandard(mu, 1) == {Tableau({(1, 1): 1})}
    assert Tableau.semistandard(mu, 2) == {
        Tableau({(1, 1): 1}),
        Tableau({(1, 1): 2}),
    }

    mu = (1, 1)
    assert Tableau.semistandard(mu, 2) == {
        Tableau({(1, 1): 1, (2, 1): 2})
    }

    mu = (2,)
    assert Tableau.semistandard(mu, 2) == {
        Tableau({(1, 1): 1, (1, 2): 1}),
        Tableau({(1, 1): 1, (1, 2): 2}),
        Tableau({(1, 1): 2, (1, 2): 2}),
    }

    mu = (1,)
    assert Tableau.semistandard(mu, 3) == {
        Tableau({(1, 1): 1}),
        Tableau({(1, 1): 2}),
        Tableau({(1, 1): 3}),
    }


def test_semistandard_setvalued():
    mu = ()
    assert Tableau.semistandard_setvalued(mu, 0) == {Tableau()}
    assert Tableau.semistandard_setvalued(mu, 1) == {Tableau()}
    assert Tableau.semistandard_setvalued(mu, 2) == {Tableau()}

    mu = (1,)
    assert Tableau.semistandard_setvalued(mu, 1) == {Tableau({(1, 1): 1})}
    assert Tableau.semistandard_setvalued(mu, 2) == {
        Tableau({(1, 1): 1}),
        Tableau({(1, 1): 2}),
        Tableau({(1, 1): (1, 2)})
    }

    mu = (1, 1)
    assert Tableau.semistandard_setvalued(mu, 2) == {
        Tableau({(1, 1): 1, (2, 1): 2})
    }

    mu = (2,)
    assert Tableau.semistandard_setvalued(mu, 2) == {
        Tableau({(1, 1): 1, (1, 2): 1}),
        Tableau({(1, 1): 1, (1, 2): 2}),
        Tableau({(1, 1): 2, (1, 2): 2}),
        Tableau({(1, 1): 1, (1, 2): (1, 2)}),
        Tableau({(1, 1): (1, 2), (1, 2): 2}),
    }


def test_semistandard_marked():
    mu = ()
    assert Tableau.semistandard_marked(mu, 0) == {Tableau()}
    assert Tableau.semistandard_marked(mu, 1) == {Tableau()}
    assert Tableau.semistandard_marked(mu, 2) == {Tableau()}

    mu = (1,)
    assert Tableau.semistandard_marked(mu, 1) == {
        Tableau({(1, 1): 1}),
        Tableau({(1, 1): -1}),
    }
    assert Tableau.semistandard_marked(mu, 2) == {
        Tableau({(1, 1): 1}),
        Tableau({(1, 1): -1}),
        Tableau({(1, 1): 2}),
        Tableau({(1, 1): -2}),
    }

    mu = (3, 1)
    print(Tableau.semistandard_marked(mu, 1))
    assert Tableau.semistandard_marked(mu, 1) == {
        Tableau({(1, 1): -1, (2, 1): -1, (1, 2): 1, (1, 3): 1}),
        Tableau({(1, 1): -1, (2, 1): 1, (1, 2): 1, (1, 3): 1}),
    }
    mu = (3, 2)
    assert Tableau.semistandard_marked(mu, 1) == set()


def test_semistandard_marked_setvalued():
    mu = ()
    assert Tableau.semistandard_marked_setvalued(mu, 0) == {Tableau()}
    assert Tableau.semistandard_marked_setvalued(mu, 1) == {Tableau()}
    assert Tableau.semistandard_marked_setvalued(mu, 2) == {Tableau()}

    mu = (1,)
    assert Tableau.semistandard_marked_setvalued(mu, 1) == {
        Tableau({(1, 1): 1}),
        Tableau({(1, 1): -1}),
        Tableau({(1, 1): (-1, 1)})
    }
    assert Tableau.semistandard_marked_setvalued(mu, 2) == {
        Tableau({(1, 1): 1}),
        Tableau({(1, 1): -1}),
        Tableau({(1, 1): (-1, 1)}),
        Tableau({(1, 1): 2}),
        Tableau({(1, 1): -2}),
        Tableau({(1, 1): (-2, 2)}),
        Tableau({(1, 1): (1, 2)}),
        Tableau({(1, 1): (1, -2)}),
        Tableau({(1, 1): (1, -2, 2)}),
        Tableau({(1, 1): (-1, 2)}),
        Tableau({(1, 1): (-1, -2)}),
        Tableau({(1, 1): (-1, -2, 2)}),
        Tableau({(1, 1): (-1, 1, 2)}),
        Tableau({(1, 1): (-1, 1, -2)}),
        Tableau({(1, 1): (-1, 1, -2, 2)})
    }


def test_semistandard_shifted_marked():
    mu = ()
    assert Tableau.semistandard_shifted_marked(mu, 0) == {Tableau()}
    assert Tableau.semistandard_shifted_marked(mu, 1) == {Tableau()}
    assert Tableau.semistandard_shifted_marked(mu, 2) == {Tableau()}

    mu = (1,)

    assert Tableau.semistandard_shifted_marked(mu, 1) == {
        Tableau({(1, 1): 1}),
        Tableau({(1, 1): -1}),
    }
    assert Tableau.semistandard_shifted_marked(mu, 1, False) == {
        Tableau({(1, 1): 1}),
    }

    assert Tableau.semistandard_shifted_marked(mu, 2) == {
        Tableau({(1, 1): 1}),
        Tableau({(1, 1): -1}),
        Tableau({(1, 1): 2}),
        Tableau({(1, 1): -2}),
    }
    assert Tableau.semistandard_shifted_marked(mu, 2, False) == {
        Tableau({(1, 1): 1}),
        Tableau({(1, 1): 2}),
    }

    mu = (3, 1)
    assert Tableau.semistandard_shifted_marked(mu, 1) == set()
    assert Tableau.semistandard_shifted_marked(mu, 1, False) == set()

    assert Tableau.semistandard_shifted_marked(mu, 2) == {
        Tableau({(1, 1): 1, (2, 2): 2, (1, 2): 1, (1, 3): 1}),
        Tableau({(1, 1): 1, (2, 2): 2, (1, 2): 1, (1, 3): 2}),
        Tableau({(1, 1): 1, (2, 2): 2, (1, 2): 1, (1, 3): -2}),
        Tableau({(1, 1): 1, (2, 2): 2, (1, 2): -2, (1, 3): 2}),
        Tableau({(1, 1): 1, (2, 2): -2, (1, 2): 1, (1, 3): 1}),
        Tableau({(1, 1): 1, (2, 2): -2, (1, 2): 1, (1, 3): 2}),
        Tableau({(1, 1): 1, (2, 2): -2, (1, 2): 1, (1, 3): -2}),
        Tableau({(1, 1): 1, (2, 2): -2, (1, 2): -2, (1, 3): 2}),
        Tableau({(1, 1): -1, (2, 2): 2, (1, 2): 1, (1, 3): 1}),
        Tableau({(1, 1): -1, (2, 2): 2, (1, 2): 1, (1, 3): 2}),
        Tableau({(1, 1): -1, (2, 2): 2, (1, 2): 1, (1, 3): -2}),
        Tableau({(1, 1): -1, (2, 2): 2, (1, 2): -2, (1, 3): 2}),
        Tableau({(1, 1): -1, (2, 2): -2, (1, 2): 1, (1, 3): 1}),
        Tableau({(1, 1): -1, (2, 2): -2, (1, 2): 1, (1, 3): 2}),
        Tableau({(1, 1): -1, (2, 2): -2, (1, 2): 1, (1, 3): -2}),
        Tableau({(1, 1): -1, (2, 2): -2, (1, 2): -2, (1, 3): 2}),
    }
    assert Tableau.semistandard_shifted_marked(mu, 2, False) == {
        Tableau({(1, 1): 1, (2, 2): 2, (1, 2): 1, (1, 3): 1}),
        Tableau({(1, 1): 1, (2, 2): 2, (1, 2): 1, (1, 3): 2}),
        Tableau({(1, 1): 1, (2, 2): 2, (1, 2): 1, (1, 3): -2}),
        Tableau({(1, 1): 1, (2, 2): 2, (1, 2): -2, (1, 3): 2}),
    }

    mu = (3, 2)
    assert Tableau.semistandard_shifted_marked(mu, 1) == set()
    assert Tableau.semistandard_shifted_marked(mu, 1, False) == set()


# def test_semistandard_shifted_marked_setvalued():
#     mu = ()
#     assert Tableau.semistandard_shifted_marked_setvalued(mu, 0) == {Tableau()}
#     assert Tableau.semistandard_shifted_marked_setvalued(mu, 1) == {Tableau()}
#     assert Tableau.semistandard_shifted_marked_setvalued(mu, 2) == {Tableau()}

#     mu = (1,)
#     assert Tableau.semistandard_shifted_marked_setvalued(mu, 1) == {
#         Tableau({(1, 1): 1}),
#         Tableau({(1, 1): -1}),
#         Tableau({(1, 1): (-1, 1)})
#     }
#     assert Tableau.semistandard_shifted_marked_setvalued(mu, 2) == {
#         Tableau({(1, 1): 1}),
#         Tableau({(1, 1): -1}),
#         Tableau({(1, 1): (-1, 1)}),
#         Tableau({(1, 1): 2}),
#         Tableau({(1, 1): -2}),
#         Tableau({(1, 1): (-2, 2)}),
#         Tableau({(1, 1): (1, 2)}),
#         Tableau({(1, 1): (1, -2)}),
#         Tableau({(1, 1): (1, -2, 2)}),
#         Tableau({(1, 1): (-1, 2)}),
#         Tableau({(1, 1): (-1, -2)}),
#         Tableau({(1, 1): (-1, -2, 2)}),
#         Tableau({(1, 1): (-1, 1, 2)}),
#         Tableau({(1, 1): (-1, 1, -2)}),
#         Tableau({(1, 1): (-1, 1, -2, 2)})
#     }
