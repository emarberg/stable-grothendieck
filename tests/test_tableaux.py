from tableaux import Tableau


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
