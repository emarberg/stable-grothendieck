from insertion import InsertionAlgorithm
from tableaux import Tableau

unprime = InsertionAlgorithm.unprime_diagonal


def test_unprime(n=4, mu=(4, 1)):
    u = list(Tableau.semistandard_shifted_marked_setvalued(n, mu, diagonal_primes=True))
    mapping = {}
    for t in u:
        tt = unprime(t)
        if tt in mapping:
            print('target:')
            print(tt)
            print('t =', t.boxes)
            nu = tt.shape()
            assert t.shape() == mu
            assert len(nu) == len(mu)
            assert all(nu[i] in [mu[i], mu[i] + 1] for i in range(len(mu)))
            assert t.is_semistandard(diagonal_primes=True)
            assert tt.is_semistandard(diagonal_primes=True)
            mapping[tt].append(t)
            for x in mapping[tt]:
                print(x)
                print(x.boxes)
        else:
            mapping[tt] = [t]
    assert not any(len(v) > 1 for v in mapping.values())


def test_instance():
    tabs = [
        Tableau({(1, 2): (1,), (1, 3): (1,), (1, 1): (1,), (2, 3): (-3, 2), (1, 4): (-3, -2, 2), (2, 2): (-2, 2), (1, 5): (3,), (1, 6): (3,), (2, 4): (3,), (3, 3): (-3,)}),
        Tableau({(1, 2): (1,), (1, 3): (1,), (1, 1): (1,), (2, 3): (2, 3), (1, 4): (2, 3), (2, 2): (-2, 2), (1, 5): (3,), (1, 6): (3,), (2, 4): (4,), (3, 3): (4,)}),
        Tableau({(1, 1): (-2, -1, 1), (1, 2): (2,), (2, 2): (3,), (1, 3): (-3,)}),
    ]
    for t in tabs:
        tt = unprime(t)
        print(t)
        print(tt)

        mu = t.shape()
        nu = tt.shape()
        assert len(nu) == len(mu)
        assert all(nu[i] in [mu[i], mu[i] + 1] for i in range(len(mu)))
        assert t.is_semistandard(diagonal_primes=True)
        assert tt.is_semistandard(diagonal_primes=True)
