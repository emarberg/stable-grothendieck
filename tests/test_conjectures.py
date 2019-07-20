from utils import GQ, GS
from symmetric import SymmetricPolynomial
from tableaux import Partition
from vectors import Vector
import pytest


@pytest.mark.slow
def test_GQ_to_GP_expansion():
    def is_binary_power(i):
        return len(list(filter(lambda x: x == '1', bin(i)))) == 1

    def sgn(mu, nu):
        boxes = sorted({(i + 1, i + j) for i in range(len(mu)) for j in range(mu[i] + 1, nu[i] + 1)})
        sgn = 1
        for i in range(len(boxes)):
            _, j = boxes[i]
            sgn *= -1 if i > 0 and boxes[i - 1][1] == j else 1
        return sgn

    def zero_one_tuples(n):
        if n == 0:
            yield ()
            return
        for i in range(2**n):
            ans = tuple(map(int, list(bin(i))[2:]))
            yield (n - len(ans)) * (0,) + ans

    for mu in Partition.all(25, strict=True):
        print('mu =', mu)
        print()
        print(Partition.printable(mu, shifted=True))
        print()
        n = len(mu)
        q = GQ(mu, n)
        expansion = SymmetricPolynomial.GP_expansion(q)
        normalized = Vector({
            tuple(nu[i] - mu[i] for i in range(len(mu))):
            c * sgn(mu, nu) * (-1)**(sum(nu) - sum(mu)) / 2**(len(mu) - sum(nu) + sum(mu))
            for nu, c in expansion.items()
        })
        unsigned = all(c > 0 for c in normalized.values())
        print('  mu =', mu, 'n =', n)
        print('  expansion =', expansion)
        print('  normalized expansion =', normalized)
        assert all(len(nu) == 0 or max(nu) <= 1 for nu in normalized)
        assert all(len(nu) == len(mu) for nu in expansion)
        assert all(Partition.contains(nu, mu) for nu in expansion)
        assert all(c % 2**(len(mu) - sum(nu) + sum(mu)) == 0 for nu, c in expansion.items())
        assert all(is_binary_power(abs(c)) for c in expansion.values())
        assert unsigned
        expected = {
            tuple(mu[i] + a[i] for i in range(len(a)))
            for a in zero_one_tuples(len(mu))
            if all(mu[i - 1] + a[i - 1] > mu[i] + a[i] for i in range(1, len(a)))
        }
        print('  expected =', expected)
        assert set(expansion) == expected
        print()
        print()


@pytest.mark.slow
def test_dewitt_conjecture(n=5):

    def mu(m, k):
        return k * (m,)

    def nu(m, k):
        ans = [m + k - 1]
        while ans[-1] - 2 >= abs(m - k) + 1:
            ans += [ans[-1] - 2]
        return tuple(ans)

    for m in range(1, 10):
        for k in range(1, 10):
            a = mu(m, k)
            b = nu(m, k)
            print(' n =', n)
            print('mu =', a)
            print('nu =', b)
            print()

            if n < len(a) and n < len(b):
                continue

            gs = GS(a, n)
            gq = GQ(b, n)
            print('GS =', gs)
            print()
            print('GQ =', gq)
            print()
            print()
            assert gs == gq
