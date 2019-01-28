from vectors import Vector
from tableaux import Tableau, Partition
from collections import defaultdict
import itertools


SCHUR_CACHE = {}
STABLE_GROTHENDIECK_CACHE = {}

SCHUR_S_CACHE = {}
STABLE_GROTHENDIECK_S_CACHE = {}

MONOMIAL_PRODUCT_CACHE = {}


class Monomial:

    def __init__(self, n, mu=()):
        assert len(mu) <= n
        self.n = n
        self.mu = Partition.sort(mu, trim=True)

    def index(self):
        return self.mu

    def order(self):
        return self.n

    def degree(self):
        return sum(self.index())

    def __lt__(self, other):
        assert type(other) == type(self)
        assert self.order() == other.order()
        return self.index() < other.index()

    def __bool__(self):
        return bool(self.mu)

    def __eq__(self, other):
        assert type(other) in [Monomial, int, Polynomial]
        if type(other) == Monomial:
            return self.order() == other.order() and self.index() == other.index()
        elif type(other) in [int, Polynomial]:
            return (self - other).is_zero()

    def __hash__(self):
        return hash((self.n, self.mu))

    def polynomial(self, coefficient=1):
        return Polynomial({self: coefficient}, printer=lambda x: str(x) if x else '')

    def __pow__(self, other):
        assert type(other) == int and other >= 0
        if other == 0:
            return Monomial(self.order())
        elif other == 1:
            return self
        else:
            f = self ** (other // 2)
            return (f * f) if other % 2 == 0 else (f * f * self)

    def __radd__(self, other):
        return self + other

    def __add__(self, other):
        assert type(other) in [int, Monomial, Polynomial]
        if type(other) == int:
            return self.polynomial() + Monomial(self.order()).polynomial(other)
        elif type(other) == Monomial:
            return self.polynomial() + other.polynomial()
        elif type(other) == Polynomial:
            return self.polynomial() + other

    def __rsub__(self, other):
        return self.polynomial() * -1 + other

    def __sub__(self, other):
        return self.polynomial() + (-1 * other)

    def __rmul__(self, other):
        return self * other

    def __mul__(self, other):
        assert type(other) in [int, Monomial, Polynomial]
        if type(other) == int:
            if other == 1:
                return self
            elif other == 0:
                return Polynomial()
            else:
                return self.polynomial(other)
        elif type(other) == Monomial:
            return Polynomial({
                Monomial(self.order(), mu): coeff
                for mu, coeff in self._multiply(self, other).items()
            })
        elif type(other) == Polynomial:
            return self.polynomial() * other

    @classmethod
    def _destandardize(cls, mu, n):
        if mu == ():
            yield n * (0,)
            return

        if n == 0:
            return

        for alpha in cls._destandardize(mu, n - 1):
            yield (0,) + alpha

        for i, part in enumerate(mu):
            if i == 0 or mu[i] != mu[i - 1]:
                nu = mu[:i] + mu[i + 1:]
                for alpha in cls._destandardize(nu, n - 1):
                    yield (part,) + alpha

    @classmethod
    def _multiply(cls, f, g):
        assert type(f) == type(g) == Monomial
        assert f.order() == g.order()
        mu = f.mu
        nu = g.mu
        n = min(f.order(), len(mu) + len(nu))
        if (mu, nu, n) not in MONOMIAL_PRODUCT_CACHE:
            ans = defaultdict(int)
            for alpha in cls._destandardize(mu, n):
                for beta in cls._destandardize(nu, n):
                    gamma = tuple(alpha[i] + beta[i] for i in range(n))
                    if Partition.is_partition(gamma):
                        while gamma and gamma[-1] == 0:
                            gamma = gamma[:-1]
                        ans[gamma] += 1
            MONOMIAL_PRODUCT_CACHE[(mu, nu, n)] = ans
        return MONOMIAL_PRODUCT_CACHE[(mu, nu, n)]

    def __repr__(self):
        if self:
            return 'm[%s;%s]' % (','.join(map(str, self.mu)), self.n)
        else:
            return '1'


class Polynomial(Vector):

    def variables(self):
        ans = set()
        for monomial in self:
            ans |= set(monomial.dictionary)
        return sorted(ans)

    def is_symmetric(self, n):
        var = ['x%i' % i for i in range(1, n + 1)]
        for monomial in self:
            for i in range(len(var) - 1):
                other = monomial.swap(var[i], var[i + 1])
                if self[other] != self[monomial]:
                    return False
        return True

    def lowest_degree_terms(self):
        degrees = {m.degree() for m in self}
        m = min(degrees) if degrees else 0
        return Polynomial({
            mon: val for mon, val in self.items()
            if mon.degree() == m
        })

    @classmethod
    def _function(cls, n, tableaux):
        dictionary = defaultdict(int)
        for partition, count in tableaux.items():
            dictionary[Monomial(n, partition)] += count
        return Polynomial(dictionary)

    @classmethod
    def schur_s(cls, mu, n,):
        cache = SCHUR_S_CACHE
        if (mu, n) not in cache:
            tableaux = Tableau.count_semistandard_marked(mu, n)
            cache[(mu, n)] = cls._function(n, tableaux)
        return cache[(mu, n)]

    @classmethod
    def stable_grothendieck_s(cls, mu, n):
        cache = STABLE_GROTHENDIECK_S_CACHE
        if (mu, n) not in cache:
            tableaux = Tableau.count_semistandard_marked_setvalued(mu, n)
            cache[(mu, n)] = cls._function(n, tableaux)
        return cache[(mu, n)]

    @classmethod
    def schur(cls, mu, n):
        cache = SCHUR_CACHE
        if (mu, n) not in cache:
            tableaux = Tableau.count_semistandard(mu, n)
            cache[(mu, n)] = cls._function(n, tableaux)
        return cache[(mu, n)]

    @classmethod
    def stable_grothendieck(cls, mu, n):
        cache = STABLE_GROTHENDIECK_CACHE
        if (mu, n) not in cache:
            tableaux = Tableau.count_semistandard_setvalued(mu, n)
            cache[(mu, n)] = cls._function(n, tableaux)
        return cache[(mu, n)]

    @classmethod
    def schur_expansion(cls, f, n):
        if f:
            t = max(f)
            c = f[t]
            mu = t.index()
            ans = cls.schur_expansion(f - c * cls.schur(mu, n), n)
            return ans + Vector({mu: c})
        else:
            return Vector()

    @classmethod
    def _slow_function(cls, n, tableaux):
        dictionary = defaultdict(int)
        for tab in tableaux:
            dictionary[tab.weight(n)] += 1
        print(dictionary)
        assert all(dictionary[Partition.sort(alpha)] == dictionary[alpha] for alpha in dictionary)
        return Polynomial({
            Monomial(n, alpha): coeff
            for alpha, coeff in dictionary.items()
            if Partition.is_partition(alpha)
        })

    @classmethod
    def slow_schur_s(cls, mu, n):
        return cls._slow_function(n, Tableau.semistandard_marked(mu, n))

    @classmethod
    def slow_stable_grothendieck_s(cls, mu, n):
        return cls._slow_function(n, Tableau.semistandard_marked_setvalued(mu, n))

    @classmethod
    def slow_schur(cls, mu, n):
        return cls._slow_function(n, Tableau.semistandard(mu, n))

    @classmethod
    def slow_stable_grothendieck(cls, mu, n):
        return cls._slow_function(n, Tableau.semistandard_setvalued(mu, n))
