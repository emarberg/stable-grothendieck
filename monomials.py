from vectors import Vector
from tableaux import Tableau
from collections import defaultdict
import itertools


class Monomial:

    def __init__(self, *variables):
        assert all(type(v) == str and len(v) > 0 for v in variables)
        self.dictionary = {v: 1 for v in variables}

    def swap(self, var1, var2):
        ans = Monomial()
        ans.dictionary = {k: v for k, v in self.dictionary.items() if k not in [var1, var2]}
        if var1 in self:
            ans.dictionary[var2] = self[var1]
        if var2 in self:
            ans.dictionary[var1] = self[var2]
        return ans

    def degree(self):
        ans = 0
        for v in self.dictionary:
            ans += self.dictionary[v]
        return ans

    @classmethod
    def from_dictionary(cls, dictionary):
        m = Monomial()
        m.dictionary = dictionary
        return m

    @classmethod
    def from_tableau(cls, tableau):
        assert type(tableau) == Tableau
        ans = Monomial()
        for i, j, values in tableau:
            for v in values:
                ans *= Monomial('x%i' % abs(v))
        return ans

    @classmethod
    def from_partition(cls, mu):
        ans = Monomial()
        for i in range(len(mu)):
            ans *= Monomial('x' + str(i + 1)) ** mu[i]
        return ans

    def __bool__(self):
        return bool(self.dictionary)

    def __getitem__(self, i):
        return self.dictionary.get(i, 0)

    def __iter__(self):
        for v in self.dictionary:
            yield v

    def __eq__(self, other):
        assert type(other) == Monomial
        return self.dictionary == other.dictionary

    def __hash__(self):
        return hash(tuple(sorted(self.dictionary.items())))

    def __rmul__(self, other):
        return self * other

    def polynomial(self, coefficient=1):
        return Polynomial({self: coefficient}, printer=lambda x: str(x) if x else '')

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
            return Monomial.from_dictionary({
                v: self[v] + other[v]
                for v in set(self.dictionary) | set(other.dictionary)
                if self[v] + other[v] != 0
            })
        elif type(other) == Polynomial:
            return self.polynomial() * other

    def __pow__(self, other):
        assert type(other) == int
        ans = Monomial()
        if other != 0:
            ans.dictionary = {v: self.dictionary[v] * other for v in self}
        return ans

    def __truediv__(self, other):
        assert (type(other) == int and other in [1, -1]) or type(other) == Monomial
        other = other ** -1 if type(other) == Monomial else other
        return self * other

    def __radd__(self, other):
        return self + other

    def __add__(self, other):
        assert type(other) in [int, Monomial, Polynomial]
        if type(other) == int:
            return self.polynomial() + Monomial().polynomial(other)
        elif type(other) == Monomial:
            return self.polynomial() + other.polynomial()
        elif type(other) == Polynomial:
            return self.polynomial() + other

    def __repr__(self):
        if self:
            return ' * '.join([
                k if v == 1 else (k + '**%i' % v)
                for k, v in sorted(self.dictionary.items())
            ])
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

    def symmetrize(self, n):
        var = ['x%i' % i for i in range(1, n + 1)]
        ans = Polynomial()
        for mon, coeff in self.items():
            for permutation in set(itertools.permutations([mon[x] for x in var])):
                ans += Monomial.from_dictionary({var[i]: a for i, a in enumerate(permutation) if a}) * coeff
        return ans

    def lowest_degree_terms(self):
        degrees = {m.degree() for m in self}
        m = min(degrees) if degrees else 0
        return Polynomial({
            mon: val for mon, val in self.items()
            if mon.degree() == m
        })

    @classmethod
    def schur_s(cls, mu, n,):
        dictionary = defaultdict(int)
        for partition, count in Tableau.count_semistandard_marked(mu, n).items():
            dictionary[Monomial.from_partition(partition)] += count
        return Polynomial(dictionary)

    @classmethod
    def stable_grothendieck_s(cls, mu, n):
        dictionary = defaultdict(int)
        for partition, count in Tableau.count_semistandard_marked_setvalued(mu, n).items():
            dictionary[Monomial.from_partition(partition)] += count
        return Polynomial(dictionary)

    @classmethod
    def schur(cls, mu, n):
        dictionary = defaultdict(int)
        for partition, count in Tableau.count_semistandard(mu, n).items():
            dictionary[Monomial.from_partition(partition)] += count
        return Polynomial(dictionary)

    @classmethod
    def stable_grothendieck(cls, mu, n):
        dictionary = defaultdict(int)
        for partition, count in Tableau.count_semistandard_setvalued(mu, n).items():
            dictionary[Monomial.from_partition(partition)] += count
        return Polynomial(dictionary)

    @classmethod
    def slow_schur_s(cls, mu, n,):
        dictionary = defaultdict(int)
        for tab in Tableau.semistandard_marked(mu, n):
            dictionary[Monomial.from_tableau(tab)] += 1
        return Polynomial(dictionary)

    @classmethod
    def slow_stable_grothendieck_s(cls, mu, n):
        dictionary = defaultdict(int)
        for tab in Tableau.semistandard_marked_setvalued(mu, n):
            dictionary[Monomial.from_tableau(tab)] += 1
        return Polynomial(dictionary)

    @classmethod
    def slow_schur(cls, mu, n):
        dictionary = defaultdict(int)
        for tab in Tableau.semistandard(mu, n):
            dictionary[Monomial.from_tableau(tab)] += 1
        return Polynomial(dictionary)

    @classmethod
    def slow_stable_grothendieck(cls, mu, n):
        dictionary = defaultdict(int)
        for tab in Tableau.semistandard_setvalued(mu, n):
            dictionary[Monomial.from_tableau(tab)] += 1
        return Polynomial(dictionary)
