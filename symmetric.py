from vectors import Vector
from tableaux import Tableau, Partition
from polynomials import Polynomial
from polynomials import beta as BETA # noqa
from cached import cached_value
from collections import defaultdict
import itertools


SCHUR_CACHE = {}
SCHUR_S_CACHE = {}
SCHUR_Q_CACHE = {}
SCHUR_P_CACHE = {}

STABLE_GROTHENDIECK_CACHE = {}
STABLE_GROTHENDIECK_S_CACHE = {}
STABLE_GROTHENDIECK_Q_CACHE = {}
STABLE_GROTHENDIECK_P_CACHE = {}

DUAL_STABLE_GROTHENDIECK_CACHE = {}
DUAL_STABLE_GROTHENDIECK_P_CACHE = {}
DUAL_STABLE_GROTHENDIECK_Q_CACHE = {}

MONOMIAL_PRODUCT_CACHE = {}


class SymmetricMonomial:

    def __init__(self, n, mu=()):
        assert len(mu) <= n
        self.n = n
        self.mu = Partition.sort(mu, trim=True)

    def serialize(self):
        return (self.mu, self.n)

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
        assert type(other) in [SymmetricMonomial, int, SymmetricPolynomial]
        if type(other) == SymmetricMonomial:
            return self.order() == other.order() and self.index() == other.index()
        elif type(other) in [int, SymmetricPolynomial]:
            return (self - other).is_zero()

    def __hash__(self):
        return hash((self.n, self.mu))

    def polynomial(self, variable='x'):
        x = Polynomial.x if variable == 'x' else Polynomial.y
        ans = Polynomial()
        for alpha in set(itertools.permutations(self.mu)):
            k = len(alpha)
            for index in itertools.combinations([i + 1 for i in range(self.n)], k):
                term = Polynomial.one()
                for i in range(k):
                    term *= x(index[i])**alpha[i]
                ans += term
        return ans

    def symmetric_polynomial(self, coefficient=1):
        return SymmetricPolynomial({self: coefficient}, printer=lambda x: str(x) if x else '')

    def __pow__(self, other):
        assert type(other) == int and other >= 0
        if other == 0:
            return SymmetricMonomial(self.order())
        elif other == 1:
            return self
        else:
            f = self ** (other // 2)
            return (f * f) if other % 2 == 0 else (f * f * self)

    def __radd__(self, other):
        return self + other

    def __add__(self, other):
        assert type(other) in [int, SymmetricMonomial, SymmetricPolynomial]
        if type(other) == int:
            return self.symmetric_polynomial() + SymmetricMonomial(self.order()).symmetric_polynomial(other)
        elif type(other) == SymmetricMonomial:
            return self.symmetric_polynomial() + other.symmetric_polynomial()
        elif type(other) == SymmetricPolynomial:
            return self.symmetric_polynomial() + other

    def __rsub__(self, other):
        return self.symmetric_polynomial() * -1 + other

    def __sub__(self, other):
        return self.symmetric_polynomial() + (-1 * other)

    def __rmul__(self, other):
        return self * other

    def __mul__(self, other):
        assert type(other) in [int, Polynomial, SymmetricMonomial, SymmetricPolynomial]
        if type(other) == int:
            if other == 1:
                return self
            elif other == 0:
                return SymmetricPolynomial()
            else:
                return self.symmetric_polynomial(other)
        elif type(other) == SymmetricMonomial:
            return SymmetricPolynomial({
                SymmetricMonomial(self.order(), mu): coeff
                for mu, coeff in self._multiply(self, other).items()
            })
        elif type(other) in [Polynomial, SymmetricPolynomial]:
            return self.symmetric_polynomial() * other

    @classmethod
    def _destandardize(cls, n, mu):
        if mu == ():
            yield n * (0,)
            return

        if n == 0:
            return

        for alpha in cls._destandardize(n - 1, mu):
            yield (0,) + alpha

        for i, part in enumerate(mu):
            if i == 0 or mu[i] != mu[i - 1]:
                nu = mu[:i] + mu[i + 1:]
                for alpha in cls._destandardize(n - 1, nu):
                    yield (part,) + alpha

    @classmethod
    def _multiply(cls, f, g):
        assert type(f) == type(g) == SymmetricMonomial
        assert f.order() == g.order()
        mu = f.mu
        nu = g.mu
        n = min(f.order(), len(mu) + len(nu))
        if (mu, nu, n) not in MONOMIAL_PRODUCT_CACHE:
            ans = defaultdict(int)
            for alpha in cls._destandardize(n, mu):
                for blpha in cls._destandardize(n, nu):
                    gamma = tuple(alpha[i] + blpha[i] for i in range(n))
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


class SymmetricPolynomial(Vector):
    def __init__(self, dictionary={}, printer=None, multiplier=None):
        self.dictionary = {key: value for key, value in dictionary.items() if value}
        self.printer = printer

        def m(a, b):
            for mu, coeff in SymmetricMonomial._multiply(a, b).items():
                yield SymmetricMonomial(a.order(), mu), coeff

        self.multiplier = m

    def __eq__(self, other):
        if type(other) == int:
            if other == 0:
                return len(self.dictionary) == 0
            if len(self.dictionary) != 1:
                return False
            mon, coeff = list(self.dictionary.items())[0]
            return len(mon.mu) == 0 and coeff == other
        return len((self - other).dictionary) == 0

    def serialize(self):
        return tuple(sorted([(key.serialize(), value) for key, value in self.items()]))

    def __hash__(self):
        return hash(self.serialize())

    def __repr__(self):
        printer = self.printer or repr
        sorted_items = sorted(
            [(key, value) for key, value in self.items()],
            key=lambda m: (sum(m[0].mu), m[0].mu)
        )
        sorted_items = [(printer(key), value) for key, value in sorted_items]
        return self._print_sorted(sorted_items)

    def polynomial(self, variable='x'):
        ans = Polynomial()
        for mon, coeff in self.items():
            ans += mon.polynomial(variable) * coeff
        return ans

    def truncate(self, degree_bound):
        return SymmetricPolynomial({
            mon: coeff
            for mon, coeff in self.items()
            if mon.degree() <= degree_bound
        })

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
        return SymmetricPolynomial({
            mon: val for mon, val in self.items()
            if mon.degree() == m
        })

    def highest_degree_terms(self):
        degrees = {m.degree() for m in self}
        m = max(degrees) if degrees else 0
        return SymmetricPolynomial({
            mon: val for mon, val in self.items()
            if mon.degree() == m
        })

    @classmethod
    def _vectorize(cls, n, tableaux, signs=None, degree_bound=None):
        dictionary = defaultdict(int)
        for partition, count in tableaux.items():
            count *= signs**sum(partition) if signs else 1
            dictionary[SymmetricMonomial(n, partition)] += count
        if degree_bound is not None:
            dictionary = {m: dictionary[m] for m in dictionary if sum(m.mu) <= degree_bound}
        return SymmetricPolynomial(dictionary)

    @classmethod
    def schur(cls, num_variables, mu, nu=()):  # noqa
        return cls._schur(num_variables, mu, nu)

    @cached_value(SCHUR_CACHE)
    def _schur(cls, num_variables, mu, nu):  # noqa
        tableaux = Tableau.count_semistandard(num_variables, mu, nu)
        return cls._vectorize(num_variables, tableaux)

    @classmethod
    def schur_s(cls, num_variables, mu, nu=()):  # noqa
        return cls._schur_s(num_variables, mu, nu)

    @cached_value(SCHUR_S_CACHE)
    def _schur_s(cls, num_variables, mu, nu):  # noqa
        tableaux = Tableau.count_semistandard_marked(num_variables, mu, nu)
        return cls._vectorize(num_variables, tableaux)

    @classmethod
    def schur_q(cls, num_variables, mu, nu=()):  # noqa
        return cls._schur_q(num_variables, mu, nu)

    @cached_value(SCHUR_Q_CACHE)
    def _schur_q(cls, num_variables, mu, nu):  # noqa
        tableaux = Tableau.count_semistandard_shifted_marked(num_variables, mu, nu)
        return cls._vectorize(num_variables, tableaux)

    @classmethod
    def schur_p(cls, num_variables, mu, nu=()):  # noqa
        return cls._schur_p(num_variables, mu, nu)

    @cached_value(SCHUR_P_CACHE)
    def _schur_p(cls, num_variables, mu, nu):  # noqa
        tableaux = Tableau.count_semistandard_shifted_marked(num_variables, mu, nu, diagonal_primes=False)
        return cls._vectorize(num_variables, tableaux)

    @classmethod
    def stable_grothendieck_s(cls, num_variables, mu, nu=(), degree_bound=None):  # noqa
        return cls._stable_grothendieck_s(num_variables, mu, nu, degree_bound)

    @cached_value(STABLE_GROTHENDIECK_S_CACHE)
    def _stable_grothendieck_s(cls, num_variables, mu, nu, degree_bound):  # noqa
        tableaux = Tableau.count_semistandard_marked_setvalued(num_variables, mu, nu)
        return BETA**(sum(nu) - sum(mu)) * cls._vectorize(num_variables, tableaux, BETA, degree_bound)

    @classmethod
    def stable_grothendieck_q_doublebar(cls, num_variables, mu, nu=(), degree_bound=None):  # noqa
        ans = SymmetricPolynomial()
        if Partition.contains(mu, nu):
            for x in Partition.remove_shifted_inner_corners(nu):
                ans += BETA**(sum(nu) - sum(x)) * cls._stable_grothendieck_q(num_variables, mu, x, degree_bound)
        return ans

    @classmethod
    def stable_grothendieck_q(cls, num_variables, mu, nu=(), degree_bound=None):  # noqa
        return cls._stable_grothendieck_q(num_variables, mu, nu, degree_bound)

    @cached_value(STABLE_GROTHENDIECK_Q_CACHE)
    def _stable_grothendieck_q(cls, num_variables, mu, nu, degree_bound):  # noqa
        tableaux = Tableau.count_semistandard_shifted_marked_setvalued(num_variables, mu, nu)
        return BETA**(sum(nu) - sum(mu)) * cls._vectorize(num_variables, tableaux, BETA, degree_bound)

    @classmethod
    def stable_grothendieck_p_doublebar(cls, num_variables, mu, nu=(), degree_bound=None):  # noqa
        ans = SymmetricPolynomial()
        if Partition.contains(mu, nu):
            for x in Partition.remove_shifted_inner_corners(nu):
                ans += BETA**(sum(nu) - sum(x)) * cls._stable_grothendieck_p(num_variables, mu, x, degree_bound)
        return ans

    @classmethod
    def stable_grothendieck_p(cls, num_variables, mu, nu=(), degree_bound=None):  # noqa
        return cls._stable_grothendieck_p(num_variables, mu, nu, degree_bound)

    @cached_value(STABLE_GROTHENDIECK_P_CACHE)
    def _stable_grothendieck_p(cls, num_variables, mu, nu, degree_bound):  # noqa
        tableaux = Tableau.count_semistandard_shifted_marked_setvalued(num_variables, mu, nu, diagonal_primes=False)
        return BETA**(sum(nu) - sum(mu)) * cls._vectorize(num_variables, tableaux, BETA, degree_bound)

    @classmethod
    def stable_grothendieck_doublebar(cls, num_variables, mu, nu=(), degree_bound=None):  # noqa
        ans = SymmetricPolynomial()
        if Partition.contains(mu, nu):
            for x in Partition.remove_inner_corners(nu):
                ans += BETA**(sum(nu) - sum(x)) * cls._stable_grothendieck(num_variables, mu, x, degree_bound)
        return ans

    @classmethod
    def stable_grothendieck(cls, num_variables, mu, nu=(), degree_bound=None):  # noqa
        return cls._stable_grothendieck(num_variables, mu, nu, degree_bound)

    @cached_value(STABLE_GROTHENDIECK_CACHE)
    def _stable_grothendieck(cls, num_variables, mu, nu, degree_bound):  # noqa
        tableaux = Tableau.count_semistandard_setvalued(num_variables, mu, nu)
        return BETA**(sum(nu) - sum(mu)) * cls._vectorize(num_variables, tableaux, BETA, degree_bound)

    @classmethod
    def mp_dual_stable_grothendieck(cls, num_variables, mu, nu=()):  # noqa
        ans = cls.dual_stable_grothendieck(num_variables, mu, nu)
        dictionary = {k: v.substitute(0, -BETA) for (k, v) in ans.dictionary.items()}
        return SymmetricPolynomial(dictionary, ans.printer, ans.multiplier)

    @classmethod
    def mp_dual_stable_grothendieck_p(cls, num_variables, mu, nu=()):  # noqa
        ans = cls.dual_stable_grothendieck_p(num_variables, mu, nu)
        dictionary = {k: v.substitute(0, -BETA) for (k, v) in ans.dictionary.items()}
        return SymmetricPolynomial(dictionary, ans.printer, ans.multiplier)

    @classmethod
    def mp_dual_stable_grothendieck_q(cls, num_variables, mu, nu=()):  # noqa
        ans = cls.dual_stable_grothendieck_q(num_variables, mu, nu)
        dictionary = {k: v.substitute(0, -BETA) for (k, v) in ans.dictionary.items()}
        return SymmetricPolynomial(dictionary, ans.printer, ans.multiplier)

    @classmethod
    def dual_stable_grothendieck(cls, num_variables, mu, nu=()):  # noqa
        return cls._dual_stable_grothendieck(num_variables, mu, nu)

    @cached_value(DUAL_STABLE_GROTHENDIECK_CACHE)
    def _dual_stable_grothendieck(cls, num_variables, mu, nu):  # noqa
        tableaux = Tableau.count_semistandard_rpp(num_variables, mu, nu)
        return (-BETA)**(sum(mu) - sum(nu)) * cls._vectorize(num_variables, tableaux, -BETA**-1)

    @classmethod
    def dual_stable_grothendieck_p(cls, num_variables, mu, nu=()):  # noqa
        return cls._dual_stable_grothendieck_p(num_variables, mu, nu)

    @cached_value(DUAL_STABLE_GROTHENDIECK_P_CACHE)
    def _dual_stable_grothendieck_p(cls, num_variables, mu, nu):  # noqa
        tableaux = Tableau.count_semistandard_marked_rpp(num_variables, mu, nu, diagonal_nonprimes=False)
        return (-BETA)**(sum(mu) - sum(nu)) * cls._vectorize(num_variables, tableaux, -BETA**-1)

    @classmethod
    def dual_stable_grothendieck_q(cls, num_variables, mu, nu=()):  # noqa
        return cls._dual_stable_grothendieck_q(num_variables, mu, nu)

    @cached_value(DUAL_STABLE_GROTHENDIECK_Q_CACHE)
    def _dual_stable_grothendieck_q(cls, num_variables, mu, nu):  # noqa
        tableaux = Tableau.count_semistandard_marked_rpp(num_variables, mu, nu, diagonal_nonprimes=True)
        return (-BETA)**(sum(mu) - sum(nu)) * cls._vectorize(num_variables, tableaux, -BETA**-1)

    @classmethod
    def schur_expansion(cls, f):
        if f:
            t = max(f)
            n = t.n
            c = f[t]
            mu = t.index()
            ans = cls.schur_expansion(f - c * cls.schur(n, mu))
            return ans + Vector({mu: c})
        else:
            return Vector()

    @classmethod
    def omega_schur_expansion(cls, f):
        return Vector({
            Partition.transpose(mu): coeff
            for mu, coeff in cls.schur_expansion(f).items()
        })

    @classmethod
    def grothendieck_expansion(cls, f):
        if f:
            t = max(f.lowest_degree_terms())
            n = t.n
            c = f[t]
            mu = t.index()
            ans = cls.grothendieck_expansion(f - c * cls.stable_grothendieck(n, mu))
            return ans + Vector({mu: c})
        else:
            return Vector()

    @classmethod
    def dual_grothendieck_expansion(cls, f):
        if f:
            t = max(f.highest_degree_terms())
            n = t.n
            c = f[t]
            mu = t.index()
            g = cls.dual_stable_grothendieck(n, mu)
            assert g[t] == 1
            ans = cls.dual_grothendieck_expansion(f - c * g)
            return ans + Vector({mu: c})
        else:
            return Vector()

    @classmethod
    def mp_dual_grothendieck_expansion(cls, f):
        if f:
            t = max(f.highest_degree_terms())
            n = t.n
            c = f[t]
            mu = t.index()
            g = cls.mp_dual_stable_grothendieck(n, mu)
            assert g[t] == 1
            ans = cls.mp_dual_grothendieck_expansion(f - c * g)
            return ans + Vector({mu: c})
        else:
            return Vector()

    @classmethod
    def transposed_dual_grothendieck_expansion(cls, f):
        if f:
            t = max(f.highest_degree_terms())
            n = t.n
            c = f[t]
            mu = t.index()
            g = cls.transposed_dual_grothendieck(n, mu)
            assert g[t] == 1
            ans = cls.transposed_dual_grothendieck(f - c * g)
            return ans + Vector({mu: c})
        else:
            return Vector()

    @classmethod
    def gp_free_expansion(cls, f):  # noqa
        if f:
            exp = cls.gp_expansion(f)
            mu = max(exp)
            n = max(f).n
            c = exp[mu]
            g = 1
            for part in mu:
                g *= cls.dual_stable_grothendieck_p(n, (part,))
            ans = cls.gp_free_expansion(f - c * g)
            return ans + Vector({mu: c})
        else:
            return Vector()

    @classmethod
    def gp_expansion(cls, f):  # noqa
        if f:
            t = max(f.highest_degree_terms())
            n = t.n
            c = f[t]
            mu = t.index()
            g = cls.dual_stable_grothendieck_p(n, mu)
            assert g[t] == 1
            ans = cls.gp_expansion(f - c * g)
            return ans + Vector({mu: c})
        else:
            return Vector()

    @classmethod
    def gq_expansion(cls, f):  # noqa
        if f:
            t = max(f.highest_degree_terms())
            n = t.n
            c = f[t]
            mu = t.index()
            g = cls.dual_stable_grothendieck_q(n, mu)
            assert g[t] == 2**len(mu)
            assert c % 2**len(mu) == 0
            c = c // 2**len(mu)
            ans = cls.gq_expansion(f - c * g)
            return ans + Vector({mu: c})
        else:
            return Vector()

    @classmethod
    def mp_gp_expansion(cls, f):  # noqa
        if f:
            t = max(f.highest_degree_terms())
            n = t.n
            c = f[t]
            mu = t.index()
            g = cls.mp_dual_stable_grothendieck_p(n, mu)
            assert g[t] == 1
            ans = cls.mp_gp_expansion(f - c * g)
            return ans + Vector({mu: c})
        else:
            return Vector()

    @classmethod
    def mp_gq_expansion(cls, f):  # noqa
        if f:
            t = max(f.highest_degree_terms())
            n = t.n
            c = f[t]
            mu = t.index()
            g = cls.mp_dual_stable_grothendieck_q(n, mu)
            assert g[t] == 2**len(mu)
            assert c % 2**len(mu) == 0
            c = c // 2**len(mu)
            ans = cls.mp_gq_expansion(f - c * g)
            return ans + Vector({mu: c})
        else:
            return Vector()


    @classmethod
    def jp_expansion(cls, f):  # noqa
        if f:
            t = max(f.highest_degree_terms())
            n = t.n
            c = f[t]
            mu = t.index()
            g = cls._slow_transposed_dual_stable_grothendieck_p(n, mu)
            assert g[t] == 1
            ans = cls.jp_expansion(f - c * g)
            return ans + Vector({mu: c})
        else:
            return Vector()

    @classmethod
    def jq_expansion(cls, f):  # noqa
        if f:
            t = max(f.highest_degree_terms())
            n = t.n
            c = f[t]
            mu = t.index()
            g = cls._slow_transposed_dual_stable_grothendieck_q(n, mu)
            assert g[t] == 2**len(mu)
            assert c % 2**len(mu) == 0
            c = c // 2**len(mu)
            ans = cls.jq_expansion(f - c * g)
            return ans + Vector({mu: c})
        else:
            return Vector()

    @classmethod
    def GP_expansion(cls, f):  # noqa
        if f:
            t = max(f.lowest_degree_terms())
            n = t.n
            c = f[t]
            mu = t.index()
            ans = cls.GP_expansion(f - c * cls.stable_grothendieck_p(n, mu))
            return ans + Vector({mu: c})
        else:
            return Vector()

    @classmethod
    def P_expansion(cls, f):  # noqa
        if f:
            t = max(f.lowest_degree_terms())
            n = t.n
            c = f[t]
            mu = t.index()
            ans = cls.P_expansion(f - c * cls.schur_p(n, mu))
            return ans + Vector({mu: c})
        else:
            return Vector()

    @classmethod
    def GQ_expansion(cls, f):  # noqa
        if f:
            t = max(f.lowest_degree_terms())
            n = t.n
            c = f[t]
            mu = t.index()
            g = cls.stable_grothendieck_q(n, mu)
            assert g[t] == 2**len(mu)
            assert c % 2**len(mu) == 0
            c = c // 2**len(mu)
            ans = cls.GQ_expansion(f - c * g)
            return ans + Vector({mu: c})
        else:
            return Vector()

    @classmethod
    def Q_expansion(cls, f):  # noqa
        if f:
            t = max(f.lowest_degree_terms())
            n = t.n
            c = f[t]
            mu = t.index()
            g = cls.schur_q(n, mu)
            assert g[t] == 2**len(mu)
            assert c % 2**len(mu) == 0
            c = c // 2**len(mu)
            ans = cls.Q_expansion(f - c * g)
            return ans + Vector({mu: c})
        else:
            return Vector()

    @classmethod
    def _slow_vectorize(cls, n, tableaux, signs=None, check=True):
        dictionary = defaultdict(int)
        for tab in tableaux:
            dictionary[tab.weight(n)] += 1
        if check:
            assert all(dictionary[Partition.sort(alpha)] == dictionary[alpha] for alpha in dictionary)
        return SymmetricPolynomial({
            SymmetricMonomial(n, alpha): coeff * (signs if signs else 1)**sum(alpha)
            for alpha, coeff in dictionary.items()
            if Partition.is_partition(alpha)
        })

    @classmethod
    def _slow_schur(cls, num_variables, mu, nu=()):
        return cls._slow_vectorize(num_variables, Tableau.semistandard(num_variables, mu, nu))

    @classmethod
    def _slow_schur_s(cls, num_variables, mu, nu=()):
        return cls._slow_vectorize(num_variables, Tableau.semistandard_marked(num_variables, mu, nu))

    @classmethod
    def _slow_schur_q(cls, num_variables, mu, nu=()):
        return cls._slow_vectorize(
            num_variables,
            Tableau.semistandard_shifted_marked(num_variables, mu, nu)
        )

    @classmethod
    def _slow_schur_p(cls, num_variables, mu, nu=()):
        return cls._slow_vectorize(
            num_variables,
            Tableau.semistandard_shifted_marked(num_variables, mu, nu, diagonal_primes=False)
        )

    @classmethod
    def _slow_stable_grothendieck_s(cls, num_variables, mu, nu=()):
        return BETA**(sum(nu) - sum(mu)) * cls._slow_vectorize(
            num_variables,
            Tableau.semistandard_marked_setvalued(num_variables, mu, nu),
            BETA
        )

    @classmethod
    def _slow_stable_grothendieck_q(cls, num_variables, mu, nu=()):
        return BETA**(sum(nu) - sum(mu)) * cls._slow_vectorize(
            num_variables,
            Tableau.semistandard_shifted_marked_setvalued(num_variables, mu, nu),
            BETA
        )

    @classmethod
    def _slow_stable_grothendieck_p(cls, num_variables, mu, nu=()):
        return BETA**(sum(nu) - sum(mu)) * cls._slow_vectorize(
            num_variables,
            Tableau.semistandard_shifted_marked_setvalued(num_variables, mu, nu, diagonal_primes=False),
            BETA
        )

    @classmethod
    def _slow_stable_grothendieck(cls, num_variables, mu, nu=()):
        return BETA**(sum(nu) - sum(mu)) * cls._slow_vectorize(
            num_variables,
            Tableau.semistandard_setvalued(num_variables, mu, nu),
            BETA
        )

    @classmethod
    def _slow_dual_stable_grothendieck(cls, num_variables, mu, nu=()):
        return (-BETA)**(sum(mu) - sum(nu)) * cls._slow_vectorize(
            num_variables,
            Tableau.semistandard_rpp(num_variables, mu, nu),
            -BETA**-1
        )

    @classmethod
    def _slow_dual_stable_grothendieck_p(cls, num_variables, mu, nu=()):
        return (-BETA)**(sum(mu) - sum(nu)) * cls._slow_vectorize(
            num_variables,
            Tableau.semistandard_marked_rpp(num_variables, mu, nu, diagonal_nonprimes=False),
            -BETA**-1
        )

    @classmethod
    def _slow_dual_stable_grothendieck_q(cls, num_variables, mu, nu=()):
        return (-BETA)**(sum(mu) - sum(nu)) * cls._slow_vectorize(
            num_variables,
            Tableau.semistandard_marked_rpp(num_variables, mu, nu, diagonal_nonprimes=True),
            -BETA**-1
        )

    @classmethod
    def _slow_transposed_dual_stable_grothendieck(cls, num_variables, mu, nu=(), beta=BETA):
        p = Polynomial.zero()
        for tab in Tableau.semistandard(num_variables, mu, nu):
            m = 1
            for i in range(1, num_variables + 1):
                r = len({x for x, y, v in tab if v[0] == i})
                a = len({(x, y) for x, y, v in tab if v[0] == i})
                x = Polynomial.x(i)
                m *= x**r * (x + 1)**(a - r)
            p += m
        dictionary = {}
        for e in p:
            tup = num_variables * [0]
            for i in e:
                tup[i - 1] = e[i]
            dictionary[tuple(tup)] = p[e]
        return SymmetricPolynomial({
            SymmetricMonomial(num_variables, alpha): coeff * (-beta**-1)**(sum(alpha))
            for alpha, coeff in dictionary.items()
            if Partition.is_partition(alpha)
        }) * (-beta)**(sum(mu) - sum(nu))

    @classmethod
    def _slow_transposed_dual_stable_grothendieck_p(cls, num_variables, mu, nu=(), beta=BETA):
        p = Polynomial.zero()
        for tab in Tableau.semistandard_shifted_marked(num_variables, mu, nu, diagonal_primes=False):
            m = 1
            for i in range(1, num_variables + 1):
                r = len({x for x, y, v in tab if v[0] == i})
                c = len({y for x, y, v in tab if v[0] == -i})
                a = len({(x, y) for x, y, v in tab if abs(v[0]) == i})
                x = Polynomial.x(i)
                m *= x**(r + c) * (x + 1)**(a - r - c)
            p += m
        dictionary = {}
        for e in p:
            tup = num_variables * [0]
            for i in e:
                tup[i - 1] = e[i]
            dictionary[tuple(tup)] = p[e]
        return SymmetricPolynomial({
            SymmetricMonomial(num_variables, alpha): coeff * (-beta**-1)**(sum(alpha))
            for alpha, coeff in dictionary.items()
            if Partition.is_partition(alpha)
        }) * (-beta)**(sum(mu) - sum(nu))

    @classmethod
    def _slow_transposed_dual_stable_grothendieck_q(cls, num_variables, mu, nu=(), beta=BETA):
        p = Polynomial.zero()
        for tab in Tableau.semistandard_shifted_marked(num_variables, mu, nu):
            m = 1
            for i in range(1, num_variables + 1):
                r = len({x for x, y, v in tab if v[0] == i})
                c = len({y for x, y, v in tab if v[0] == -i})
                a = len({(x, y) for x, y, v in tab if abs(v[0]) == i})
                x = Polynomial.x(i)
                m *= x**(r + c) * (x + 1)**(a - r - c)
            p += m
        dictionary = {}
        for e in p:
            tup = num_variables * [0]
            for i in e:
                tup[i - 1] = e[i]
            dictionary[tuple(tup)] = p[e]
        return SymmetricPolynomial({
            SymmetricMonomial(num_variables, alpha): coeff * (-beta**-1)**(sum(alpha))
            for alpha, coeff in dictionary.items()
            if Partition.is_partition(alpha)
        }) * (-beta)**(sum(mu) - sum(nu))

    @classmethod
    def _diag_vector(cls, rpp):
        ans = []
        i, j = (1, 1)
        while (i, j) in rpp:
            ans.append(abs(rpp.get(i, j, unpack=True)))
            i, j = (i + 1, j + 1)
        return tuple(ans)

    @classmethod
    def _slow_refined_dual_stable_grothendieck_p(cls, num_variables, mu, nu=()):
        dictionary = defaultdict(list)
        for t in Tableau.semistandard_marked_rpp(num_variables, mu, nu, diagonal_nonprimes=False):
            dictionary[cls._diag_vector(t)].append(t)
        mapping = {}
        for t in dictionary:
            mapping[t] = (-BETA)**(sum(mu) - sum(nu)) * cls._slow_vectorize(
                num_variables,
                dictionary[t],
                -BETA**-1,
                check=False,
            )
        return Vector(mapping)

    @classmethod
    def _slow_refined_dual_stable_grothendieck_q(cls, num_variables, mu, nu=()):
        dictionary = defaultdict(list)
        for t in Tableau.semistandard_marked_rpp(num_variables, mu, nu, diagonal_nonprimes=True):
            dictionary[cls._diag_vector(t)].append(t)
        mapping = {}
        for t in dictionary:
            mapping[t] = (-BETA)**(sum(mu) - sum(nu)) * cls._slow_vectorize(
                num_variables,
                dictionary[t],
                -BETA**-1,
                check=False,
            )
        return Vector(mapping)
