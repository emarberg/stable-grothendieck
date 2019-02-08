from cached import cached_value
from collections import defaultdict

COUNT_SEMISTANDARD_CACHE = {}
COUNT_SEMISTANDARD_MARKED_CACHE = {}
COUNT_SEMISTANDARD_SHIFTED_MARKED_CACHE = {}

SEMISTANDARD_CACHE = {}
SEMISTANDARD_MARKED_CACHE = {}
SEMISTANDARD_SHIFTED_MARKED_CACHE = {}

HORIZONTAL_STRIPS_CACHE = {}
SHIFTED_HORIZONTAL_STRIPS_CACHE = {}
SHIFTED_VERTICAL_STRIPS_CACHE = {}


def nchoosek(m, k):
    ans = 1
    for i in range(k):
        ans *= m - i
    for i in range(k):
        ans //= i + 1
    return ans


class Partition:

    @classmethod
    def sort(cls, mu, trim=False):
        return tuple(reversed(sorted(filter(bool, mu)))) if trim else tuple(reversed(sorted(mu)))

    @classmethod
    def is_partition(cls, mu):
        return all(mu[i - 1] >= mu[i] for i in range(1, len(mu))) and (mu == () or mu[-1] >= 0)

    @classmethod
    def is_strict_partition(cls, mu):
        return all(mu[i - 1] > mu[i] for i in range(1, len(mu)))

    @classmethod
    def transpose(cls, mu):
        if mu:
            return tuple(len([i for i in range(len(mu)) if mu[i] > j]) for j in range(mu[0]))
        else:
            return mu

    @classmethod
    def generate(cls, n, max_part=None, strict=False):
        if n == 0:
            yield ()
        else:
            max_part = n if (max_part is None or max_part > n) else max_part
            for i in range(1, max_part + 1):
                for mu in cls.generate(n - i, i):
                    nu = (i,) + mu
                    if not strict or Partition.is_strict_partition(nu):
                        yield (i,) + mu


class Tableau:

    def __init__(self, mapping={}):
        def tuplize(i, j):
            v = mapping[(i, j)]
            if type(v) == tuple:
                assert all(type(i) == int for i in v)
                return tuple(sorted(v))
            else:
                assert type(v) == int
                return (v,)

        self.boxes = {(i, j): tuplize(i, j) for i, j in mapping}

    def __eq__(self, other):
        assert type(other) == Tableau
        return self.boxes == other.boxes

    def __hash__(self):
        return hash(tuple(sorted(self.boxes.items())))

    def __bool__(self):
        return bool(self.boxes)

    def __contains__(self, box):
        return box in self.boxes

    def __iter__(self):
        for i, j in sorted(self.boxes):
            yield i, j, self.boxes[(i, j)]

    def find(self, v):
        return [(i, j) for i, j, values in self if v in values]

    def get(self, i, j, default=None):
        return self.boxes.get((i, j), default)

    def add(self, i, j, v):
        mapping = self.boxes.copy()
        if (i, j) not in self:
            mapping[(i, j)] = v
            return Tableau(mapping)
        mapping[(i, j)] = tuple(sorted(mapping[(i, j)] + (v,)))
        return Tableau(mapping)

    def remove(self, i, j, v):
        assert (i, j) in self
        assert v in self.get(i, j)
        mapping = self.boxes.copy()
        k = 0
        while self.get(i, j)[k] != v:
            k += 1
        mapping[(i, j)] = self.get(i, j)[:k] + self.get(i, j)[k + 1:]
        if len(mapping[(i, j)]) == 0:
            del mapping[(i, j)]
        return Tableau(mapping)

    def clear(self, i, j):
        assert (i, j) in self
        mapping = self.boxes.copy()
        del mapping[(i, j)]
        return Tableau(mapping)

    def max_row(self):
        return max([i for i, j in self.boxes]) if self.boxes else 0

    def max_column(self):
        return max([j for i, j in self.boxes]) if self.boxes else 0

    def _string_array(self):
        boxes = {k: ','.join(map(str, v)) for k, v in self.boxes.items()}
        maximum = max([len(str(v)) for v in boxes.values()]) if self.boxes else 0

        def pad(x):
            return (maximum - len(str(x))) * ' ' + str(x)

        array = []
        for i in range(1, self.max_row() + 1):
            row = []
            for j in range(1, self.max_column() + 1):
                row += [pad(boxes.get((i, j), '.'))]
            array += [row]
        return array

    def __repr__(self):
        array = self._string_array()
        return '\n' + '\n'.join([' '.join(line) for line in reversed(array)]) + '\n'

    def weight(self, n):
        ans = n * [0]
        for i, j, v in self:
            for x in v:
                ans[abs(x) - 1] += 1
        return tuple(ans)

    @cached_value(HORIZONTAL_STRIPS_CACHE)
    def _horizontal_strips(cls, mu):  # noqa
        core = tuple(mu[i + 1] if i + 1 < len(mu) else 0 for i in range(len(mu)))
        ans = []
        level = {core}
        while level:
            for nu in level:
                diff = {(i + 1, j + 1) for i in range(len(mu)) for j in range(nu[i], mu[i])}
                nu = nu if nu and nu[-1] > 0 else nu[:-1]
                corners = [(i + 1, nu[i]) for i in range(len(nu)) if core[i] < nu[i]]
                ans.append((nu, diff, corners))
            level = {
                nu[:i] + (nu[i] + 1,) + nu[i + 1:]
                for i in range(len(mu))
                for nu in level
                if nu[i] < mu[i]
            }
        return ans

    @classmethod
    def _vertical_strips(cls, mu):
        ans = []
        for nu, diff, corners in cls._horizontal_strips(Partition.transpose(mu)):
            nu = Partition.transpose(nu)
            diff = {(j, i) for (i, j) in diff}
            corners = [(j, i) for (i, j) in corners]
            ans.append((nu, diff, corners))
        return ans

    @cached_value(SHIFTED_HORIZONTAL_STRIPS_CACHE)
    def _shifted_horizontal_strips(cls, mu):  # noqa
        assert Partition.is_strict_partition(mu)
        core = tuple(mu[i + 1] + 1 if i + 1 < len(mu) else 0 for i in range(len(mu)))
        ans = []
        level = {core}
        while level:
            for nu in level:
                diff = {(i + 1, i + j + 1) for i in range(len(mu)) for j in range(nu[i], mu[i])}
                nu = nu if nu and nu[-1] > 0 else nu[:-1]
                corners = [(i + 1, i + nu[i]) for i in range(len(nu)) if core[i] < nu[i]]
                ans.append((nu, diff, corners))
            level = {
                nu[:i] + (nu[i] + 1,) + nu[i + 1:]
                for i in range(len(mu))
                for nu in level
                if nu[i] < mu[i]
            }
        return ans

    @cached_value(SHIFTED_VERTICAL_STRIPS_CACHE)
    def _shifted_vertical_strips(cls, mu):  # noqa
        assert Partition.is_strict_partition(mu)
        core = tuple(a - 1 for a in mu)
        ans = []
        level = {core}
        while level:
            for nu in level:
                diff = {(i + 1, i + j + 1) for i in range(len(mu)) for j in range(nu[i], mu[i])}
                nu = nu if nu and nu[-1] > 0 else nu[:-1]
                corners = [
                    (i + 1, i + nu[i]) for i in range(len(nu))
                    if core[i] < nu[i] and (i == len(nu) - 1 or 1 + nu[i + 1] < nu[i])
                ]
                ans.append((nu, diff, corners))
            level = {
                nu[:i] + (nu[i] + 1,) + nu[i + 1:]
                for i in range(len(mu))
                for nu in level
                if nu[i] < mu[i] and (i == 0 or 1 + nu[i] < nu[i - 1])
            }
        return ans

    @cached_value(COUNT_SEMISTANDARD_CACHE)
    def count_semistandard(cls, mu, n, setvalued=False):  # noqa
        ans = defaultdict(int)
        if mu == tuple():
            ans[()] = 1
        elif n > 0:
            for nu, diff, corners in cls._horizontal_strips(mu):
                for partition, count in cls.count_semistandard(nu, n - 1, setvalued).items():
                    for i in range(len(corners) + 1 if setvalued else 1):
                        m = len(diff) + i
                        if m == 0:
                            ans[partition] += count
                        elif len(partition) < n - 1 or (partition and m > partition[-1]):
                            break
                        else:
                            updated_partition = partition + (m,)
                            ans[updated_partition] += count * nchoosek(len(corners), i)
        return ans

    @classmethod
    def count_semistandard_setvalued(cls, mu, n):
        return cls.count_semistandard(mu, n, True)

    @cached_value(COUNT_SEMISTANDARD_MARKED_CACHE)
    def count_semistandard_marked(cls, mu, n, setvalued=False):  # noqa
        ans = defaultdict(int)
        if mu == tuple():
            ans[()] = 1
        elif n > 0:
            for nu1, diff1, corners1 in cls._horizontal_strips(mu):
                for nu2, diff2, corners2 in cls._vertical_strips(nu1):
                    for partition, count in cls.count_semistandard_marked(nu2, n - 1, setvalued).items():
                        for i in range(len(corners1) + 1 if setvalued else 1):
                            for j in range(len(corners2) + 1 if setvalued else 1):
                                m = len(diff1) + len(diff2) + i + j
                                if m == 0:
                                    ans[partition] += count
                                elif len(partition) < n - 1 or (partition and m > partition[-1]):
                                    break
                                else:
                                    updated_partition = partition + (m,)
                                    ans[updated_partition] += count * nchoosek(len(corners1), i) * nchoosek(len(corners2), j)
        return ans

    @classmethod
    def count_semistandard_marked_setvalued(cls, mu, n):
        return cls.count_semistandard_marked(mu, n, True)

    @cached_value(COUNT_SEMISTANDARD_SHIFTED_MARKED_CACHE)
    def count_semistandard_shifted_marked(cls, mu, n, diagonalprimes=True, setvalued=False):  # noqa
        assert Partition.is_strict_partition(mu)
        ans = defaultdict(int)
        if mu == tuple():
            ans[()] = 1
        elif n > 0:
            for nu1, diff1, corners1 in cls._shifted_horizontal_strips(mu):
                for nu2, diff2, corners2 in cls._shifted_vertical_strips(nu1):
                    if not diagonalprimes:
                        if any(i == j for (i, j) in diff2):
                            continue
                        corners2 = [(i, j) for (i, j) in corners2 if i != j]
                    for partition, count in cls.count_semistandard_shifted_marked(nu2, n - 1, diagonalprimes, setvalued).items():
                        for i in range(len(corners1) + 1 if setvalued else 1):
                            for j in range(len(corners2) + 1 if setvalued else 1):
                                m = len(diff1) + len(diff2) + i + j
                                if m == 0:
                                    ans[partition] += count
                                elif len(partition) < n - 1 or (partition and m > partition[-1]):
                                    break
                                else:
                                    updated_partition = partition + (m,)
                                    ans[updated_partition] += count * nchoosek(len(corners1), i) * nchoosek(len(corners2), j)
        return ans

    @classmethod
    def count_semistandard_shifted_marked_setvalued(cls, mu, n, diagonalprimes=True):
        return cls.count_semistandard_shifted_marked(mu, n, diagonalprimes, True)

    @cached_value(SEMISTANDARD_CACHE)
    def semistandard(cls, mu, n, setvalued=False):  # noqa
        ans = set()
        if mu == tuple():
            ans = {Tableau()}
        elif n > 0:
            for nu, diff, corners in cls._horizontal_strips(mu):
                    for aug in cls._subsets(diff, corners, setvalued):
                        for tab in cls.semistandard(nu, n - 1, setvalued):
                            for (i, j) in aug:
                                tab = tab.add(i, j, n)
                            ans.add(tab)
        return ans

    @classmethod
    def semistandard_setvalued(cls, mu, n):
        return cls.semistandard(mu, n, True)

    @cached_value(SEMISTANDARD_MARKED_CACHE)
    def semistandard_marked(cls, mu, n, setvalued=False):  # noqa
        ans = set()
        if mu == tuple():
            ans = {Tableau()}
        elif n > 0:
            for nu1, diff1, corners1 in cls._horizontal_strips(mu):
                for nu2, diff2, corners2 in cls._vertical_strips(nu1):
                    for aug1 in cls._subsets(diff1, corners1, setvalued):
                        for aug2 in cls._subsets(diff2, corners2, setvalued):
                            for tab in cls.semistandard_marked(nu2, n - 1, setvalued):
                                for (i, j) in aug1:
                                    tab = tab.add(i, j, n)
                                for (i, j) in aug2:
                                    tab = tab.add(i, j, -n)
                                ans.add(tab)
        return ans

    @classmethod
    def semistandard_marked_setvalued(cls, mu, n):
        return cls.semistandard_marked(mu, n, True)

    @cached_value(SEMISTANDARD_SHIFTED_MARKED_CACHE)
    def semistandard_shifted_marked(cls, mu, n, diagonalprimes=True, setvalued=False):  # noqa
        assert Partition.is_strict_partition(mu)
        ans = set()
        if mu == tuple():
            ans = {Tableau()}
        elif n > 0:
            for nu1, diff1, corners1 in cls._shifted_horizontal_strips(mu):
                for nu2, diff2, corners2 in cls._shifted_vertical_strips(nu1):
                    if not diagonalprimes:
                        if any(i == j for (i, j) in diff2):
                            continue
                        corners2 = [(i, j) for (i, j) in corners2 if i != j]
                    for aug1 in cls._subsets(diff1, corners1, setvalued):
                        for aug2 in cls._subsets(diff2, corners2, setvalued):
                            for tab in cls.semistandard_shifted_marked(nu2, n - 1, diagonalprimes, setvalued):
                                for (i, j) in aug1:
                                    tab = tab.add(i, j, n)
                                for (i, j) in aug2:
                                    tab = tab.add(i, j, -n)
                                ans.add(tab)
        return ans

    @classmethod
    def semistandard_shifted_marked_setvalued(cls, mu, n, diagonalprimes=True):
        return cls.semistandard_shifted_marked(mu, n, diagonalprimes, True)

    @classmethod
    def _subsets(cls, diff, corners, setvalued):
        if setvalued:
            for v in range(2**len(corners)):
                thisdiff = diff
                for i in range(len(corners)):
                    thisdiff = thisdiff if v % 2 == 0 else thisdiff | {corners[i]}
                    v = v // 2
                yield thisdiff
        else:
            yield diff
