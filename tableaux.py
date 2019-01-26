from collections import defaultdict

COUNT_SEMISTANDARD_CACHE = {}
COUNT_SEMISTANDARD_SETVALUED_CACHE = {}
COUNT_SEMISTANDARD_MARKED_CACHE = {}
COUNT_SEMISTANDARD_MARKED_SETVALUED_CACHE = {}

SEMISTANDARD_CACHE = {}
SEMISTANDARD_SETVALUED_CACHE = {}
SEMISTANDARD_MARKED_CACHE = {}
SEMISTANDARD_MARKED_SETVALUED_CACHE = {}

HORIZONTAL_STRIPS_CACHE = {}


def nchoosek(m, k):
    ans = 1
    for i in range(k):
        ans *= m - i
    for i in range(k):
        ans //= i + 1
    return ans


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

    @classmethod
    def transpose_partition(cls, mu):
        if mu:
            return tuple(len([i for i in range(len(mu)) if mu[i] > j]) for j in range(mu[0]))
        else:
            return mu

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

    @classmethod
    def _horizontal_strips(cls, mu):
        core = tuple(mu[i + 1] if i + 1 < len(mu) else 0 for i in range(len(mu)))
        if mu not in HORIZONTAL_STRIPS_CACHE:
            HORIZONTAL_STRIPS_CACHE[mu] = []
            level = {core}
            while level:
                for nu in level:
                    diff = {(i + 1, j + 1) for i in range(len(mu)) for j in range(nu[i], mu[i])}
                    nu = nu if nu and nu[-1] > 0 else nu[:-1]
                    corners = [(i + 1, nu[i]) for i in range(len(nu)) if core[i] < nu[i]]
                    HORIZONTAL_STRIPS_CACHE[mu].append((nu, diff, corners))
                level = {
                    nu[:i] + (nu[i] + 1,) + nu[i + 1:]
                    for i in range(len(mu))
                    for nu in level
                    if nu[i] < mu[i]
                }
        return HORIZONTAL_STRIPS_CACHE[mu]

    @classmethod
    def _vertical_strips(cls, mu):
        ans = []
        for nu, diff, corners in cls._horizontal_strips(cls.transpose_partition(mu)):
            nu = cls.transpose_partition(nu)
            diff = {(j, i) for (i, j) in diff}
            corners = [(j, i) for (i, j) in corners]
            ans.append((nu, diff, corners))
        return ans

    @classmethod
    def semistandard(cls, mu, n, setvalued=False):
        cache = SEMISTANDARD_SETVALUED_CACHE if setvalued else SEMISTANDARD_CACHE
        if (mu, n) not in cache:
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
            cache[(mu, n)] = ans
        return cache[(mu, n)]

    @classmethod
    def semistandard_setvalued(cls, mu, n):
        return cls.semistandard(mu, n, True)

    @classmethod
    def semistandard_marked(cls, mu, n, setvalued=False):
        cache = SEMISTANDARD_MARKED_SETVALUED_CACHE if setvalued else SEMISTANDARD_MARKED_CACHE
        if (mu, n) not in cache:
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
            cache[(mu, n)] = ans
        return cache[(mu, n)]

    @classmethod
    def semistandard_marked_setvalued(cls, mu, n):
        return cls.semistandard_marked(mu, n, True)

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
