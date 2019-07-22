from cached import cached_value
from words import Word
from collections import defaultdict

COUNT_SEMISTANDARD_CACHE = {}
COUNT_SEMISTANDARD_MARKED_CACHE = {}
COUNT_SEMISTANDARD_SHIFTED_MARKED_CACHE = {}
COUNT_SEMISTANDARD_RPP_CACHE = {}
COUNT_SEMISTANDARD_MARKED_RPP_CACHE = {}

SEMISTANDARD_CACHE = {}
SEMISTANDARD_RPP_CACHE = {}
SEMISTANDARD_MARKED_RPP_CACHE = {}
SEMISTANDARD_MARKED_CACHE = {}
SEMISTANDARD_SHIFTED_MARKED_CACHE = {}

HORIZONTAL_STRIPS_CACHE = {}
SHIFTED_HORIZONTAL_STRIPS_CACHE = {}
SHIFTED_VERTICAL_STRIPS_CACHE = {}

RPP_HORIZONTAL_STRIPS_CACHE = {}
SHIFTED_RPP_HORIZONTAL_STRIPS_CACHE = {}
SHIFTED_RPP_VERTICAL_STRIPS_CACHE = {}

PARTITIONS = {}


def nchoosek(m, k):
    ans = 1
    for i in range(k):
        ans *= m - i
    for i in range(k):
        ans //= i + 1
    return ans


class Partition:

    @classmethod
    def complement(cls, n, mu):
        assert all(mu[i] <= n - i for i in range(len(mu)))
        dictionary = defaultdict(int)
        for i in range(n):
            start = i + (mu[i] if i < len(mu) else 0) + 1
            for j in range(start, n + 1):
                dictionary[j] += 1
        return Partition.sort(dictionary.values(), trim=True)

    @classmethod
    def printable(cls, mu, shifted=False):
        s = []
        for i, a in enumerate(mu):
            s = [(i * '  ' if shifted else '') + a * '* '] + s
        return '\n'.join(s)

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
    def skew(cls, mu, nu, shifted=False):
        ans = set()
        for i, part in enumerate(mu):
            subpart = nu[i] if i < len(nu) else 0
            for j in range(subpart, part):
                ans.add((i + 1, j + 1 + (i if shifted else 0)))
        return ans

    @classmethod
    def contains(cls, bigger, smaller):
        """Returns true if mu subseteq nu as partitions."""
        if len(smaller) > len(bigger):
            return False
        return all(0 <= smaller[i] <= bigger[i] for i in range(len(smaller)))

    @classmethod
    def all(cls, n, strict=False):
        for i in range(n + 1):
            for mu in cls.generate(i, strict=strict):
                yield mu

    @classmethod
    def generate(cls, n, max_part=None, strict=False):
        if n == 0:
            yield ()
        else:
            if (n, max_part, strict) not in PARTITIONS:
                ans = []
                max_part = n if (max_part is None or max_part > n) else max_part
                for i in range(1, max_part + 1):
                    for mu in cls.generate(n - i, i):
                        nu = (i,) + mu
                        if not strict or Partition.is_strict_partition(nu):
                            ku = (i,) + mu
                            ans.append(ku)
                            yield ku
                PARTITIONS[(n, max_part, strict)] = ans
            else:
                for mu in PARTITIONS[(n, max_part, strict)]:
                    yield mu


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
        self._sorting_word = tuple(
            (2 * v if v > 0 else -1 - 2 * v)
            for b in sorted(self.boxes) for v in self.boxes[b]
        )
        self._string_array = None

    def __eq__(self, other):
        assert type(other) == type(self)
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

    def __len__(self):
        return len(self._sorting_word)

    def __lt__(self, other):
        return len(self) < len(other) or (len(self) == len(other) and self._sorting_word < other._sorting_word)

    def transpose(self):
        return Tableau({(j, i): v for i, j, v in self})

    def standardize(self):
        entries = sorted(
            [(v, i, j) for i, j, value in self for v in value],
            key=lambda x: (-2 * x[0] - 1, x[1]) if x[0] < 0 else (2 * x[0], x[2])
        )
        ans = Tableau()
        for a, (v, i, j) in enumerate(entries):
            a = a + 1 if v > 0 else -a - 1
            ans = ans.add(i, j, a)
        return ans

    def is_semistandard(self):
        for i, j, values in self:
            if (i, j + 1) in self and max(values) > min(self.get(i, j + 1, unpack=False)):
                return False
            if (i + 1, j) in self and max(values) >= min(self.get(i + 1, j, unpack=False)):
                return False
        return True

    def crystal_reading_word(self):
        columns = defaultdict(dict)
        for i, j, value in self:
            columns[j][i] = value
        columns = [[columns[j][i] for i in sorted(columns[j])] for j in sorted(columns)]
        ans = ()
        for col in columns:
            pre, post = (), ()
            for row in col:
                pre = row[:1] + pre
                post = post + row[1:]
            ans += (pre + post,)
        return ans

    @classmethod
    def from_crystal_reading_word(cls, crystal_reading_word):
        if crystal_reading_word is None:
            return None
        ans = Tableau()
        for j, col in enumerate(crystal_reading_word):
            rows = []
            for i in range(len(col)):
                if i == 0 or col[i - 1] > col[i]:
                    rows = [col[i]] + rows
                else:
                    r = 1 + [a for a in range(len(rows)) if rows[a] <= col[i]][-1]
                    ans = ans.add(r, j + 1, col[i])
            for i, a in enumerate(rows):
                ans = ans.add(i + 1, j + 1, a)
        return ans

    def e_crystal_operator(self, i, multisetvalued=True):
        if multisetvalued:
            crystal_word = self.crystal_reading_word()
            operated = Word.e_crystal_operator(i, *crystal_word, unpack=False)
            return self.from_crystal_reading_word(operated)
        else:
            s = sorted(
                [(x, y, v) for x, y, value in self for v in value if v in [i, i + 1]],
                key=lambda t: (t[1], -t[0], -t[2])
            )
            ell = len(s) + 1
            while len(s) < ell:
                ell = len(s)
                for x in range(len(s) - 1):
                    if s[x][-1] > s[x + 1][-1]:
                        s = s[:x] + s[x + 2:]
                        break
            s = [p for p in s if p[-1] == i + 1]
            if len(s) == 0:
                return None
            x, y, _ = s[0]
            if (x, y - 1) in self and i + 1 in self.get(x, y - 1, unpack=False):
                return self.remove(x, y - 1, i + 1).add(x, y, i)
            else:
                return self.remove(x, y, i + 1).add(x, y, i)

    def f_crystal_operator(self, i, multisetvalued=True):
        if multisetvalued:
            crystal_word = self.crystal_reading_word()
            operated = Word.f_crystal_operator(i, *crystal_word, unpack=False)
            return self.from_crystal_reading_word(operated)
        else:
            s = sorted(
                [(x, y, v) for x, y, value in self for v in value if v in [i, i + 1]],
                key=lambda t: (t[1], -t[0], -t[2])
            )
            ell = len(s) + 1
            while len(s) < ell:
                ell = len(s)
                for x in range(len(s) - 1):
                    if s[x][-1] > s[x + 1][-1]:
                        s = s[:x] + s[x + 2:]
                        break
            s = [p for p in s if p[-1] == i]
            if len(s) == 0:
                return None
            x, y, _ = s[-1]
            if (x, y + 1) in self and i in self.get(x, y + 1, unpack=False):
                return self.remove(x, y + 1, i).add(x, y, i + 1)
            else:
                return self.remove(x, y, i).add(x, y, i + 1)

    def apply(self, function):
        return Tableau(mapping={
            (i, j): tuple(function(v) for v in value)
            for i, j, value in self
        })

    def double(self):
        return self.apply(lambda x: 2 * x)

    def halve(self):
        assert all(v % 2 == 0 for _, _, value in self for v in value)
        return self.apply(lambda x: x // 2)

    def negate(self, skipdiagonal=True):
        if skipdiagonal:
            return Tableau(mapping={
                (i, j): tuple(-v if i != j else v for v in value)
                for i, j, value in self
            })
        else:
            return self.apply(lambda x: -x)

    def serialize(self):
        return self.boxes

    def find(self, v):
        return [(i, j) for i, j, values in self if v in values]

    def get(self, i, j, default=None, unpack=True):
        ans = self.boxes.get((i, j), default)
        if unpack and type(ans) == tuple and len(ans) == 1:
            ans = ans[0]
        return ans

    def add(self, i, j, v):
        mapping = self.boxes.copy()
        if (i, j) not in self:
            mapping[(i, j)] = v
            return self.__class__(mapping)
        mapping[(i, j)] = tuple(sorted(mapping[(i, j)] + (v,)))
        return self.__class__(mapping)

    def remove(self, i, j, v=None):
        assert (i, j) in self
        mapping = self.boxes.copy()
        if v is None:
            del mapping[(i, j)]
            return self.__class__(mapping)

        assert v in self.get(i, j, unpack=False)
        k = 0
        while self.get(i, j, unpack=False)[k] != v:
            k += 1
        mapping[(i, j)] = self.get(i, j, unpack=False)[:k] + self.get(i, j, unpack=False)[k + 1:]
        if len(mapping[(i, j)]) == 0:
            del mapping[(i, j)]
        return self.__class__(mapping)

    def replace(self, i, j, v):
        assert (i, j) in self
        mapping = self.boxes.copy()
        mapping[(i, j)] = v
        return self.__class__(mapping)

    def clear(self, i, j):
        assert (i, j) in self
        mapping = self.boxes.copy()
        del mapping[(i, j)]
        return self.__class__(mapping)

    def values(self):
        ans = set()
        for i, j, v in self:
            ans |= set(v) if type(v) != int else {v}
        return ans

    def last(self):
        v = self.values()
        n = max(abs(min(v)), abs(max(v))) if v else 0
        assert not (n in v and -n in v)
        return n if n in v else -n

    def pop(self):
        n = self.last()
        boxes = self.find(n)
        assert len(boxes) == 1
        i, j = list(boxes)[0]
        record = i, j, self.get(i, j)
        v = self.get(i, j)
        if n == v or -n == v or (n,) == v or (-n,) == v:
            return self.remove(i, j), record
        else:
            v = tuple(sorted(set(self.get(i, j)) - {n, -n}))
            return self.remove(i, j).add(i, j, v), record

    def max_row(self):
        return max([i for i, j in self.boxes]) if self.boxes else 0

    def max_column(self):
        return max([j for i, j in self.boxes]) if self.boxes else 0

    COLUMN_SPACING_OFFSET = 2

    @property
    def string_array(self):
        if self._string_array is None:
            boxes = {
                k: ''.join([
                    str(-x) + '\'' if x < 0 else str(x)
                    for x in sorted(v, key=lambda x:(abs(x), x))
                ]) for k, v in self.boxes.items()}

            def maximum(j):
                column = [len(str(boxes[b])) for b in boxes if b[1] == j]
                return max(column) if column else 0

            def pad(x, j):
                return str(x) + (self.COLUMN_SPACING_OFFSET + maximum(j) - len(str(x))) * ' '

            array = []
            for i in range(1, self.max_row() + 1):
                row = []
                for j in range(1, self.max_column() + 1):
                    row += [pad(boxes.get((i, j), ''), j)]
                array += [row]
            self._string_array = array
        return self._string_array

    def __repr__(self):
        # array = self.string_array  # english
        array = reversed(self.string_array)  # french
        return '\n\n' + '\n'.join([' '.join(line) for line in array]) + '\n\n'

    def weight(self, n):
        ans = n * [0]
        for i, j, v in self:
            for x in v:
                ans[abs(x) - 1] += 1
        return tuple(ans)

    @cached_value(HORIZONTAL_STRIPS_CACHE)
    def _horizontal_strips(cls, mu, lam):  # noqa
        if not Partition.contains(mu, lam):
            return []

        core = [mu[i + 1] if i + 1 < len(mu) else 0 for i in range(len(mu))]
        for i in range(len(lam)):
            core[i] = max(core[i], lam[i])
        core = tuple(core)

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
    def _vertical_strips(cls, mu, lam):
        ans = []
        strips = cls._horizontal_strips(Partition.transpose(mu), Partition.transpose(lam))
        for nu, diff, corners in strips:
            nu = Partition.transpose(nu)
            diff = {(j, i) for (i, j) in diff}
            corners = [(j, i) for (i, j) in corners]
            ans.append((nu, diff, corners))
        return ans

    @cached_value(SHIFTED_HORIZONTAL_STRIPS_CACHE)
    def _shifted_horizontal_strips(cls, mu, lam):  # noqa
        assert Partition.is_strict_partition(mu)
        if not Partition.contains(mu, lam):
            return []

        core = [mu[i + 1] + 1 if i + 1 < len(mu) else 0 for i in range(len(mu))]
        for i in range(len(lam)):
            core[i] = max(core[i], lam[i])
        core = tuple(core)

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
    def _shifted_vertical_strips(cls, mu, lam):  # noqa
        assert Partition.is_strict_partition(mu)
        if not Partition.contains(mu, lam):
            return []

        core = [a - 1 for a in mu]
        for i in range(len(lam)):
            core[i] = max(core[i], lam[i])
        core = tuple(core)

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

    @classmethod
    def _rpp_vertical_strips(cls, mu, lam):  # noqa
        strips = cls._rpp_horizontal_strips(Partition.transpose(mu), Partition.transpose(lam))
        return [
            (Partition.transpose(nu), {(j, i) for i, j in diff})
            for nu, diff in strips
        ]

    @cached_value(RPP_HORIZONTAL_STRIPS_CACHE)
    def _rpp_horizontal_strips(cls, mu, lam):  # noqa
        if mu == lam:
            return [(mu, set())]

        if not Partition.contains(mu, lam):
            return []

        def remove_box(nu, i):
            if i < len(nu) and nu[i] > 0:
                nu = nu[:i] + (nu[i] - 1,) + nu[i + 1:]
                while nu and nu[-1] == 0:
                    nu = nu[:-1]
                if all(nu[j] >= nu[j + 1] for j in range(len(nu) - 1)):
                    if Partition.contains(nu, lam):
                        yield nu

        def remove_all_boxes(nu, i):
            queue = [nu]
            while queue:
                nu, queue = queue[0], queue[1:]
                yield nu
                for x in remove_box(nu, i):
                    queue.append(x)

        ans = set()
        queue = [(mu, len(mu) - 1)]
        while queue:
            nu, i = queue[0]
            queue = queue[1:]
            if i >= 0:
                for nu in remove_all_boxes(nu, i):
                    ans.add(nu)
                    queue.append((nu, i - 1))

        return [(nu, Partition.skew(mu, nu)) for nu in ans]

    @cached_value(SHIFTED_RPP_HORIZONTAL_STRIPS_CACHE)
    def _shifted_rpp_horizontal_strips(cls, mu):  # noqa
        assert Partition.is_strict_partition(mu)
        if mu == ():
            return [(mu, set())]

        def remove_box(nu, i):
            if i < len(nu) and nu[i] > 0:
                nu = nu[:i] + (nu[i] - 1,) + nu[i + 1:]
                while nu and nu[-1] == 0:
                    nu = nu[:-1]
                if all(nu[j] > nu[j + 1] for j in range(len(nu) - 1)):
                    yield nu

        def remove_all_boxes(nu, i):
            queue = [nu]
            while queue:
                nu, queue = queue[0], queue[1:]
                yield nu
                for x in remove_box(nu, i):
                    queue.append(x)

        ans = set()
        queue = [(mu, len(mu) - 1)]
        while queue:
            nu, i = queue[0]
            queue = queue[1:]
            if i >= 0:
                for nu in remove_all_boxes(nu, i):
                    ans.add(nu)
                    queue.append((nu, i - 1))

        return [(nu, Partition.skew(mu, nu, shifted=True)) for nu in ans]

    @cached_value(SHIFTED_RPP_VERTICAL_STRIPS_CACHE)
    def _shifted_rpp_vertical_strips(cls, mu):  # noqa
        assert Partition.is_strict_partition(mu)
        if mu == ():
            return [(mu, set())]

        def remove_box(nu, i):
            for j in range(len(nu) - 1, -1, -1):
                if j + nu[j] == i + 1:
                    nu = nu[:j] + (nu[j] - 1,) + nu[j + 1:]
                    while nu and nu[-1] == 0:
                        nu = nu[:-1]
                    yield nu
                    return

        def remove_all_boxes(nu, i):
            queue = [nu]
            while queue:
                nu, queue = queue[0], queue[1:]
                yield nu
                for x in remove_box(nu, i):
                    queue.append(x)

        ans = set()
        queue = [(mu, (mu[0] if mu else 0) - 1)]
        while queue:
            nu, i = queue[0]
            queue = queue[1:]
            if i >= 0:
                for nu in remove_all_boxes(nu, i):
                    ans.add(nu)
                    queue.append((nu, i - 1))

        return [(nu, Partition.skew(mu, nu, shifted=True)) for nu in ans]

    @classmethod
    def count_semistandard(cls, max_entry, mu, nu=(), setvalued=False):  # noqa
        return cls._count_semistandard(max_entry, mu, nu, setvalued)

    @cached_value(COUNT_SEMISTANDARD_CACHE)
    def _count_semistandard(cls, max_entry, mu, lam, setvalued):  # noqa    
        ans = defaultdict(int)
        if mu == lam:
            ans[()] = 1
        elif Partition.contains(mu, lam) and max_entry > 0:
            for nu, diff, corners in cls._horizontal_strips(mu, lam):
                for partition, count in cls._count_semistandard(max_entry - 1, nu, lam, setvalued).items():
                    for i in range(len(corners) + 1 if setvalued else 1):
                        m = len(diff) + i
                        if m == 0:
                            ans[partition] += count
                        elif len(partition) < max_entry - 1 or (partition and m > partition[-1]):
                            break
                        else:
                            updated_partition = partition + (m,)
                            ans[updated_partition] += count * nchoosek(len(corners), i)
        return ans

    @classmethod
    def count_semistandard_setvalued(cls, max_entry, mu, nu=()):
        return cls.count_semistandard(max_entry, mu, nu, True)

    @classmethod
    def count_semistandard_marked(cls, max_entry, mu, nu=(), setvalued=False):  # noqa
        return cls._count_semistandard_marked(max_entry, mu, nu, setvalued)

    @cached_value(COUNT_SEMISTANDARD_MARKED_CACHE)
    def _count_semistandard_marked(cls, max_entry, mu, lam, setvalued):  # noqa
        ans = defaultdict(int)
        if mu == lam:
            ans[()] = 1
        elif Partition.contains(mu, lam) and max_entry > 0:
            for nu1, diff1, corners1 in cls._horizontal_strips(mu, lam):
                for nu2, diff2, corners2 in cls._vertical_strips(nu1, lam):
                    for partition, count in cls._count_semistandard_marked(max_entry - 1, nu2, lam, setvalued).items():
                        for i in range(len(corners1) + 1 if setvalued else 1):
                            for j in range(len(corners2) + 1 if setvalued else 1):
                                m = len(diff1) + len(diff2) + i + j
                                if m == 0:
                                    ans[partition] += count
                                elif len(partition) < max_entry - 1 or (partition and m > partition[-1]):
                                    break
                                else:
                                    updated_partition = partition + (m,)
                                    ans[updated_partition] += count * nchoosek(len(corners1), i) * nchoosek(len(corners2), j)
        return ans

    @classmethod
    def count_semistandard_marked_setvalued(cls, max_entry, mu, nu=()):
        return cls.count_semistandard_marked(max_entry, mu, nu, True)

    @classmethod
    def count_semistandard_shifted_marked(cls, max_entry, mu, nu=(), diagonal_primes=True, setvalued=False):  # noqa
        return cls._count_semistandard_shifted_marked(max_entry, mu, nu, diagonal_primes, setvalued)

    @cached_value(COUNT_SEMISTANDARD_SHIFTED_MARKED_CACHE)
    def _count_semistandard_shifted_marked(cls, max_entry, mu, lam, diagonal_primes, setvalued):  # noqa
        assert Partition.is_strict_partition(mu)
        ans = defaultdict(int)
        if mu == lam:
            ans[()] = 1
        elif Partition.contains(mu, lam) and max_entry > 0:
            for nu1, diff1, corners1 in cls._shifted_horizontal_strips(mu, lam):
                for nu2, diff2, corners2 in cls._shifted_vertical_strips(nu1, lam):
                    if not diagonal_primes:
                        if any(i == j for (i, j) in diff2):
                            continue
                        corners2 = [(i, j) for (i, j) in corners2 if i != j]
                    for partition, count in cls._count_semistandard_shifted_marked(max_entry - 1, nu2, lam, diagonal_primes, setvalued).items():
                        for i in range(len(corners1) + 1 if setvalued else 1):
                            for j in range(len(corners2) + 1 if setvalued else 1):
                                m = len(diff1) + len(diff2) + i + j
                                if m == 0:
                                    ans[partition] += count
                                elif len(partition) < max_entry - 1 or (partition and m > partition[-1]):
                                    break
                                else:
                                    updated_partition = partition + (m,)
                                    ans[updated_partition] += count * nchoosek(len(corners1), i) * nchoosek(len(corners2), j)
        return ans

    @classmethod
    def count_semistandard_rpp(cls, max_entry, mu, nu=()):  # noqa
        return cls._count_semistandard_rpp(max_entry, mu, nu)

    @cached_value(COUNT_SEMISTANDARD_RPP_CACHE)
    def _count_semistandard_rpp(cls, max_entry, mu, lam):  # noqa
        ans = defaultdict(int)
        if mu == lam:
            ans[()] = 1
        elif Partition.contains(mu, lam) and max_entry > 0:
            for nu, diff in cls._rpp_vertical_strips(mu, lam):
                if mu == nu:
                    continue
                m = len({j for _, j in diff})
                for partition, count in cls._count_semistandard_rpp(max_entry - 1, nu, lam).items():
                    if not (partition and m > partition[-1]):
                        ans[partition + (m,)] += count
        return ans

    @classmethod
    def count_semistandard_marked_rpp(cls, max_entry, mu, nu=(), diagonal_nonprimes=True):  # noqa
        return cls._count_semistandard_marked_rpp(max_entry, mu, nu, diagonal_nonprimes)

    @cached_value(COUNT_SEMISTANDARD_MARKED_RPP_CACHE)
    def _count_semistandard_marked_rpp(cls, max_entry, mu, lam, diagonal_nonprimes):  # noqa
        assert Partition.is_strict_partition(mu)
        ans = defaultdict(int)
        if mu == lam:
            ans[()] = 1
        elif Partition.contains(mu, lam) and max_entry > 0:
            for nu1, diff1 in cls._shifted_rpp_horizontal_strips(mu):
                if not diagonal_nonprimes and any(i == j for i, j in diff1):
                    continue
                for nu2, diff2 in cls._shifted_rpp_vertical_strips(nu1):
                    if mu == nu2:
                        continue
                    m = len({j for _, j in diff1}) + len({i for i, _ in diff2})
                    for partition, count in cls._count_semistandard_marked_rpp(max_entry - 1, nu2, lam, diagonal_nonprimes).items():
                        if not (partition and m > partition[-1]):
                            ans[partition + (m,)] += count
        return ans

    @classmethod
    def count_semistandard_shifted_marked_setvalued(cls, max_entry, mu, nu=(), diagonal_primes=True):
        return cls.count_semistandard_shifted_marked(max_entry, mu, nu, diagonal_primes, True)

    @classmethod
    def semistandard_rpp(cls, max_entry, mu, nu=()):  # noqa
        return cls._semistandard_rpp(max_entry, mu, nu)

    @cached_value(SEMISTANDARD_RPP_CACHE)
    def _semistandard_rpp(cls, max_entry, mu, lam):  # noqa
        ans = set()
        if mu == lam:
            ans = {ReversePlanePartition()}
        elif Partition.contains(mu, lam) and max_entry > 0:
            for nu, diff in cls._rpp_vertical_strips(mu, lam):
                for tab in cls._semistandard_rpp(max_entry - 1, nu, lam):
                    for (i, j) in diff:
                        tab = tab.add(i, j, max_entry)
                    ans.add(tab)
        return ans

    @classmethod
    def semistandard_shifted_rpp(cls, max_entry, mu, nu=(), diagonal_nonprimes=True):
        return cls.semistandard_marked_rpp(max_entry, mu, nu, diagonal_nonprimes)

    @classmethod
    def semistandard_marked_rpp(cls, max_entry, mu, nu=(), diagonal_nonprimes=True):  # noqa
        return cls._semistandard_marked_rpp(max_entry, mu, nu, diagonal_nonprimes)

    @cached_value(SEMISTANDARD_MARKED_RPP_CACHE)
    def _semistandard_marked_rpp(cls, max_entry, mu, lam, diagonal_nonprimes):  # noqa
        assert Partition.is_strict_partition(mu)
        ans = set()
        if mu == lam:
            ans = {MarkedReversePlanePartition()}
        elif Partition.contains(mu, lam) and max_entry > 0:
            for nu1, diff1 in cls._shifted_rpp_horizontal_strips(mu):
                for nu2, diff2 in cls._shifted_rpp_vertical_strips(nu1):
                    if diagonal_nonprimes or not any(i == j for i, j in diff1):
                        for tab in cls._semistandard_marked_rpp(max_entry - 1, nu2, lam, diagonal_nonprimes):
                            for (i, j) in diff1:
                                tab = tab.add(i, j, max_entry)
                            for (i, j) in diff2:
                                tab = tab.add(i, j, -max_entry)
                            ans.add(tab)
        return ans

    @classmethod
    def semistandard(cls, max_entry, mu, nu=(), setvalued=False):  # noqa
        return cls._semistandard(max_entry, mu, nu, setvalued)

    @cached_value(SEMISTANDARD_CACHE)
    def _semistandard(cls, max_entry, mu, lam, setvalued):  # noqa
        ans = set()
        if mu == lam:
            ans = {Tableau()}
        elif Partition.contains(mu, lam) and max_entry > 0:
            for nu, diff, corners in cls._horizontal_strips(mu, lam):
                    for aug in cls._subsets(diff, corners, setvalued):
                        for tab in cls._semistandard(max_entry - 1, nu, lam, setvalued):
                            for (i, j) in aug:
                                tab = tab.add(i, j, max_entry)
                            ans.add(tab)
        return ans

    @classmethod
    def semistandard_setvalued(cls, max_entry, mu, nu=()):
        return cls.semistandard(max_entry, mu, nu, setvalued=True)

    @classmethod
    def semistandard_marked(cls, max_entry, mu, nu=(), setvalued=False):  # noqa
        return cls._semistandard_marked(max_entry, mu, nu, setvalued)

    @cached_value(SEMISTANDARD_MARKED_CACHE)
    def _semistandard_marked(cls, max_entry, mu, lam, setvalued):  # noqa
        ans = set()
        if mu == lam:
            ans = {Tableau()}
        elif Partition.contains(mu, lam) and max_entry > 0:
            for nu1, diff1, corners1 in cls._horizontal_strips(mu, lam):
                for nu2, diff2, corners2 in cls._vertical_strips(nu1, lam):
                    for aug1 in cls._subsets(diff1, corners1, setvalued):
                        for aug2 in cls._subsets(diff2, corners2, setvalued):
                            for tab in cls._semistandard_marked(max_entry - 1, nu2, lam, setvalued):
                                for (i, j) in aug1:
                                    tab = tab.add(i, j, max_entry)
                                for (i, j) in aug2:
                                    tab = tab.add(i, j, -max_entry)
                                ans.add(tab)
        return ans

    @classmethod
    def semistandard_marked_setvalued(cls, max_entry, mu, nu=()):
        return cls.semistandard_marked(max_entry, mu, nu, True)

    @classmethod
    def semistandard_shifted(cls, max_entry, mu, nu=(), diagonal_primes=True, setvalued=False):
        return cls.semistandard_shifted_marked(max_entry, mu, nu, diagonal_primes, setvalued)

    @classmethod
    def semistandard_shifted_marked(cls, max_entry, mu, nu=(), diagonal_primes=True, setvalued=False):  # noqa
        return cls._semistandard_shifted_marked(max_entry, mu, nu, diagonal_primes, setvalued)

    @cached_value(SEMISTANDARD_SHIFTED_MARKED_CACHE)
    def _semistandard_shifted_marked(cls, max_entry, mu, lam, diagonal_primes, setvalued):  # noqa
        assert Partition.is_strict_partition(mu)
        ans = set()
        if mu == lam:
            ans = {Tableau()}
        elif Partition.contains(mu, lam) and max_entry > 0:
            for nu1, diff1, corners1 in cls._shifted_horizontal_strips(mu, lam):
                for nu2, diff2, corners2 in cls._shifted_vertical_strips(nu1, lam):
                    if not diagonal_primes:
                        if any(i == j for (i, j) in diff2):
                            continue
                        corners2 = [(i, j) for (i, j) in corners2 if i != j]
                    for aug1 in cls._subsets(diff1, corners1, setvalued):
                        for aug2 in cls._subsets(diff2, corners2, setvalued):
                            for tab in cls._semistandard_shifted_marked(max_entry - 1, nu2, lam, diagonal_primes, setvalued):
                                for (i, j) in aug1:
                                    tab = tab.add(i, j, max_entry)
                                for (i, j) in aug2:
                                    tab = tab.add(i, j, -max_entry)
                                ans.add(tab)
        return ans

    @classmethod
    def semistandard_shifted_setvalued(cls, max_entry, mu, nu=(), diagonal_primes=True):
        return cls.semistandard_shifted_marked_setvalued(max_entry, mu, nu, diagonal_primes)

    @classmethod
    def semistandard_shifted_marked_setvalued(cls, max_entry, mu, nu=(), diagonal_primes=True):
        return cls.semistandard_shifted_marked(max_entry, mu, nu, diagonal_primes, True)

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


class ReversePlanePartition(Tableau):

    def weight(self, n):
        ans = [set() for i in range(n)]
        for i, j, v in self:
            for x in v:
                ans[abs(x) - 1].add(j)
        return tuple(len(s) for s in ans)


class MarkedReversePlanePartition(Tableau):

    def weight(self, n):
        ans = [set() for i in range(n)]
        bns = [set() for i in range(n)]
        for i, j, v in self:
            for x in v:
                assert x != 0
                if x > 0:
                    ans[x - 1].add(j)
                if x < 0:
                    bns[-x - 1].add(i)

        return tuple(len(ans[i]) + len(bns[i]) for i in range(n))
