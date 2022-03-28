from cached import cached_value
from words import Word
from partitions import Partition
from collections import defaultdict
from operator import itemgetter
import itertools

FRENCH = True

COUNT_SEMISTANDARD_CACHE = {}
COUNT_SEMISTANDARD_MARKED_CACHE = {}
COUNT_SEMISTANDARD_SHIFTED_MARKED_CACHE = {}
COUNT_SEMISTANDARD_RPP_CACHE = {}
COUNT_SEMISTANDARD_MARKED_RPP_CACHE = {}

STANDARD_CACHE = {}
STANDARD_SHIFTED_MARKED_CACHE = {}

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

KOG_CACHE = {}
KLG_CACHE = {}

KOG_COUNTS = {}
KOG_MAPS = {}

KLG_COUNTS = {}
KLG_MAPS = {}

KOG_COUNTS_HELPER = {}
KLG_COUNTS_HELPER = {}

VALUED_SET_CACHE = {}


def nchoosek(m, k):
    ans = 1
    for i in range(k):
        ans *= m - i
    for i in range(k):
        ans //= i + 1
    return ans


class Tableau:

    def __init__(self, mapping=None):
        mapping = {} if mapping is None else mapping
        if type(mapping) == str:
            mapping = self.decode_string_input(mapping)

        def tuplize(i, j):
            v = mapping[(i, j)]
            if type(v) == tuple:
                assert all(type(i) == int for i in v)
                return tuple(sorted(v))
            else:
                assert type(v) in [int, bool]
                return (int(v),)

        assert all(i > 0 and j > 0 for i, j in mapping)
        self.boxes = {(i, j): tuplize(i, j) for i, j in mapping}
        self._sorting_word = tuple(
            (2 * v if v > 0 else -1 - 2 * v)
            for b in sorted(self.boxes) for v in self.boxes[b]
        )
        self._string_array = None

    @classmethod
    def decode_string_input(cls, s):
        def decode(cell):
            for x in cell.split(','):
                if x:
                    yield -int(x[:-1]) if x[-1] == '\'' else int(x)

        mapping = {}
        s = s.strip().replace('\n', ';')
        for i, row in enumerate(s.split(';')):
            j = 0
            for cell in row.strip().split(' '):
                cell = cell.strip()
                if cell:
                    if cell != '.':
                        mapping[(i + 1, j + 1)] = tuple(decode(cell))
                    j += 1
        return mapping

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

    def tex(self):
        s = defaultdict(list)
        for box in sorted(self.boxes):
            v = ','.join(map(str, self.boxes[box]))
            s[box[0]] += [v]
        s = '\n\\\\\n'.join([' & '.join(s[row]) for row in reversed(sorted(s))])
        return '\\begin{ytableau}\n' + s + '\n\\end{ytableau}'

    def row_reading_word(self):
        return tuple(
            a for key in sorted(self.boxes, key=lambda x: (-x[0], x[1]))
            for a in self.boxes[key]
        )

    @classmethod
    def from_row_reading_word(cls, word):
        subwords = [[]]
        for i, a in enumerate(word):
            if i == 0 or word[i - 1] <= a:
                subwords[-1].append(a)
            else:
                subwords.append([a])
        mapping = {}
        for i, row in enumerate(reversed(subwords)):
            for j, entry in enumerate(row):
                mapping[(i + 1, j + 1)] = entry
        return Tableau(mapping)

    def column_reading_word(self):
        return tuple(
            a for key in sorted(self.boxes, key=lambda x: (x[1], -x[0]))
            for a in self.boxes[key]
        )

    def shifted_reading_word(self):
        x = tuple(-a for a in reversed(self.column_reading_word()) if a < 0)
        y = tuple(a for a in self.row_reading_word() if a > 0)
        return x + y

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

    def is_shifted_column_superstandard(self):
        if any(len(v) > 1 for i, j, v in self):
            return False
        c = 1
        for i, j in sorted(self.boxes, key=itemgetter(1, 0)):
            if i == j and self.get(i, j) != c:
                return False
            if i != j and self.get(i, j) != c:
                return False
            c += 1
        return True

    def is_semistandard(self, diagonal_primes=False):
        def dbl(v):
            return (2 * v) if v > 0 else (-1 - 2 * v)

        for i, j, values in self:
            if i == j and not diagonal_primes and min(values) < 0:
                return False

            maxval = max({dbl(v) for v in values})
            if (i, j + 1) in self:
                minval = min({dbl(v) for v in self.get(i, j + 1, unpack=False)})
                if minval < maxval or (maxval == minval and minval % 2 != 0):
                    return False
            if (i + 1, j) in self:
                minval = min({dbl(v) for v in self.get(i + 1, j, unpack=False)})
                if minval < maxval or (maxval == minval and minval % 2 == 0):
                    return False
        return True

    def is_standard(self):
        if not self.is_semistandard():
            return False
        val = [abs(x) for i, j, value in self for x in value]
        return tuple(sorted(val)) == tuple(range(1, len(val) + 1))

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

    def last_box_in_row(self, row):
        cols = [c for (i, c) in self.boxes if i == row]
        return max(cols) if cols else None

    def last_box_in_column(self, col):
        rows = [r for (r, j) in self.boxes if j == col]
        return max(rows) if rows else None

    def find(self, v):
        return [(i, j) for i, j, values in self if v in values]

    def get(self, i, j, default=None, unpack=True):
        ans = self.boxes.get((i, j), default)
        if unpack and type(ans) == tuple and len(ans) == 1:
            ans = ans[0]
        return ans

    def add(self, i, j, v):
        assert i > 0 and j > 0
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

    def contains_box(self, i, j):
        return (i, j) in self.boxes

    def set(self, i, j, v):
        if (i, j) in self:
            return self.replace(i, j, v)
        else:
            return self.add(i, j, v)

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

    def max_row(self, col=None):
        return max([0] + [i for i, j in self.boxes if (col is None or j == col)])

    def max_column(self, row=None):
        return max([0] + [j for i, j in self.boxes if (row is None or i == row)])

    def shape(self):
        ans = []
        for i in range(1, self.max_row() + 1):
            ans += [len([j for i_, j, v in self if i == i_])]
        return tuple(ans)

    def setvalued_excess(self):
        return sum([len(v) - 1 for _, _, v in self])

    def rpp_excess(self):
        return sum([(i, j) for i, j, v in self if v is self.get(i + 1, j)])

    COLUMN_SPACING_OFFSET = 2

    @property
    def string_array(self):
        if self._string_array is None:
            boxes = {
                k: ','.join([
                    str(-x) + '\'' if x < 0 else str(x)
                    for x in sorted(v, key=lambda x:(abs(x), x))
                ]) for k, v in self.boxes.items()}

            allmax = max([2] + [len(str(boxes[b])) for b in boxes])

            def maximum(j):
                if allmax <= 2:
                    return allmax
                column = [len(str(boxes[b])) for b in boxes if b[1] == j]
                return max(column) if column else 0

            def pad(x, j):
                return str(x) + (self.COLUMN_SPACING_OFFSET + maximum(j) - len(str(x))) * ' '

            array = []
            for i in range(1, self.max_row() + 1):
                row = []
                for j in range(1, self.max_column() + 1):
                    row += [pad(boxes.get((i, j), '.'), j)]
                array += [row]
            self._string_array = array if not FRENCH else list(reversed(array))
        return self._string_array

    def __repr__(self):
        array = self.string_array
        return '\n\n' + '\n'.join([' '.join(line) for line in array]) + '\n\n'

    def weight(self, n):
        ans = n * [0]
        for i, j, v in self:
            for x in v:
                ans[abs(x) - 1] += 1
        return tuple(ans)

    def descent_set(self):
        ans = set()
        rows = {}
        cols = {}
        for i, j, values in self:
            for v in values:
                assert abs(v) not in rows and abs(v) not in cols
                if v > 0:
                    rows[v] = i
                else:
                    cols[-v] = j
        n = len(rows) + len(cols)
        assert set(rows) | set(cols) == set(range(1, n + 1))
        for i in range(1, n):
            if i in rows and i + 1 in cols:
                ans.add(i)
            elif i in rows and i + 1 in rows and rows[i + 1] > rows[i]:
                ans.add(i)
            elif i in cols and i + 1 in cols and cols[i + 1] > cols[i]:
                ans.add(i)
        return ans

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
    def standard(cls, mu, nu=()):  # noqa
        return cls._standard(mu, nu)

    @cached_value(STANDARD_CACHE)
    def _standard(cls, mu, lam):  # noqa
        ans = set()
        if mu == lam:
            ans = {Tableau()}
        elif Partition.contains(mu, lam):
            n = sum(mu) - sum(lam)
            for i in range(len(mu)):
                row, col = (i + 1), mu[i]
                nu = list(mu)
                nu[i] -= 1
                nu = Partition.trim(nu)
                if Partition.is_partition(nu):
                    for tab in cls._standard(nu, lam):
                        ans.add(tab.add(row, col, n))
        return ans

    @classmethod
    def standard_shifted_marked(cls, mu, nu=(), diagonal_primes=False):  # noqa
        return cls._standard_shifted_marked(mu, nu, diagonal_primes)

    @cached_value(STANDARD_SHIFTED_MARKED_CACHE)
    def _standard_shifted_marked(cls, mu, lam, diagonal_primes):  # noqa
        assert Partition.is_strict_partition(mu)
        ans = set()
        if mu == lam:
            ans = {Tableau()}
        elif Partition.contains(mu, lam):
            n = sum(mu) - sum(lam)
            for i in range(len(mu)):
                row, col = (i + 1), (i + mu[i])
                nu = list(mu)
                nu[i] -= 1
                nu = Partition.trim(nu)
                if Partition.is_strict_partition(nu):
                    for tab in cls._standard_shifted_marked(nu, lam, diagonal_primes):
                        ans.add(tab.add(row, col, n))
                        if diagonal_primes or row != col:
                            ans.add(tab.add(row, col, -n))
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

    @cached_value(KOG_CACHE)
    def _KOG_helper(cls, n, cells_filled, cells_left):  # noqa
        assert n >= 0
        if len(cells_left) == 0 and n > 0:
            return []
        elif len(cells_left) > 0 and n == 0:
            return []
        elif len(cells_left) == n == 0:
            return [cells_filled]
        else:
            corners = {
                (i, j) for (i, j) in cells_left
                if (i, j + 1) not in cells_left and (i + 1, j) not in cells_left and
                not any(
                    k1 <= i and j <= l1 and k1 <= k2 and l2 <= l1 and a1 < a2
                    for (k1, l1, a1) in cells_filled
                    for (k2, l2, a2) in cells_filled
                )
            }
            ans = []
            for k in range(1, len(corners) + 1):
                for subset in itertools.combinations(corners, k):
                    new_cells = cells_filled + tuple((i, j, n) for (i, j) in subset)
                    new_left = tuple(sorted(set(cells_left) - set(subset)))
                    ans.extend(cls._KOG_helper(n - 1, new_cells, new_left))
            return ans

    @classmethod
    def KOG(cls, content_max, mu):  # noqa
        for shape in Partition.rims(mu, content_max):
            for tab in cls._KOG_helper(content_max, (), shape):
                yield Tableau({(i, j): a for (i, j, a) in tab})

    @classmethod
    def _add_rim(cls, mu, tab):
        nu = list(mu) + [0]
        for (i, j) in (tab.boxes if type(tab) == Tableau else tab):
            nu[i - 1] += 1
        return Partition.trim(nu)

    @classmethod
    def KOG_by_shape(cls, content_max, mu):  # noqa
        ans = defaultdict(list)
        for tab in cls.KOG(content_max, mu):
            nu = cls._add_rim(mu, tab)
            ans[nu].append(tab)
        return ans

    @cached_value(KOG_COUNTS)
    def KOG_counts(cls, nu, mu, p):  # noqa
        rim = tuple(sorted(Partition.skew(nu, mu, shifted=True)))
        return cls._KOG_count_helper(p, rim)

    @cached_value(KOG_COUNTS_HELPER)
    def _KOG_count_helper(cls, p, rim): # noqa
        if len(rim) == 0:
            return 1 if p == 0 else 0
        if p <= 0:
            return 0
        term = Partition.rim_terminal_part(rim)
        if len(term) == 2 and p == 1:
            return 0
        elif len(term) == 2:
            left = tuple(a for a in rim if a not in term)
            x = cls._KOG_count_helper(p, left)
            y = cls._KOG_count_helper(p - 1, left)
            z = cls._KOG_count_helper(p - 2, left)
            c = 2 if p > 2 else 1
            return x + (1 + c) * y + c * z
        else:
            left = tuple(a for a in rim if a != term[-1])
            x = cls._KOG_count_helper(p, left)
            y = cls._KOG_count_helper(p - 1, left)
            c = 2 if p > 1 else 1
            if len(term) == 1:
                return c * (x + y)
            if len(term) == 3 and term[0][0] != term[-1][0] and term[0][1] != term[-1][1]:
                return x + y
            if len(term) == 3:
                return y

    @cached_value(KOG_MAPS)
    def KOG_counts_by_shape(cls, content_max, mu): # noqa
        ans = {}
        for rim in Partition.rims(mu, content_max):
            count = cls._KOG_count_helper(content_max, rim)
            if count:
                nu = cls._add_rim(mu, rim)
                ans[nu] = count
        return ans

    @cached_value(KLG_CACHE)
    def _KLG_helper(cls, n, cells_filled, cells_left):  # noqa
        if len(cells_left) == 0 and n > 0:
            return []
        elif len(cells_left) > 0 and n == 0:
            return []
        elif len(cells_left) == n == 0:
            return [cells_filled]
        else:
            corners = {
                (i, j) for (i, j) in cells_left
                if (i, j + 1) not in cells_left and (i + 1, j) not in cells_left and
                (i != j or n > 0) and
                not any(
                    k <= i and j <= l and a < 0
                    for (k, l, a) in cells_filled
                ) and (n < 0 or not any(
                    i <= k and l <= j
                    for (k, l, _) in cells_filled
                ))
            }
            ans = []
            s = 1 if (n < 0 and not any(a == -n for (_, _, a) in cells_filled)) else 0
            for k in range(s, len(corners) + 1):
                for subset in itertools.combinations(corners, k):
                    new_cells = cells_filled + tuple((i, j, n) for (i, j) in subset)
                    new_left = tuple(sorted(set(cells_left) - set(subset)))
                    new_n = -n if n > 0 else -n - 1
                    ans.extend(cls._KLG_helper(new_n, new_cells, new_left))
            return ans

    @classmethod
    def KLG(cls, content_max, mu):  # noqa
        for shape in Partition.rims(mu, content_max + 1):
            for tab in cls._KLG_helper(content_max, (), shape):
                yield Tableau({(i, j): a for (i, j, a) in tab})

    @classmethod
    def KLG_by_shape(cls, content_max, mu):  # noqa
        ans = defaultdict(list)
        for tab in cls.KLG(content_max, mu):
            nu = cls._add_rim(mu, tab)
            ans[nu].append(tab)
        return ans

    @cached_value(KLG_COUNTS_HELPER)
    def _KLG_count_helper(cls, p, rim): # noqa
        if len(rim) == 0:
            return 1 if p == 0 else 0
        if p <= 0:
            return 0
        term = Partition.rim_terminal_part(rim)
        if len(term) == 2:
            left = tuple(a for a in rim if a not in term)
            x = cls._KLG_count_helper(p, left)
            y = cls._KLG_count_helper(p - 1, left)
            z = cls._KLG_count_helper(p - 2, left)
            if p > 2:
                return x + 3 * y + 2 * z
            elif p == 2 and term[0][0] == term[0][1]:
                return 1
            elif p == 2:
                return x + 3 * y + 2 * z
            elif p == 1:
                if term[0][0] == term[0][1] == term[1][0]:
                    return 0
                elif term[0][0] == term[0][1] == term[1][1]:
                    return 1
                else:
                    return x + y
        else:
            left = tuple(a for a in rim if a != term[-1])
            x = cls._KLG_count_helper(p, left)
            y = cls._KLG_count_helper(p - 1, left)
            if len(term) == 1:
                if p == 1 and term[0][0] == term[0][1]:
                    return 1
                else:
                    return 2 * x + 2 * y
            if len(term) == 3 and term[0][0] != term[-1][0] and term[0][1] != term[-1][1]:
                return x + y
            if len(term) == 3:
                return y

    @cached_value(KLG_MAPS)
    def KLG_counts_by_shape(cls, content_max, mu): # noqa
        ans = {}
        for rim in Partition.rims(mu, content_max + 1):
            count = cls._KLG_count_helper(content_max, rim)
            if count:
                nu = cls._add_rim(mu, rim)
                ans[nu] = count
        return ans

    @cached_value(KLG_COUNTS)
    def KLG_counts(cls, nu, mu, p):  # noqa
        rim = tuple(sorted(Partition.skew(nu, mu, shifted=True)))
        return cls._KLG_count_helper(p, rim)

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

    @classmethod
    def from_shifted_growth_diagram(cls, growth, edges, corners):
        def shdiff(nu, lam):
            return next(iter(Partition.skew(nu, lam, shifted=True)))

        p, q = Tableau(), Tableau()
        n, m = len(growth) - 1, len(growth[0]) - 1 if growth else 0

        for i in range(1, n + 1):
            mu, nu = growth[i][m], growth[i - 1][m]
            for a, b in Partition.skew(mu, nu, shifted=True):
                p = p.add(a, b, i)

        for i in range(1, m + 1):
            mu, nu = growth[n][i], growth[n][i - 1]
            v = -i if edges[n][i] else i
            j = corners[n][i]
            if mu != nu:
                a, b = shdiff(mu, nu)
                q = q.add(a, b, v)
            elif v < 0:
                x = q.last_box_in_column(j)
                q = q.add(x, j, v)
            else:
                x = q.last_box_in_row(j)
                q = q.add(j, x, v)

        return p, q


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


class ValuedSetTableau:

    def __hash__(self):
        return hash((self.tableau, self.grouping))

    def is_semistandard(self, order=None):
        if order is None:
            ofn = (lambda x: -2 - 2 * x if x < 0 else 2 * x - 1)
        if type(order) in [list, tuple]:
            order = {a: i for (i, a) in enumerate(order)}
        if type(order) == dict:
            ofn = (lambda x: order[x])
        rowint = {(i, j) for (i, j) in self.tableau.boxes if (i + 1, j) in self.tableau.boxes}
        colint = {(i, j) for (i, j) in self.tableau.boxes if (i, j + 1) in self.tableau.boxes}
        for (i, j) in rowint:
            p = ofn(self.tableau.get(i, j))
            q = ofn(self.tableau.get(i + 1, j))
            if p > q or (p == q and p % 2 != 0):
                return False
        for (i, j) in colint:
            p = ofn(self.tableau.get(i, j))
            q = ofn(self.tableau.get(i, j + 1))
            if p > q or (p == q and p % 2 == 0):
                return False
        return True

    @classmethod
    def is_valid_position(cls, tab, grp, i, j):
        assert (i, j) in tab.boxes
        grp = grp.get(i, j) if type(grp) != bool else grp
        if not grp:
            return True
        if tab.get(i, j) > 0:
            return tab.get(i, j) == tab.get(i, j + 1)
        else:
            return tab.get(i, j) == tab.get(i + 1, j)

    def __init__(self, entry_mapping=None, group_mapping=None):
        self.tableau = entry_mapping if type(entry_mapping) == Tableau else Tableau(entry_mapping)
        self.grouping = group_mapping if type(entry_mapping) == Tableau else Tableau(group_mapping)
        assert all(self.is_valid_position(self.tableau, self.grouping, i, j) for (i, j) in self.tableau.boxes)
        self._strval = None

    def get_verticals(self, v):
        assert v < 0
        vertical_ends = sorted([(i, j) for (i, j) in self.tableau.boxes if self.tableau.get(i, j) == v and not self.grouping.get(i, j)], key=lambda ij: (-ij[1], ij[0]))
        vertical_starts = []
        for i, j in vertical_ends:
            while self.grouping.get(i - 1, j) and self.tableau.get(i - 1, j) == self.tableau.get(i, j) < 0:
                i -= 1
            vertical_starts.append((i, j))
        return vertical_starts, vertical_ends

    def get_horizontals(self, v):
        assert v > 0
        horizontal_ends = sorted([(i, j) for (i, j) in self.tableau.boxes if self.tableau.get(i, j) == v and not self.grouping.get(i, j)], key=lambda ij: (ij[0], -ij[1]))
        horizontal_starts = []
        for i, j in horizontal_ends:
            while self.grouping.get(i, j - 1) and self.tableau.get(i, j - 1) == self.tableau.get(i, j) > 0:
                j -= 1
            horizontal_starts.append((i, j))
        return horizontal_starts, horizontal_ends

    def forward_transition(self, value):
        tab, grp = self.tableau.boxes.copy(), self.grouping.boxes.copy()
        vertical_starts, vertical_ends = self.get_verticals(-value - 1)
        horizontal_starts, horizontal_ends = self.get_horizontals(value)

        for i in range(len(vertical_starts)):
            (vx1, vy1), (vx2, vy2) = vertical_starts[i], vertical_ends[i]
            for j in range(len(horizontal_starts)):
                (hx1, hy1), (hx2, hy2) = horizontal_starts[j], horizontal_ends[j]
                if (vx1, vy1) == (hx1 + 1, hy1):
                    if hy1 < hy2:
                        hy1 += 1
                        vx1 -= 1
                    else:
                        hx1 = vx2
                        hx2 = vx2
                        vx1 -= 1
                        vx2 -= 1
                elif (vx2, vy2) == (hx2, hy2 + 1):
                    if vx1 < vx2:
                        hy2 += 1
                        vx2 -= 1
                    else:
                        vy1 = hy1
                        vy2 = hy1
                        hy1 += 1
                        hy2 += 1
                horizontal_starts[j], horizontal_ends[j] = (hx1, hy1), (hx2, hy2)
            vertical_starts[i], vertical_ends[i] = (vx1, vy1), (vx2, vy2)

        for (vx1, vy1), (vx2, vy2) in zip(vertical_starts, vertical_ends):
            for a in range(vx1, vx2 + 1):
                for b in range(vy1, vy2 + 1):
                    tab[a, b] = -value - 1
                    grp[a, b] = 1
            grp[vx2, vy2] = 0

        for (hx1, hy1), (hx2, hy2) in zip(horizontal_starts, horizontal_ends):
            for a in range(hx1, hx2 + 1):
                for b in range(hy1, hy2 + 1):
                    tab[a, b] = value
                    grp[a, b] = 1
            grp[hx2, hy2] = 0

        return ValuedSetTableau(Tableau(tab), Tableau(grp))

    def middle_transition(self, value):
        tab, grp = self.tableau.boxes.copy(), self.grouping.boxes.copy()

        h1 = [(hx1, hy1, hy2, 0) for ((hx1, hy1), (hx2, hy2)) in zip(*self.get_horizontals(value))]
        h2 = [(hx1, hy1, hy2, 1) for ((hx1, hy1), (hx2, hy2)) in zip(*self.get_horizontals(value + 1))]
        hh = sorted(h1 + h2, key=lambda t: (-t[0], t[1]))
        hh_by_rows = []
        for tup in hh:
            if not hh_by_rows or hh_by_rows[-1][-1][0] != tup[0]:
                hh_by_rows.append([])
            hh_by_rows[-1].append(tup)

        one_row_groups = []
        two_row_groups = []

        for i in range(len(hh_by_rows)):
            currow = hh_by_rows[i]
            nextrow = hh_by_rows[i + 1] if i + 1 < len(hh_by_rows) else None
            for (x, y1, y2, f) in currow:
                if nextrow and nextrow[0][0] == x - 1 and nextrow[0][1] <= y2:
                    row2, start2, stop2 = x, y1, currow[-1][2]
                    row1, start1, stop1 = nextrow[0][:3]
                    while hh_by_rows[i + 1]:
                        _, z1, z2, _ = hh_by_rows[i + 1][0]
                        if z1 <= stop2:
                            stop1 = z2
                            hh_by_rows[i + 1] = hh_by_rows[i + 1][1:]
                        else:
                            break
                    two_row_groups.append((row1, start1, stop1, row2, start2, stop2))
                    break
                else:
                    if not one_row_groups or one_row_groups[-1][-1][0][0] != x:
                        one_row_groups.append([0, 0, []])
                    one_row_groups[-1][f] += 1
                    one_row_groups[-1][-1].append((x, y1, y2))

        v1 = [(vy1, vx1, vx2, 0) for ((vx1, vy1), (vx2, vy2)) in zip(*self.get_verticals(-value))]
        v2 = [(vy1, vx1, vx2, 1) for ((vx1, vy1), (vx2, vy2)) in zip(*self.get_verticals(-value - 1))]
        vv = sorted(v1 + v2, key=lambda t: (-t[0], t[1]))
        vv_by_cols = []
        for tup in vv:
            if not vv_by_cols or vv_by_cols[-1][-1][0] != tup[0]:
                vv_by_cols.append([])
            vv_by_cols[-1].append(tup)

        one_col_groups = []
        two_col_groups = []

        for i in range(len(vv_by_cols)):
            currcol = vv_by_cols[i]
            nextcol = vv_by_cols[i + 1] if i + 1 < len(vv_by_cols) else None
            for (y, x1, x2, f) in currcol:
                if nextcol and nextcol[0][0] == y - 1 and nextcol[0][1] <= x2:
                    col2, start2, stop2 = y, x1, currcol[-1][2]
                    col1, start1, stop1 = nextcol[0][:3]
                    while vv_by_cols[i + 1]:
                        _, z1, z2, _ = vv_by_cols[i + 1][0]
                        if z1 <= stop2:
                            stop1 = z2
                            vv_by_cols[i + 1] = vv_by_cols[i + 1][1:]
                        else:
                            break
                    two_col_groups.append((col1, start1, stop1, col2, start2, stop2))
                    break
                else:
                    if not one_col_groups or one_col_groups[-1][-1][0][0] != y:
                        one_col_groups.append([0, 0, []])
                    one_col_groups[-1][f] += 1
                    one_col_groups[-1][-1].append((y, x1, x2))

        for p, q, g in one_row_groups:
            assert len(g) == p + q
            for i in range(len(g)):
                row, col1, col2 = g[i]
                for y in range(col1, col2 + 1):
                    tab[row, y] = (value + 1) if i < q else value
                    grp[row, y] = 1
                grp[row, col2] = 0

        for row1, start1, stop1, row2, start2, stop2 in two_row_groups:
            for y in range(start1, stop1 + 1):
                tab[row1, y] = value + 1
            for y in range(start2, stop2 + 1):
                tab[row2, y] = value
            for y in range(max(start1, start2), min(stop1, stop2)):
                grp[row1, y], grp[row2, y] = grp[row2, y], grp[row1, y]

        for p, q, g in one_col_groups:
            assert len(g) == p + q
            for i in range(len(g)):
                col, row1, row2 = g[i]
                for x in range(row1, row2 + 1):
                    tab[x, col] = (-value - 1) if i < q else -value
                    grp[x, col] = 1
                grp[row2, col] = 0

        for col1, start1, stop1, col2, start2, stop2 in two_col_groups:
            for x in range(start1, stop1 + 1):
                tab[x, col1] = -value - 1
            for x in range(start2, stop2 + 1):
                tab[x, col2] = -value
            for x in range(max(start1, start2), min(stop1, stop2)):
                grp[x, col1], grp[x, col2] = grp[x, col2], grp[x, col1]

        return ValuedSetTableau(Tableau(tab), Tableau(grp))

    def backward_transition(self, value):
        tab, grp = self.tableau.boxes.copy(), self.grouping.boxes.copy()
        vertical_starts, vertical_ends = self.get_verticals(-value)
        horizontal_starts, horizontal_ends = self.get_horizontals(value + 1)

        for i in range(len(vertical_starts) - 1, -1, -1):
            (vx1, vy1), (vx2, vy2) = vertical_starts[i], vertical_ends[i]
            for j in range(len(horizontal_starts) - 1, -1, -1):
                (hx1, hy1), (hx2, hy2) = horizontal_starts[j], horizontal_ends[j]
                if (vx1, vy1) == (hx1, hy1 - 1):
                    if vx1 < vx2:
                        vx1 += 1
                        hy1 -= 1
                    else:
                        print('here?', (vx1, vy1), (vx2, vy2), (hx1, hy1), (hx2, hy2))
                        print(self, self.tableau, self.grouping, horizontal_ends)
                        vy1 = hy2
                        vy2 = hy2
                        hy1 -= 1
                        hy2 -= 1
                elif (vx2, vy2) == (hx2 - 1, hy2):
                    if hy1 < hy2:
                        vx2 += 1
                        hy2 -= 1
                    else:
                        hx1 = vx1
                        hx2 = vx1
                        vx1 += 1
                        vx2 += 1
                horizontal_starts[j], horizontal_ends[j] = (hx1, hy1), (hx2, hy2)
            vertical_starts[i], vertical_ends[i] = (vx1, vy1), (vx2, vy2)

        for (vx1, vy1), (vx2, vy2) in zip(vertical_starts, vertical_ends):
            for a in range(vx1, vx2 + 1):
                for b in range(vy1, vy2 + 1):
                    tab[a, b] = -value
                    grp[a, b] = 1
            grp[vx2, vy2] = 0

        for (hx1, hy1), (hx2, hy2) in zip(horizontal_starts, horizontal_ends):
            for a in range(hx1, hx2 + 1):
                for b in range(hy1, hy2 + 1):
                    tab[a, b] = value + 1
                    grp[a, b] = 1
            grp[hx2, hy2] = 0

        return ValuedSetTableau(Tableau(tab), Tableau(grp))

    def transition(self, index):
        ans = self.forward_transition(index).middle_transition(index).backward_transition(index)
        mapping = {index: index + 1, -index: -index - 1, index + 1: index, -index - 1: -index}
        boxes = {}
        for i, j in ans.tableau.boxes:
            v = ans.tableau.get(i, j)
            boxes[i, j] = mapping.get(v, v)
        return ValuedSetTableau(Tableau(boxes), ans.grouping)

    @classmethod
    def all(cls, max_entry, mu, nu=(), diagonal_nonprimes=True):
        return cls._all(max_entry, mu, nu, diagonal_nonprimes)

    @cached_value(VALUED_SET_CACHE)
    def _all(cls, max_entry, mu, nu, diagonal_nonprimes):  # noqa
        ans = set()
        tabs = Tableau.semistandard_shifted(max_entry, mu, nu, diagonal_nonprimes)
        for tab in tabs:
            grp = {(i, j): int(cls.is_valid_position(tab, True, i, j)) for (i, j) in tab.boxes}
            g = sorted(grp)
            e = sum(grp.values())
            for v in range(2**e):
                next_grp = {}
                for i, j in g:
                    if grp[i, j]:
                        next_grp[i, j] = v % 2
                        v = v // 2
                    else:
                        next_grp[i, j] = 0
                next_grp = Tableau(next_grp)
                ans.add(cls(tab, next_grp))
        return ans

    def compute_string_array(self):
        boxes = {
            k: ','.join([
                str(-x) + '\'' if x < 0 else str(x)
                for x in sorted(v, key=lambda x:(abs(x), x))
            ]) for k, v in self.tableau.boxes.items()}

        allmax = max([2] + [len(str(boxes[b])) for b in boxes])
        column_spacing_offset = 0

        def pad(x, j):
            return str(x) + (column_spacing_offset + allmax - len(str(x))) * ' '

        array = []
        for i in range(1, self.tableau.max_row() + 1):
            row = []
            for j in range(1, self.tableau.max_column() + 1):
                row += [pad(boxes.get((i, j), '.'), j)]
            array += [row]
        return array, allmax

    def in_same_group(self, i, j, x, y):
        if i == x and y == j + 1:
            return bool(self.grouping.get(i, j, default=False)) and self.tableau.get(i, j) > 0
        if i == x and y + 1 == j:
            return self.in_same_group(x, y, i, j)
        if i + 1 == x and y == j:
            return bool(self.grouping.get(i, j, default=False)) and self.tableau.get(i, j) < 0
        if i == x + 1 and y == j:
            return self.in_same_group(x, y, i, j)
        return False

    def __lt__(self, other):
        return self.tableau < other.tableau and self.grouping < other.grouping

    def __eq__(self, other):
        return type(other) == type(self) and self.tableau == other.tableau and self.grouping == other.grouping

    def __repr__(self):
        if self._strval is None:
            array, column_width = self.compute_string_array()
            m, n = len(array), (len(array[0]) if array else 0)
            boxes = self.tableau.boxes
            newarray = [(2 * n + 1) * [" "] for _ in range(2 * m + 1)]

            for i in range(2 * m + 1):
                for j in range(2 * n + 1):
                    a, b = i // 2 + 1, j // 2 + 1
                    if i % 2 != 0 and j % 2 != 0:
                        newarray[i][j] = array[i // 2][j // 2]
                    elif i % 2 == 0 and j % 2 != 0:
                        bordering = (a, b) in boxes or (a - 1, b) in boxes
                        newarray[i][j] = column_width * ('─' if bordering else ' ')
                    elif i % 2 != 0 and j % 2 == 0:
                        bordering = (a, b) in boxes or (a, b - 1) in boxes
                        newarray[i][j] = '│' if bordering else ' '
                    else:
                        if ((a - 1, b) in boxes and (a, b - 1) in boxes) or ((a, b) in boxes and (a - 1, b - 1) in boxes):
                            newarray[i][j] = '┼'
                        elif (a, b) in boxes and (a, b - 1) in boxes:
                            newarray[i][j] = '┴'
                        elif (a - 1, b) in boxes and (a - 1, b - 1) in boxes:
                            newarray[i][j] = '┬'
                        elif (a, b) in boxes and (a - 1, b) in boxes:
                            newarray[i][j] = '├'
                        elif (a, b - 1) in boxes and (a - 1, b - 1) in boxes:
                            newarray[i][j] = '┤'
                        elif (a, b) in boxes:
                            newarray[i][j] = '└'
                        elif (a - 1, b) in boxes:
                            newarray[i][j] = '┌'
                        elif (a, b - 1) in boxes:
                            newarray[i][j] = '┘'
                        elif (a - 1, b - 1) in boxes:
                            newarray[i][j] = '┐'

            array = newarray
            for (a, b) in boxes:
                i, j = 2 * a - 1, 2 * b - 1
                if not self.grouping.get(a, b):
                    continue
                if self.tableau.get(a, b) > 0:
                    array[i][j + 1] = ' '
                    c = array[i + 1][j + 1]
                    array[i + 1][j + 1] = '┴' if c == '┼' else '─'
                    c = array[i - 1][j + 1]
                    array[i - 1][j + 1] = '┬' if c == '┼' else '─'
                else:
                    array[i + 1][j] = ' ' * column_width
                    c = array[i + 1][j + 1]
                    array[i + 1][j + 1] = '├' if c == '┼' else '│'
                    c = array[i + 1][j - 1]
                    array[i + 1][j - 1] = '┤' if c == '┼' else '│'

            array = array if not FRENCH else list(reversed(array))
            self._strval = '\n\n' + '\n'.join([''.join(line) for line in array]) + '\n\n'
        return self._strval

