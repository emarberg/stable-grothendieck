from tableaux import Tableau
from partitions import Partition
from vectors import Vector
import pytest
import utils
import time
from collections import defaultdict
import numpy

LSHAPES = {}
RSHAPES = {}

LCACHE = {}
RCACHE = {}

SPLIT_CACHE = {}


class RowVector:

    def __init__(self, row=()):
        self.row = list(row)
        self.size = len(row)

    def __iter__(self):
        return iter(self.row)

    def __getitem__(self, i):
        return self.row[i]

    def __len__(self):
        return self.size

    def __add__(self, other):
        assert self.size == other.size
        return RowVector([self[i] + other[i] for i in range(self.size)])

    def __sub__(self, other):
        assert self.size == other.size
        return RowVector([self[i] - other[i] for i in range(self.size)])

    def __mul__(self, other):
        assert type(other) == int
        return RowVector([self[i] * other for i in range(self.size)])

    def __rmul__(self, other):
        return self * other

    def __rshift__(self, n):
        assert type(n) == int and n >= 0
        return RowVector(n * [0] + self.row[:self.size - n])

    def __eq__(self, other):
        if type(other) == int and other == 0:
            return all(self[i] == 0 for i in range(self.size))
        assert self.size == other.size
        return all(self[i] == other[i] for i in range(self.size))

    def append(self, x):
        self.row.append(x)
        self.size += 1

    def subvector(self, n):
        assert n <= self.size
        return RowVector(self.row[:n])

    def __hash__(self):
        return hash(tuple(self.row))

    def __repr__(self):
        m = max(1 + max([0] + [len(str(i)) for i in self.row]), 5)

        def pad(i):
            s = str(i)
            return s + (m - len(s)) * ' '

        return '[ ' + (''.join([pad(i) for i in self.row])).strip() + ' ]'

    def copy(self):
        return RowVector(self.row.copy())

    def __or__(self, other):
        return RowVector(self.row + other.row)

    @classmethod
    def linear_combination(cls, coefficients, vectors):
        assert len(vectors) > 0
        ans = coefficients[0] * vectors[0]
        for i in range(1, len(vectors)):
            ans += coefficients[i] * vectors[i]
        return ans

    @classmethod
    def print_matrix(cls, vectors):
        print()
        for v in vectors:
            print(v)
        print()

    @classmethod
    def solve(cls, vectors, target, integral=False):
        def find_nonzero(arr, j):
            for i in range(j, len(arr[j])):
                if arr[j][i] != 0:
                    return i

        def cancel(arr, i, j, k):
            a = arr[j][i]
            b = arr[j][k]
            if b != 0:
                lcm = int(numpy.lcm(a, b))
                for t in range(len(arr)):
                    arr[t].row[k] = arr[t][k] * (lcm // b) - arr[t][i] * (lcm // a)

        def kernel(arr):
            return []

        m, n = len(target), len(vectors)
        arr = [v.copy() for v in vectors] + [target.copy()]

        for j in range(n):
            i = find_nonzero(arr, j)
            if i is None:
                continue
            for k in range(m):
                if i != k:
                    cancel(arr, i, j, k)
            for v in range(n + 1):
                arr[v].row[j], arr[v].row[i] = arr[v][i], arr[v][j]

        solution = {}
        for i in range(m - 1, -1, -1):
            for j in range(n + 1):
                if arr[j][i] != 0:
                    if j == n:
                        return []
                    v = arr[n][i]
                    for k in range(j + 1, n):
                        v -= arr[k][i] * solution.get(k, 0)
                    if integral and v % arr[j][i] != 0:
                        return[]
                    solution[j] = v // arr[j][i] if v % arr[j][i] == 0 else v / arr[j][i]
                    break
        s = RowVector([solution.get(i, 0) for i in range(n)])
        return [s] + [s + vec for vec in kernel(arr[:-1])]


def split(nu, mu):
    if (nu, mu) not in SPLIT_CACHE:
        i, j = None, None
        for lam in lshapes(nu, mu):
            term = Partition.rim_terminal_part(left_rim(nu, mu, lam))
            i = max(([] if i is None else [i]) + [a for (a, _) in term])
            j = min(([] if j is None else [j]) + [b for (_, b) in term])
        for lam in rshapes(nu, mu):
            term = Partition.rim_terminal_part(right_rim(nu, mu, lam))
            i = max(([] if i is None else [i]) + [a for (a, _) in term])
            j = min(([] if j is None else [j]) + [b for (_, b) in term])
        skew = Partition.skew(nu, mu, shifted=True)
        # new_nu = Partition.trim(tuple(min(nu[t], j - t - 1) for t in range(i)) + nu[i:])
        # new_mu = Partition.trim(tuple(min(mu[t], j - t - 1) for t in range(i)) + mu[i:])
        while (i + 1, j) in skew:
            i += 1
        rest = {(a, b) for (a, b) in skew if (a <= i or b >= j)}
        split_nu, split_mu = Partition.from_skew(rest, shifted=False)
        diagonal = any(a == b for (a, b) in skew)
        SPLIT_CACHE[(nu, mu)] = (split_nu, split_mu, diagonal)
    return SPLIT_CACHE[(nu, mu)]


def left_rim(nu, mu, lam):
    return Partition.skew(nu, lam, shifted=True)


def right_rim(nu, mu, lam):
    return Partition.skew(lam, mu, shifted=True)


def print_lshape(nu, mu, lam):
    lines = []
    for i in range(1, len(nu) + 1):
        line = (i - 1) * [' '] + Partition.get(mu, i) * ['.']
        line += (Partition.get(lam, i) - Partition.get(mu, i)) * ['*']
        line += (Partition.get(nu, i) - Partition.get(lam, i)) * ['R']
        lines += [' '.join(line)]
    print('\n'.join(reversed(lines)) + '\n')
    print('term:', Partition.rim_terminal_part(left_rim(nu, mu, lam)))
    print()


def print_rshape(nu, mu, lam):
    lines = []
    for i in range(1, len(nu) + 1):
        line = (i - 1) * [' '] + Partition.get(mu, i) * ['.']
        line += (Partition.get(lam, i) - Partition.get(mu, i)) * ['R']
        line += (Partition.get(nu, i) - Partition.get(lam, i)) * ['*']
        lines += [' '.join(line)]
    print('\n'.join(reversed(lines)) + '\n')
    print('term:', Partition.rim_terminal_part(right_rim(nu, mu, lam)))
    print()


def columns(lam, mu):
    return len({j for (i, j) in Partition.skew(lam, mu, True)})


def lvector(nu, mu, p_max):
    ans = LCACHE.get((nu, mu), RowVector())
    if len(ans) >= p_max:
        return ans.subvector(p_max)
    for p in range(len(ans), p_max):
        v = 0
        for lam in lshapes(nu, mu):
            c = sum(lam) - sum(mu)
            d = Tableau.KOG_counts(nu, lam, p) if p == 0 else 2 * Tableau.KOG_counts(nu, lam, p) + Tableau.KOG_counts(nu, lam, p + 1)
            v += 2 ** (len(mu) - c) * (-1) ** (columns(lam, mu) + c) * d
        ans.append(v)
    LCACHE[(nu, mu)] = ans
    return ans


def rvector(nu, mu, p_max):
    ans = RCACHE.get((nu, mu), RowVector())
    if len(ans) >= p_max:
        return ans.subvector(p_max)
    for p in range(len(ans), p_max):
        v = 0
        for lam in rshapes(nu, mu):
            c = sum(nu) - sum(lam)
            v += 2 ** (len(lam) - c) * (-1) ** (columns(nu, lam) + c) * Tableau.KLG_counts(lam, mu, p)
        ans.append(v)
    RCACHE[(nu, mu)] = ans
    return ans


def lshapes(nu, mu=None):
    if nu not in LSHAPES:
        LSHAPES[nu] = {}
        for lam in Partition.decompose_shifted_shape_by_rims(nu):
            for x in Partition.decompose_shifted_shape_by_vertical_strips(lam):
                if len(x) == len(lam):
                    if x not in LSHAPES[nu]:
                        LSHAPES[nu][x] = []
                    LSHAPES[nu][x].append(lam)
    if mu is None:
        return LSHAPES[nu]
    else:
        return LSHAPES[nu].get(mu, [])


def rshapes(nu, mu=None):
    if nu not in RSHAPES:
        RSHAPES[nu] = {}
        for lam in Partition.decompose_shifted_shape_by_vertical_strips(nu):
            if len(lam) == len(nu):
                for x in Partition.decompose_shifted_shape_by_rims(lam):
                    if x not in RSHAPES[nu]:
                        RSHAPES[nu][x] = []
                    RSHAPES[nu][x].append(lam)
    if mu is None:
        return RSHAPES[nu]
    else:
        return RSHAPES[nu].get(mu, [])


def is_reducible(nu, mu):
    if len(nu) == 0:
        return True
    if any(nu[i] == mu[i] for i in range(len(mu))):
        return True
    if any(nu[i] + 1 < mu[i - 1] for i in range(1, len(nu))):
        return True
    if sum(nu) - sum(mu) <= 5:
        return True
    # if Partition.get(nu, 1) <= Partition.get(mu, 1) + 1:
    #    return True

    skew = Partition.skew(nu, mu, shifted=True)
    (i, j) = min(skew, key=lambda ij: (-ij[0], ij[1]))
    while True:
        if (i - 1, j) in skew:
            (i, j) = (i - 1, j)
        elif (i, j + 1) in skew:
            (i, j) = (i, j + 1)
        else:
            break
    if (i, j) != max(skew, key=lambda ij: (-ij[0], ij[1])):
        return True


def deep_dive(nu, mu):
    split_nu, split_mu, diagonal = split(nu, mu)
    print()
    print('TYPE:')
    print()
    print(split_nu, split_mu)
    print()
    print(Partition.printable(split_nu, split_mu, shifted=False))
    print()

    print('nu =', nu, 'mu =', mu)
    print()
    print(Partition.printable(nu, mu, shifted=True))
    print()

    m = sum(nu) - sum(mu) + 2
    target = lvector(nu, mu, m)

    print(target)
    print()

    vectors = []
    for (nu2, mu2) in [
        ((4, 2, 1), (2, 1)),
        ((3, 2, 1), (2, 1)),
    ]:
        v = lvector(nu2, mu2, m)
        w = lvector(nu2, mu2, m) >> 1
        print(v, nu2, mu2)
        print(w)
        print()
        print(Partition.printable(nu2, mu2, shifted=True))
        print()
        vectors += [v, w]
    # for mu2 in Partition.subpartitions(mu, strict=True):
    #     for nu2 in Partition.subpartitions(nu, strict=True):
    #         if (nu, mu) != (nu2, mu2) and Partition.contains(nu2, mu2):
    #             if not Partition.skew(nu2, mu2).issubset(Partition.skew(nu, mu)):
    #                 continue
    #             v = lvector(nu2, mu2, m)
    #             w = lvector(nu2, mu2, m) >> 1
    #             print(v, nu2, mu2)
    #             print(w)
    #             print()
    #             print(Partition.printable(nu2, mu2, shifted=True))
    #             print()
    #             vectors += [v, w]
    print('solutions:', RowVector.solve(vectors, target))


def test_all(n=6):
    types = {}
    success = [0, 0]
    unique = set()
    for nu in Partition.all(n, strict=True):
        for mu in set(lshapes(nu)) | set(rshapes(nu)):
            if is_reducible(nu, mu):
                continue
            triple = split(nu, mu)
            if triple not in types:
                types[triple] = []
            types[triple].append((nu, mu))

    for (split_nu, split_mu, diagonal) in types:
        target = RowVector()
        w = [RowVector() for i in range(11)]

        # print(10 * '\n')
        domain = types[(split_nu, split_mu, diagonal)]
        for (nu, mu) in domain:
            if mu == ():
                continue

            nu1, mu1 = Partition.decrement_one(nu), mu
            nu2, mu2 = Partition.decrement_one(nu), Partition.decrement_one(mu)
            nu3, mu3 = Partition.decrement_one(Partition.decrement_one(nu)), Partition.decrement_one(mu)
            nu4, mu4 = Partition.trim(nu[1:]), Partition.trim(mu[1:])
            nu5, mu5 = Partition.trim(nu[2:]), Partition.trim(mu[2:])

            m = sum(nu) - sum(mu) + 2

            target |= lvector(nu, mu, m)

            w[0] |= lvector(nu1, mu1, m)
            w[1] |= lvector(nu1, mu1, m) >> 1

            w[2] |= lvector(nu2, mu2, m)
            w[3] |= lvector(nu2, mu2, m) >> 1

            w[4] |= lvector(nu3, mu3, m)
            w[5] |= lvector(nu3, mu3, m) >> 1

            w[6] |= lvector(nu4, mu4, m)
            w[7] |= lvector(nu4, mu4, m) >> 1

            w[8] |= lvector(nu5, mu5, m)
            w[9] |= lvector(nu5, mu5, m) >> 1
            w[10] |= lvector(nu5, mu5, m) >> 2
            assert lvector(nu, mu, m) == rvector(nu, mu, m)

        solved = RowVector.solve(w, target)
        # if solved:
        #     solved = solved[0]
        #     while not all(type(t) == int and abs(t) <= 2 for t in solved):
        #         keep = {i for i in range(len(w)) if type(solved[i]) == int and abs(solved[i]) <= 2}
        #         w = [w[i] for i in keep]
        #         test = RowVector.solve(w, target)
        #         if test:
        #             solved = test[0]
        #         else:
        #             break
        if solved: # and all(type(t) == int and abs(t) <= 2 for t in solved):
            solved = solved[0]
            success[0] += 1
            unique.add(solved)
            if len(domain) > 1:
                print(len(domain), ':', solved)
            continue
        else:
            success[1] += 1
            for (nu, mu) in domain:
                if mu == ():
                    continue

                nu1, mu1 = Partition.decrement_one(nu), mu
                nu2, mu2 = Partition.decrement_one(nu), Partition.decrement_one(mu)
                nu3, mu3 = Partition.decrement_one(Partition.decrement_one(nu)), Partition.decrement_one(mu)
                nu4, mu4 = Partition.trim(nu[1:]), Partition.trim(mu[1:])
                nu5, mu5 = Partition.trim(nu[2:]), Partition.trim(mu[2:])

                m = sum(nu) - sum(mu) + 2

                print('nu =', nu, '; mu =', mu)
                print()
                print(Partition.printable(nu, mu, shifted=True))
                print()

                v = [
                    lvector(nu1, mu1, m),
                    lvector(nu1, mu1, m) >> 1,
                    lvector(nu2, mu2, m),
                    lvector(nu2, mu2, m) >> 1,
                    lvector(nu3, mu3, m),
                    lvector(nu3, mu3, m) >> 1,
                    lvector(nu4, mu4, m),
                    lvector(nu4, mu4, m) >> 1,
                    lvector(nu5, mu5, m),
                    lvector(nu5, mu5, m) >> 1,
                    lvector(nu5, mu5, m) >> 2,
                ]
                print('target =')
                print(lvector(nu, mu, m))
                print()
                RowVector.print_matrix(v)
                print()
                print('partial solutions =', RowVector.solve(v, lvector(nu, mu, m), integral=False))
                print()
                # print('solutions so far:')
                # for u in unique:
                #     print('*', u)

        print()
        print('TYPE:')
        print()
        print(split_nu, split_mu, diagonal)
        print()
        print(Partition.printable(split_nu, split_mu, shifted=False))
        print()
        print()
        print(solved)
        print()
        input('')
        print(10 * '\n')
    print(success)


def seconds(t):
    return str(int(1000 * t) / 1000.0)


# def pairs_io(n=5):
#     assert 0 <= n < 256

#     def tobytes(nu, lam, mu):
#         return bytes(nu + (0,) + lam + (0,) + mu + (0,))

#     def frombytes(b, i):
#         nu = []
#         while i < len(b) and b[i] != 0:
#             nu.append(b[i])
#             i += 1
#         i += 1
#         lam = []
#         while i < len(b) and b[i] != 0:
#             lam.append(b[i])
#             i += 1
#         i += 1
#         mu = []
#         while i < len(b) and b[i] != 0:
#             mu.append(b[i])
#             i += 1
#         i += 1
#         return tuple(nu), tuple(lam), tuple(mu), i

#     def read(filename):
#         print()
#         print('Trying to read from `%s`' % filename)
#         t0 = time.time()
#         with open(filename, 'rb') as file:
#             b = bytearray(file.read())
#         print('* Opened file in %s seconds' % seconds(time.time() - t0))
#         ans = defaultdict(list)
#         i = 0
#         while i < len(b):
#             nu, lam, mu, i = frombytes(b, i)
#             ans[(nu, mu)].append(lam)
#         print('* Succeeded in %s seconds' % seconds(time.time() - t0))
#         print()
#         return ans

#     def write(file, ans):
#         t0 = time.time()
#         b = bytearray()
#         for nu, mu in ans:
#             for lam in ans[(nu, mu)]:
#                 b += tobytes(nu, lam, mu)
#         print('* Writing to file `%s`' % file)
#         with open(file, 'wb') as f:
#             f.write(b)
#         print('* Succeeded in %s seconds' % seconds(time.time() - t0))
#         print()

#     directory = '/Users/emarberg/examples/gpq/'
#     lfile = directory + 'rims_then_strips_%s.b' % n
#     rfile = directory + 'strips_then_rims_%s.b' % n

#     try:
#         lhs = read(lfile)
#     except FileNotFoundError:
#         print('* Failed, computing instead')
#         lhs = pairs_lhs(n)
#         write(lfile, lhs)

#     try:
#         rhs = read(rfile)
#     except FileNotFoundError:
#         print('* Failed, computing instead')
#         rhs = pairs_rhs(n)
#         write(rfile, rhs)

#     return lhs, rhs


def test_rims():
    assert set(Partition.rims((), 0)) == {()}
    assert set(Partition.rims((), 1)) == {(), ((1, 1),)}
    assert set(Partition.rims((), 2)) == {(), ((1, 1),), ((1, 1), (1, 2))}
    assert set(Partition.rims((1,), 1)) == {(), ((1, 2),), ((1, 2), (2, 2))}
    assert set(Partition.rims((3, 2))) == {(), ((1, 4),), ((1, 4), (2, 4)), ((3, 3),), ((1, 4), (3, 3)), ((1, 4), (2, 4), (3, 3)),  ((1, 4), (2, 4), (3, 3), (3, 4))} # noqa


def test_KOG():  # noqa
    mu = (6, 4, 1)
    assert set(Tableau.KOG(0, mu)) == {Tableau()}

    assert set(Tableau.KOG(1, ())) == {Tableau({(1, 1): 1})}

    nu = (7, 6, 3, 1)
    assert set(Tableau.KOG_by_shape(5, mu)[nu]) == {
        Tableau({
            (1, 7): 1,
            (2, 6): 1,
            (2, 7): 5,
            (3, 4): 2,
            (3, 5): 4,
            (4, 4): 3,
        }),
        Tableau({
            (1, 7): 1,
            (2, 6): 2,
            (2, 7): 5,
            (3, 4): 2,
            (3, 5): 4,
            (4, 4): 3,
        }),
        Tableau({
            (1, 7): 1,
            (2, 6): 2,
            (2, 7): 5,
            (3, 4): 3,
            (3, 5): 5,
            (4, 4): 4,
        }),
        Tableau({
            (1, 7): 1,
            (2, 6): 2,
            (2, 7): 5,
            (3, 4): 3,
            (3, 5): 4,
            (4, 4): 4,
        }),
        Tableau({
            (1, 7): 1,
            (2, 6): 4,
            (2, 7): 5,
            (3, 4): 1,
            (3, 5): 3,
            (4, 4): 2,
        }),
        Tableau({
            (1, 7): 1,
            (2, 6): 4,
            (2, 7): 5,
            (3, 4): 2,
            (3, 5): 4,
            (4, 4): 3,
        }),
        Tableau({
            (1, 7): 1,
            (2, 6): 4,
            (2, 7): 5,
            (3, 4): 2,
            (3, 5): 3,
            (4, 4): 3,
        })
    }


def test_KLG():  # noqa
    mu = (6, 4, 1)
    assert set(Tableau.KLG(0, mu)) == {Tableau()}

    assert set(Tableau.KLG(1, ())) == {Tableau({(1, 1): 1})}

    nu = (7, 6, 3, 1)
    assert set(Tableau.KLG_by_shape(5, mu)[nu]) == {
        Tableau({
            (1, 7): -1,
            (2, 6): -1,
            (2, 7): 5,
            (3, 4): -2,
            (3, 5): 4,
            (4, 4): 3,
        }),
        Tableau({
            (1, 7): -1,
            (2, 6): -2,
            (2, 7): 5,
            (3, 4): -2,
            (3, 5): 4,
            (4, 4): 3,
        }),
        Tableau({
            (1, 7): -1,
            (2, 6): -2,
            (2, 7): 5,
            (3, 4): -3,
            (3, 5): 5,
            (4, 4): 4,
        }),
        Tableau({
            (1, 7): -1,
            (2, 6): -2,
            (2, 7): 5,
            (3, 4): -3,
            (3, 5): 4,
            (4, 4): 4,
        }),
        Tableau({
            (1, 7): -1,
            (2, 6): 4,
            (2, 7): 5,
            (3, 4): -1,
            (3, 5): 3,
            (4, 4): 2,
        }),
        Tableau({
            (1, 7): -1,
            (2, 6): 4,
            (2, 7): 5,
            (3, 4): -2,
            (3, 5): 4,
            (4, 4): 3,
        }),
        Tableau({
            (1, 7): -1,
            (2, 6): 4,
            (2, 7): 5,
            (3, 4): -2,
            (3, 5): 3,
            (4, 4): 3,
        }),
        Tableau({
            (1, 7): -1,
            (2, 6): 4,
            (2, 7): 5,
            (3, 4): -2,
            (3, 5): 3,
            (4, 4): 2,
        }),
        Tableau({
            (1, 7): -1,
            (2, 6): -2,
            (2, 7): 5,
            (3, 4): -3,
            (3, 5): 4,
            (4, 4): 3,
        })
    }


def test_fast_KOG(): # noqa
    for mu in Partition.all(7, strict=True):
        for p in range(7):
            counts = Tableau.KOG_counts_by_shape(p, mu)
            tabs = Tableau.KOG_by_shape(p, mu)
            print(mu, p)
            print(counts)
            print(tabs)
            print()
            assert set(counts) == set(tabs)
            for nu in counts:
                assert counts[nu] == len(tabs[nu])


def test_fast_KLG(): # noqa
    for mu in Partition.all(7, strict=True):
        for p in range(7):
            counts = Tableau.KLG_counts_by_shape(p, mu)
            tabs = Tableau.KLG_by_shape(p, mu)
            print(mu, p)
            print(counts)
            print(tabs)
            print()
            assert set(counts) == set(tabs)
            for nu in counts:
                print('*', nu)
                assert counts[nu] == len(tabs[nu])


def test_GP_pieri(): # noqa
    for mu in Partition.all(6, strict=True):
        for p in [0, 1, 2]:
            ans = Vector()
            for nu, tabs in Tableau.KOG_by_shape(p, mu).items():
                ans += Vector({nu: utils.beta**(sum(nu) - sum(mu) - p) * len(tabs)})
            f = utils.GP(len(mu) + 1, mu) * utils.GP(len(mu) + 1, (p,))
            assert ans == utils.GP_expansion(f)


def test_GQ_pieri(): # noqa
    for mu in Partition.all(5, strict=True):
        for p in [0, 1, 2]:
            ans = Vector()
            for nu, tabs in Tableau.KLG_by_shape(p, mu).items():
                ans += Vector({nu: utils.beta**(sum(nu) - sum(mu) - p) * len(tabs)})
            f = utils.GQ(len(mu) + 1, mu) * utils.GQ(len(mu) + 1, (p,))
            assert ans == utils.GQ_expansion(f)


@pytest.mark.slow
def test_GP_pieri_slow(): # noqa
    for mu in Partition.all(10, strict=True):
        for p in [1, 2, 3]:
            ans = Vector()
            for nu, tabs in Tableau.KOG_by_shape(p, mu).items():
                ans += Vector({nu: utils.beta**(sum(nu) - sum(mu) - p) * len(tabs)})
            f = utils.GP(len(mu) + 1, mu) * utils.GP(len(mu) + 1, (p,))
            assert ans == utils.GP_expansion(f)


@pytest.mark.slow
def test_GQ_pieri_slow(): # noqa
    for mu in Partition.all(10, strict=True):
        for p in [1, 2, 3]:
            ans = Vector()
            for nu, tabs in Tableau.KLG_by_shape(p, mu).items():
                ans += Vector({nu: utils.beta**(sum(nu) - sum(mu) - p) * len(tabs)})
            f = utils.GQ(len(mu) + 1, mu) * utils.GQ(len(mu) + 1, (p,))
            assert ans == utils.GQ_expansion(f)
