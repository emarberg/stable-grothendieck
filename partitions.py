from collections import defaultdict
import itertools

PARTITIONS = {}


class Partition:

    @classmethod
    def shape(cls, mu):
        ans = set()
        for i, a in enumerate(mu):
            for j in range(a):
                ans.add((i + 1, j + 1))
        return ans

    @classmethod
    def remove_inner_corners(cls, mu):
        rows = []
        for i in range(1, len(mu) + 1):
            try:
                cls.remove_inner_corner(mu, i)
                rows += [i]
            except:
                continue
        for k in range(len(rows) + 1):
            for subset in itertools.combinations(rows, k):
                ans = mu
                for row in subset:
                    ans = cls.remove_inner_corner(ans, row)
                yield ans

    @classmethod
    def remove_inner_corner(cls, mu, row):
        if row == 1 and len(mu) > 1 and mu[1] < mu[0]:
            return (mu[0] - 1,) + mu[1:]
        elif 1 < row < len(mu) and mu[row] < mu[row - 1]:
            ans = list(mu)
            ans[row - 1] -= 1
            return tuple(ans)
        elif row == len(mu) and mu[-1] > 1:
            return mu[:-1] + (mu[-1] - 1,)
        elif row == len(mu) and mu[-1] == 1:
            return mu[:-1]
        raise Exception

    @classmethod
    def remove_shifted_inner_corners(cls, mu):
        rows = []
        for i in range(1, len(mu) + 1):
            try:
                cls.remove_shifted_inner_corner(mu, i)
                rows += [i]
            except:
                continue
        for k in range(len(rows) + 1):
            for subset in itertools.combinations(rows, k):
                ans = mu
                for row in subset:
                    ans = cls.remove_shifted_inner_corner(ans, row)
                yield ans

    @classmethod
    def remove_shifted_inner_corner(cls, mu, row):
        if row == 1 and len(mu) > 1 and mu[1] + 1 < mu[0]:
            return (mu[0] - 1,) + mu[1:]
        elif 1 < row < len(mu) and mu[row] + 1 < mu[row - 1]:
            ans = list(mu)
            ans[row - 1] -= 1
            return tuple(ans)
        elif row == len(mu) and mu[-1] > 1:
            return mu[:-1] + (mu[-1] - 1,)
        elif row == len(mu) and mu[-1] == 1:
            return mu[:-1]
        raise Exception

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
    def trim(cls, mu):
        while mu and mu[-1] == 0:
            mu = mu[:-1]
        return tuple(mu)

    @classmethod
    def sort(cls, mu, trim=False):
        ans = tuple(reversed(sorted(mu)))
        if trim:
            ans = cls.trim(ans)
        return ans

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
    def all(cls, n, max_part=None, max_row=None, strict=False, even_parts=False):
        for i in range(n + 1):
            for mu in cls.generate(i, max_part=max_part, max_row=max_row, strict=strict, even_parts=even_parts):
                yield mu

    @classmethod
    def generate(cls, n, max_part=None, max_row=None, strict=False, even_parts=False):
        if n == 0:
            yield ()
        else:
            if (n, max_part, max_row, strict, even_parts) not in PARTITIONS:
                ans = []
                max_part = n if (max_part is None or max_part > n) else max_part
                max_row = n if (max_row is None or max_row > n) else max_row
                if max_row > 0:
                    for i in range(1, max_part + 1):
                        if even_parts and i % 2 != 0:
                            continue
                        for mu in cls.generate(n - i, i, max_row - 1, strict=strict, even_parts=even_parts):
                            nu = (i,) + mu
                            if not strict or Partition.is_strict_partition(nu):
                                ans.append(nu)
                                yield nu
                PARTITIONS[(n, max_part, strict, even_parts)] = ans
            else:
                for mu in PARTITIONS[(n, max_part, strict, even_parts)]:
                    yield mu

    @classmethod
    def subpartitions(cls, mu, strict=False):

        def _subpartitions(mu, strict):
            if mu:
                for nu in _subpartitions(mu[1:], strict):
                    lb = (nu[0] + (1 if strict else 0)) if (nu and nu[0] > 0) else 0
                    ub = mu[0]
                    for a in range(lb, ub + 1):
                        yield (a,) + nu
            else:
                yield ()

        for nu in _subpartitions(mu, strict):
            yield cls.trim(nu)
