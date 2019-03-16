import itertools
from symmetric import SymmetricPolynomial, SymmetricMonomial
from words import Word
from vectors import Vector


INVOLUTION_WORDS = {(): {()}}
SIGNED_REDUCED_WORDS = {(): {()}}
EVEN_SIGNED_REDUCED_WORDS = {(): [()]}
EVEN_SIGNED_REDUCED_COUNTS = {(): 1}
EVEN_SIGNED_INVOLUTION_WORDS = {(): [()]}

atoms_b_cache = {}
atoms_d_cache = {}


class Permutation:

    """
    Simple class for (signed) permutatations
    """

    def __init__(self, *oneline):
        if len(oneline) == 1 and type(oneline[0]) != int:
            oneline = tuple(oneline[0])
        while oneline and oneline[-1] == len(oneline):
            oneline = oneline[:-1]
        assert set(range(1, len(oneline) + 1)) == set(abs(i) for i in oneline)

        self._oneline = tuple(oneline)
        self._rank = len(oneline)
        # cached fields
        self._rdes = None
        self._ldes = None
        self._len = None

    @property
    def oneline(self):
        return self._oneline

    @property
    def rank(self):
        return self._rank

    def __repr__(self):
        # return 'Permutation(' + ', '.join([repr(i) for i in self.oneline]) + ')'
        return str(self)

    def __str__(self):
        s = []
        for i in self.oneline:
            s += [str(abs(i))]
            if i < 0:
                s += ['\u0305']
        if s:
            return ''.join(s)
        else:
            return '1'

    @classmethod
    def all(cls, n, signed=False):
        for args in itertools.permutations(range(1, n + 1)):
            if not signed:
                yield Permutation(*args)
            else:
                for v in range(2**n):
                    oneline = []
                    for i in range(n):
                        oneline.append(args[i] * (-1) ** (v % 2))
                        v = v // 2
                    yield Permutation(*oneline)

    @classmethod
    def involutions(cls, n, signed=False, fpf=False):
        if n % 2 != 0 and fpf:
            return

        def create(delta, fixed):
            w = Permutation()
            numbers = [i for i in range(1, n + 1) if i not in fixed]
            for i in range(len(delta)):
                w *= cls.transposition(numbers[0], numbers[delta[i]])
                numbers = numbers[1:delta[i]] + numbers[delta[i] + 1:]
                print(numbers)
            return w

        if not signed:
            for k in [0] if fpf else range(n, -1, -2):
                for fixed in itertools.combinations(range(1, n + 1), k):
                    queue = [tuple(range(n - k - 1, 0, -2))]
                    while queue:
                        delta, queue = queue[0], queue[1:]
                        yield create(delta, fixed)
                        queue += [
                            delta[:i] + (delta[i] - 1,) + delta[i + 1:]
                            for i in range(len(delta))
                            if delta[i] > 1
                        ]
        else:
            for w in cls.involutions(n, False):
                a = [i for i in range(1, n + 1) if i <= w(i)]
                for k in range(len(a) + 1):
                    for conjugators in itertools.combinations(a, k):
                        v = w
                        for i in conjugators:
                            s = cls.reflection_s(i, i)
                            v = v * s if v(i) == i else s * v * s
                        yield v

    @classmethod
    def fpf_involutions(cls, n):
        return cls.involutions(n, False, True)

    def fpf_involution_shape(self):
        cycles = [(a, self(a)) for a in range(1, self.rank + 1) if a < self(a)]
        filtered = sorted([(a, b) for a, b in cycles if any(i < b and a < j < b for i, j in cycles)])
        r = len(filtered)
        if r == 0:
            return ()
        phi = [a for a, b in filtered]
        top = [b for a, b in filtered]
        n = top[0] - 1
        if phi[-1] >= n + 1 or top != [n + i + 1 for i in range(r)]:
            return None
        return tuple(n - a for a in phi)

    def is_fpf_grassmannian(self):
        if any(self(i) < 0 or self(i) == i or self(self(i)) != i for i in range(1, 1 + self.rank)):
            return False
        return self.fpf_involution_shape() is not None

    def is_involution(self):
        return self == self.inverse()

    def is_identity(self):
        return self == self.identity()

    @classmethod
    def symplectic_hecke_words(cls, n, length_bound=-1):
        assert n % 2 == 0
        for level in cls.symplectic_hecke_levels(n, length_bound):
            for pi, w in level:
                yield w

    @classmethod
    def symplectic_hecke_levels(cls, n, length_bound=-1):
        assert n % 2 == 0
        a = cls.identity()
        for i in range(1, n, 2):
            a *= cls.s_i(i)
        start = (a, ())
        level = {start}
        while level:
            next_level = set()
            yield level
            for pi, w in level:
                for i in range(n):
                    s = Permutation.s_i(i)
                    if pi(i) != i + 1:
                        sigma = s % pi % s
                        next_level.add((sigma, w + (i,)))
            level = next_level
            if length_bound == 0:
                break
            length_bound -= 1

    def get_symplectic_hecke_words(self, length_bound):
        for level in self.symplectic_hecke_levels(self.rank, length_bound):
            for pi, w in level:
                if self == pi:
                    yield w

    @classmethod
    def involution_hecke_words(cls, n, length_bound=-1):
        for level in cls.involution_hecke_levels(n, length_bound):
            for pi, w in level:
                yield w

    @classmethod
    def involution_hecke_levels(cls, n, length_bound=-1):
        start = (cls.identity(), ())
        level = {start}
        while level:
            next_level = set()
            yield level
            for pi, w in level:
                for i in range(n):
                    s = Permutation.s_i(i)
                    sigma = s % pi % s
                    next_level.add((sigma, w + (i,)))
            level = next_level
            if length_bound == 0:
                break
            length_bound -= 1

    def get_involution_hecke_words(self, length_bound):
        for level in self.involution_hecke_levels(self.rank, length_bound):
            for pi, w in level:
                if self == pi:
                    yield w

    def filter_involution_hecke_words(self, length_bound):
        for level in self.involution_hecke_levels(self.rank, length_bound):
            for w in [w for pi, w in level if self == pi]:
                print(w, ':', Word.quasisymmetrize(w, Word.decreasing_zeta))
            print()

    @classmethod
    def _symmetrize(cls, vector, n, check=True):
        def sort(t):
            return tuple(reversed(sorted(t)))

        ans = SymmetricPolynomial({
            SymmetricMonomial(n, alpha): max(vector[a] for a in vector if sort(a) == sort(alpha))
            for alpha in vector if sort(alpha) == alpha
        })

        if check:
            assert all(
                vector.dictionary[sort(alpha)] == vector.dictionary[alpha]
                for alpha in vector.dictionary
            )
        return ans

    def signed_involution_stable_grothendieck(self, degree_bound):
        ans = Vector()
        ell = self.involution_length
        for w in self.get_involution_hecke_words(degree_bound):
            ans += (-1)**(len(w) - ell) * Word.quasisymmetrize(w, Word.unimodal_zeta)
        return self._symmetrize(ans, degree_bound)

    def symplectic_stable_grothendieck(self, degree_bound):
        ans = Vector()
        for w in self.get_symplectic_hecke_words(degree_bound):
            ans += Word.quasisymmetrize(w, Word.decreasing_zeta)
        return self._symmetrize(ans, degree_bound)

    def involution_stable_grothendieck(self, degree_bound):
        ans = Vector()
        for w in self.get_involution_hecke_words(degree_bound):
            ans += Word.quasisymmetrize(w, Word.decreasing_zeta)
        return self._symmetrize(ans, degree_bound)

    def marked_stable_grothendieck(self, degree_bound):
        ans = Vector()
        ell = self.involution_length
        for w in self.get_marked_hecke_words(degree_bound):
            ans += (-1)**(len(w) - ell) * Word.quasisymmetrize(w, Word.decreasing_zeta)

        def sort(t):
            return tuple(reversed(sorted(t)))

        base = Vector({alpha: ans[alpha] for alpha in ans if sort(alpha) == alpha})
        left = ans - base
        rows = [base]
        while left:
            r = {}
            for alpha in left:
                if sort(alpha) not in r:
                    r[sort(alpha)] = (alpha, left[alpha])
            base = Vector({alpha: c for (alpha, c) in r.values()})
            left = left - base
            rows.append(base)

        return self._symmetrize(ans, degree_bound, check=False), rows

    @property
    def reduced_word(self):
        return self.get_reduced_word()

    def get_reduced_word(self):
        if self.left_descent_set:
            i = min(self.left_descent_set)
            s = Permutation.s_i(i)
            return (i,) + (s * self).get_reduced_word()
        else:
            return ()

    @property
    def reduced_words(self):
        return self.get_reduced_words()

    def get_reduced_words(self):
        w = self
        oneline = w.oneline
        if oneline not in SIGNED_REDUCED_WORDS:
            words = set()
            for i in w.right_descent_set:
                s = Permutation.s_i(i)
                words |= {e + (i,) for e in (w * s).get_reduced_words()}
            SIGNED_REDUCED_WORDS[oneline] = words
        return SIGNED_REDUCED_WORDS[oneline]

    @property
    def involution_words(self):
        return self.get_involution_words()

    def get_involution_words(self):
        w = self
        assert w.inverse() == w
        oneline = w.oneline
        if oneline not in INVOLUTION_WORDS:
            INVOLUTION_WORDS[oneline] = {word for a in w.get_atoms() for word in a.get_reduced_words()}
        return INVOLUTION_WORDS[oneline]

    def __call__(self, i):
        if 0 < i <= self.rank:
            return self.oneline[i - 1]
        elif -self.rank <= i < 0:
            return -self.oneline[-i - 1]
        else:
            return i

    def __hash__(self):
        return hash(self.oneline)

    def reduce(self):
        newline = self.oneline
        while newline and newline[-1] == len(newline):
            newline = newline[:-1]
        return Permutation(*newline)

    def pair(self):
        n = self.rank
        return [
            (a, self(a))
            for a in range(-n, n + 1)
            if 0 < abs(a) < self(a)
        ]

    def neg(self):
        n = self.rank
        return [(-a, -a) for a in range(1, n + 1) if self(a) == -a]

    def fix(self):
        n = self.rank
        return [(a, a) for a in range(1, n + 1) if self(a) == a]

    def cyc(self):
        return sorted(self.pair() + self.neg() + self.fix())

    @classmethod
    def identity(cls):
        return Permutation()

    @classmethod
    def longest_element(cls, n):
        return Permutation(*[-i for i in range(1, n + 1)])

    @classmethod
    def longest_unsigned(cls, n):
        return Permutation(*[i for i in range(n, 0, -1)])

    @classmethod
    def s_i(cls, i):
        assert 0 <= i
        if i == 0:
            oneline = [-1]
        else:
            oneline = list(range(1, i)) + [i + 1, i]
        return Permutation(*oneline)

    @property
    def right_descent_set(self):
        if self._rdes is None:
            self._rdes = set()
            if self.rank >= 1 and self(1) < 0:
                self._rdes.add(0)
            for i in range(1, self.rank):
                if self(i) > self(i + 1):
                    self._rdes.add(i)
        return self._rdes

    @property
    def left_descent_set(self):
        if self._ldes is None:
            self._ldes = self.inverse().right_descent_set
        return self._ldes

    def __len__(self):
        if self._len is None:
            biline = [-i for i in reversed(self.oneline)] + list(self.oneline)
            n = self.rank * 2
            inv = len([(i, j) for i in range(n) for j in range(i + 1, n) if biline[i] > biline[j]])
            inv_zero = len([i for i in self.oneline if i < 0])
            assert inv % 2 == inv_zero % 2
            self._len = (inv + inv_zero) // 2
        return self._len

    @property
    def length(self):
        return self._len

    @property
    def involution_length(self):
        return (len(self.neg()) + len(self.pair()) + len(self)) // 2

    def __mod__(self, other):
        assert type(other) == Permutation
        if other.left_descent_set:
            i = next(iter(other.left_descent_set))
            s = Permutation.s_i(i)
            if i in self.right_descent_set:
                return self * (s * other)
            else:
                return (self * s) % (s * other)
        else:
            return self

    def __mul__(self, other):
        assert type(other) == Permutation
        newline = [self(other(i)) for i in range(1, max(self.rank, other.rank) + 1)]
        return Permutation(*newline)

    def inverse(self):
        newline = self.rank * [0]
        for i in range(1, self.rank + 1):
            j = self(i)
            if j > 0:
                newline[j - 1] = i
            else:
                newline[-j - 1] = -i
        return Permutation(*newline)

    def __lt__(self, other):
        assert type(other) == Permutation
        return self.oneline < other.oneline

    def __eq__(self, other):
        assert type(other) == Permutation
        return self.oneline == other.oneline

    def __pow__(self, n):
        if n < 0:
            return self.inverse().__pow__(-n)
        elif n == 0:
            return Permutation.identity(self.rank)
        elif n == 1:
            return Permutation(*self.oneline)
        else:
            p = n // 2
            q = n - p
            return self.__pow__(p) * self.__pow__(q)

    def is_even_signed(self):
        return len([i for i in self.oneline if i < 0]) % 2 == 0

    def last_descent(self):
        n = self.rank
        descents = [i for i in range(1, n) if self(i) > self(i + 1)]
        if descents:
            return max(descents)

    @classmethod
    def reflection_s(cls, i, j):
        caller = list(range(1, j + 1))
        caller[i - 1] = -j
        caller[j - 1] = -i
        return cls(*caller)

    @classmethod
    def transposition(cls, i, j):
        return cls.reflection_t(i, j)

    @classmethod
    def reflection_t(cls, i, j):
        assert i != j
        caller = list(range(1, j + 1))
        caller[i - 1] = j
        caller[j - 1] = i
        return cls(*caller)

    def tex(self):
        s = '$'
        for i in self.oneline:
            if i > 0:
                s += str(i) + '\\hs '
            else:
                s += '\\bar' + str(-i) + '\\hs '
        s = s[:-len('\\hs ')]
        s += '$'
        return s

    def _min_inv_atom_oneline(self):
        tup = tuple(i for p in self.cyc() for i in reversed(p))
        minimum = []
        for i in tup:
            if minimum and minimum[-1] == i:
                continue
            minimum += [i]
        return tuple(minimum)

    def get_min_atom(self):
        assert self == self.inverse()
        return Permutation(*self._min_inv_atom_oneline())

    def get_atoms(self):
        assert self == self.inverse()
        w = self.reduce()
        if w not in atoms_b_cache:
            atoms_b_cache[w] = list(w._get_atoms())
        return atoms_b_cache[w]

    @classmethod
    def from_involution_word(cls, word, strict=True):
        w = cls.identity()
        for i in word:
            s = Permutation.s_i(i)
            if i in w.right_descent_set:
                if strict:
                    return None
            elif s * w == w * s:
                w = w * s
            else:
                w = s * w * s
        return w

    def _get_atoms(self):
        def involution(oneline):
            word = reversed(Permutation(*oneline).get_reduced_word())
            return self.from_involution_word(word)

        def next(oneline):
            y = involution(oneline)
            assert y is not None
            for i in range(len(oneline) - 2):
                c, a, b = oneline[i:i + 3]
                if a < b < c:
                    newline = oneline[:i] + (b, c, a) + oneline[i + 3:]
                    yield newline
            for i in range(len(oneline) - 1):
                b_, a_ = oneline[i:i + 2]
                a, b = -a_, -b_
                if 0 < a < b == min(map(abs, oneline[:i + 1])):
                    newline = oneline[:i] + (a, -b) + oneline[i + 2:]
                    z = involution(newline)
                    if z and y == z:
                        yield newline

        minimum = self._min_inv_atom_oneline()
        add = {minimum}
        while add:
            for w in add:
                yield Permutation(*w).inverse()
            add = {new for w in add for new in next(w)}

    def get_atoms_d(self):
        assert self.is_even_signed()
        assert self == self.inverse()
        w = self.reduce()
        if w not in atoms_d_cache:
            atoms_d_cache[w] = list(w._get_atoms_d())
        return atoms_d_cache[w]

    def _get_atoms_d(self):
        def length(w):
            ans = 0
            for i in range(1, w.rank + 1):
                for j in range(i + 1, w.rank + 1):
                    if w(i) > w(j):
                        ans += 1
                    if -w(i) > w(j):
                        ans += 1
            return ans

        if length(self) == 0:
            yield self
            return

        def s_i(i):
            return self.s_i(i) if i != 0 else self.s_i(0) * self.s_i(1) * self.s_i(0)

        for i in range(self.rank):
            s = s_i(i)
            w = self * s
            if length(w) < length(self):
                if w == s * self:
                    for a in w.get_atoms_d():
                        yield a * s
                else:
                    for a in (s * w).get_atoms_d():
                        yield a * s

    @classmethod
    def marked_hecke_words(cls, n, length_bound=-1):
        for level in cls.marked_hecke_levels(n, length_bound):
            for pi, w in level:
                yield w

    @classmethod
    def marked_hecke_levels(cls, n, length_bound=-1):
        def valid(w):
            # if len(w) > 2 and (w[0] + 1) // 2 == (w[1] + 1) // 2 and w[0] != w[1]:
            #     return False
            # if w[:4] in [
            #     (6, 5, 1, 4), (6, 5, 2, 4),
            #     (2, 1, 6, 4), (2, 1, 5, 4),
            #     (1, 2, 5, 4), (1, 2, 6, 4),
            #     (5, 6, 1, 4), (5, 6, 2, 4),
            # ]:
            #     return False
            # if w[:5] in [
            #     (6, 5, 2, 1, 4), (2, 1, 6, 5, 4),
            #     (6, 2, 6, 5, 4), (2, 6, 2, 1, 4),
            # ]:
            #     return False
            return True
        start = (cls.identity(), ())
        level = {start}
        while level:
            next_level = set()
            yield level
            for pi, w in level:
                for i in range(n):
                    s = Permutation.s_i(i)
                    sigma = s % pi % s
                    ww = w + (2 * i,)
                    if valid(ww):
                        next_level.add((sigma, ww))
                    if s * sigma == sigma * s:
                        ww = w + (2 * i - 1,)
                        if valid(ww):
                            next_level.add((sigma, ww))
            level = next_level
            if length_bound == 0:
                break
            length_bound -= 1

    def get_marked_hecke_words(self, length_bound):
        for level in self.marked_hecke_levels(self.rank, length_bound):
            for pi, w in level:
                if self == pi:
                    yield w

    def filter_marked_hecke_words(self, length_bound):
        for level in self.marked_hecke_levels(self.rank, length_bound):
            qq = []
            for w in [w for pi, w in level if self == pi]:
                qq += [(w, Word.quasisymmetrize(w, Word.decreasing_zeta))]
            for w, q in sorted(qq, key=lambda m: tuple(reversed(sorted(max(m[1]))))):
                # letters = tuple(sorted([i // 2 if i % 2 == 0 else -(i + 1) // 2 for i in w]))
                ww = '(' + ', '.join([
                    '+' + str(i // 2) if i % 2 == 0 else str(-(i + 1) // 2)
                    for i in w
                ]) + ')'
                print(ww, ':', q)
            print()
        ans, rows = self.marked_stable_grothendieck(length_bound)

        print(ans)
        print()

        def sort(t):
            return tuple(reversed(sorted(t)))

        first = rows[0]
        sorted_keys = sorted([key for key in first])
        split = str(first).split(' + ')
        for i in range(len(split) - 1):
            split[i] += ' + '
        print(first)
        for row in rows[1:]:
            line = ''
            for i in range(len(sorted_keys)):
                k = [k for k in row if sort(k) == sorted_keys[i]]
                if k:
                    s = str(Vector({k[0]: row[k[0]]}))
                    line += s + (len(split[i]) - len(s)) * ' '
                else:
                    line += len(split[i]) * ' '
            print(line)
        print()
