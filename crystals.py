from tableaux import Tableau, Partition
from permutations import Permutation
from insertion import InsertionAlgorithm
import subprocess


class TableauCrystalGenerator:

    DIRECTORY = '/Users/emarberg/Desktop/crystals/tableaux/'

    def _filename(self, j):
        mu = ''.join([str(i) for i in self.mu])
        return 'size%s_%s_%s' % (str(len(self._components[j])), str(self.max_entry), mu)

    def dot_filename(self, i):
        return self.DIRECTORY + 'dot/' + '%s_component%s.dot' % (self._filename(i), i)

    def png_filename(self, i):
        return self.DIRECTORY + 'png/' + '%s_component_%s.png' % (self._filename(i), i)

    def node_label(self, i):
        return str(self[i])

    def write_dotfile(self, i):
        s = []
        s += ['digraph G {']
        s += ['    overlap=false;']
        s += ['    splines=line;']
        s += ['    node [shape=box; fontname="courier"; style=filled];']
        s += [
            '    "%s" -> "%s" [label="%s"];' % (self.node_label(x), self.node_label(y), j)
            for (x, y, j) in self.edges if x in self._components[i]
        ]
        s += ['}']
        s = '\n'.join(s)

        with open(self.dot_filename(i), 'w') as f:
            f.write(s)

    def generate(self):
        self._components = dict(enumerate(self.components))
        self._components['all'] = list(range(len(self)))
        for i in self._components:
            self.write_dotfile(i)
            subprocess.run(["dot", "-Tpng", self.dot_filename(i), "-o", self.png_filename(i)])

    @classmethod
    def all(cls, n, max_entry):
        for mu in Partition.generate(n, max_entry):
            cg = cls(mu, max_entry)
            if cg.edges:
                cg.generate()

    def __init__(self, mu, max_entry):
        self.mu = mu
        self.max_entry = max_entry
        self.tableaux = list(Tableau.semistandard_setvalued(mu, max_entry))
        self._edges = None
        self._components = None

    def __len__(self):
        return len(self.tableaux)

    def __getitem__(self, i):
        return self.tableaux[i]

    @property
    def components(self):
        def dfs(seed):
            seen = {seed}
            added = {seed}
            while added:
                connected = set()
                for s in added:
                    connected |= \
                        {y for (x, y, e) in self.edges if s == x} | \
                        {x for (x, y, e) in self.edges if s == y}
                added = (connected - seen)
                seen |= added
            for x in seen:
                yield x

        rest = set(range(len(self)))
        components = []
        while rest:
            seed = rest.pop()
            comp = set(dfs(seed))
            rest -= comp
            components.append(comp)
        return components

    @property
    def edges(self):
        if self._edges is None:
            self._edges = []
            for x in range(len(self)):
                for i in range(1, self.max_entry):
                    y = self.f_crystal_operator(x, i)
                    if y is not None:
                        self._edges.append((x, y, i))
        return self._edges

    def f_crystal_operator(self, element_index, operator_index):
        t = self[element_index].f_crystal_operator(operator_index, multisetvalued=False)
        if t is None:
            return None
        else:
            element_indices = [x for x in range(len(self)) if self[x] == t]
            assert len(element_indices) == 1
            return element_indices[0]


# class ShiftedCrystalGenerator:

#     DIRECTORY = '/Users/emarberg/Desktop/crystals/shifted/'

#     def _filename(self, j):
#         w = ''.join([str(i) for i in self.oneline])
#         return 'size%s_%s_%s_%s' % (str(self.num_factors), str(len(self._components[j])), w, str(self.max_length))

#     def dot_filename(self, i):
#         return self.DIRECTORY + 'dot/' + '%s_component%s.dot' % (self._filename(i), i)

#     def png_filename(self, i):
#         return self.DIRECTORY + 'png/' + '%s_component%s.png' % (self._filename(i), i)

#     def node_label(self, i):
#         return str(self[i])

#     def write_dotfile(self, i):
#         s = []
#         s += ['digraph G {']
#         s += ['    overlap=false;']
#         s += ['    splines=line;']
#         s += ['    node [shape=box; fontname="courier"; style=filled];']
#         s += [
#             '    "%s" -> "%s" [label="%s"];' % (self.node_label(x), self.node_label(y), j)
#             for (x, y, j) in self.edges if x in self._components[i]
#         ]
#         s += ['}']
#         s = '\n'.join(s)

#         with open(self.dot_filename(i), 'w') as f:
#             f.write(s)

#     def generate(self):
#         self._components = dict(enumerate(self.components))
#         for i in self._components:
#             self.write_dotfile(i)
#             subprocess.run(["dot", "-Tpng", self.dot_filename(i), "-o", self.png_filename(i)])

#     @classmethod
#     def all(cls, n, k, l):
#         for w in Permutation.all(n):
#             if w(n) == n:
#                 continue
#             cg = cls(w.oneline, k, l)
#             if cg.edges:
#                 cg.generate()

#     def __init__(self, oneline, k, length):
#         self.oneline = oneline
#         self.words = [w for w in Permutation(*oneline).get_hecke_words(length_bound=length) if len(w) == length]
#         self.num_factors = k
#         self.max_length = length
#         self._edges = None
#         self._components = None

#     def __len__(self):
#         return len(self.factorizations)

#     def __getitem__(self, i):
#         return self.factorizations[i]

#     @property
#     def components(self):
#         def dfs(seed):
#             seen = {seed}
#             added = {seed}
#             while added:
#                 connected = set()
#                 for s in added:
#                     connected |= \
#                         {y for (x, y, e) in self.edges if s == x} | \
#                         {x for (x, y, e) in self.edges if s == y}
#                 added = (connected - seen)
#                 seen |= added
#             for x in seen:
#                 yield x

#         rest = set(range(len(self)))
#         components = []
#         while rest:
#             seed = rest.pop()
#             comp = set(dfs(seed))
#             rest -= comp
#             components.append(comp)
#         return components

#     @property
#     def edges(self):
#         if self._edges is None:
#             self._edges = []
#             for x in range(len(self)):
#                 for i in range(self.num_factors):
#                     y = self.f_crystal_operator(x, i)
#                     if y is not None:
#                         self._edges.append((x, y, i))
#         return self._edges

#     def f_crystal_operator(self, element_index, operator_index):
#         w, i = self.factorizations[element_index]
#         p, q = InsertionAlgorithm.hecke(w, i)

#         q = q.f_crystal_operator(operator_index)
#         if q is None:
#             return None
#         else:
#             ans = InsertionAlgorithm.inverse_hecke(p, q)
#             element_indices = [x for x in range(len(self)) if self[x] == ans]
#             assert len(element_indices) == 1
#             return element_indices[0]


class WordCrystalGenerator:

    DIRECTORY = '/Users/emarberg/Desktop/crystals/words/'

    def _filename(self, j):
        w = ''.join([str(i) for i in self.oneline])
        return 'size%s_%s_%s_%s' % (str(len(self._components[j])), str(self.num_factors), str(self.max_length), w)

    @property
    def symbol(self):
        return '>' if self._decreasing else '<='

    def dot_filename(self, i):
        return self.DIRECTORY + 'dot/' + '%s %s_component_%s.dot' % (self.symbol, self._filename(i), i)

    def png_filename(self, i):
        return self.DIRECTORY + 'png/' + '%s %s_component_%s.png' % (self.symbol, self._filename(i), i)

    def node_label(self, i):
        word, record = self[i]
        w = tuple(tuple(word[x] for x in range(len(record)) if record[x] == y) for y in range(1, self.num_factors + 1))
        pre = ''.join(str(part) for part in w)
        # top = ' '.join(map(lambda x: str(x) + (2 - len(str(x))) * ' ', self[i][0]))
        # bot = ' '.join(map(lambda x: str(x) + (2 - len(str(x))) * ' ', self[i][1]))
        return pre # + '\n\n' + top + '\n' + bot

    def write_dotfile(self, i):
        s = []
        s += ['digraph G {']
        s += ['    overlap=false;']
        s += ['    splines=line;']
        s += ['    node [shape=box; fontname="courier"; style=filled];']
        s += [
            '    "%s" -> "%s" [label="%s"];' % (self.node_label(x), self.node_label(y), j)
            for (x, y, j) in self.edges if x in self._components[i]
        ]
        s += ['}']
        s = '\n'.join(s)

        with open(self.dot_filename(i), 'w') as f:
            f.write(s)

    def generate(self):
        self._components = dict(enumerate(self.components))
        self._components['all'] = list(range(len(self)))
        for i in self._components:
            self.write_dotfile(i)
            subprocess.run(["dot", "-Tpng", self.dot_filename(i), "-o", self.png_filename(i)])

    @classmethod
    def all(cls, n, k, l, decreasing=False):
        for w in Permutation.all(n):
            if w(n) == n:
                continue
            cg = cls(w.oneline, k, l, decreasing)
            if cg.edges:
                cg.generate()

    def __init__(self, oneline, k, length, decreasing=False):
        self.oneline = oneline
        self.words = [w for w in Permutation(*oneline).get_hecke_words(length_bound=length) if len(w) == length]
        self.num_factors = k
        self.max_length = length
        self._factorizations = None
        self._edges = None
        self._components = None
        self._decreasing = decreasing

    @classmethod
    def factorization_tuple_to_array(cls, tup):
        w, i = (), ()
        index = 1
        for part in tup:
            w += part
            i += len(part) * (index,)
            index += 1
        return w, i

    @property
    def factorizations(self):
        if self._factorizations is None:
            fac = [
                self.factorization_tuple_to_array(tup)
                for w in self.words
                for tup in self.get_increasing_factorizations(w, self.num_factors, self._decreasing)
            ]
            self._factorizations = sorted(fac, key=lambda f: tuple(len(a) for a in f), reverse=True)
        return self._factorizations

    def __len__(self):
        return len(self.factorizations)

    def __getitem__(self, i):
        return self.factorizations[i]

    @property
    def components(self):
        def dfs(seed):
            seen = {seed}
            added = {seed}
            while added:
                connected = set()
                for s in added:
                    connected |= \
                        {y for (x, y, e) in self.edges if s == x} | \
                        {x for (x, y, e) in self.edges if s == y}
                added = (connected - seen)
                seen |= added
            for x in seen:
                yield x

        rest = set(range(len(self)))
        components = []
        while rest:
            seed = rest.pop()
            comp = set(dfs(seed))
            rest -= comp
            components.append(comp)
        return components

    @property
    def edges(self):
        if self._edges is None:
            self._edges = []
            for x in range(len(self)):
                for i in range(self.num_factors):
                    y = self.f_crystal_operator(x, i)
                    if y is not None:
                        self._edges.append((x, y, i))
        return self._edges

    def f_crystal_operator(self, element_index, operator_index):
        multisetvalued = not self._decreasing
        w, i = self.factorizations[element_index]
        p, q = InsertionAlgorithm.hecke(w, i, multisetvalued=multisetvalued)

        print('?')
        print(operator_index)
        print(p)
        print(q)
        q = q.f_crystal_operator(operator_index, multisetvalued=multisetvalued)
        print(q)
        if q is None:
            return None
        else:
            print(w, i)
            ans = InsertionAlgorithm.inverse_hecke(p, q, multisetvalued=multisetvalued)
            print(ans)
            element_indices = [x for x in range(len(self)) if self[x] == ans]
            assert len(element_indices) == 1
            return element_indices[0]

    @classmethod
    def get_increasing_factorizations(cls, w, k, decreasing=False):
        def is_increasing(x):
            if decreasing:
                return all(x[i] < x[i - 1] for i in range(1, len(x)))
            else:
                return all(x[i] >= x[i - 1] for i in range(1, len(x)))
        #
        if k == 0 and len(w) == 0:
            yield tuple()
        elif k == 0:
            return
        elif len(w) == 0:
            yield tuple(k * [tuple()])
        else:
            for i in range(len(w) + 1):
                if is_increasing(w[:i]):
                    for tup in cls.get_increasing_factorizations(w[i:], k - 1, decreasing):
                        yield (tuple(w[:i]),) + tup
                else:
                    break
