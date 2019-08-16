from permutations import Permutation
from partitions import Partition
from insertion import InsertionAlgorithm
from tableaux import Tableau
from words import Word


def test_grassmannian():
    w = Permutation.get_grassmannian(4, 4, 3, 2, 2, 2, 1)
    assert w.shape() == (4, 4, 3, 2, 2, 2, 1)

    shapes = set(Partition.subpartitions([4, 3, 2, 1]))
    gr = list(Permutation.grassmannians(5))
    assert len(gr) == len(shapes)
    assert {w.shape() for w in gr} == shapes


def test_inv_grassmannian():
    w = Permutation.get_inv_grassmannian(4, 3, 1)
    t = Permutation.transposition
    assert w == t(1, 5) * t(2, 6) * t(4, 7)
    assert w.involution_shape() == (4, 3, 1)

    shapes = set(Partition.subpartitions([5, 3, 1], True))
    gr = list(Permutation.inv_grassmannians(6))
    assert len(gr) == len(shapes)
    assert {w.involution_shape() for w in gr} == shapes


def test_fpf_grassmannian():
    w = Permutation.get_fpf_grassmannian(4, 3, 1)
    t = Permutation.transposition

    assert w == t(1, 6) * t(2, 7) * t(4, 8) * t(3, 5)
    assert w.fpf_involution_shape() == (4, 3, 1)

    shapes = set(Partition.subpartitions([6, 4, 2], True))
    gr = list(Permutation.fpf_grassmannians(8))
    assert len(gr) == len(shapes)
    assert {w.fpf_involution_shape() for w in gr} == shapes


def test_eg_insertion():
    rank = 4
    for w in Permutation.grassmannians(rank):
        print(w.shape(), '=', 'shape(', w, ') = shape(', w.inverse(), '^-1 )')
        print()

        n = 0 if len(w) == 0 else list(w.left_descent_set)[0]
        a = tuple(w.inverse()(i) for i in range(n, 0, -1))
        b = tuple(w.inverse()(i) for i in range(n + 1, w.rank + 1))
        mapping = {(x, y): (i + 1, j + 1) for i, x in enumerate(a) for j, y in enumerate(b)}
        print(mapping)
        for word in w.get_reduced_words():

            tab = Tableau()
            partial = [Permutation()]
            for i in reversed(word):
                partial = [partial[0] * Permutation.s_i(i)] + partial
            for i, e in enumerate(word):
                x, y = partial[i](e), partial[i](e + 1)
                x, y = mapping[(x,y)]
                tab = tab.add(x, y, i + 1)

            print(w)
            print(word)
            print(Word.wiring_diagram(word))
            p, q = InsertionAlgorithm.hecke(word)
            print(p)
            print(q)
            print()

            assert q == tab
        print()
        print()
        print()


def double_inv_word(word):
    w = Permutation()
    ans = []
    for i in word:
        s = Permutation.s_i(i)
        if w * s == s * w:
            ans = [None] + ans + [i]
            w = w * s
        else:
            ans = [i] + ans + [i]
            w = s * w * s
    return tuple(ans)


def test_inv_eg_insertion():
    rank = 4
    for w in Permutation.inv_grassmannians(rank):
        w = w.star()
        print(w.involution_shape(), '=', 'shape(', w, ')')
        print()

#        n = 0 if len(w) == 0 else list(w.visible_descent_set)[0]
#        a = tuple(w.inverse()(i) for i in range(n, 0, -1))
#        b = tuple(w.inverse()(i) for i in range(n + 1, w.rank + 1))
#        mapping = {(x, y): (i + 1, j + 1) for i, x in enumerate(a) for j, y in enumerate(b)}

        for word in w.get_involution_words():
            p, q = InsertionAlgorithm.orthogonal_hecke(word)

            # if not q.is_shifted_column_superstandard():
            #     continue

            # tab = Tableau()
            # partial = [Permutation()]
            # for i in reversed(word):
            #     partial = [partial[0] * Permutation.s_i(i)] + partial
            # for i, e in enumerate(word):
            #     x, y = partial[i](e), partial[i](e + 1)
            #     x, y = mapping[(x,y)]
            #     tab = tab.add(x, y, i + 1)

            print()
            print()
            print()
            print(w.cycle_repr())
            print()
            print(w.get_min_atom())
            print()
            print(word)

            labels = {}
            for i in range(1, len(word) + 1):
                j = len(word) + i
                k = len(word) + 1 - i
                if q.find(i):
                    a = q.find(i)[0]
                else:
                    a = tuple(reversed(q.find(-i)[0]))
                labels[j] = a
                labels[k] = ''

            double = double_inv_word(word)
            #print(double)
            #print()
            print(q)
            print(Word.wiring_diagram(double, labels))
            print()
            print(p)
            #print()

        print()
        print()
        print()
#    assert False


def test_inv_grassmannian_braids():

    def standardize(t):
        seen = set()
        offset = 0
        offsets = {}
        for i, j in sorted(t.boxes, key=lambda box: (box[1], box[0])):
            v = t.get(i, j)
            while v + offset in seen:
                offset += 1
            seen.add(v + offset)
            offsets[(i, j)] = offset
        for i, j in offsets:
            t = t.replace(i, j, t.get(i, j) + offsets[(i, j)])
        return t

    rank = 6
    for w in Permutation.inv_grassmannians(rank):
        w = w.star()
        for v in w.get_involution_words():
            print()
            p, q = InsertionAlgorithm.orthogonal_hecke(v)
            print(p)
            t = standardize(p)
            print(t)
            u, _ = InsertionAlgorithm.inverse_orthogonal_hecke(t, q)
            print(v)
            print(u)
            d = tuple(u[i] - v[i] for i in range(len(u)))
            print(d)
            print(Word.wiring_diagram(v))
    assert False

