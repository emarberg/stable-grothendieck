from permutations import Permutation
from partitions import Partition
from insertion import InsertionAlgorithm
from tableaux import Tableau
from words import Word
from tests.test_little import get_inv_ck_moves


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


def rectify_print(v, printw=True):
    p, q = InsertionAlgorithm.orthogonal_hecke(v)
    t = standardize(p)
    u, _ = InsertionAlgorithm.inverse_orthogonal_hecke(t, q)

    def st(word):
        w = len(word) * [0]
        for i, e in enumerate(sorted(enumerate(word), key=lambda p: (p[1], p[0]))):
            w[e[0]] = i + 1
        return tuple(w)

    if printw:
        labels = get_crossings(v)
        print(Word.involution_wiring_diagram(v, [''] + labels))
        # print(Word.wiring_diagram(u))
        print()
        # print(', '.join(tuple(str(i) + (' ' if i < 10 else '') for i in v)))
        # print(', '.join(tuple(str(i) + (' ' if i < 10 else '') for i in u)))
        # print()
    return u


def minimal_word(mu):
    if len(mu) == 0:
        return ()
    columns = [sum([1 for i in range(len(mu)) if mu[i] + i >= j]) for j in range(len(mu) + 1, mu[0] + 1)]
    word = [i for j in range(len(mu)) for i in range(2 * j + 1, j, -1)]
    m = max(word) + 1
    for c in columns:
        word += [i for i in range(m, m - c, -1)]
        m += 1
    return tuple(word)


def get_crossings(word):
    w = Permutation.from_involution_word(word)

    def endpoint(arrows, index, i):
        key = (index, i)
        while key in arrows:
            key = arrows[key]
        _, a = key
        if w(a) < a:
            a = w(a)
        return a

    n = (max(word) + 1) if word else 0
    arrows = {}
    for index, i in enumerate(word):
        arrows[(index, i)] = (index + 1, i + 1)
        arrows[(index, i + 1)] = (index + 1, i)
        for j in range(1, n + 1):
            if j not in [i, i + 1]:
                arrows[(index, j)] = (index + 1, j)
    ans = []
    for index, i in enumerate(word):
        a = endpoint(arrows, index, i)
        b = endpoint(arrows, index, i + 1)
        ans += [(a, b)]
    return ans


def test_inv_grassmannian_braids():
    rank = 7
    delta = tuple(range(rank - 1, 0, -2))
    rect = {}
    for mu in Partition.subpartitions(delta, strict=True):
        # w = (1, 3, 2, 7, 9, 8, 4, 6, 5, 6, 4, 3, 7, 4, 5, 6, 5, 4)
        # z = Permutation.from_involution_word(w)
        # mu = z.involution_shape()

        minword = minimal_word(mu)
        rect[minword] = rectify_print(minword, False)

        def getmap(word):
            return {
                tuple(sorted(x)): rect[word][i]
                for i, x in enumerate(get_crossings(word))
            }

        mmap = getmap(minword)

        w = minword
        # w = (1, 3, 2, 7, 9, 8, 4, 6, 5, 6, 4, 3, 7, 4, 5, 6, 5, 4)

        seen = set()
        add = {w}
        while True:
            nextadd = set()
            for v in add:
                rect[v] = rectify_print(v, False)
                vmap = getmap(v)
                for i, u in get_inv_ck_moves(v):

                    r = rect[v]
                    if i >= 0 and v[i] == v[i + 2] < v[i + 1]:
                        assert r[i + 2] < r[i] < r[i + 1]
                    if i >= 0 and v[i] == v[i + 2] > v[i + 1]:
                        assert r[i + 1] < r[i] < r[i + 2]
                    if i >= 0 and v[i + 1] < v[i] < v[i + 2]:
                        assert r[i + 1] < r[i] < r[i + 2]
                    if i >= 0 and v[i + 2] < v[i] < v[i + 1]:
                        assert r[i + 2] < r[i] < r[i + 1]

                    if u not in seen:
                        rect[u] = rectify_print(u, False)
                        umap = getmap(u)
                        nextadd.add(u)

                        if len([key for key in mmap if vmap[key] != umap[key]]) >= 3:
                            rectify_print(minword, True)
                            rectify_print(v, True)
                            rectify_print(u, True)

                            print(' ' + i * '   ' + '*')
                            print(v, '          ', rect[v])
                            print(u, '          ', rect[u])
                            print()

                            for key in mmap:
                                if vmap[key] != umap[key]:
                                    print(key, ':', mmap[key], '=>', vmap[key], '=>', umap[key])

                            print(10 * '\n')
                            input('pause')

            if len(nextadd) == 0:
                break
            seen |= add
            add = nextadd

    # for w in Permutation.inv_grassmannians(rank):
    #     w = w.star()
    #     rect = {}
    #     print(w.cycle_repr())
    #     # print(10 * '*')
    #     print()
    #     for v in w.get_involution_words():
    #         rect[v] = rectify_print(v)
        # for v in sorted(rect, key=lambda v: tuple(rect[v][i] - v[i] for i in range(len(v)))):
        #     r = {}
        #     for i, u in get_inv_ck_moves(rect[v]):
        #         r[i] = u
        #     for i, u in get_inv_ck_moves(v):
        #         d = tuple(rect[v][i] - v[i] for i in range(len(u)))
        #         e = tuple(rect[u][i] - u[i] for i in range(len(u)))
        #         f = tuple(d[i] - e[i] for i in range(len(u)))
        #         f = tuple(i for i in range(len(u)) if f[i] != 0)
        #         if f == () or (len(f) == 2 and d[f[0]] == e[f[1]] and d[f[1]] == e[f[0]]):
        #             continue
        #         # print(Word.wiring_diagram(v))
        #         print()
        #         print(' ' + i * '   ' + '*')
        #         print(v, '          ', rect[v])
        #         print(u, '          ', rect[u])
        #         print()
        #         print(d)
        #         print(e)
        #         print()
        #         # print(Word.wiring_diagram(u))
        #         print()
        #         assert v[i] == v[i + 2]
        #         assert rect[u][i + 1] == rect[v][i + 2]
        #         assert rect[v][i + 1] == rect[u][i + 2]

    print('\ndone')
    # assert False


#    assert False
# w = (1,3,2,7,9,8,4,5,6,5,4,3,7,4,6,5,6,4)
# w = (1 , 3 , 2 , 7 , 9 , 8 , 4 , 6 , 5 , 4 , 3 , 7 , 6 , 7 , 4 , 5 , 6 , 4)
