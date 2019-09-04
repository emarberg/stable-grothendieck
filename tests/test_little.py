from permutations import Permutation
from insertion import InsertionAlgorithm
from tableaux import Tableau


def test_push():
    word = (1, 3, 4, 3, 2, 3, 1)
    i = 6

    new, j = (1, 3, 4, 3, 2, 3, 2), 3
    assert Permutation.little_push(word, i) == (new, j)
    word, i = new, j

    new, j = (1, 3, 4, 4, 2, 3, 2), 2
    assert Permutation.little_push(word, i) == (new, j)
    word, i = new, j

    new, j = (1, 3, 5, 4, 2, 3, 2), None
    assert Permutation.little_push(word, i) == (new, j)


def test_bump():
    word = (1, 3, 4, 3, 2, 3, 1)
    end = (1, 3, 5, 4, 2, 3, 2)
    j, k = 1, 2
    assert Permutation.little_bump(word, j, k) == end


def get_ck_type_two(word):
    for i in range(len(word) - 2):
        a, c, b = word[i:i + 3]
        if a < b < c:
            yield i, word[:i] + (c, a, b) + word[i + 3:]
        c, a, b = word[i:i + 3]
        if a < b < c:
            yield i, word[:i] + (a, c, b) + word[i + 3:]


def get_ck_type_one(word):
    for i in range(len(word) - 2):
        b, a, c = word[i:i + 3]
        if a < b < c:
            yield i, word[:i] + (b, c, a) + word[i + 3:]
        b, c, a = word[i:i + 3]
        if a < b < c:
            yield i, word[:i] + (b, a, c) + word[i + 3:]


def get_ck_type_three(word):
    for i in range(len(word) - 2):
        b, a, bb = word[i:i + 3]
        if a < b == bb:
            yield i, word[:i] + (a, b, a) + word[i + 3:]
        a, b, aa = word[i:i + 3]
        if a == aa < b:
            yield i, word[:i] + (b, a, b) + word[i + 3:]


def get_ck_type_inv(word):
    if len(word) >= 2 and word[0] != word[1]:
        yield -1, (word[1], word[0]) + word[2:]


def get_ck_type_fpf(word):
    if len(word) >= 2 and word[0] != word[1]:
        if word[0] % 2 == word[1] % 2:
            yield -1, (word[1], word[0]) + word[2:]
        elif word[1] == word[0] - 1:
            yield -1, (word[0], word[1] + 2) + word[2:]
        elif word[1] == word[0] + 1:
            yield -1, (word[0], word[1] - 2) + word[2:]


def get_inv_ck_moves(word):
    for i, v in get_ck_moves(word):
        yield i, v
    for i, v in get_ck_type_inv(word):
        yield i, v


def get_fpf_ck_moves(word):
    for i, v in get_ck_moves(word):
        yield i, v
    for i, v in get_ck_type_fpf(word):
        yield i, v


def get_ck_moves(word):
    for i, v in get_ck_type_one(word):
        yield i, v
    for i, v in get_ck_type_two(word):
        yield i, v
    for i, v in get_ck_type_three(word):
        yield i, v


def get_bruhat_covers(w):
    n = w.rank
    for j in range(1, n + 1):
        for k in range(j + 1, n + 1):
            t = Permutation.transposition(j, k)
            if len(w * t) == len(w) - 1:
                yield j, k


def get_inv_bruhat_covers(w):
    n = w.rank
    for j in range(1, n + 1):
        for k in range(j + 1, n + 1):
            t = Permutation.transposition(j, k)
            subatoms = [x * t for x in w.get_atoms() if len(x * t) == len(x) - 1]
            subatoms = [x for x in subatoms if (x.inverse() % x).involution_length == len(x)]
            if subatoms:
                yield j, k


def adjacent_unmatched_primes(q, k):
    if k in q.values() and -k - 1 in q.values():
        (x1, y1) = q.find(k)[0]
        (x2, y2) = q.find(-k - 1)[0]
        return x1 == x2 or y1 == y2
    elif -k in q.values() and k + 1 in q.values():
        (x1, y1) = q.find(-k)[0]
        (x2, y2) = q.find(k + 1)[0]
        return x1 == x2 or y1 == y2
    elif k in q.values() and k + 1 in q.values():
        (x1, y1) = q.find(k)[0]
        (x2, y2) = q.find(k + 1)[0]
        return x1 == x2 and x1 == y1
    return False


def on_diagonal(q, k):
    if k in q.values():
        i, j = q.find(k)[0]
        return i == j
    return False


def test_ck_moves():
    for w in Permutation.all(4):
        covers = list(get_bruhat_covers(w))

        for word in w.get_reduced_words():
            alphas = {}
            p, q = InsertionAlgorithm.hecke(word)
            for i, other in get_ck_type_one(word):
                pp, qq = InsertionAlgorithm.hecke(other)
                print('  type I:', word, '<->', other, 'at i =', i + 1)
                print(q)
                print(qq)
                assert p == pp
                s = Permutation.s_i(i + 2)
                assert q == s * qq
                alphas[i] = s

            for i, other in get_ck_type_three(word):
                pp, qq = InsertionAlgorithm.hecke(other)
                print('type III:', word, '<->', other, 'at i =', i + 1)
                print(q)
                print(qq)
                assert p == pp
                s = Permutation.s_i(i + 2)
                assert q == s * qq
                alphas[i] = s

            for i, other in get_ck_type_two(word):
                pp, qq = InsertionAlgorithm.hecke(other)
                print(' type II:', word, '<->', other, 'at i =', i + 1)
                print(q)
                print(qq)
                assert p == pp
                s = Permutation.s_i(i + 2)
                t = Permutation.s_i(i + 1)
                assert q == s * qq or q == t * qq
                alphas[i] = s if q == s * qq else t

            for j, k in covers:
                bumped = Permutation.little_bump(word, j, k)
                p, q = InsertionAlgorithm.hecke(bumped)
                label = 'bump(%s, %s, %s) =' % (str(word), j, k)

                moves = []
                for i, other in get_ck_type_one(bumped):
                    pp, qq = InsertionAlgorithm.hecke(other)
                    print('  type I:', label, bumped, '<->', other, 'at i =', i + 1)
                    print(q)
                    print(qq)
                    assert i in alphas and q == alphas[i] * qq
                    moves.append(i)

                for i, other in get_ck_type_three(bumped):
                    pp, qq = InsertionAlgorithm.hecke(other)
                    print('type III:', label, bumped, '<->', other, 'at i =', i + 1)
                    print(q)
                    print(qq)
                    assert i in alphas and q == alphas[i] * qq
                    moves.append(i)

                for i, other in get_ck_type_two(bumped):
                    pp, qq = InsertionAlgorithm.hecke(other)
                    print(' type II:', label, bumped, '<->', other, 'at i =', i + 1)
                    print(q)
                    print(qq)
                    assert i in alphas and q == alphas[i] * qq
                    moves.append(i)

                assert sorted(moves) == sorted(alphas)


def test_unique_transposition():
    for w in Permutation.all(3):
        for word in w.get_reduced_words():
            p, q = InsertionAlgorithm.hecke(word)
            for i, other in get_ck_moves(word):
                pp, qq = InsertionAlgorithm.hecke(other)

                print(word, '<->', other, 'at i =', i)
                print(q)
                print(qq)

                s = Permutation.s_i(i + 2)
                t = Permutation.s_i(i + 1)
                pi = [s, t, s * t * s]

                a = [x for x in pi if qq == x * q]
                b = [x for x in pi if qq.descent_set() == (x * q).descent_set() and (x * q).is_semistandard()]
                print('a =', a)
                print('b =', b)
                print()
                assert a == b
                assert len(a) == 1


def test_unique_transposition_invol():
    for w in Permutation.involutions(5):
        for word in w.get_involution_words():
            p, q = InsertionAlgorithm.orthogonal_hecke(word)
            for i, other in get_inv_ck_moves(word):
                pp, qq = InsertionAlgorithm.orthogonal_hecke(other)

                if i == -1 and other == (word[1], word[0]) + word[2:]:
                    r = Permutation(1, -2)
                    pi = [r]
                else:
                    s = Permutation.s_i(i + 1)
                    t = Permutation.s_i(i + 2)

                    f = Permutation.reflection_s(i + 1)
                    g = Permutation.reflection_s(i + 2)
                    h = Permutation.reflection_s(i + 3)

                    pi = []
                    if adjacent_unmatched_primes(q, i + 1):
                        pi += [g] if on_diagonal(q, i + 1) else [f * g]
                    else:
                        pi.append(s)
                    if adjacent_unmatched_primes(q, i + 2):
                        pi += [h] if on_diagonal(q, i + 2) else [g * h]
                    else:
                        pi.append(t)

                a = set(x for x in pi if qq == x * q)
                b = set(x for x in pi if qq.descent_set() == (x * q).descent_set() and (x * q).is_semistandard())

                if len(b) != 1:
                    print(pi)
                    print(word, '<->', other, 'at i =', i + 1)
                    print(q)
                    print(qq)
                    print('a =', a)
                    print('b =', b)
                    print()

                assert a == b
                assert len(b) == 1


def test_inv_ck_moves():
    for w in Permutation.involutions(4):
        covers = list(get_inv_bruhat_covers(w))

        for word in w.get_involution_words():
            alphas = {}
            p, q = InsertionAlgorithm.orthogonal_hecke(word)
            for i, other in get_ck_type_one(word):
                pp, qq = InsertionAlgorithm.orthogonal_hecke(other)
                print('  type I:', word, '<->', other, 'at i =', i + 1)
                print(q)
                print(qq)
                assert p == pp
                if adjacent_unmatched_primes(q, i + 2):
                    s = Permutation.reflection_s(i + 2) * Permutation.reflection_s(i + 3)
                else:
                    s = Permutation.s_i(i + 2)
                assert q == s * qq
                alphas[i] = s

            for i, other in get_ck_type_three(word):
                pp, qq = InsertionAlgorithm.orthogonal_hecke(other)
                print('type III:', word, '<->', other, 'at i =', i + 1)
                print(q)
                print(qq)
                assert p == pp
                if adjacent_unmatched_primes(q, i + 2):
                    s = Permutation.reflection_s(i + 2) * Permutation.reflection_s(i + 3)
                else:
                    s = Permutation.s_i(i + 2)
                assert q == s * qq
                alphas[i] = s

            for i, other in get_ck_type_two(word):
                pp, qq = InsertionAlgorithm.orthogonal_hecke(other)
                print(' type II:', word, '<->', other, 'at i =', i + 1)
                print(q)
                print(qq)
                assert p == pp
                if adjacent_unmatched_primes(q, i + 2):
                    s = Permutation.reflection_s(i + 2) * Permutation.reflection_s(i + 3)
                else:
                    s = Permutation.s_i(i + 2)
                if adjacent_unmatched_primes(q, i + 1):
                    if on_diagonal(q, i + 1):
                        t = Permutation.reflection_s(i + 2)
                    else:
                        t = Permutation.reflection_s(i + 1) * Permutation.reflection_s(i + 2)
                else:
                    t = Permutation.s_i(i + 1)
                assert q == s * qq or q == t * qq
                alphas[i] = s if q == s * qq else t

            for i, other in get_ck_type_inv(word):
                pp, qq = InsertionAlgorithm.orthogonal_hecke(other)
                print(' type IV:', word, '<->', other, 'at i =', i + 1)
                print(q)
                print(qq)
                assert p == pp
                r = Permutation(1, -2)
                assert q == r * qq

            for j, k in covers:
                bumped = Permutation.involution_little_bump(word, j, k)
                p, q = InsertionAlgorithm.orthogonal_hecke(bumped)
                label = 'ibump(%s, %s, %s) =' % (str(word), j, k)

                moves = []
                for i, other in get_ck_type_one(bumped):
                    pp, qq = InsertionAlgorithm.orthogonal_hecke(other)
                    print('  type I:', label, bumped, '<->', other, 'at i =', i + 1)
                    print(q)
                    print(qq)
                    print('s =', alphas[i])
                    print()
                    assert i in alphas and q == alphas[i] * qq
                    moves.append(i)

                for i, other in get_ck_type_three(bumped):
                    pp, qq = InsertionAlgorithm.orthogonal_hecke(other)
                    print('type III:', label, bumped, '<->', other, 'at i =', i + 1)
                    print(q)
                    print(qq)
                    print('s =', alphas[i])
                    print()
                    assert i in alphas and q == alphas[i] * qq
                    moves.append(i)

                for i, other in get_ck_type_two(bumped):
                    pp, qq = InsertionAlgorithm.orthogonal_hecke(other)
                    print(' type II:', label, bumped, '<->', other, 'at i =', i + 1)
                    print(q)
                    print(qq)
                    print('s =', alphas[i])
                    print()
                    assert i in alphas and q == alphas[i] * qq
                    moves.append(i)

                for i, other in get_ck_type_inv(bumped):
                    pp, qq = InsertionAlgorithm.orthogonal_hecke(other)
                    print(' type IV:', label, bumped, '<->', other, 'at i =', i + 1)
                    print(q)
                    print(qq)
                    s = Permutation(1, -2)
                    print('s =', s)
                    print()
                    assert q == s * qq

                assert sorted(moves) == sorted(alphas)


def test_ck_moves_type_two():
    for w in Permutation.all(3):
        for word in w.get_reduced_words():
            p, q = InsertionAlgorithm.hecke(word)
            for i, other in get_ck_type_two(word):
                pp, qq = InsertionAlgorithm.hecke(other)
                s = Permutation.s_i(i + 2)
                t = Permutation.s_i(i + 1)
                a = Tableau({(x, y): v for x, y, v in q if i + 3 in v or i + 1 in v or i + 2 in v})
                aa = Tableau({(x, y): v for x, y, v in qq if i + 3 in v or i + 1 in v or i + 2 in v})
                assert q == s * qq or q == t * qq
                if q == t * qq:
                    print('type IIb:', word, '<->', other, 'at i =', i + 1)
                    print(p)
                    print(a)
                    print(aa)
                else:
                    print('type IIa:', word, '<->', other, 'at i =', i + 1)
                    print(p)
                    print(a)
                    print(aa)
#    assert False


def test_tab():
    for w in Permutation.involutions(3):
        if w.is_inv_grassmannian():
            for word in w.get_involution_words():
                p, q = InsertionAlgorithm.orthogonal_hecke(word)
                print(w)
                print(word)
                print(q)
                print()
#    assert False


def test_q_tab():
    n = 4
    for pi in Permutation.all(n):
        for w in pi.get_reduced_words():
            _, q = InsertionAlgorithm.hecke(w)
            for i in range(1, n):
                for j in range(i + 1, n + 1):
                    bumped = Permutation.little_bump(w, i, j)
                    if w != bumped:
                        _, qq = InsertionAlgorithm.hecke(bumped)
                        print(pi)
                        print(i, j, ':', w, '->', bumped)
                        print(q)
                        print(qq)
                        assert q == qq


def test_inv_q_tab():
    n = 5
    countmax = 5
    for pi in Permutation.involutions(n):
        if len(pi) == 0:
            continue
        for count, w in enumerate(pi.get_involution_words()):
            if count >= countmax:
                break
            _, q = InsertionAlgorithm.orthogonal_hecke(w)
            cases = 0
            for i in range(1, n):
                for j in range(i + 1, n + 1):
                    bumped = Permutation.involution_little_bump(w, i, j)
                    if w != bumped:
                        cases += 1
                        _, qq = InsertionAlgorithm.orthogonal_hecke(bumped)
                        print(pi)
                        print(i, j, ':', w, '->', bumped)
                        print(q)
                        print(qq)
                        assert q == qq
            assert cases != 0


def test_inv_bumping_multiplicity():
    n = 5
    countmax = 100
    for pi in Permutation.involutions(n):
        for count, w in enumerate(pi.get_involution_words()):
            if count >= countmax:
                break
            for i in range(1, n):
                for j in range(i + 1, n + 1):
                    bumped = Permutation.involution_little_bump(w, i, j)
                    if w != bumped:
                        assert {bumped[e] - w[e] for e in range(len(w))}.issubset({0, 1})


def test_fpf_bumping_multiplicity():
    n = 6
    countmax = 100
    for pi in Permutation.fpf_involutions(n):
        for count, w in enumerate(pi.get_fpf_involution_words()):
            if count >= countmax:
                break
            for i in range(1, n):
                for j in range(i + 1, n + 1):
                    bumped = Permutation.fpf_involution_little_bump(w, i, j)
                    if w != bumped:
                        assert {bumped[e] - w[e] for e in range(len(w))}.issubset({0, 1, 2})


def test_fpf_q_tab():
    n = 6
    countmax = None
    for pi in Permutation.fpf_involutions(n):
        if len(pi) == 0:
            continue
        for count, w in enumerate(pi.get_fpf_involution_words()):
            if len(w) == 0:
                continue
            if countmax and count >= countmax:
                break
            _, q = InsertionAlgorithm.symplectic_hecke(w)
            cases = 0
            for i in range(1, n):
                for j in range(i + 1, n + 1):
                    bumped = Permutation.fpf_involution_little_bump(w, i, j)
                    print(w, bumped, w == bumped)
                    if w != bumped:
                        cases += 1
                        _, qq = InsertionAlgorithm.symplectic_hecke(bumped)
                        if q != qq:
                            print(pi)
                            print(i, j, ':', w, '->', bumped)
                            print(q)
                            print(qq)
                        assert q == qq
            assert cases != 0
