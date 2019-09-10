from tableaux import Tableau
from insertion import InsertionAlgorithm
from permutations import Permutation


def ck(i, w):
    if i == -1:
        assert len(w) >= 2
        a, b = w[:2]
        if a % 2 == b % 2 == 0:
            return (b, a) + w[2:]
        if a % 2 == 0 and b in [a - 1, a + 1]:
            pre = (a, a + 1) if b < a else (a, a - 1)
            return pre + w[2:]
        raise Exception
    else:
        assert i + 2 < len(w)
        pre, mid, pos = w[:i], w[i:i + 3], w[i + 3:]

        while True:
            c, a, b = mid
            if a < b < c:
                mid = (a, c, b)
                break
            a, c, b = mid
            if a < b < c:
                mid = (c, a, b)
                break

            b, c, a = mid
            if a < b < c:
                mid = (b, a, c)
                break
            b, a, c = mid
            if a < b < c:
                mid = (b, c, a)
                break

            a, b, c = mid
            if a == c:
                mid = (b, a, b)
                break

            raise Exception

        return pre + mid + pos


def descent(i, w):
    assert 0 <= i < len(w) - 1
    return w[i] > w[i + 1]


def get_rows(word):
    rows = []
    for x in word:
        if len(rows) == 0 or rows[-1][-1] > x:
            rows += [(x,)]
        else:
            rows[-1] += (x,)
    return rows


def get_columns(word):
    cols = []
    for x in word:
        if len(cols) == 0 or cols[-1][-1] < x:
            cols += [(x,)]
        else:
            cols[-1] += (x,)
    return cols


def to_row_reading(word):
    columns = get_columns(word)
    t = Tableau()
    for j, col in enumerate(columns):
        for i, v in enumerate(reversed(col)):
            t = t.add(i + 1, j + 1, v)
    row = t.row_reading_word()
    ans, col = to_column_reading(row)
    assert col == word
    return list(reversed(ans)), row


def to_column_reading(word):
    rows = get_rows(word)
    if len(rows) <= 1:
        return [], word

    starts = []
    for r in rows[:-1]:
        starts += [len(r) + (starts[-1] if starts else -1)]
    m = len(rows[-1])

    starts = [starts[i] for i in range(len(rows) - 1) if len(rows[i]) + len(rows) - i - 1 == m]
    m = len(starts) + 1

    ans = []
    for i, a in enumerate(starts):
        for b in range(i + 1):
            c = (starts[i + 1] - 1) if i + 1 < len(starts) else len(word) - 2
            for x in range(a - b, c - b):
                ans += [x]
                word = ck(x, word)
    subword = word[:-m]
    bns, newword = to_column_reading(subword)

    ans += bns
    word = newword + word[-m:]
    return ans, word


def sp_ck_compute(tab, letter):
    ans = []
    rows = get_rows(tab.row_reading_word())

    _print('rows:', rows, '<-', letter)
    _print()

    if len(rows) == 0:
        return 1, 1, True, ans

    left = ()

    for i in range(1, len(rows) + 1):
        o = sum([len(r) for r in rows[:-i]])
        row = rows[-i]

        _print('*', i)
        _print('*', row, ':', letter, ':', left)
        _print()

        working = row + (letter,)
        if not descent(len(row) - 1, working):
            j = len(row) + i
            return i, j, True, ans

        if i == len(rows):
            testing = working
            for a in range(-1, len(row) - 2):
                testing = ck(a, testing)
            if not descent(len(row) - 1, testing):
                ans += [a + o for a in range(len(row) - 2, -1, -1)]
                return i + 1, i + 1, True, ans

        for a in range(len(row) - 2, -1, -1):
            ans += [a + o]
            working = ck(a, working)

        letter = working[0]
        left = working[1:] + left

    working = (letter,) + left

    _print('* collecting diagonal')
    _print('*', working, '|', ans)
    _print()

    for i in range(1, len(rows)):
        for b in range(i):
            s = sum([len(r) for r in rows[-i:]])
            n = len(working) - s - 2 + b
            m = len(working) - s - len(rows[-i - 1]) - 1 + b
            for a in range(n, m, -1):
                ans += [a]
                working = ck(a, working)

    bit = working[:len(rows) + 1]
    left = working[len(rows) + 1:]

    _print('* reversing bit')
    _print('*', bit, ':', left, '|', ans)
    _print()

    for b in range(len(bit) - 2, -1, -1):
        for a in range(-1, b):
            ans += [a]
            bit = ck(a, bit)

    _print('* left to column word')
    _print('*', bit, ':', left, '|', ans)
    _print()

    seq, left = to_column_reading(left)
    ans += [a + len(bit) for a in seq]

    _print('* integration')
    _print('*', bit, ':', left, '|', ans)
    _print()

    m = len(rows) - 1
    working = bit + left
    start = m
    for a in range(m):
        for _ in range(a + 1):
            for b in range(start, start - m + a, -1):
                ans += [b]
                working = ck(b, working)
            start += 1

    i = m * (m + 1) // 2 + len(rows)
    bit, letter, left = working[:i], working[i], working[i + 1:]

    _print('* column insertion')
    _print('*', bit, ':', letter, ':', left, '|', ans)
    _print()

    offset = len(rows)
    columns = get_columns(left)

    if len(columns) == 0:
        ans += to_row_reading(bit + (letter,))[0]
        return 1, offset + 1, False, ans

    for j in range(1, 1 + len(columns)):
        o = i + sum([len(c) for c in columns[:j - 1]])
        col = columns[j - 1]
        working = (letter,) + col

        rest = ()
        for k in range(j, len(columns)):
            rest += columns[k]

        _print('* column', j + offset, j, len(columns))
        _print('*', working)
        _print()

        if descent(0, working):
            i = len(col) + 1
            ans += to_row_reading(bit + working + rest)[0]
            return i, j + offset, False, ans

        if j == len(columns):
            for a in range(len(working) - 2):
                ans += [a + o]
                working = ck(a, working)
            ans += to_row_reading(bit + working + rest)[0]
            return 1, j + offset + 1, False, ans

        for a in range(len(working) - 2):
            ans += [a + o]
            working = ck(a, working)

        letter = working[-1]
        bit += working[:-1]

    raise Exception


def _print(*args):
    print(*args)


def test_fpf(n=10, maxcount=10):
    for a, pi in enumerate(Permutation.fpf_involutions(n)):
        print(a, pi.fpf_involution_length, pi)
        for b, w in enumerate(pi.get_fpf_involution_words()):
            if len(w) == 0:
                continue
            if b > maxcount > 0:
                break

            v = w[:-1]
            p, q = InsertionAlgorithm.symplectic_hecke(v)
            t, r = InsertionAlgorithm.symplectic_hecke(w)

            if r.max_row() > 3:
                continue

            _print()
            _print('word =', w)
            _print(p)
            _print(q)
            _print(r)

            newboxes = [b for b in r.boxes if b not in q.boxes]
            assert len(newboxes) == 1

            i, j, sgn, seq = sp_ck_compute(p, w[-1])
            assert (i, j) == newboxes[0]
            assert sgn == (r.get(i, j) > 0)

            _print(seq)

            word = p.row_reading_word() + (w[-1],)
            for a in seq:
                word = ck(a, word)
            assert word == t.row_reading_word()
