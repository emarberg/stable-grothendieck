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
            if  a < b < c:
                mid = (a, c, b)
                break
            a, c, b = mid
            if a < b < c:
                mid = (c, a, b)
                break

            b, c, a = mid
            if  a < b < c:
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
    # _print()
    # _print((1 + 3 * i) * ' ' + '*')
    # _print(w, 'descent at index', i)
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


def sp_ck_compute(tab, letter):
    rows = get_rows(tab.row_reading_word())

    _print('rows:', rows, '<-', letter)
    _print()

    if len(rows) == 0:
        return (1, 1, True)

    bit = (letter,)
    left = ()

    for i in range(1, len(rows) + 1):
        row = rows[-i]

        _print('*', i)
        _print('*', row, ':', bit, ':', left)
        _print()

        assert len(bit) == i

        working = row + bit
        if not descent(len(row) - 1, working):
            j = len(row) + i
            return i, j, True

        if i == len(rows):
            testing = working
            for a in range(-1, len(row) - 2):
                testing = ck(a, testing)
            if not descent(len(row) - 1, testing):
                return i + 1, i + 1, True
            elif i == 1:
                return 1, len(row) + 1, False

        for a in range(len(row) - 2, -1, -1):
            for b in range(i):
                working = ck(a + b, working)

        bit = working[:i + 1]
        left = working[i + 1:] + left

    _print('* reversing bit')
    _print('*', bit, ':', left)
    _print()

    for b in range(len(bit) - 2, -1, -1):
        for a in range(-1, b):
            bit = ck(a, bit)

    _print('* left to column word')
    _print('*', bit, ':', left)
    _print()

    leftrows = get_rows(left)
    lefttab = Tableau()
    for i, row in enumerate(reversed(leftrows)):
        for j, a in enumerate(row):
            lefttab = lefttab.add(i + 1, i + j + 1, a)
    left = lefttab.column_reading_word()

    _print('* integration')
    _print('*', bit, ':', left)
    _print()

    m  = len(rows) - 1
    working = bit + left
    start = m
    for a in range(m):
        for _ in range(a + 1):
            for b in range(start, start - m + a, -1):
                working = ck(b, working)
            start += 1

    i = m * (m + 1) // 2 + len(rows) + 1
    bit, left = working[:i], working[i:]
    letter = bit[-1]

    _print('* ready for column insertion')
    _print('*', bit[:-1], ':', letter, ':', left)
    _print()

    offset = len(rows)
    columns = get_columns(left)

    if len(columns) == 0:
        return (1, offset + 1, False)

    for j in range(1, 1 + len(columns)):
        col = columns[j - 1]
        working = (letter,) + col

        _print('* column', j + offset, j, len(columns))
        _print('*', working)
        _print()

        if descent(0, working):
            i = len(col) + 1
            return i, j + offset, False

        if j == len(columns):
            return 1, j + offset + 1, False

        for a in range(len(working) - 2):
            working = ck(a, working)

        letter = working[-1]

    raise Exception



def _print(*args):
    # print(*args)
    pass


def test_fpf(n=8, maxcount=20):
    for a, pi in enumerate(Permutation.fpf_involutions(n)):
        print(a, pi.fpf_involution_length, pi)
        for b, w in enumerate(pi.get_fpf_involution_words()):
            if len(w) == 0:
                continue
            if b > maxcount:
                break

            v = w[:-1]
            p, q = InsertionAlgorithm.symplectic_hecke(v)
            _, r = InsertionAlgorithm.symplectic_hecke(w)

            if r.max_row() > 3:
                continue

            _print()
            _print('word =', w)
            _print(p)
            _print(q)
            _print(r)

            newboxes = [b for b in r.boxes if b not in q.boxes]
            assert len(newboxes) == 1
            i, j = newboxes[0]
            sgn = r.get(i, j) > 0

            assert (i, j, sgn) == sp_ck_compute(p, w[-1])
