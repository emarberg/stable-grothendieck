from vectors import Vector
from polynomials import Polynomial


def r(i):
    return BoxOperator(0, i)


def c(i):
    return BoxOperator(0, -i)


def box_operator_printer(tup):
    b, tup = tup[0], tup[1:]
    beta = '' if b == 0 else 'β' if b == 1 else 'β^%i' % b
    ans = [beta] + ['r%i' % i if i >= 0 else 'c%i' % -i for i in tup]
    return (' '.join(ans)).strip()


def box_operator_multipler(a, b):
    beta = (a[0] + b[0],)
    for tup, coeff in box_operator_relations(a[1:] + b[1:]):
        yield beta + tup, coeff


def box_operator_relations(tup):
    if any(tup[i] == tup[i + 1] for i in range(len(tup) - 1)):
        return
    if any(tup[i] == -tup[i + 1] and abs(tup[i]) >= 0 for i in range(len(tup) - 1)):
        return
    yield tup, 1


class BoxOperator:
    def __init__(self, *args):
        if len(args) == 1 and type(args[0]) in [list, tuple]:
            args = args[0]
        assert len(args) == 0 or args[0] >= 0
        dictionary = {tuple(args): 1} if args else {}
        self.vector = Vector(dictionary, printer=box_operator_printer, multiplier=box_operator_multipler)

    def __repr__(self):
        return repr(self.vector)

    def __add__(self, other):
        assert type(other) == BoxOperator
        ans = BoxOperator()
        ans.vector = self.vector + other.vector
        return ans

    def __sub__(self, other):
        assert type(other) == BoxOperator
        ans = BoxOperator()
        ans.vector = self.vector - other.vector
        return ans

    def __mul__(self, other):
        ans = BoxOperator()
        if type(other) in [int, Polynomial]:
            ans.vector = self.vector * other
            return ans
        if type(other) == str:
            other = other.lower()
            assert other and other[0] in ['r', 'c']
            if other[0] == 'r':
                return self * BoxOperator(0, int(other[1:]))
            else:
                return self * BoxOperator(0, -int(other[1:]))
        if type(other) == BoxOperator:
            ans.vector = self.vector * other.vector
            return ans
        raise Exception

    def __rmul__(self, other):
        if type(other) == BoxOperator:
            return other * self
        if type(other) == str:
            other = other.lower()
            assert other and other[0] in ['r', 'c']
            if other[0] == 'r':
                return BoxOperator(0, int(other[1:])) * self
            else:
                return BoxOperator(0, -int(other[1:])) * self
        return self * other

    def __eq__(self, other):
        assert type(other) == BoxOperator
        return self.vector == other.vector
