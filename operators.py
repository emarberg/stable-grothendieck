from vectors import Vector
from polynomials import Polynomial, X, Y
from partitions import Partition


beta = X(0)
xvar = X()
yvar = Y()


def operator_decorator(func):
    def op(vector):
        if type(vector) == tuple:
            vector = Vector({vector: 1})
        return func(vector)
    return op


def operator_shifted_add(index):
    @operator_decorator
    def op(vector):
        ans = 0
        for mu, coefficient in vector.items():
            row = Partition.find_shifted_inner_corner(mu, abs(index))
            if row is not None:
                ans += Vector({mu: beta * coefficient})
                continue
            row = Partition.find_shifted_outer_corner(mu, abs(index))
            if row is not None:
                nu = Partition.add(mu, row)
                ans += Vector({nu: coefficient})
        return ans
    return op


def operator_AL(n, var=yvar):  # noqa
    @operator_decorator
    def op(vector):
        ans = vector
        for i in range(n + 1):
            ans += var * operator_shifted_add(i)(ans)
        return ans
    return op


def operator_AR(n, var=yvar):  # noqa
    @operator_decorator
    def op(vector):
        ans = vector
        for i in range(n, -1, -1):
            ans += var * operator_shifted_add(i)(ans)
        return ans
    return op


def operator_shifted_row_remove(index):
    @operator_decorator
    def op(vector):
        ans = 0
        for mu, coefficient in vector.items():
            e = 0
            i = abs(index)
            row = Partition.find_shifted_inner_corner(mu, i)
            while i >= 0 and row is not None:
                i = i - 1
                mu = Partition.remove_shifted_inner_corner(mu, row)
                if mu is None:
                    break
                ans += Vector({mu: (-beta) ** e})
                e += 1
        return ans
    return op


def operator_shifted_column_remove(index):
    @operator_decorator
    def op(vector):
        ans = 0
        for mu, coefficient in vector.items():
            e = 0
            i = abs(index)
            row = Partition.find_shifted_inner_corner(mu, i)
            while row is not None:
                i = i + 1
                mu = Partition.remove_inner_corner(mu, row)
                ans += Vector({mu: (-beta) ** e})
                row = Partition.find_shifted_inner_corner(mu, i)
                e += 1
        return ans
    return op


def operator_R(n, var=xvar):  # noqa
    @operator_decorator
    def op(vector):
        ans = vector
        for i in range(n + 1):
            ans += var * operator_shifted_row_remove(i)(ans)
        return ans
    return op


def operator_C(n, var=xvar):  # noqa
    @operator_decorator
    def op(vector):
        ans = vector
        for i in range(n, 0, -1):
            ans += var * operator_shifted_column_remove(i)(ans)
        return ans
    return op


AR = operator_AR
AL = operator_AL
R = operator_R
C = operator_C
