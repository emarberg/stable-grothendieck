from polynomials import Polynomial


class Vector:

    def __init__(self, dictionary={}, printer=None, multiplier=None):
        self.dictionary = {key: value for key, value in dictionary.items() if value}
        self.printer = printer
        self.multiplier = multiplier

    @classmethod
    def base(cls, key, printer=None, multiplier=None):
        return cls({key: 1}, printer, multiplier)

    def keys(self):
        return self.dictionary.keys()

    def values(self):
        return self.dictionary.values()

    def items(self):
        return self.dictionary.items()

    def is_singleton(self):
        if len(self.dictionary) != 1:
            return False
        return list(self.dictionary.values())[0] == 1

    def is_nonnegative(self):
        return all(v > 0 for v in self.values())

    def is_positive(self):
        return not self.is_zero() and self.is_nonnegative()

    def __len__(self):
        return len(self.dictionary)

    def __eq__(self, other):
        return len((self - other).dictionary) == 0

    def __iter__(self):
        return self.dictionary.__iter__()

    def __getitem__(self, item):
        return self.dictionary.get(item, 0)

    def __radd__(self, other):
        return self.__add__(other)

    def __add__(self, other):
        if other == 0:
            return self
        if type(other) == type(self):
            keys = self.keys() | other.keys()
            return self.__class__(
                {key: self[key] + other[key] for key in keys},
                self.printer or other.printer,
                self.multiplier or other.multiplier
            )
        else:
            return other.__radd__(self)

    def __sub__(self, other):
        if other == 0:
            return self
        if type(other) == type(self):
            keys = self.keys() | other.keys()
            return self.__class__(
                {key: self[key] - other[key] for key in keys},
                self.printer or other.printer,
                self.multiplier or other.multiplier
            )
        else:
            return other.__rsub__(self)

    def __mul__(self, other):
        if type(other) in [int, Polynomial]:
            return self.__class__(
                {key: self[key] * other for key in self.keys()},
                self.printer,
                self.multiplier
            )
        elif type(other) == type(self):
            ans = {}
            for a, x in self.items():
                for b, y in other.items():
                    for m, coeff in self.multiplier(a, b) if self.multiplier else [(a * b, 1)]:
                        ans[m] = ans.get(m, 0) + x * y * coeff
            return self.__class__(
                ans,
                self.printer or other.printer,
                self.multiplier or other.multiplier
            )
        else:
            return self * self.base(other)

    def __rmul__(self, other):
        return self.__mul__(other)

    def __neg__(self):
        return self.__mul__(-1)

    def __pow__(self, other):
        assert type(other) == int and other > 0
        if other == 1:
            return self
        x = self ** (other // 2)
        if other % 2 == 0:
            return x * x
        else:
            return x * x * self

    def is_zero(self):
        return len(self) == 0

    def _repr_coeff(self, coeff, key):
        if type(coeff) == Polynomial and len(coeff) > 1:
            return ' + (%s)*' % str(coeff)
        try:
            if key == '' and coeff > 0:
                return ' + %s' % coeff
            if key == '' and coeff < 0:
                return ' - %s' % -coeff
            if coeff == 1:
                return ' + '
            elif coeff == -1:
                return ' - '
            elif coeff > 0:
                return ' + %s*' % coeff
            else:
                return ' - %s*' % -coeff
        except TypeError:
            return ' + (%s)*' % str(coeff)

    def _print_sorted(self, sorted_items):
        base = ''.join(self._repr_coeff(value, key) + key for key, value in sorted_items)
        if base.startswith(' + '):
            return base[3:]
        elif base.startswith(' - '):
            return '-' + base[3:]
        else:
            return '0'

    def __repr__(self):
        printer = self.printer or repr
        sorted_items = sorted([(printer(key), value) for key, value in self.items()])
        return self._print_sorted(sorted_items)
