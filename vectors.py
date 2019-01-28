
class Vector:
    def __init__(self, dictionary={}, printer=None):
        self.dictionary = {key: value for key, value in dictionary.items() if value}
        self.printer = printer

    @classmethod
    def base(cls, key, printer=None):
        return cls({key: 1}, printer)

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

    def __add__(self, other):
        if type(other) == type(self):
            keys = self.keys() | other.keys()
            return self.__class__({key: self[key] + other[key] for key in keys}, self.printer or other.printer)
        else:
            return other.__radd__(self)

    def __sub__(self, other):
        if type(other) == type(self):
            keys = self.keys() | other.keys()
            return self.__class__({key: self[key] - other[key] for key in keys}, self.printer or other.printer)
        else:
            return other.__rsub__(self)

    def __mul__(self, other):
        if type(other) == int:
            return self.__class__({key: self[key] * other for key in self.keys()}, self.printer)
        elif type(other) == type(self):
            ans = self.__class__(printer=self.printer or other.printer)
            for a, x in self.items():
                for b, y in other.items():
                    ans += (x * y) * (a * b)
            return ans
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

    def __repr__(self):
        printer = self.printer or repr
        sorted_items = sorted([(printer(key), value) for key, value in self.items()])
        base = ''.join(self._repr_coeff(value, key) + key for key, value in sorted_items)
        if base.startswith(' + '):
            return base[3:]
        elif base.startswith(' - '):
            return '-' + base[3:]
        else:
            return '0'
