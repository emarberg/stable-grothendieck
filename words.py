from vectors import Vector
from collections import defaultdict


class Word:

    @classmethod
    def ascents(cls, v):
        return len([i for i in range(len(v) - 1) if v[i] <= v[i + 1]])

    @classmethod
    def peaks(cls, v):
        s = [i for i in range(1, len(v) - 1) if v[i - 1] <= v[i] >= v[i + 1]]
        return (1 + cls.peaks(v[s[0] + 1:])) if s else (1 if v else 0)

    @classmethod
    def increasing_zeta(cls, word):
        return 1 if all(word[i] < word[i + 1] for i in range(len(word) - 1)) else 0

    @classmethod
    def decreasing_zeta(cls, word):
        return 1 if all(word[i] > word[i + 1] for i in range(len(word) - 1)) else 0

    @classmethod
    def weakly_increasing_zeta(cls, word):
        return 1 if all(word[i] <= word[i + 1] for i in range(len(word) - 1)) else 0

    @classmethod
    def weakly_decreasing_zeta(cls, word):
        return 1 if all(word[i] >= word[i + 1] for i in range(len(word) - 1)) else 0

    @classmethod
    def unimodal_zeta(cls, word):
        return sum([cls.decreasing_zeta(word[:i]) * cls.increasing_zeta(word[i:]) for i in range(len(word) + 1)])

    @classmethod
    def right_weakly_unimodal_zeta(cls, word):
        return sum([cls.decreasing_zeta(word[:i]) * cls.weakly_increasing_zeta(word[i:]) for i in range(len(word) + 1)])

    @classmethod
    def left_weakly_unimodal_zeta(cls, word):
        return sum([cls.weakly_decreasing_zeta(word[:i]) * cls.increasing_zeta(word[i:]) for i in range(len(word) + 1)])

    @classmethod
    def weakly_unimodal_zeta(cls, word):
        return sum([cls.weakly_decreasing_zeta(word[:i]) * cls.weakly_increasing_zeta(word[i:]) for i in range(len(word) + 1)])

    @classmethod
    def quasisymmetrize(cls, word, zeta):
        if len(word) == 0:
            return Vector({(): 1})
        ans = defaultdict(int)
        for i in range(1, len(word) + 1):
            pre = word[:i]
            sub = word[i:]
            z = zeta(pre)
            if z != 0:
                for k, v in cls.quasisymmetrize(sub, zeta).items():
                    ans[(i,) + k] += z * v
        ans = Vector({k: v for k, v in ans.items() if v != 0})
        return ans
