from vectors import Vector
from tableaux import Tableau
from states import FState
from collections import defaultdict


class Word:

    @classmethod
    def _check_record(cls, word, record):
        record = tuple(range(1, len(word) + 1)) if record is None else record
        assert len(word) == len(record)
        assert all(record[i] <= record[i + 1] for i in range(len(record) - 1))
        assert all(record[i] < record[i + 1] for i in range(len(record) - 1) if word[i] > word[i + 1])
        mapping = {}
        for i, a in enumerate(record):
            mapping[i + 1] = a
            mapping[-i - 1] = -a
        return mapping

    @classmethod
    def _check_recording_tableau(cls, recording_tableau):
        standard_tableau = recording_tableau.standardize()
        record = len(recording_tableau) * [0]
        for i, j, value in standard_tableau:
            for k, v in enumerate(value):
                record[abs(v) - 1] = abs(recording_tableau.get(i, j, unpack=False)[k])
        return standard_tableau, tuple(record)

    @classmethod
    def hecke_insertion(cls, word, record=None):
        record = cls._check_record(word, record)

    @classmethod
    def orthogonal_hecke_insertion(cls, word, record=None):
        record = cls._check_record(word, record)

    @classmethod
    def symplectic_hecke_insertion(cls, word, record=None):
        record = cls._check_record(word, record)
        insertion_tableau, recording_tableau = FState.insertion_tableaux(*word)
        recording_tableau = Tableau(mapping={
            (i, j): tuple(record[v] for v in value)
            for i, j, value in recording_tableau
        })
        return insertion_tableau, recording_tableau

    @classmethod
    def inverse_hecke_insertion(cls, insertion_tableau, recording_tableau):
        pass

    @classmethod
    def inverse_orthogonal_hecke_insertion(cls, insertion_tableau, recording_tableau):
        pass

    @classmethod
    def inverse_symplectic_hecke_insertion(cls, insertion_tableau, recording_tableau):
        standard_tableau, record = cls._check_recording_tableau(recording_tableau)
        word = FState.inverse_insertion(insertion_tableau, standard_tableau)
        return word, record

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
