from permutations import Permutation


def test_transitions_finite():
    n = 4
    for w in Permutation.all(n):
        for r in range(1, n + 1):
            phi_plus = {w * Permutation.transposition(r, j) for j in w.upper_transitions(r)}
            phi_minus = {w * Permutation.transposition(i, r) for i in w.lower_transitions(r)}

            assert all(w in v.bruhat_covers() for v in phi_plus)
            assert all(w in v.bruhat_covers() for v in phi_minus)


def test_transitions_all():
    n = 4
    for w in Permutation.all(n):
        for r in range(1, n + 1):
            phi_plus = {w * Permutation.transposition(r, j, n) for j in w.upper_transitions(r)}
            phi_minus = {w * Permutation.transposition(i, r, n) for i in w.lower_transitions(r)}

            assert all(w in v.bruhat_covers() for v in phi_plus)
            assert all(w in v.bruhat_covers() for v in phi_minus)

            print (w.get_reduced_word())
            print(r)
            print([x.get_reduced_word() for x in phi_plus])
            print([x.get_reduced_word() for x in phi_minus])


def test_fpf_transitions_simple():
    y = Permutation(2, 1, 4, 3)
    p = 1
    q = 2

    def t(i, j):
        return Permutation.transposition(i, j)

    phi_plus = {t(q, j) * y * t(q, j) for j in y.upper_fpf_involution_transitions(q)}
    assert phi_plus == {Permutation(3, 4, 1, 2)}
