from tableaux import Tableau, ValuedSetTableau
from partitions import Partition


def combine_str(a, *b):
    if not b:
        return str(a)
    elif len(b) > 1:
        return combine_str(a, combine_str(b[0], *b[1:]))
    else:
        b = b[0]
    ziplines = zip(str(a).split("\n"), str(b).split("\n"))
    lines = ["  ->  ".join(pair) if "".join(pair).strip() else "" for pair in ziplines]
    toprint = "\n".join(lines)
    return toprint


def test_simple():
    test = []

    lhs = ValuedSetTableau(Tableau(". 1' 2';. 1' 2'"), Tableau(". 1 1;. 0 0"))
    forward = ValuedSetTableau(Tableau(". 1' 2';. 1' 2'"), Tableau(". 1 1;. 0 0"))
    middle = ValuedSetTableau(Tableau(". 2' 1';. 2' 1'"), Tableau(". 1 1;. 0 0"))
    backward = ValuedSetTableau(Tableau(". 2' 1';. 2' 1'"), Tableau(". 1 1;. 0 0"))
    test += [(lhs, forward, middle, backward)]

    lhs = ValuedSetTableau(Tableau(". 1' 2';. 1' 2'"), Tableau(". 1 0;. 0 0"))
    forward = ValuedSetTableau(Tableau(". 1' 2';. 1' 2'"), Tableau(". 1 0;. 0 0"))
    middle = ValuedSetTableau(Tableau(". 2' 1';. 2' 1'"), Tableau(". 0 1;. 0 0"))
    backward = ValuedSetTableau(Tableau(". 2' 1';. 2' 1'"), Tableau(". 0 1;. 0 0"))
    test += [(lhs, forward, middle, backward)]

    lhs = ValuedSetTableau(Tableau(". 1' 2';. 1' 2"), Tableau(". 1 0;. 0 0"))
    forward = ValuedSetTableau(Tableau(". 1' 2';. 1' 2"), Tableau(". 1 0;. 0 0"))
    middle = ValuedSetTableau(Tableau(". 2' 1';. 2' 2"), Tableau(". 1 0;. 0 0"))
    backward = ValuedSetTableau(Tableau(". 2' 2;. 2' 1'"), Tableau(". 1 0;. 0 0"))
    test += [(lhs, forward, middle, backward)]

    lhs = ValuedSetTableau(Tableau(". 1' 1;. 1' 2'"), Tableau(". 1 0;. 0 0"))
    forward = ValuedSetTableau(Tableau(". 1' 2';. 1' 1"), Tableau(". 1 0;. 0 0"))
    middle = ValuedSetTableau(Tableau(". 2' 1';. 2' 1"), Tableau(". 1 0;. 0 0"))
    backward = ValuedSetTableau(Tableau(". 2' 1';. 2' 1"), Tableau(". 1 0;. 0 0"))
    test += [(lhs, forward, middle, backward)]

    #

    lhs = ValuedSetTableau(Tableau(". 1' 1;. 1' 2"), Tableau(". 1 0;. 0 0"))
    forward = ValuedSetTableau(Tableau(". 1' 1;. 1' 2"), Tableau(". 1 0;. 0 0"))
    middle = ValuedSetTableau(Tableau(". 1' 2;. 1' 1"), Tableau(". 1 0;. 0 0"))
    backward = ValuedSetTableau(Tableau(". 2 2;. 1' 1"), Tableau(". 1 0;. 0 0"))
    test += [(lhs, forward, middle, backward)]

    lhs = ValuedSetTableau(Tableau(". 1' 2';. 1' 2'"), Tableau(". 0 0;. 0 0"))
    forward = ValuedSetTableau(Tableau(". 1' 2';. 1' 2'"), Tableau(". 0 0;. 0 0"))
    middle = ValuedSetTableau(Tableau(". 2' 1';. 2' 1'"), Tableau(". 0 0;. 0 0"))
    backward = ValuedSetTableau(Tableau(". 2' 1';. 2' 1'"), Tableau(". 0 0;. 0 0"))
    test += [(lhs, forward, middle, backward)]

    lhs = ValuedSetTableau(Tableau(". 1' 1;. 1' 2'"), Tableau(". 0 0;. 0 0"))
    forward = ValuedSetTableau(Tableau(". 1' 2';. 1' 1"), Tableau(". 0 0;. 0 0"))
    middle = ValuedSetTableau(Tableau(". 2' 1';. 1' 1"), Tableau(". 0 0;. 0 0"))
    backward = ValuedSetTableau(Tableau(". 2' 1';. 1' 1"), Tableau(". 0 0;. 0 0"))
    test += [(lhs, forward, middle, backward)]

    lhs = ValuedSetTableau(Tableau(". 1' 1;. 1' 2"), Tableau(". 0 0;. 0 0"))
    forward = ValuedSetTableau(Tableau(". 1' 1;. 1' 2"), Tableau(". 0 0;. 0 0"))
    middle = ValuedSetTableau(Tableau(". 1' 2;. 1' 1"), Tableau(". 0 0;. 0 0"))
    backward = ValuedSetTableau(Tableau(". 2 1';. 1' 1"), Tableau(". 0 0;. 0 0"))
    test += [(lhs, forward, middle, backward)]

    #

    lhs = ValuedSetTableau(Tableau(". 1' 2';. 1' 2"), Tableau(". 0 0;. 0 0"))
    forward = ValuedSetTableau(Tableau(". 1' 2';. 1' 2"), Tableau(". 0 0;. 0 0"))
    middle = ValuedSetTableau(Tableau(". 2' 1';. 1' 2"), Tableau(". 0 0;. 0 0"))
    backward = ValuedSetTableau(Tableau(". 2' 1';. 2 1'"), Tableau(". 0 0;. 0 0"))
    test += [(lhs, forward, middle, backward)]

    lhs = ValuedSetTableau(Tableau(". 1' 1;. 1 2'"), Tableau(". 0 0;. 0 0"))
    forward = ValuedSetTableau(Tableau(". 1' 2';. 1 1"), Tableau(". 0 0;. 0 0"))
    middle = ValuedSetTableau(Tableau(". 2' 1';. 1 1"), Tableau(". 0 0;. 0 0"))
    backward = ValuedSetTableau(Tableau(". 2' 1';. 1 1"), Tableau(". 0 0;. 0 0"))
    test += [(lhs, forward, middle, backward)]

    lhs = ValuedSetTableau(Tableau(". 1' 1;. 1 2"), Tableau(". 0 0;. 0 0"))
    forward = ValuedSetTableau(Tableau(". 1' 1;. 1 2"), Tableau(". 0 0;. 0 0"))
    middle = ValuedSetTableau(Tableau(". 1' 2;. 1 1"), Tableau(". 0 0;. 0 0"))
    backward = ValuedSetTableau(Tableau(". 2 1';. 1 1"), Tableau(". 0 0;. 0 0"))
    test += [(lhs, forward, middle, backward)]

    lhs = ValuedSetTableau(Tableau(". 1' 2';. 1 2'"), Tableau(". 0 1;. 0 0"))
    forward = ValuedSetTableau(Tableau(". 1' 2';. 1 1"), Tableau(". 0 0;. 1 0"))
    middle = ValuedSetTableau(Tableau(". 2' 1';. 1 1"), Tableau(". 0 0;. 1 0"))
    backward = ValuedSetTableau(Tableau(". 2' 1';. 1 1"), Tableau(". 0 0;. 1 0"))
    test += [(lhs, forward, middle, backward)]

    #

    lhs = ValuedSetTableau(Tableau(". 1' 2';. 1 2"), Tableau(". 0 0;. 0 0"))
    forward = ValuedSetTableau(Tableau(". 1' 2';. 1 2"), Tableau(". 0 0;. 0 0"))
    middle = ValuedSetTableau(Tableau(". 2' 1';. 2 1"), Tableau(". 0 0;. 0 0"))
    backward = ValuedSetTableau(Tableau(". 2' 1';. 2 1"), Tableau(". 0 0;. 0 0"))
    test += [(lhs, forward, middle, backward)]

    lhs = ValuedSetTableau(Tableau(". 1' 1;. 2' 2"), Tableau(". 0 0;. 0 0"))
    forward = ValuedSetTableau(Tableau(". 1' 1;. 2' 2"), Tableau(". 0 0;. 0 0"))
    middle = ValuedSetTableau(Tableau(". 2' 2;. 1' 1"), Tableau(". 0 0;. 0 0"))
    backward = ValuedSetTableau(Tableau(". 2' 2;. 1' 1"), Tableau(". 0 0;. 0 0"))
    test += [(lhs, forward, middle, backward)]

    lhs = ValuedSetTableau(Tableau(". 1' 1;. 2 2"), Tableau(". 0 0;. 1 0"))
    forward = ValuedSetTableau(Tableau(". 1' 1;. 2 2"), Tableau(". 0 0;. 1 0"))
    middle = ValuedSetTableau(Tableau(". 1' 2;. 1 1"), Tableau(". 0 0;. 1 0"))
    backward = ValuedSetTableau(Tableau(". 2 1';. 1 1"), Tableau(". 0 0;. 1 0"))
    test += [(lhs, forward, middle, backward)]

    lhs = ValuedSetTableau(Tableau(". 1' 1;. 2 2"), Tableau(". 0 0;. 0 0"))
    forward = ValuedSetTableau(Tableau(". 1' 1;. 2 2"), Tableau(". 0 0;. 0 0"))
    middle = ValuedSetTableau(Tableau(". 1' 2;. 2 1"), Tableau(". 0 0;. 0 0"))
    backward = ValuedSetTableau(Tableau(". 2 2;. 1' 1"), Tableau(". 0 0;. 0 0"))
    test += [(lhs, forward, middle, backward)]

    for lhs, expected_forward, expected_middle, expected_backward in test:
        forward = lhs.forward_transition(1)
        middle = forward.middle_transition(1)
        backward = middle.backward_transition(1)
        print(combine_str(lhs, forward, middle, backward))
        print(combine_str(lhs, expected_forward, expected_middle, expected_backward))
        print()
        print()
        assert expected_forward == forward
        assert expected_middle == middle
        assert expected_backward == backward
        assert lhs.is_semistandard([-1, 1, -2, 2])
        assert forward.is_semistandard([-1, None, -2, 1, None, 2])
        assert middle.is_semistandard([-2, None, -1, 2, None, 1])
        assert backward.is_semistandard([-2, 2, -1, 1])


def test_small():
    for mu in Partition.all(10, strict=True):
        for nu in Partition.subpartitions(mu, strict=True):
            test = sorted(ValuedSetTableau.all(2, mu, nu, diagonal_nonprimes=False))
            seen = set()
            for i, vst in enumerate(test):
                print('mu =', mu, 'nu =', nu, 'case:', i)
                f = vst.forward_transition(1)
                m = f.middle_transition(1)
                b = m.backward_transition(1)
                image = vst.transition(1)
                print(combine_str(vst, f, m, b, image))
                print()
                print()
                assert vst.is_semistandard([-1, 1, -2, 2])
                assert f.is_semistandard([-1, None, -2, 1, None, 2])
                assert m.is_semistandard([-2, None, -1, 2, None, 1])
                assert b.is_semistandard([-2, 2, -1, 1])

                assert image.is_semistandard()
                assert image not in seen
                seen.add(image)
