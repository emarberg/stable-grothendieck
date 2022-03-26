from tableaux import Tableau, ValuedSetTableau


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
    test += [(lhs, forward, middle)]

    lhs = ValuedSetTableau(Tableau(". 1' 2';. 1' 2'"), Tableau(". 1 0;. 0 0"))
    forward = ValuedSetTableau(Tableau(". 1' 2';. 1' 2'"), Tableau(". 1 0;. 0 0"))
    middle = ValuedSetTableau(Tableau(". 2' 1';. 2' 1'"), Tableau(". 0 1;. 0 0"))
    test += [(lhs, forward, middle)]

    lhs = ValuedSetTableau(Tableau(". 1' 2';. 1' 2"), Tableau(". 1 0;. 0 0"))
    forward = ValuedSetTableau(Tableau(". 1' 2';. 1' 2"), Tableau(". 1 0;. 0 0"))
    middle = ValuedSetTableau(Tableau(". 2' 1';. 2' 2"), Tableau(". 1 0;. 0 0"))
    test += [(lhs, forward, middle)]

    lhs = ValuedSetTableau(Tableau(". 1' 1;. 1' 2'"), Tableau(". 1 0;. 0 0"))
    forward = ValuedSetTableau(Tableau(". 1' 2';. 1' 1"), Tableau(". 1 0;. 0 0"))
    middle = ValuedSetTableau(Tableau(". 2' 1';. 2' 1"), Tableau(". 1 0;. 0 0"))
    test += [(lhs, forward, middle)]

    #

    lhs = ValuedSetTableau(Tableau(". 1' 1;. 1' 2"), Tableau(". 1 0;. 0 0"))
    forward = ValuedSetTableau(Tableau(". 1' 1;. 1' 2"), Tableau(". 1 0;. 0 0"))
    middle = ValuedSetTableau(Tableau(". 1' 2;. 1' 1"), Tableau(". 1 0;. 0 0"))
    test += [(lhs, forward, middle)]

    lhs = ValuedSetTableau(Tableau(". 1' 2';. 1' 2'"), Tableau(". 0 0;. 0 0"))
    forward = ValuedSetTableau(Tableau(". 1' 2';. 1' 2'"), Tableau(". 0 0;. 0 0"))
    middle = ValuedSetTableau(Tableau(". 2' 1';. 2' 1'"), Tableau(". 0 0;. 0 0"))
    test += [(lhs, forward, middle)]

    lhs = ValuedSetTableau(Tableau(". 1' 1;. 1' 2'"), Tableau(". 0 0;. 0 0"))
    forward = ValuedSetTableau(Tableau(". 1' 2';. 1' 1"), Tableau(". 0 0;. 0 0"))
    middle = ValuedSetTableau(Tableau(". 2' 1';. 1' 1"), Tableau(". 0 0;. 0 0"))
    test += [(lhs, forward, middle)]

    lhs = ValuedSetTableau(Tableau(". 1' 1;. 1' 2"), Tableau(". 0 0;. 0 0"))
    forward = ValuedSetTableau(Tableau(". 1' 1;. 1' 2"), Tableau(". 0 0;. 0 0"))
    middle = ValuedSetTableau(Tableau(". 1' 2;. 1' 1"), Tableau(". 0 0;. 0 0"))
    test += [(lhs, forward, middle)]

    #

    lhs = ValuedSetTableau(Tableau(". 1' 2';. 1' 2"), Tableau(". 0 0;. 0 0"))
    forward = ValuedSetTableau(Tableau(". 1' 2';. 1' 2"), Tableau(". 0 0;. 0 0"))
    middle = ValuedSetTableau(Tableau(". 2' 1';. 1' 2"), Tableau(". 0 0;. 0 0"))
    test += [(lhs, forward, middle)]

    lhs = ValuedSetTableau(Tableau(". 1' 1;. 1 2'"), Tableau(". 0 0;. 0 0"))
    forward = ValuedSetTableau(Tableau(". 1' 2';. 1 1"), Tableau(". 0 0;. 0 0"))
    middle = ValuedSetTableau(Tableau(". 2' 1';. 1 1"), Tableau(". 0 0;. 0 0"))
    test += [(lhs, forward, middle)]

    lhs = ValuedSetTableau(Tableau(". 1' 1;. 1 2"), Tableau(". 0 0;. 0 0"))
    forward = ValuedSetTableau(Tableau(". 1' 1;. 1 2"), Tableau(". 0 0;. 0 0"))
    middle = ValuedSetTableau(Tableau(". 1' 2;. 1 1"), Tableau(". 0 0;. 0 0"))
    test += [(lhs, forward, middle)]

    lhs = ValuedSetTableau(Tableau(". 1' 2';. 1 2'"), Tableau(". 0 1;. 0 0"))
    forward = ValuedSetTableau(Tableau(". 1' 2';. 1 1"), Tableau(". 0 0;. 1 0"))
    middle = ValuedSetTableau(Tableau(". 2' 1';. 1 1"), Tableau(". 0 0;. 1 0"))
    test += [(lhs, forward, middle)]

    #

    lhs = ValuedSetTableau(Tableau(". 1' 2';. 1 2"), Tableau(". 0 0;. 0 0"))
    forward = ValuedSetTableau(Tableau(". 1' 2';. 1 2"), Tableau(". 0 0;. 0 0"))
    middle = ValuedSetTableau(Tableau(". 2' 1';. 2 1"), Tableau(". 0 0;. 0 0"))
    test += [(lhs, forward, middle)]

    lhs = ValuedSetTableau(Tableau(". 1' 1;. 2' 2"), Tableau(". 0 0;. 0 0"))
    forward = ValuedSetTableau(Tableau(". 1' 1;. 2' 2"), Tableau(". 0 0;. 0 0"))
    middle = ValuedSetTableau(Tableau(". 2' 2;. 1' 1"), Tableau(". 0 0;. 0 0"))
    test += [(lhs, forward, middle)]

    lhs = ValuedSetTableau(Tableau(". 1' 1;. 2 2"), Tableau(". 0 0;. 1 0"))
    forward = ValuedSetTableau(Tableau(". 1' 1;. 2 2"), Tableau(". 0 0;. 1 0"))
    middle = ValuedSetTableau(Tableau(". 1' 2;. 1 1"), Tableau(". 0 0;. 1 0"))
    test += [(lhs, forward, middle)]

    lhs = ValuedSetTableau(Tableau(". 1' 1;. 2 2"), Tableau(". 0 0;. 0 0"))
    forward = ValuedSetTableau(Tableau(". 1' 1;. 2 2"), Tableau(". 0 0;. 0 0"))
    middle = ValuedSetTableau(Tableau(". 1' 2;. 2 1"), Tableau(". 0 0;. 0 0"))
    test += [(lhs, forward, middle)]

    for lhs, expected_forward, expected_middle in test:
        forward = lhs.forward_transition(1)
        middle = forward.middle_transition(1)
        print(combine_str(lhs, forward, middle))
        print(combine_str(lhs, expected_forward, expected_middle))
        print()
        print()
        assert expected_forward == forward
        assert expected_middle == middle


def test_forward():
    mu=(3, 2)
    nu=(1,)
    test = sorted(ValuedSetTableau.all(2, mu, nu))
    for i, vst in enumerate(test):
        forward_vst = vst.forward_transition(1)
        ziplines = zip(str(vst).split("\n"), str(forward_vst).split("\n"))
        lines = ["  ->  ".join(pair) for pair in ziplines]
        toprint = "\n".join(lines[2:-2])
        print('case:', i)
        print(combine_str(vst, forward_vst))
        print()
        print()
