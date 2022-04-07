from tableaux import Tableau
from vst import ValuedSetTableau
from partitions import Partition
import traceback


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


def print_transition(vst, i):
    f = vst.forward_transition(1)
    m = f.middle_transition(1)
    b = m.backward_transition(1)
    image = vst.transition(1)
    post = image.transition(1)
    print(combine_str(vst, f, m, b, image, post))


def test_small(dnp=True):
    n = 2
    for mu in Partition.all(18, strict=True):
        for nu in Partition.subpartitions(mu, strict=True):
            print('mu =', mu, 'nu =', nu)
            test = sorted(ValuedSetTableau.all(n, mu, nu, diagonal_nonprimes=dnp))

            # test = [vst for vst in test if not vst.is_highest_weight()]

            _unseen = set(test)
            images = {}
            multiplicities = {}
            for vst in test:
                image = vst.transition(1)
                multiplicities[image] = multiplicities.get(image, 0) + 1
                key = image.unprime_diagonal()
                images[key] = images.get(key, []) + [vst]
            unseen = {}
            for vst in _unseen:
                key = (tuple(sorted(vst.tableau.boxes)), vst.weight(n))
                unseen[key] = unseen.get(key, []) + [vst]

            seen = {}
            for vst in test:
                try:
                    f = vst.forward_transition(1)
                    m = f.middle_transition(1)
                    b = m.backward_transition(1)
                    image = vst.transition(1)
                    post = image.transition(1)
                    key = image
                    seen[key] = seen.get(key, []) + [vst]

                    p = vst.primed_groups()
                    q = f.primed_groups()

                    assert vst.is_semistandard([-1, 1, -2, 2])
                    assert f.is_semistandard([-1, None, -2, 1, None, 2])
                    assert m.is_semistandard([-2, None, -1, 2, None, 1])
                    assert b.is_semistandard([-2, 2, -1, 1])
                    assert image.is_semistandard()
                    assert tuple(reversed(image.weight(n))) == vst.weight(n)
                    assert all(preimage.unprime_diagonal() == vst.unprime_diagonal() for preimage in seen[key])
                    assert dnp or not image.diagonal_primes()

                    assert len(seen[key]) == 1
                except:
                    print(5 * '\n')

                    for preimage in seen[key]:
                        if preimage.transition(1).transition(1) == preimage:
                            continue
                        print_transition(preimage, 1)
                    print('preimages:', len(seen[key]))
                    print()
                    print()

                    print('alternatives:')
                    alts = images.get(key.unprime_diagonal(), [])
                    for u in alts:
                        if u not in seen[key] or u.transition(1).transition(1) == u:
                            print_transition(u, 1)
                    print()
                    print()

                    print('unseen:')
                    ukey = (tuple(sorted(image.tableau.boxes)), image.weight(n))
                    uns = unseen.get(ukey, [])
                    for u in uns:
                        print_transition(u, 1)

                    print(5 * '\n')

                    traceback.print_exc()
                    assert tuple(reversed(image.weight(n))) == vst.weight(n)

                    input('')
