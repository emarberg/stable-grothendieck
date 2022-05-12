from tableaux import Tableau
from vst import combine_str, ValuedSetTableau
from partitions import Partition
from utils import jp_expansion, jp, jq_expansion, jq
import traceback


def is_first_in_group(boxes, i):
    x, y = boxes[i]
    return i == 0 or (boxes[i - 1][0] != x and boxes[i - 1][1] != y)


def is_second_in_group(boxes, i):
    return not is_first_in_group(boxes, i) and is_first_in_group(boxes, i - 1)


def is_last_in_group(boxes, i):
    x, y = boxes[i]
    return i + 1 == len(boxes) or (boxes[i + 1][0] != x and boxes[i + 1][1] != y)


def is_only_in_group(boxes, i):
    return is_first_in_group(boxes, i) and is_last_in_group(boxes, i)


def is_p_corner_box(boxes, i):
    x, y = boxes[i]
    return ((x, y - 1) in boxes and (x + 1, y - 1) in boxes) or ((x + 1, y) in boxes and (x + 1, y - 1) in boxes) or (x + 1 == y and (x + 1, y) in boxes)


def is_q_corner_box(boxes, i):
    x, y = boxes[i]
    return ((x, y - 1) in boxes and (x + 1, y - 1) in boxes) or ((x + 1, y) in boxes and (x + 1, y - 1) in boxes)


def is_first_component(boxes, i):
    u, v = boxes[i]
    for j in range(i - 1, -1, -1):
        a, b = boxes[j]
        if a != u and b != v:
            return False
        u, v = a, b
    return True


def is_first_component_with_multiple_boxes(boxes, i):
    if is_only_in_group(boxes, i):
        return False
    u, v = boxes[i]
    for j in range(i - 1, -1, -1):
        a, b = boxes[j]
        if a != u and b != v:
            return all(is_only_in_group(boxes, k) for k in range(j + 1))
        u, v = a, b
    return True


def is_last_component(boxes, i):
    u, v = boxes[i]
    for j in range(i + 1, len(boxes)):
        a, b = boxes[j]
        if a != u and b != v:
            return False
        u, v = a, b
    return True


def comarked_q_ribbons(nu, lam):
    def generate(boxes, i=0):
        if i == len(boxes):
            yield ()
            return

        first_component_with_multiple_boxes = is_first_component_with_multiple_boxes(boxes, i)
        first_in_group = is_first_in_group(boxes, i)
        second_in_group = is_second_in_group(boxes, i)
        only_in_group = is_only_in_group(boxes, i)
        corner_box = is_q_corner_box(boxes, i)

        choices = [1, 2] if corner_box else [2]
        if first_component_with_multiple_boxes:
            pass
        else:
            if only_in_group:
                choices = [1, -1, 2, -2]
            elif first_in_group:
                choices = [1, 2, -2]
            elif second_in_group:
                choices = [1, 2]

        x, y = boxes[i]
        for c in choices:
            for rest in generate(boxes, i + 1):
                yield (((x, y), c),) + rest

    ans = {}
    boxes = sorted(Partition.shifted_shape(nu, lam), key=lambda xy: (-xy[0], xy[1]))
    for tup in generate(boxes):
        tab = Tableau({box: val for box, val in tup})
        weight = tab.abs_sum() - tab.size()
        ans[weight] = ans.get(weight, []) + [tab]

    decomp = {key: len(val) for key, val in ans.items()}
    expected = {sum(key): val.substitute(0, 1) for key, val in jq_expansion(jq(1, nu, lam)).items()}
    return ans, decomp, expected


def test_comarked_q_ribbons(n=7):
    for nu, lam in Partition.shifted_ribbons(n):
        boxes = sorted(Partition.shifted_shape(nu, lam), key=lambda xy: (-xy[0], xy[1]))
        if all(is_only_in_group(boxes, i) for i in range(len(boxes))):
            continue
        ans, decomp, expected = comarked_q_ribbons(nu, lam)
        print('    nu =', nu, 'lambda =', lam)
        if decomp != expected:
            Partition.print_shifted(nu, lam)
            print()
            print('  decomp:', decomp)
            print('expected:', expected)
            print()
            print(ans)
            input('\n')
        # assert decomp == expected


def comarked_p_ribbons(nu, lam):
    def generate(boxes, i=0):
        if i == len(boxes):
            yield ()
            return

        first_component = is_first_component(boxes, i)
        first_component_on_diagonal = boxes[0][0] == boxes[0][1]

        first_in_group = is_first_in_group(boxes, i)
        second_in_group = is_second_in_group(boxes, i)
        only_in_group = is_only_in_group(boxes, i)
        corner_box = is_p_corner_box(boxes, i)

        choices = [1, 2] if corner_box else [2]
        if first_component and first_component_on_diagonal:
            pass
        elif first_component and not first_component_on_diagonal:
            if only_in_group:
                choices = [2, -2]
            elif first_in_group:
                choices = [1, 2, -2]
        else:
            if only_in_group:
                choices = [1, -1, 2, -2]
            elif first_in_group:
                choices = [1, 2, -2]
            elif second_in_group:
                choices = [1, 2]

        x, y = boxes[i]
        for c in choices:
            for rest in generate(boxes, i + 1):
                yield (((x, y), c),) + rest

    ans = {}
    boxes = sorted(Partition.shifted_shape(nu, lam), key=lambda xy: (-xy[0], xy[1]))
    for tup in generate(boxes):
        tab = Tableau({box: val for box, val in tup})
        weight = tab.abs_sum() - tab.size()
        ans[weight] = ans.get(weight, []) + [tab]

    decomp = {key: len(val) for key, val in ans.items()}
    expected = {sum(key): val.substitute(0, 1) for key, val in jp_expansion(jp(1, nu, lam)).items()}
    return ans, decomp, expected


def test_comarked_p_ribbons(n=7):
    for nu, lam in Partition.shifted_ribbons(n):
        ans, decomp, expected = comarked_p_ribbons(nu, lam)
        print('    nu =', nu, 'lambda =', lam)
        if decomp != expected:
            Partition.print_shifted(nu, lam)
            print()
            print('  decomp:', decomp)
            print('expected:', expected)
            print()
            print(ans)
            input('\n')
        # assert decomp == expected


# def test_simple():
#     test = []

#     lhs =      ValuedSetTableau(Tableau(". 1' 2';. 1' 2'"), Tableau(". 0 0;. 1 1"))
#     forward =  ValuedSetTableau(Tableau(". 1' 2';. 1' 2'"), Tableau(". 0 0;. 1 1"))
#     middle =   ValuedSetTableau(Tableau(". 2' 1';. 2' 1'"), Tableau(". 0 0;. 1 1"))
#     backward = ValuedSetTableau(Tableau(". 2' 1';. 2' 1'"), Tableau(". 0 0;. 1 1"))
#     test += [(lhs, forward, middle, backward)]

#     lhs =      ValuedSetTableau(Tableau(". 1' 2';. 1' 2'"), Tableau(". 0 0;. 1 0"))
#     forward =  ValuedSetTableau(Tableau(". 1' 2';. 1' 2'"), Tableau(". 0 0;. 1 0"))
#     middle =   ValuedSetTableau(Tableau(". 2' 1';. 2' 1'"), Tableau(". 0 0;. 0 1"))
#     backward = ValuedSetTableau(Tableau(". 2' 1';. 2' 1'"), Tableau(". 0 0;. 0 1"))
#     test += [(lhs, forward, middle, backward)]

#     lhs =      ValuedSetTableau(Tableau(". 1' 2';. 1' 2"), Tableau(". 0 0;. 1 0"))
#     forward =  ValuedSetTableau(Tableau(". 1' 2';. 1' 2"), Tableau(". 0 0;. 1 0"))
#     middle =   ValuedSetTableau(Tableau(". 2' 1';. 2' 2"), Tableau(". 0 0;. 1 0"))
#     backward = ValuedSetTableau(Tableau(". 2' 2;. 2' 1'"), Tableau(". 0 0;. 1 0"))
#     test += [(lhs, forward, middle, backward)]

#     lhs =      ValuedSetTableau(Tableau(". 1' 1;. 1' 2'"), Tableau(". 0 0;. 1 0"))
#     forward =  ValuedSetTableau(Tableau(". 1' 2';. 1' 1"), Tableau(". 0 0;. 1 0"))
#     middle =   ValuedSetTableau(Tableau(". 2' 1';. 2' 1"), Tableau(". 0 0;. 1 0"))
#     backward = ValuedSetTableau(Tableau(". 2' 1';. 2' 1"), Tableau(". 0 0;. 1 0"))
#     test += [(lhs, forward, middle, backward)]

#     #

#     lhs =      ValuedSetTableau(Tableau(". 1' 1;. 1' 2"), Tableau(". 0 0;. 1 0"))
#     forward =  ValuedSetTableau(Tableau(". 1' 1;. 1' 2"), Tableau(". 0 0;. 1 0"))
#     middle =   ValuedSetTableau(Tableau(". 1' 2;. 1' 1"), Tableau(". 0 0;. 1 0"))
#     backward = ValuedSetTableau(Tableau(". 2 2;. 1' 1"), Tableau(". 1 0;. 0 0"))
#     test += [(lhs, forward, middle, backward)]

#     lhs =      ValuedSetTableau(Tableau(". 1' 2';. 1' 2'"), Tableau(". 0 0;. 0 0"))
#     forward =  ValuedSetTableau(Tableau(". 1' 2';. 1' 2'"), Tableau(". 0 0;. 0 0"))
#     middle =   ValuedSetTableau(Tableau(". 2' 1';. 2' 1'"), Tableau(". 0 0;. 0 0"))
#     backward = ValuedSetTableau(Tableau(". 2' 1';. 2' 1'"), Tableau(". 0 0;. 0 0"))
#     test += [(lhs, forward, middle, backward)]

#     lhs =      ValuedSetTableau(Tableau(". 1' 1;. 1' 2'"), Tableau(". 0 0;. 0 0"))
#     forward =  ValuedSetTableau(Tableau(". 1' 2';. 1' 1"), Tableau(". 0 0;. 0 0"))
#     middle =   ValuedSetTableau(Tableau(". 2' 1';. 1' 1"), Tableau(". 0 0;. 0 0"))
#     backward = ValuedSetTableau(Tableau(". 2' 1';. 1' 1"), Tableau(". 0 0;. 0 0"))
#     test += [(lhs, forward, middle, backward)]

#     lhs =      ValuedSetTableau(Tableau(". 1' 1;. 1' 2"), Tableau(". 0 0;. 0 0"))
#     forward =  ValuedSetTableau(Tableau(". 1' 1;. 1' 2"), Tableau(". 0 0;. 0 0"))
#     middle =   ValuedSetTableau(Tableau(". 1' 2;. 1' 1"), Tableau(". 0 0;. 0 0"))
#     backward = ValuedSetTableau(Tableau(". 2 1';. 1' 1"), Tableau(". 0 0;. 0 0"))
#     test += [(lhs, forward, middle, backward)]

#     #

#     lhs =      ValuedSetTableau(Tableau(". 1' 2';. 1' 2"), Tableau(". 0 0;. 0 0"))
#     forward =  ValuedSetTableau(Tableau(". 1' 2';. 1' 2"), Tableau(". 0 0;. 0 0"))
#     middle =   ValuedSetTableau(Tableau(". 2' 1';. 1' 2"), Tableau(". 0 0;. 0 0"))
#     backward = ValuedSetTableau(Tableau(". 2' 1';. 2 1'"), Tableau(". 0 0;. 0 0"))
#     test += [(lhs, forward, middle, backward)]

#     lhs =      ValuedSetTableau(Tableau(". 1' 1;. 1 2'"), Tableau(". 0 0;. 0 0"))
#     forward =  ValuedSetTableau(Tableau(". 1' 2';. 1 1"), Tableau(". 0 0;. 0 0"))
#     middle =   ValuedSetTableau(Tableau(". 2' 1';. 1 1"), Tableau(". 0 0;. 0 0"))
#     backward = ValuedSetTableau(Tableau(". 2' 1';. 1 1"), Tableau(". 0 0;. 0 0"))
#     test += [(lhs, forward, middle, backward)]

#     lhs =      ValuedSetTableau(Tableau(". 1' 1;. 1 2"), Tableau(". 0 0;. 0 0"))
#     forward =  ValuedSetTableau(Tableau(". 1' 1;. 1 2"), Tableau(". 0 0;. 0 0"))
#     middle =   ValuedSetTableau(Tableau(". 1' 2;. 1 1"), Tableau(". 0 0;. 0 0"))
#     backward = ValuedSetTableau(Tableau(". 2 1';. 1 1"), Tableau(". 0 0;. 0 0"))
#     test += [(lhs, forward, middle, backward)]

#     lhs =      ValuedSetTableau(Tableau(". 1' 2';. 1 2'"), Tableau(". 0 0;. 0 1"))
#     forward =  ValuedSetTableau(Tableau(". 1' 2';. 1 1"), Tableau(". 0 0;. 1 0"))
#     middle =   ValuedSetTableau(Tableau(". 2' 1';. 1 1"), Tableau(". 0 0;. 1 0"))
#     backward = ValuedSetTableau(Tableau(". 2' 1';. 1 1"), Tableau(". 0 0;. 1 0"))
#     test += [(lhs, forward, middle, backward)]

#     #

#     lhs =      ValuedSetTableau(Tableau(". 1' 2';. 1 2"), Tableau(". 0 0;. 0 0"))
#     forward =  ValuedSetTableau(Tableau(". 1' 2';. 1 2"), Tableau(". 0 0;. 0 0"))
#     middle =   ValuedSetTableau(Tableau(". 2' 1';. 2 1"), Tableau(". 0 0;. 0 0"))
#     backward = ValuedSetTableau(Tableau(". 2' 1';. 2 1"), Tableau(". 0 0;. 0 0"))
#     test += [(lhs, forward, middle, backward)]

#     lhs =      ValuedSetTableau(Tableau(". 1' 1;. 2' 2"), Tableau(". 0 0;. 0 0"))
#     forward =  ValuedSetTableau(Tableau(". 1' 1;. 2' 2"), Tableau(". 0 0;. 0 0"))
#     middle =   ValuedSetTableau(Tableau(". 2' 2;. 1' 1"), Tableau(". 0 0;. 0 0"))
#     backward = ValuedSetTableau(Tableau(". 2' 2;. 1' 1"), Tableau(". 0 0;. 0 0"))
#     test += [(lhs, forward, middle, backward)]

#     lhs =      ValuedSetTableau(Tableau(". 1' 1;. 2 2"), Tableau(". 0 0;. 1 0"))
#     forward =  ValuedSetTableau(Tableau(". 1' 1;. 2 2"), Tableau(". 0 0;. 1 0"))
#     middle =   ValuedSetTableau(Tableau(". 1' 2;. 1 1"), Tableau(". 0 0;. 1 0"))
#     backward = ValuedSetTableau(Tableau(". 2 1';. 1 1"), Tableau(". 0 0;. 1 0"))
#     test += [(lhs, forward, middle, backward)]

#     lhs =      ValuedSetTableau(Tableau(". 1' 1;. 2 2"), Tableau(". 0 0;. 0 0"))
#     forward =  ValuedSetTableau(Tableau(". 1' 1;. 2 2"), Tableau(". 0 0;. 0 0"))
#     middle =   ValuedSetTableau(Tableau(". 1' 2;. 2 1"), Tableau(". 0 0;. 0 0"))
#     backward = ValuedSetTableau(Tableau(". 2 2;. 1' 1"), Tableau(". 0 0;. 0 0"))
#     test += [(lhs, forward, middle, backward)]

#     for lhs, expected_forward, expected_middle, expected_backward in test:
#         forward = lhs.forward_transition(1)
#         middle, case = forward.middle_transition(1)
#         backward = middle.backward_transition(1)
#         image = lhs.transition(1)
#         post = image.transition(1)
#         print(combine_str(lhs, forward, middle, backward, image, post))
#         print(combine_str(lhs, expected_forward, expected_middle, expected_backward))
#         print()
#         print()
#         assert expected_forward == forward
#         assert expected_middle == middle
#         assert expected_backward == backward
#         assert lhs.is_semistandard()
#         assert forward.is_interstandard()
#         assert middle.is_interstandard()
#         assert backward.is_semistandard()
#         assert lhs == post


def print_transition(vst, i):
    altered = vst.is_altered(i)
    f = vst.forward_transition(i)
    print(combine_str(vst, f))
    m, case = f.middle_transition(i)
    b = m.backward_transition(i)
    image = vst.transition(i)
    post = image.transition(i)
    print(combine_str(vst, f, m, b, image, post))
    print('is altered:', altered, '->', image.is_altered(i), '| middle case:', case)


def _test_small(k=9, dnp=True, verbose=False):
    n = 2
    for mu in Partition.all(k, strict=True):
        for nu in Partition.subpartitions(mu, strict=True):
            print('mu =', mu, 'nu =', nu)
            test = sorted(ValuedSetTableau.all(n, mu, nu, diagonal_primes=dnp))

            if verbose:
                _unseen = set(test)
                images = {}
                multiplicities = {}
                for vst in test:
                    image = vst.transition(1)
                    _unseen.discard(image)
                    multiplicities[image] = multiplicities.get(image, 0) + 1
                    key = image.unprime()
                    images[key] = images.get(key, []) + [vst]
                unseen = {}
                for vst in _unseen:
                    key = vst.unprime_diagonal()
                    unseen[key] = unseen.get(key, []) + [vst]

            seen = {}
            for vst in test:
                try:
                    f = vst.forward_transition(1)
                    m, case = f.middle_transition(1)
                    b = m.backward_transition(1)
                    image = vst.transition(1)
                    post = image.transition(1)

                    assert vst.is_semistandard()
                    assert f.is_interstandard()
                    assert m.is_interstandard()
                    assert b.is_semistandard()
                    assert image.is_semistandard()
                    assert tuple(reversed(image.weight(n))) == vst.weight(n)
                    assert dnp or not image.diagonal_primes()
                    assert vst.diagonal_primes() == image.diagonal_primes()
                    assert vst == post

                    if verbose:
                        key = image
                        seen[key] = seen.get(key, []) + [vst]
                        assert len(seen[key]) == 1
                except:
                    print('\nvst:')
                    print_transition(vst, 1)
                    print_transition(image, 1)
                    print_transition(post, 1)
                    print_transition(post.transition(1), 1)
                    print(5 * '\n')
                    input('')

                    if verbose:
                        for preimage in seen[key]:
                            print_transition(preimage, 1)
                        print('preimages:', len(seen[key]))
                        print()
                        print('alternatives:')
                        alts = images.get(key.unprime(), [])
                        for u in alts:
                            if u not in seen[key]:
                                print_transition(u, 1)
                        print()
                        print('unseen:')
                        ukey = image.unprime_diagonal()
                        uns = unseen.get(ukey, [])
                        for u in uns:
                            print_transition(u, 1)

                    print(5 * '\n')

                    traceback.print_exc()
                    assert tuple(reversed(image.weight(n))) == vst.weight(n)
                    assert vst == post


def test_small_p(k=8, verbose=False):
    _test_small(k, False, verbose)


def test_small_q(k=8, verbose=False):
    _test_small(k, True, verbose)


def test_interstandard(k=8):
    for mu in Partition.all(k, strict=True):
        for nu in Partition.subpartitions(mu, strict=True):
            print('mu =', mu, 'nu =', nu)
            source = sorted(ValuedSetTableau.all(2, mu, nu, diagonal_primes=True))
            target = sorted(ValuedSetTableau.all_interstandard(mu, nu))

            for vst in source:
                try:
                    image = vst.forward_transition(1)
                    post = image.backward_transition(1)
                    assert vst.is_semistandard()
                    assert image.is_interstandard()
                    assert post.is_semistandard()
                    assert vst == post
                except:
                    print('\nswap direction:\n')
                    print(combine_str(vst, image, post))
                    traceback.print_exc()
                    input('\n?\n')

            for vst in target:
                try:
                    image = vst.backward_transition(1)
                    post = image.forward_transition(1)
                    assert vst.is_interstandard()
                    assert image.is_semistandard()
                    assert post.is_interstandard()
                    assert vst == post
                except:
                    print('\nunswap direction:\n')
                    print(combine_str(vst, image, post))
                    traceback.print_exc()
                    input('\n?\n')
