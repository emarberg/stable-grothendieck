from vectors import Vector
from partitions import Partition 


def apply(*args):
    assert len(args) > 0
    if len(args) == 1:
        return args[0] if type(args[0]) == Vector else Vector({args[0]: 1})
    else:
        vec = args[-1] if type(args[-1]) == Vector else Vector({args[-1]: 1})
        op = args[-2]
        ans = Vector(printer=vec.printer)
        for mu in vec:
            ans += vec[mu] * op(mu)
        return apply(*(list(args)[:-2] + [ans]))


def u(i):
    def op(mu):
        for j in range(len(mu)):
            if mu[j] == i - 1:
                nu = mu[:j] + (i,) + mu[j + 1:]
                return Vector({nu: 1})
        if i == 1:
            nu = mu + (1,)
            return Vector({nu: 1})
        return Vector()
    return op


def d(i):
    def op(mu):
        for j in reversed(range(len(mu))):
            if mu[j] == i:
                nu = (mu[:j] + (i - 1,) + mu[j + 1:]) if i > 1 else mu[:j]
                return Vector({nu: 1})
        return Vector()
    return op


def small_schur_expressions(n):
    yield 'id', lambda x: x
    # for i in range(1, n + 1):
    #     yield 'u_%i' % i, u(i)
    #     yield 'd_%i' % i, d(i)

    for i in range(1, n + 1):
        for j in range(1, n + 1):
            if i != j:
                yield 'u_%i u_%i' % (i, j), lambda x: apply(u(i), u(j), x)
                yield 'd_%i d_%i' % (i, j), lambda x: apply(d(i), d(j), x)
            yield 'u_%i d_%i' % (i, j), lambda x: apply(u(i), d(j), x)
            yield 'd_%i u_%i' % (i, j), lambda x: apply(d(i), u(j), x)

    # for i in range(1, n + 1):
    #     for j in range(1, n + 1):
    #         for k in range(1, n + 1):
    #             yield (i, j, k), lambda x: u(i)(u(j)(u(k)(x)))
    #             yield (i, j, k), lambda x: d(i)(d(j)(d(k)(x)))


def test_find_relations(n=5):
    square = (n + 1) * (n + 1,)
    partitions = list(Partition.subpartitions(square))
    small = {k: v for k, v in small_schur_expressions(n)}
    for k1 in small:
        print()
        for k2 in small:
            if k1 != k2:
                op1 = small[k1]
                op2 = small[k2]
                if all(op1(mu) == op2(mu) for mu in partitions) and not all(op1(mu).is_zero() for mu in partitions):
                    print('\n\n\n\n')
                    for mu in partitions:
                        print(k1, mu, '=', op1(mu), '==', op2(mu), '=', k2, mu, '??', apply(d(2), d(1), mu))
                        print()
                    print('*', k1, '=', k2)
    print()

