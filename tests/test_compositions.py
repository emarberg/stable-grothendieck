from partitions import Partition
from polynomials import beta
from vectors import Vector


def is_peak_composition(alpha):
    return all(a > 0 for a in alpha) and all(a >= 2 for a in alpha[:-1])


def I(alpha):
    des = [0]
    for a in alpha[:-1]:
        des.append(des[-1] + a)
    des = des[1:]
    return des


def Lambda(alpha):
    if sum(alpha) == 0:
        return ()
    des = I(alpha)
    des = [0] + [i for i in des if 0 < i - 1 and i - 1 not in des] + [sum(alpha)]
    return tuple(des[i] - des[i - 1] for i in range(1, len(des)))


def ribbon_multiplier(alpha, gamma):
    if sum(alpha) == 0:
        yield gamma, 1
    elif sum(gamma) == 0:
        yield alpha, 1
    else:
        a = alpha + gamma
            
        b = list(alpha)
        b[-1] += gamma[0] 
        b += list(gamma[1:])
        b = tuple(b)

        c = list(alpha)
        c[-1] += gamma[0] - 1
        c += list(gamma[1:])
        c = tuple(c)

        yield a, 1
        yield b, 1
        yield c, beta


def shribbon_multiplier(alpha, gamma):
    if sum(alpha) == 0:
        yield gamma, 1
    elif sum(gamma) == 0:
        yield alpha, 1
    else:
        a = alpha + gamma
            
        b = list(alpha)
        b[-1] += gamma[0] 
        b += list(gamma[1:])
        b = tuple(b)

        c = list(alpha)
        c[-1] += gamma[0] - 1
        c += list(gamma[1:])
        c = tuple(c)

        g = list(gamma)
        g[0] -= 1

        d = list(alpha)
        d[-1] += 1
        d += g
        d = tuple(d)

        e = list(alpha)
        e += g
        e = tuple(e)

        if is_peak_composition(a):
            yield a, 1
        if is_peak_composition(b):
            yield b, 1
        if is_peak_composition(c):
            yield c, beta
        if is_peak_composition(d):
            yield d, 1
        if is_peak_composition(e):
            yield e, beta


def Ribbon(alpha, multiplier=ribbon_multiplier):
    return Vector({alpha: 1}, multiplier=multiplier)


def ShRibbon(alpha):
    ans = 0
    n = sum(alpha)
    for gamma in Partition.compositions(n):
        if Lambda(gamma) == alpha:
            ans += Ribbon(gamma)
    return ans


def ShRibbon_expansion(f):
    ans = 0
    while f != 0:
        alpha = max(f)
        coeff = f[alpha]
        ans += Ribbon(alpha) * coeff
        f -= ShRibbon(alpha) * coeff
    return ans


def test_shribbon_product(n=5):
    for p in range(1, n + 1):
        for q in range(1, n + 1):
            for alpha in Partition.peak_compositions(p):
                for gamma in Partition.peak_compositions(q):
                    expected = Ribbon(alpha, shribbon_multiplier) * Ribbon(gamma, shribbon_multiplier)
                    computed = ShRibbon_expansion(ShRibbon(alpha) * ShRibbon(gamma))
                    if expected != computed:
                        print('alpha =', alpha, 'gamma =', gamma)
                        print('expected =', expected)
                        print('computed =', computed)
                        input('\n')
                    else:
                        print('alpha =', alpha, 'gamma =', gamma)
                        print('expected =', expected)
                        print('.')
                        print('\n')
                    # assert expected == computed

