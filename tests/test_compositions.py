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

        c = list(alpha)
        c[-1] += gamma[0] - 1
        c += list(gamma[1:])

        g = list(gamma)
        g[0] -= 1

        d = list(alpha)
        d[-1] += 1
        d += g

        e = list(alpha)
        e += g

        if is_peak_composition(a):
            yield tuple(a), 1
        if is_peak_composition(b):
            yield tuple(b), 1
        if is_peak_composition(c):
            yield tuple(c), beta
        if is_peak_composition(d):
            yield tuple(d), 1
        if is_peak_composition(e):
            yield tuple(e), beta


def bar_shribbon_multiplier(alpha, gamma):
    if sum(alpha) == 0:
        yield gamma, 1
    elif sum(gamma) == 0:
        yield alpha, 1
    else:
        a = alpha + gamma
            
        b = list(alpha)
        b[-1] += gamma[0] 
        b += list(gamma[1:])

        c = list(alpha)
        c[-1] += gamma[0] - 1
        c += list(gamma[1:])

        g = list(gamma)
        g[0] -= 1

        d = list(alpha)
        d[-1] += 1
        d += g

        e = list(alpha)
        e += g

        f = list(alpha)
        f[-1] += gamma[0] - 2
        f += list(gamma[1:])

        p = int(alpha[-1] > 1)
        q = int(gamma[0] > 2 or gamma == (2,))

        if is_peak_composition(a):
            yield tuple(a), 1
        if is_peak_composition(b):
            yield tuple(b), 2
        if is_peak_composition(c):
            yield tuple(c), (1 + p + q) * beta
        if is_peak_composition(d):
            yield tuple(d), 1
        if is_peak_composition(e):
            yield tuple(e), beta
        if is_peak_composition(f):
            yield tuple(f), p * q * beta**2


def Ribbon(alpha, multiplier=ribbon_multiplier):
    return Vector({alpha: 1}, multiplier=multiplier)


def barShRibbon(alpha):
    ans = 0
    n = len(alpha)
    for v in range(2**n):
        delta = list(alpha)
        for i in range(n):
            delta[i] -= v % 2
            v = v // 2
        if any(d < 2 for d in delta[:-1]):
            continue
        ans += ShRibbon(tuple(delta)) * 2**(n + sum(delta) - sum(alpha)) * beta**(sum(alpha) - sum(delta))
    return ans


def ShRibbon(alpha):
    ans = 0
    n = sum(alpha)
    for gamma in Partition.compositions(n):
        if Lambda(gamma) == alpha:
            ans += Ribbon(gamma)
    return ans



def barShRibbon_expansion(f):
    ans = 0
    while f != 0:
        alpha = max(f)
        two = 2**len(alpha)
        coeff = f[alpha]
        assert coeff % two == 0
        ans += Ribbon(alpha) * (coeff // two)
        f -= barShRibbon(alpha) * (coeff // two)
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


def test_bar_shribbon_product(n=5):
    for p in range(1, n + 1):
        for q in range(1, n + 1):
            for alpha in Partition.peak_compositions(p):
                for gamma in Partition.peak_compositions(q):
                    expected = Ribbon(alpha, bar_shribbon_multiplier) * Ribbon(gamma, bar_shribbon_multiplier)
                    computed = barShRibbon_expansion(barShRibbon(alpha) * barShRibbon(gamma))
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
