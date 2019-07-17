from utils import GQ, GS


def test_dewitt_conjecture(n=3):

    def mu(m, k):
        return k * (m,)

    def nu(m, k):
        ans = [m + k - 1]
        while ans[-1] - 2 >= abs(m - k) + 1:
            ans += [ans[-1] - 2]
        return tuple(ans)

    for m in range(1, 10):
        for k in range(1, 10):
            a = mu(m, k)
            b = nu(m, k)
            print(' n =', n)
            print('mu =', a)
            print('nu =', b)
            print()

            if n < len(a) and n < len(b):
                continue

            gs = GS(a, n)
            gq = GQ(b, n)
            print('GS =', gs)
            print()
            print('GQ =', gq)
            print()
            print()
            assert gs == gq
