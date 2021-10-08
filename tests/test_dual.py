from symmetric import SymmetricPolynomial
from partitions import Partition


def test_simple():
    f = SymmetricPolynomial._slow_transposed_dual_stable_grothendieck_q(4, (4,))
    g = SymmetricPolynomial.dual_stable_grothendieck_q(4, (4,))
    assert f.omega_schur_expansion(f) == g.schur_expansion(g)

    f = SymmetricPolynomial._slow_transposed_dual_stable_grothendieck_p(4, (4,))
    g = SymmetricPolynomial.dual_stable_grothendieck_p(4, (4,))
    assert f.omega_schur_expansion(f) == g.schur_expansion(g)


def test_p(n=4):
    lam = tuple([n - i for i in range(n)])
    for mu in Partition.subpartitions(lam, strict=True):
        for nu in Partition.subpartitions(mu, strict=True):
            f = SymmetricPolynomial._slow_transposed_dual_stable_grothendieck_p(n, mu, nu)
            g = SymmetricPolynomial.dual_stable_grothendieck_p(n, mu, nu)
            ef, eg = f.omega_schur_expansion(f), g.schur_expansion(g)
            print('n =', n, 'mu =', mu)
            print(ef)
            print(eg)
            print()
            assert ef == eg


def test_q(n=4):
    lam = tuple([n - i for i in range(n)])
    for mu in Partition.subpartitions(lam, strict=True):
        for nu in Partition.subpartitions(mu, strict=True):
            f = SymmetricPolynomial._slow_transposed_dual_stable_grothendieck_q(n, mu, nu)
            g = SymmetricPolynomial.dual_stable_grothendieck_q(n, mu, nu)
            ef, eg = f.omega_schur_expansion(f), g.schur_expansion(g)
            print('n =', n, 'mu =', mu)
            print(ef)
            print(eg)
            print()
            assert ef == eg
