"""Check the spacing representation and evaluate the tied-rank approximation."""
from math import comb, exp, fsum

import numpy as np

from rank_sweep import neighboring_keys


def check_spacing_tests():
    rng = np.random.default_rng(7393620)
    checked = 0
    for n in (32, 128, 512):
        for _ in range(100):
            spacings = rng.exponential(size=n)
            arc = np.r_[0, np.cumsum(spacings[:-1])] / spacings.sum()
            ranks, _, _ = neighboring_keys(arc)
            left = np.zeros(n)
            right = np.zeros(n)
            for m in range(1, 13):
                left += np.roll(spacings, m)
                right += np.roll(spacings, -m)
                np.testing.assert_array_equal(
                    (spacings > left) & (spacings > right), ranks['OR_kNN'] > m)
                np.testing.assert_array_equal(
                    (spacings > left) | (spacings > right), ranks['mutual_kNN'] > m)
                checked += 2 * n
    print(f'{checked} spacing-threshold comparisons agreed.')


def symmetric_tail(m):
    return fsum(comb(2*m + u + v, m + u) / 3**(2*m + u + v + 1)
                for u in range(80) for v in range(80))


def tie_approximations(sample_sizes):
    # The omitted rank levels contribute at most n*pi_65.
    symmetric = [symmetric_tail(m) for m in range(66)]
    mutual = [2.**(1-m) - tail for m, tail in enumerate(symmetric)]
    return {
        rule: tuple(1 - fsum(n * (a - b) * exp(-n * a)
                             for a, b in zip(tails, tails[1:]))
                    for n in sample_sizes)
        for rule, tails in [('symmetric', symmetric), ('mutual', mutual)]
    }


if __name__ == '__main__':
    check_spacing_tests()
    sizes = (128, 512, 4096, 32768)
    for rule, rates in tie_approximations(sizes).items():
        print(rule, ', '.join(f'n={n}: {rate:.1%}' for n, rate in zip(sizes, rates)))
