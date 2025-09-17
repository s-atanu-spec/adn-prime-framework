# prime_sums_parity_check.py
# Purpose: Empirically verify parity laws for sums of consecutive primes:
#   1) If P1 = P + 2m, then P + P1 = 2(P + m) is always even.
#   2) For twins (m=1), P + P1 is always divisible by 4.
# Outputs: audit-style summary + a few sample lines (referee friendly)

import argparse
import sys

# Deterministic Miller–Rabin good for very large ints with these bases in practice.
# (For our audit it's more than enough; if you test 64-bit, it's provably correct.)
MR_BASES = (2, 3, 5, 7, 11, 13, 17)

def is_probable_prime(n: int) -> bool:
    if n < 2: return False
    # small prime quick path
    small_primes = (2, 3, 5, 7, 11, 13, 17)
    for p in small_primes:
        if n == p: return True
        if n % p == 0: return n == p
    # write n-1 = d * 2^s
    d = n - 1
    s = 0
    while d % 2 == 0:
        d //= 2
        s += 1
    def check(a):
        x = pow(a, d, n)
        if x == 1 or x == n - 1:
            return True
        for _ in range(s - 1):
            x = (x * x) % n
            if x == n - 1:
                return True
        return False
    return all(check(a) for a in MR_BASES)

def next_prime_at_or_above(n: int) -> int:
    # start at an odd >=3
    x = n if n % 2 == 1 else n + 1
    if x < 3:
        x = 3
    while not is_probable_prime(x):
        x += 2
    return x

def next_prime_after(n: int) -> int:
    x = n + 2 if n % 2 == 1 else n + 1
    while not is_probable_prime(x):
        x += 2
    return x

def run(seed: int, pairs: int, show: int):
    p = next_prime_at_or_above(seed)
    results = []
    counts = {
        "pairs": 0,
        "even_gap_pairs": 0,   # sanity (should equal pairs)
        "sum_even": 0,         # sanity (should equal pairs)
        "sum_mod4_0": 0,
        "sum_mod4_2": 0,
        "twins_total": 0,
        "twins_sum_mod4_0": 0
    }

    # collect pairs
    for _ in range(pairs):
        q = next_prime_after(p)
        gap = q - p
        m = gap // 2  # half-gap
        s = p + q
        sum_mod4 = s % 4

        # tallies
        counts["pairs"] += 1
        if gap % 2 == 0:
            counts["even_gap_pairs"] += 1
        if s % 2 == 0:
            counts["sum_even"] += 1
        if sum_mod4 == 0:
            counts["sum_mod4_0"] += 1
        elif sum_mod4 == 2:
            counts["sum_mod4_2"] += 1
        # twins
        if gap == 2:
            counts["twins_total"] += 1
            if sum_mod4 == 0:
                counts["twins_sum_mod4_0"] += 1

        # keep a few samples for print
        if len(results) < show:
            results.append((p, q, gap, m, sum_mod4))
        p = q

    # Print audit summary
    print("Prime Sums & Parity Audit (strict ADN parity lens)")
    print(f"seed={seed}; pairs={pairs}; samples_shown={len(results)}")
    print(f"pairs_total={counts['pairs']}")
    print(f"even_gap_pairs={counts['even_gap_pairs']}   (expected = pairs_total)")
    print(f"sum_even={counts['sum_even']}               (expected = pairs_total)")
    print(f"sum_mod4: 0->{counts['sum_mod4_0']} , 2->{counts['sum_mod4_2']}")
    print(f"twins_total={counts['twins_total']}")
    print(f"twins_sum_mod4_0={counts['twins_sum_mod4_0']}   (expected = twins_total)")
    print("\nSamples (P, P1, gap, m=gap/2, (P+P1) mod 4):")
    for (P, Q, gap, m, mod4) in results:
        print(f"P={P}, P1={Q}, gap={gap}, m={m}, sum_mod4={mod4}")

if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="Audit parity of sums for consecutive primes.")
    ap.add_argument("--seed", type=int, default=10_000, help="Start search at or above this integer (default 10000)")
    ap.add_argument("--pairs", type=int, default=1000, help="Number of consecutive prime pairs to test (default 1000)")
    ap.add_argument("--show", type=int, default=12, help="How many sample rows to print (default 12)")
    # Add the following line to parse known args and ignore others
    args, unknown = ap.parse_known_args()
    run(args.seed, args.pairs, args.show)