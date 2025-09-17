# ======================================================================
# ADN — Additional Script: Gap Decomposition Audit (strict ADN, no I/O)
# ----------------------------------------------------------------------
# Purpose: Empirically validate and characterize the lane-based gap split
#          for consecutive primes P2 - P1 = (a2 - a1) + (n2*d2 - n1*d1),
#          under a *canonical* strict-ADN mapping P -> (a,d,n).
#
# Design:
#   - Online-compiler friendly: no input(), no files/CSVs. Pure stdout.
#   - Strict ADN everywhere: a,d odd; n even; gcd(a,d)=gcd(a,n)=gcd(d,n)=1.
#   - Deterministic Miller–Rabin for n < 2^64 using bases {2,3,5,7,11,13,17}.
#   - Two canonicalization modes; choose via CONFIG:
#       (A) 'sband': pick admissible (a,d,n) minimizing |n - P/d|
#                    within small odd A_ALLOWED and D_CANDIDATES.
#       (B) 'backbone_n2': pick the smallest odd d with n=2 and gcd(P-2d,d)=1.
#   - Verifies invariants and reports summaries that address referee critiques:
#       * Evenness of components (Δa, Δ(nd)) and identity Δa + Δ(nd) = gap
#       * Same-lane structure: gap is multiple of 2d; min same-lane gap = 2d
#       * Twin gap (2) cannot be same-lane for d>1 (asserted)
#       * Breakdown by gap size and lane relation (same/cross; d=1 vs d>1)
#
# Tuning:
#   - Adjust PRIME_SEED and PRIME_COUNT below. Default is modest; scale up
#     cautiously in online environments.
# ======================================================================

from typing import List, Tuple, Dict

# ---------------- Shared utilities ----------------

def _gcd(a: int, b: int) -> int:
    import math
    return math.gcd(a, b)

_DEF_MR_BASES = (2,3,5,7,11,13,17)

def is_probable_prime(n: int) -> bool:
    if n < 2:
        return False
    # small trial division first
    for p in _DEF_MR_BASES:
        if n == p:
            return True
        if n % p == 0:
            return n == p
    # write n-1 = d*2^s
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
    for a in _DEF_MR_BASES:
        if not check(a):
            return False
    return True

# Square-free kernel (radical)

def rad(m: int) -> int:
    if m == 0:
        return 0
    r = 1
    x = m
    p = 2
    while p * p <= x:
        if x % p == 0:
            r *= p
            while x % p == 0:
                x //= p
        p = 3 if p == 2 else p + 2
    if x > 1:
        r *= x
    return r

# Strict ADN gate check

def strict_adn_ok(a: int, d: int, n: int) -> bool:
    return (
        a % 2 == 1 and d % 2 == 1 and n % 2 == 0 and
        _gcd(a, d) == 1 and _gcd(a, n) == 1 and _gcd(d, n) == 1
    )

# --------------- Canonicalization modes ---------------

A_ALLOWED_DEFAULT = (1,3,5,7,9)
D_CANDIDATES_DEFAULT = (1,3,5,7,9,11,13,17)

def canonical_map_sband(P: int,
                        A_ALLOWED=A_ALLOWED_DEFAULT,
                        D_CANDS=D_CANDIDATES_DEFAULT,
                        include_11_lane: bool = True) -> Tuple[int,int,int]:
    """Pick (a,d,n) minimizing |n - P/d| within small A_ALLOWED and D_CANDS.
    If none admissible, fall back to backbone_n2.
    """
    best = None
    for d in D_CANDS:
        if d % 2 == 0:
            continue
        for a in A_ALLOWED:
            if not include_11_lane and (a == 1 and d == 1):
                continue
            if _gcd(a, d) != 1:
                continue
            num = P - a
            if num % d != 0:
                continue
            n = num // d
            if n % 2 != 0:
                continue
            if _gcd(a, n) != 1 or _gcd(d, n) != 1:
                continue
            # admissible
            delta = abs(n - (P / d))
            key = (delta, d, a)  # tie-breakers: smaller delta, then d, then a
            if best is None or key < best[0]:
                best = (key, (a, d, n))
    if best is not None:
        return best[1]
    # fallback
    return canonical_map_backbone_n2(P, include_11_lane=include_11_lane)

def canonical_map_backbone_n2(P: int, *, include_11_lane: bool = True,
                              D_MAX: int | None = None) -> Tuple[int,int,int]:
    """Pick the smallest odd d giving an admissible n=2 triple (a=P-2d).
    Existence is guaranteed for odd P; counts match ceil((P-1)/4) at scale.
    """
    if D_MAX is None:
        D_MAX = (P - 1) // 2
    for d in range(1, D_MAX + 1, 2):
        a = P - 2 * d
        if a <= 0:
            break
        if not include_11_lane and (a == 1 and d == 1):
            continue
        if _gcd(a, d) != 1:
            continue
        n = 2
        # gcd(a,n)=1 and gcd(d,n)=1 hold automatically with n=2 and odd a,d
        if strict_adn_ok(a, d, n):
            return a, d, n
    # As a safety net (should not happen), return the trivial d=1 representation if admissible
    a = P - 2
    if _gcd(a, 1) == 1 and strict_adn_ok(a, 1, 2):
        return a, 1, 2
    raise RuntimeError("No admissible canonical mapping found (unexpected)")

# --------------- Prime generation ---------------

def first_prime_at_or_above(n: int) -> int:
    x = n if n % 2 == 1 else n + 1
    if x < 3:
        x = 3
    while not is_probable_prime(x):
        x += 2
    return x

def next_prime_after(n: int) -> int:
    x = n + 2 if n % 2 == 1 else n + 1
    if x % 2 == 0:
        x += 1
    while not is_probable_prime(x):
        x += 2
    return x

def consecutive_primes(start_at: int, count_pairs: int) -> List[int]:
    """Return a list of length count_pairs+1 of consecutive primes starting at >= start_at."""
    p = first_prime_at_or_above(start_at)
    primes = [p]
    while len(primes) < count_pairs + 1:
        p = next_prime_after(p)
        primes.append(p)
    return primes

# --------------- Analysis core ---------------

def analyze_gaps(primes: List[int], mode: str = 'sband', include_11_lane: bool = True,
                 A_ALLOWED=A_ALLOWED_DEFAULT, D_CANDS=D_CANDIDATES_DEFAULT) -> None:
    maps: Dict[int, Tuple[int,int,int]] = {}
    def cmap(P: int) -> Tuple[int,int,int]:
        if P in maps:
            return maps[P]
        if mode == 'sband':
            triple = canonical_map_sband(P, A_ALLOWED=A_ALLOWED, D_CANDS=D_CANDS, include_11_lane=include_11_lane)
        elif mode == 'backbone_n2':
            triple = canonical_map_backbone_n2(P, include_11_lane=include_11_lane)
        else:
            raise ValueError("mode must be 'sband' or 'backbone_n2'")
        maps[P] = triple
        return triple

    # Counters
    total_pairs = 0
    id_fail = 0
    even_fail = 0
    twin_same_d_gt1 = 0
    same_lane_pairs = 0
    same_lane_mult_2d_fail = 0

    gap_hist: Dict[int, int] = {}
    gap_breakdown: Dict[int, Dict[str,int]] = {}  # per gap: {'same_d1':..,'same_dgt1':..,'cross':..}

    samples_same_lane: List[Tuple[int,Tuple[int,int,int],int,Tuple[int,int,int]]] = []
    samples_twin: List[Tuple[int,Tuple[int,int,int],int,Tuple[int,int,int]]] = []

    for i in range(len(primes) - 1):
        P1 = primes[i]
        P2 = primes[i+1]
        total_pairs += 1
        a1, d1, n1 = cmap(P1)
        a2, d2, n2 = cmap(P2)
        gap = P2 - P1
        da = a2 - a1
        dnd = n2 * d2 - n1 * d1
        if da + dnd != gap:
            id_fail += 1
        if (da % 2 != 0) or (dnd % 2 != 0) or (gap % 2 != 0):
            even_fail += 1

        same_lane = (a1 == a2 and d1 == d2)
        if same_lane:
            same_lane_pairs += 1
            if gap % (2 * d1) != 0:
                same_lane_mult_2d_fail += 1
            if len(samples_same_lane) < 8:
                samples_same_lane.append((P1, (a1,d1,n1), gap, (a2,d2,n2)))
        if gap == 2:
            if same_lane and d1 > 1:
                twin_same_d_gt1 += 1
            if len(samples_twin) < 8:
                samples_twin.append((P1, (a1,d1,n1), gap, (a2,d2,n2)))

        gap_hist[gap] = gap_hist.get(gap, 0) + 1
        b = gap_breakdown.setdefault(gap, {'same_d1':0,'same_dgt1':0,'cross':0})
        if same_lane and d1 == 1:
            b['same_d1'] += 1
        elif same_lane and d1 > 1:
            b['same_dgt1'] += 1
        else:
            b['cross'] += 1

    # Reports
    print("ADN Gap Decomposition Audit (strict ADN; mode=", mode, ")\n", sep='')
    print(f"pairs_analyzed: {total_pairs}")
    print(f"identity_failures (da+dnd!=gap): {id_fail}")
    print(f"evenness_failures (da or dnd or gap odd): {even_fail}")
    print(f"same-lane pairs: {same_lane_pairs}")
    print(f"same-lane multiple-of-2d failures: {same_lane_mult_2d_fail}")
    print(f"co-lane twins with d>1 (should be 0): {twin_same_d_gt1}")

    # Small gap table (up to 20)
    print("\n[SMALL_GAP_BREAKDOWN]")
    print("gap,count,same_d1,same_d>1,cross_lane")
    for g in sorted(k for k in gap_hist if k <= 20):
        d = gap_breakdown[g]
        print(f"{g},{gap_hist[g]},{d['same_d1']},{d['same_dgt1']},{d['cross']}")

    # Overall gap histogram (top 10 by frequency)
    print("\n[TOP_GAPS]")
    top = sorted(gap_hist.items(), key=lambda kv: (-kv[1], kv[0]))[:10]
    for g,c in top:
        d = gap_breakdown[g]
        print(f"gap={g:>3}  count={c:>7}  same_d1={d['same_d1']:>7}  same_d>1={d['same_dgt1']:>7}  cross={d['cross']:>7}")

    # Sample lines
    if samples_twin:
        print("\n[SAMPLES: TWIN GAP PAIRS]")
        print("P1,(a1,d1,n1),gap,(a2,d2,n2)")
        for rec in samples_twin:
            print(rec)
    if samples_same_lane:
        print("\n[SAMPLES: SAME-LANE PAIRS]")
        print("P1,(a1,d1,n1),gap,(a2,d2,n2)")
        for rec in samples_same_lane:
            print(rec)

# ------------------- CONFIG & RUN -------------------
if __name__ == "__main__":
    # Prime window configuration
    PRIME_SEED = 10_000_000     # start at/above this integer
    PRIME_PAIRS = 2000          # analyze this many consecutive prime gaps (adjust as needed)

    # Canonicalization configuration
    MODE = 'sband'              # 'sband' or 'backbone_n2'
    INCLUDE_11_LANE = True      # allow (1,1) in mapping (predictor excludes it; analysis may include)
    A_ALLOWED = (1,3,5,7,9)
    D_CANDS   = (1,3,5,7,9,11,13,17)

    primes = consecutive_primes(PRIME_SEED, PRIME_PAIRS)
    analyze_gaps(primes, mode=MODE, include_11_lane=INCLUDE_11_LANE,
                 A_ALLOWED=A_ALLOWED, D_CANDS=D_CANDS)
