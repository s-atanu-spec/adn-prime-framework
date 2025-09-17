# ======================================================================
# ADN — Benchmark Script (odd-step, wheel-30, segmented sieve, ADN predictors)
# Fixed: AP-slice & lane-stride starts now avoid all-even / all-composite traps.
# Online-compiler friendly: no input(), no file I/O; prints as it goes.
# Deterministic Miller–Rabin bases {2,3,5,7,11,13,17} for n < 2^64.
# Metrics per method: unique candidates tested; first-hit index (hit - P1); wall-clock time.
# ======================================================================

from math import gcd, isqrt
from time import perf_counter

# --------------------------- Config ---------------------------
P1 = 100003              # anchor (search strictly > P1)
WINDOW_W = 5000          # segmented sieve lookahead window (shrink/grow as needed)

AP_A_VALUES = (3, 5, 7, 11)        # AP-slice lanes (d=1)
LANE_STRIDE = ((3, 7), (5, 11))    # lane-stride (a,d) choices (d>1)

WHEEL_MOD = 30
WHEEL_RES = (1, 7, 11, 13, 17, 19, 23, 29)

# Time/robustness controls
ENABLE_SIEVE = True
TIMEOUT_BASE_SEC   = 3.0    # per-method timeout for odd-step & wheel
TIMEOUT_SIEVE_SEC  = 5.0    # timeout for segmented sieve
TIMEOUT_PRED_SEC   = 5.0    # timeout for each predictor
BLOCK_SIZE = 4096           # sieve block size (smaller = more responsive)

# ------------------ Number theory helpers --------------------
DEF_MR_BASES = (2, 3, 5, 7, 11, 13, 17)

def is_probable_prime(n: int) -> bool:
    if n < 2:
        return False
    for p in DEF_MR_BASES:
        if n == p:
            return True
        if n % p == 0:
            return n == p
    # n-1 = d*2^s
    d = n - 1
    s = 0
    while d % 2 == 0:
        d //= 2
        s += 1
    def check(a: int) -> bool:
        x = pow(a, d, n)
        if x == 1 or x == n - 1:
            return True
        for _ in range(s - 1):
            x = (x * x) % n
            if x == n - 1:
                return True
        return False
    for a in DEF_MR_BASES:
        if not check(a):
            return False
    return True

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

def next_odd(n: int) -> int:
    return n if n % 2 == 1 else n + 1

# ------------------ Baseline methods -------------------------
def method_odd_step(p1: int):
    t0 = perf_counter()
    tested = 0
    n = next_odd(p1 + 1)
    while True:
        if perf_counter() - t0 > TIMEOUT_BASE_SEC:
            return {"method":"odd-step","hit":None,"first_hit_index":None,
                    "tested":tested,"time":perf_counter()-t0,"timeout":True}
        tested += 1
        if is_probable_prime(n):
            return {"method":"odd-step","hit":n,"first_hit_index":n - p1,
                    "tested":tested,"time":perf_counter()-t0}
        n += 2

def method_wheel30(p1: int):
    t0 = perf_counter()
    tested = 0
    n = next_odd(p1 + 1)
    while n % WHEEL_MOD not in WHEEL_RES:
        if perf_counter() - t0 > TIMEOUT_BASE_SEC:
            return {"method":"wheel-30","hit":None,"first_hit_index":None,
                    "tested":tested,"time":perf_counter()-t0,"timeout":True}
        n += 2
    while True:
        if perf_counter() - t0 > TIMEOUT_BASE_SEC:
            return {"method":"wheel-30","hit":None,"first_hit_index":None,
                    "tested":tested,"time":perf_counter()-t0,"timeout":True}
        tested += 1
        if is_probable_prime(n):
            return {"method":"wheel-30","hit":n,"first_hit_index":n - p1,
                    "tested":tested,"time":perf_counter()-t0}
        # hop to next allowed residue by +2
        while True:
            n += 2
            if n % WHEEL_MOD in WHEEL_RES:
                break

def segmented_sieve_first(p1: int, W: int):
    hi = p1 + W
    R = isqrt(hi)

    # base primes up to sqrt(hi)
    base = [True]*(R+1)
    base[0:2] = [False, False]
    for i in range(2, isqrt(R)+1):
        if base[i]:
            start = i*i
            base[start:R+1:i] = [False]*(((R - start)//i)+1)
    bases = [i for i,pr in enumerate(base) if pr]

    counted_odds = 0
    t0 = perf_counter()
    for lo in range(p1 + 1, hi + 1, BLOCK_SIZE):
        if perf_counter() - t0 > TIMEOUT_SIEVE_SEC:
            return {"method":"segmented-sieve","hit":None,"first_hit_index":None,
                    "tested":counted_odds,"time":perf_counter()-t0,"timeout":True}
        hi2 = min(lo + BLOCK_SIZE - 1, hi)
        seg = [True]*(hi2 - lo + 1)
        # mark composites in [lo..hi2]
        for p in bases:
            start = max(p*p, ((lo + p - 1)//p)*p)
            for x in range(start, hi2 + 1, p):
                seg[x - lo] = False
        # scan for first prime > p1
        s = max(lo, 2)
        for x in range(s, hi2 + 1):
            if x % 2 == 1:  # count odd candidates considered
                counted_odds += 1
            if seg[x - lo] and x > p1:
                return {"method":"segmented-sieve","hit":x,"first_hit_index":x - p1,
                        "tested":counted_odds,"time":perf_counter()-t0}
    # no hit in window
    return {"method":"segmented-sieve","hit":None,"first_hit_index":None,
            "tested":counted_odds,"time":perf_counter()-t0}

# ------------------ ADN predictors (fixed) -------------------
def method_ap_slice(p1: int, a: int):
    """
    G = 2*rad(a). We must ensure candidates are odd and not locked to a class sharing a factor with G.
    For a prime/odd anchor p1, using p1 + k*G keeps parity odd. Start at p1 + G (>p1).
    """
    G = 2*rad(a)
    t0 = perf_counter()
    tested = 0

    # Prefer anchor-based residue: start at p1 + G (keeps parity of p1; if p1 is odd -> odd).
    n = p1 + G
    # If somehow even (e.g., even p1), shift to an odd residue by adding G//2 once.
    if n % 2 == 0 and (G % 2 == 0):
        n += (G // 2)

    # Safety: time out if something goes wrong
    while True:
        if perf_counter() - t0 > TIMEOUT_PRED_SEC:
            return {"method":f"ADN AP-slice (a={a})","G":G,"hit":None,"first_hit_index":None,
                    "tested":tested,"time":perf_counter()-t0,"timeout":True}
        tested += 1
        if is_probable_prime(n):
            return {"method":f"ADN AP-slice (a={a})","G":G,"hit":n,
                    "first_hit_index":n - p1,"tested":tested,"time":perf_counter()-t0}
        n += G  # even step; parity preserved

def method_lane_stride(p1: int, a: int, d: int):
    """
    Step M = 2*rad(a)*d. Start P0 is the smallest >= p1+1 with P0 ≡ a (mod d),
    then adjust by +d until P0 is odd and gcd(P0, rad(a)) == 1, so gcd(P0, M) == 1.
    This guarantees the arithmetic progression can contain primes (Dirichlet condition).
    """
    G = 2*rad(a)
    M = G * d
    t0 = perf_counter()
    tested = 0

    # base congruence P0 ≡ a (mod d)
    P0 = (p1 + 1) + ((a - (p1 + 1)) % d)
    # adjust to avoid even and shared rad(a) factors; keep congruence by stepping +d
    ra = rad(a)
    while (P0 % 2 == 0) or (gcd(P0, ra) != 1):
        P0 += d

    # safety: if after a few steps we're still not coprime to M, give up
    guard = 0
    while gcd(P0, M) != 1 and guard < 10:
        P0 += d
        guard += 1

    n = P0
    while True:
        if perf_counter() - t0 > TIMEOUT_PRED_SEC:
            return {"method":f"ADN lane-stride (a={a},d={d})","M":M,"hit":None,
                    "first_hit_index":None,"tested":tested,"time":perf_counter()-t0,"timeout":True}
        tested += 1
        if is_probable_prime(n):
            return {"method":f"ADN lane-stride (a={a},d={d})","M":M,"hit":n,
                    "first_hit_index":n - p1,"tested":tested,"time":perf_counter()-t0}
        n += M

# --------------------------- Runner ---------------------------
def print_header():
    print(f"Benchmark window: (P1, P1+W] where P1={P1}, W={WINDOW_W}")
    print("\n| method                     | step/sieve            | unique candidates | first-hit index | wall-clock (s) | notes")
    print("|---------------------------|----------------------|-------------------|-----------------|----------------|------")

def print_row(m):
    meth = m.get("method")
    step = ("+2" if meth=="odd-step" else
            ("wheel-30 residues" if meth=="wheel-30" else
             ("segmented sieve" if meth=="segmented-sieve" else
              (f"G={m.get('G')}" if 'G' in m else (f"M={m.get('M')}" if 'M' in m else "-")))))
    tested = m.get("tested")
    fh = m.get("first_hit_index")
    t = m.get("time")
    note = "timeout" if m.get("timeout") else ""
    print(f"| {meth:<25} | {step:<20} | {tested:>17} | {str(fh):>15} | {t:.6f}       | {note}", flush=True)

if __name__ == "__main__":
    print_header()

    r = method_odd_step(P1);    print_row(r)
    r = method_wheel30(P1);     print_row(r)

    if ENABLE_SIEVE:
        r = segmented_sieve_first(P1, WINDOW_W); print_row(r)

    for a in AP_A_VALUES:
        r = method_ap_slice(P1, a); print_row(r)

    for (a, d) in LANE_STRIDE:
        r = method_lane_stride(P1, a, d); print_row(r)
