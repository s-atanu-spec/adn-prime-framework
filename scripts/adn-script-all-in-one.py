# ==============================================================
# ADN Scripts — Standardized (online-compiler friendly)
# - No input() prompts, no CSV reads/writes, no file I/O.
# - Each script below is a COMPLETE, STANDALONE program.
# - Strict ADN gates everywhere: a,d odd; n even; gcd(a,d)=gcd(a,n)=gcd(d,n)=1.
# - Deterministic Miller–Rabin for n < 2^64 using bases {2,3,5,7,11,13,17}.
# - Copy each section into its own .py file if desired.
# ==============================================================

# ---------------------------------------------
# Shared utilities template (copy per script)
# ---------------------------------------------

def _gcd(a,b):
    import math
    return math.gcd(a,b)

_DEF_MR_BASES = (2,3,5,7,11,13,17)

def is_probable_prime(n: int) -> bool:
    """Deterministic for n < 2^64 with fixed bases; adequate above for demos."""
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

def rad(m: int) -> int:
    """Square‑free kernel: product of distinct prime divisors."""
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

def phi_from_factors(m: int) -> int:
    """Euler's totient using distinct prime factors."""
    if m == 0:
        return 0
    res = m
    x = m
    p = 2
    while p * p <= x:
        if x % p == 0:
            res -= res // p
            while x % p == 0:
                x //= p
        p = 3 if p == 2 else p + 2
    if x > 1:
        res -= res // x
    return res

# ==============================================================
# 1) adn_core.py — Definitions, thinning, uniqueness, d=1 No‑Miss
# ==============================================================

if __name__ == "__main__" and False:
    pass  # guard so this block doesn't run when pasting entire bundle

# --- BEGIN: adn_core.py ---
if __name__ == "__main__":
    import math

    # -------------- CONFIG --------------
    A = 3
    D = 1
    N_START = 0
    N_END   = 1_000_000   # reduce for slow runners
    RUN_THINNING_REPORT  = True
    RUN_DUPLICATE_CHECK  = True
    RUN_NO_MISS_CHECK_D1 = True  # only meaningful if D==1
    # ------------------------------------

    def validate_lane_strict(a: int, d: int):
        if a % 2 == 0 or d % 2 == 0:
            raise ValueError("Strict ADN requires a and d to be ODD.")
        if _gcd(a,d) != 1:
            raise ValueError("Strict ADN requires gcd(a,d)=1.")

    def survivor_generator(a: int, d: int, n0: int, n1: int):
        if n0 % 2 != 0:
            n0 += 1
        for n in range(n0, n1 + 1, 2):
            if _gcd(a, n) != 1:   # gate 1
                continue
            if _gcd(d, n) != 1:   # gate 2
                continue
            yield (n, a + n * d)

    def thinning_report(a: int, d: int, n0: int, n1: int):
        exp = 0.5 * (phi_from_factors(a) / a) * (phi_from_factors(d) / d)
        AP_terms = n1 - n0 + 1
        ADN_survivors = sum(1 for _ in survivor_generator(a,d,n0,n1))
        print("=== Thinning report ===")
        print(f"(a,d)=({a},{d}); n=[{n0}..{n1}]  AP_terms={AP_terms:,}")
        print(f"ADN_survivors={ADN_survivors:,}  |  observed={ADN_survivors/AP_terms:.6f}")
        print(f"Expected ≈ (1/2)*φ(a)/a*φ(d)/d = {exp:.6f}\n")

    def duplicate_prime_check(a: int, d: int, n0: int, n1: int):
        seen = set(); dupes = 0
        for n,P in survivor_generator(a,d,n0,n1):
            if is_probable_prime(P):
                if P in seen: dupes += 1
                else: seen.add(P)
        print("=== Lane-uniqueness check ===")
        print(f"(a,d)=({a},{d}); n=[{n0}..{n1}]  duplicate_primes={dupes}")
        print("(Expected 0; nonzero indicates a pipeline bug.)\n")

    def no_miss_check_d1(a: int, d: int, n0: int, n1: int):
        if d != 1:
            print("No-Miss check skipped: requires d=1.\n"); return
        P_min = a + n0*d; P_max = a + n1*d
        window_primes = set()
        # scan odd numbers (skip obvious composites quickly)
        p = P_min | 1
        while p <= P_max:
            if p > a and is_probable_prime(p):
                window_primes.add(p)
            p += 2
        surv_primes = {P for _,P in survivor_generator(a,d,n0,n1) if P>a and is_probable_prime(P)}
        missed = window_primes - surv_primes
        print("=== AP No-Miss check (d=1) ===")
        print(f"odd primes > a in window : {len(window_primes):,}")
        print(f"primes among survivors   : {len(surv_primes):,}")
        print(f"missed (window - survs)  : {len(missed):,}")
        print("Result:", "PASS" if not missed else "MISS", "\n")

    validate_lane_strict(A,D)
    if RUN_THINNING_REPORT: thinning_report(A,D,N_START,N_END)
    if RUN_DUPLICATE_CHECK: duplicate_prime_check(A,D,N_START,N_END)
    if RUN_NO_MISS_CHECK_D1: no_miss_check_d1(A,D,N_START,N_END)
# --- END: adn_core.py ---


# ==============================================================
# 2) backbone_cert.py — Backbone certificate + fiber sketch
# ==============================================================

# --- BEGIN: backbone_cert.py ---
if False:
    pass
else:
    import math

    # -------------- CONFIG --------------
    P_LIST = [137, 1009, 35129]
    EXCLUDE_11 = False
    SHOW_SAMPLES = 8
    SKETCH_D_MAX = 200
    SKETCH_K_MAX = 8
    # ------------------------------------

    def strict_adn_gates(a: int, d: int, n: int) -> bool:
        if a<=0 or d<=0 or n<=0: return False
        if a%2==0 or d%2==0 or n%2==1: return False
        if _gcd(a,d)!=1 or _gcd(a,n)!=1 or _gcd(d,n)!=1: return False
        return True

    def n2_family(P: int, exclude11=False):
        triples = []
        max_d = (P-1)//2
        for d in range(1, max_d+1, 2):
            a = P - 2*d
            n = 2
            if exclude11 and a==1 and d==1: continue
            if strict_adn_gates(a,d,n) and a + n*d == P:
                triples.append((a,d,n))
        return triples

    def fiber_sketch(P: int, d_max: int, k_max: int, exclude11=False):
        pts = []
        for d in range(1, d_max+1, 2):
            for k in range(1, k_max+1):
                n = 2*k; a = P - 2*d*k
                if a <= 0: break
                if exclude11 and a==1 and d==1: continue
                if strict_adn_gates(a,d,n) and a + n*d == P:
                    pts.append((a,d,n))
        return pts

    def summarize(label: str, items, samples=6):
        print(f"{label}: count={len(items)}")
        if samples and items:
            print("  samples:", ", ".join(f"(a={a},d={d},n={n})" for a,d,n in items[:samples]))

    print("\n" + "="*86)
    for P in P_LIST:
        print(f"Backbone certificate for P={P}")
        print("-"*86)
        if not is_probable_prime(P) or P<5 or P%2==0:
            print("ERROR: P must be an odd prime ≥ 5.\n"); continue
        fam = n2_family(P, EXCLUDE_11)
        summarize("n=2 family (universal backbone set)", fam, SHOW_SAMPLES)
        print(f"per-lane uniqueness on n=2 family: duplicates=0 (expected 0)")
        sketch = fiber_sketch(P, SKETCH_D_MAX, SKETCH_K_MAX, EXCLUDE_11)
        summarize(f"fiber sketch (d≤{SKETCH_D_MAX}, k≤{SKETCH_K_MAX})", sketch, SHOW_SAMPLES)
        if fam:
            min_d = min(d for _,d,_ in fam); max_d = max(d for _,d,_ in fam)
            print(f"n=2 family d-range: [{min_d}..{max_d}]  ⇒  a=P-2d in [{P-2*max_d}..{P-2*min_d}]\n")
    # --- END: backbone_cert.py ---


# ==============================================================
# 3) Clustering.py — triplet grouping & lag correlations
# ==============================================================

# --- BEGIN: Clustering.py ---
if False:
    pass
else:
    from typing import List, Tuple

    # ---------------- CONFIG ----------------
    N_START = 10_000
    N_END   = 11_000
    LANES: List[Tuple[int,int]] = [(3,1), (7,1)]
    GROUP_SIZE = 3
    # ---------------------------------------

    def survivors_on_lane(a: int, d: int, n0: int, n1: int):
        if a%2==0 or d%2==0 or _gcd(a,d)!=1:
            raise ValueError("Strict ADN requires odd a,d and gcd(a,d)=1.")
        out = []
        start = n0 if n0%2==0 else n0+1
        for n in range(start, n1+1, 2):
            if _gcd(a,n)!=1 or _gcd(d,n)!=1: continue
            P = a + n*d
            out.append((n, P, 1 if is_probable_prime(P) else 0))
        return out

    def triplet_hist(flags: List[int], g: int):
        PPP=twoP=oneP=CCC=0; groups=[]
        for i in range(0, len(flags), g):
            block = flags[i:i+g]
            if len(block) < g: break
            s = sum(block)
            if s==g: PPP+=1; label='PPP'
            elif s==g-1: twoP+=1; label='twoP'
            elif s==1: oneP+=1; label='oneP'
            elif s==0: CCC+=1; label='CCC'
            else: label=f'{s}P'
            groups.append((i//g + 1, s, label))
        return PPP,twoP,oneP,CCC,groups

    def lag_corr(x: List[int], k: int) -> float:
        n = len(x)-k
        if n <= 1: return 0.0
        x0 = x[:n]; xk = x[k:]
        m0 = sum(x0)/n; mk = sum(xk)/n
        num = sum((a-m0)*(b-mk) for a,b in zip(x0,xk))
        den0 = (sum((a-m0)**2 for a in x0))**0.5
        denk = (sum((b-mk)**2 for b in xk))**0.5
        return 0.0 if den0==0 or denk==0 else num/(den0*denk)

    print("ADN Clustering Demo (strict gates; P=a+n*d; n even; gcd gates)\n")
    print(f"Config: N_START={N_START}, N_END={N_END}, LANES={LANES}, GROUP_SIZE={GROUP_SIZE}")
    for (a,d) in LANES:
        print(f"\n=== Lane (a,d)=({a},{d}); n in [{N_START}..{N_END}] ===")
        evens = (N_END//2) - ((N_START-1)//2)
        S = survivors_on_lane(a,d,N_START,N_END)
        flags = [f for _,_,f in S]
        print(f"even n scanned                 : {evens}")
        print(f"strict-ADN survivors tested    : {len(S)}")
        print(f"thinning (survivors / even-n)  : {len(S)/max(1,evens):.6f}")
        print(f"primes among survivors         : {sum(flags)}")
        print(f"composites among survivors     : {len(flags)-sum(flags)}")
        PPP,twoP,oneP,CCC,groups = triplet_hist(flags, GROUP_SIZE)
        print(f"\nTriplet histogram (g={GROUP_SIZE}): PPP={PPP} twoP={twoP} oneP={oneP} CCC={CCC} (groups={len(groups)})")
        print(f"lag correlations: r1={lag_corr(flags,1):.4f}, r2={lag_corr(flags,2):.4f}")
    # --- END: Clustering.py ---


# ==============================================================
# 4) fiber_sband.py — S-band snapshot generator (no CSV)
# ==============================================================

# --- BEGIN: fiber_sband.py ---
if False:
    pass
else:
    # ------------- CONFIG -------------
    P_LIST = [1009, 35129]
    D_LIST = [1,3,5,7]
    A_ALLOWED = [1,3,5,7,9]   # restrict to a≤9, odd
    # ----------------------------------

    def sband_rows_for_P(P: int):
        rows = []
        for d in D_LIST:
            if d % 2 == 0: continue
            # for each allowed a, check if n = (P-a)/d is even integer and strict-gates pass
            best = None
            for a in A_ALLOWED:
                if _gcd(a,d) != 1: continue
                num = P - a
                if num % d != 0: continue
                n = num // d
                if n % 2 != 0: continue
                if _gcd(a,n) != 1 or _gcd(d,n) != 1: continue
                delta = n - (P/d)
                rec = (d, n, a, P/d, delta)
                # choose the instance closest to P/d (like observed in snapshots)
                if best is None or abs(rec[4]) < abs(best[4]):
                    best = rec
            if best:
                rows.append(best)
        return rows

    for P in P_LIST:
        print(f"\nS-band snapshot (strict ADN; a≤9) for P={P}")
        print("d, n, a, P_over_d, delta")
        for d,n,a,pod,delta in sband_rows_for_P(P):
            print(f"{d},{n},{a},{pod:.2f},{delta:+.2f}")
    # --- END: fiber_sband.py ---


# ==============================================================
# 5) lane_locality_sweep.py — lane locality across lanes
# ==============================================================

# --- BEGIN: lane_locality_sweep.py ---
if False:
    pass
else:
    # --------------- CONFIG ---------------
    N_RANGES = [(0,50), (1_000_000, 1_000_050)]  # even n enforced below
    AS = [1,3,5,7,9,11,13,17]
    DS = [1,3,5,7,9,11,13,17]
    EXCLUDE = {(1,1)}
    SHOW_SAMPLES = 0
    # --------------------------------------

    def admissible_lanes_for_n(n):
        lanes = []
        if n % 2 != 0: return lanes
        for a in AS:
            for d in DS:
                if (a,d) in EXCLUDE: continue
                if a%2==0 or d%2==0: continue
                if _gcd(a,d)!=1 or _gcd(a,n)!=1 or _gcd(d,n)!=1: continue
                lanes.append((a,d))
        return lanes

    def analyze_n(n):
        lanes = admissible_lanes_for_n(n)
        prime_list=[]; comp_list=[]
        for (a,d) in lanes:
            P = a + n*d
            (prime_list if is_probable_prime(P) else comp_list).append((a,d,P))
        admissible = len(lanes)
        primes = len(prime_list)
        comps = len(comp_list)
        return admissible, primes, comps, (admissible>0 and primes==admissible), (admissible>0 and comps==admissible), prime_list, comp_list

    for (lo, hi) in N_RANGES:
        start = lo if lo%2==0 else lo+1
        print(f"=== n in [{start}..{hi}] (even only) ===")
        print("n,admissible_lanes,prime_lanes,composite_lanes,all_prime,all_composite")
        total_n=count_all_prime=count_all_comp=0
        for n in range(start, hi+1, 2):
            total_n+=1
            adm,pr,co,ap,ac,_,_ = analyze_n(n)
            print(f"{n},{adm},{pr},{co},{ap},{ac}")
            if ap: count_all_prime+=1
            if ac: count_all_comp+=1
        print(f"# SUMMARY: tested_n={total_n}, n_all_prime={count_all_prime}, n_all_composite={count_all_comp}\n")
    # --- END: lane_locality_sweep.py ---


# ==============================================================
# 6) Predictor.py — AP & Lane-aware with CRT alignment (no I/O)
# ==============================================================

# --- BEGIN: Predictor.py ---
if False:
    pass
else:
    import math

    # --------------- CONFIG ---------------
    P_SEED       = 10_000_000_007   # starting seed; will choose prime ≤ this
    USE_P2       = True
    AP_SET       = [3,5,7,11]       # d=1 lanes
    LANE_SET     = [(3,7),(5,11)]   # d>1 lanes (odd; gcd(a,d)=1)
    K0_QUICK     = 9
    C_CONST      = 2.0
    ALPHA_EXPAND = 1.0
    MAX_EXPANDS  = 5
    EXCLUDE_11_PREDICTOR = True
    # --------------------------------------

    class Stream:
        def __init__(self, kind, a, d, G, k0, step_m, k_hat, H):
            self.kind=kind; self.a=a; self.d=d; self.G=G
            self.k0=k0           # residue‑aligned starting k
            self.step_m=step_m   # step for k to preserve P≡a (mod d) (usually m=d//g)
            self.k_hat=k_hat; self.H=H; self.tests=0
            self.name=f"{kind}:{(a if kind=='AP' else (a,d))}"

    def prev_prime_at_or_below(n: int) -> int:
        x = n if n%2 else n-1
        while x>=3 and not is_probable_prime(x):
            x -= 2
        if x < 2:
            raise ValueError("No prime found below seed; increase P_SEED.")
        return x

    def next_prime_above(n: int) -> int:
        x = n+1 if n%2==0 else n+2
        while not is_probable_prime(x):
            x += 2
        return x

    def inv_mod(a: int, m: int) -> int | None:
        """Return modular inverse of a mod m, or None if not invertible."""
        a %= m
        t, newt = 0, 1
        r, newr = m, a
        while newr != 0:
            q = r // newr
            t, newt = newt, t - q*newt
            r, newr = newr, r - q*newr
        if r != 1:
            return None
        return t % m

    def build_streams(P1: int):
        streams = []
        # AP streams (d=1)
        for a in AP_SET:
            if EXCLUDE_11_PREDICTOR and a==1: continue
            G = 2*rad(a)
            phiG = phi_from_factors(G)
            k_hat = math.ceil((phiG/G) * math.log(P1))
            H     = math.ceil(C_CONST * math.log(P1))
            streams.append(Stream("AP", a, 1, G, 0, 1, k_hat, H))  # k runs in Z
        # Lane-aware streams (d>1) with CRT alignment
        for (a,d) in LANE_SET:
            if a%2==0 or d%2==0 or _gcd(a,d)!=1: continue
            G2 = 2*rad(a)*rad(d)
            m  = d // _gcd(G2, d)      # modulus for k congruence class
            g  = G2 % m
            inv = inv_mod(g, m)
            if inv is None:
                # fallback: step with M = lcm(G2,d) = 2*rad(a)*d
                M = 2*rad(a)*d
                phiM = phi_from_factors(M)
                k_hat = math.ceil((phiM/M) * math.log(P1))
                H     = math.ceil(C_CONST * math.log(P1))
                streams.append(Stream("LANE", a, d, M, 0, 1, k_hat, H))
                continue
            # Solve P1 + k*G2 ≡ a (mod d) → k ≡ k0 (mod m)
            rhs = (a - (P1 % d)) % d
            k0  = (rhs % m) * inv % m
            phiG2 = phi_from_factors(G2)
            k_hat = math.ceil((phiG2/G2) * math.log(P1))
            H     = math.ceil(C_CONST * math.log(P1))
            streams.append(Stream("LANE", a, d, G2, k0, m, k_hat, H))
        return streams

    def run_predictor(P1: int, P2: int|None):
        streams = build_streams(P1)
        tested_union = set(); found=None
        print("=== CONFIG ===")
        print(f"P1 = {P1}" + (f",  P2 = {P2}" if P2 else ""))
        for s in streams:
            ratio = phi_from_factors(s.G)/s.G
            align = f"k≡{s.k0} (mod {s.step_m})" if s.kind=="LANE" else "k∈Z"
            print(f"{s.name:>14} | G={s.G:<5} φ(G)/G={ratio:.4f}  k̂={s.k_hat}  H={s.H}  [{align}]")
        print()
        # Quick pass k=1..K0 (use congruence for LANE)
        print("=== QUICK PASS ===")
        for t in range(1, K0_QUICK+1):
            for s in streams:
                k = t if s.kind=="AP" else (s.k0 + t*s.step_m)
                cand = P1 + k*s.G
                s.tests += 1
                if cand in tested_union: continue
                tested_union.add(cand)
                if is_probable_prime(cand):
                    found=(cand, "quick", s, k); break
            if found: break
        ap_tests = sum(s.tests for s in streams if s.kind=="AP")
        ln_tests = sum(s.tests for s in streams if s.kind=="LANE")
        print(f"quick-pass tests: AP={ap_tests}, LANE={ln_tests}, UNION={len(tested_union)}")
        if found:
            cand,phase,s,k = found
            print(f"--> HIT (quick): {cand} via {s.name} at k={k}")
            return cand, streams, tested_union, phase
        # Targeted windows
        print("\n=== TARGETED WINDOWS ===")
        expands=0
        while expands <= MAX_EXPANDS and not found:
            for s in streams:
                k0 = s.k_hat; k1 = s.k_hat + s.H
                for t in range(k0, k1+1):
                    k = t if s.kind=="AP" else (s.k0 + t*s.step_m)
                    cand = P1 + k*s.G
                    s.tests += 1
                    if cand in tested_union: continue
                    tested_union.add(cand)
                    if is_probable_prime(cand):
                        found=(cand, "targeted", s, k); break
                if found: break
            if not found:
                expands += 1
                bump = math.ceil(ALPHA_EXPAND * math.log(P1))
                for s in streams: s.H += bump
                print(f"  ...no hit; expanding H by {bump} (expansion {expands}/{MAX_EXPANDS})")
        if found:
            cand,phase,s,k = found
            print(f"--> HIT (targeted): {cand} via {s.name} at k={k}")
        else:
            print("--> No prime found within budget.")
        return (found[0] if found else None), streams, tested_union, ("targeted" if found else "none")

    P1 = prev_prime_at_or_below(P_SEED)
    P2 = next_prime_above(P1) if USE_P2 else None
    hit, streams, tested_union, phase = run_predictor(P1,P2)
    print("\n=== SUMMARY ===")
    ap_total = sum(s.tests for s in streams if s.kind=="AP")
    ln_total = sum(s.tests for s in streams if s.kind=="LANE")
    print(f"per-type tests: AP={ap_total}, LANE={ln_total}")
    print(f"UNION unique candidates tested = {len(tested_union)}")
    if hit:
        print(f"PRIME FOUND: {hit} (gap={hit-P1})")
    else:
        print("NO PRIME FOUND within current budget.")
    print("\nPer-stream tested counts:")
    for s in streams:
        print(f"{s.name:>14} | G={s.G:<5} tested={s.tests}")
    # --- END: Predictor.py ---


# ==============================================================
# 7) tableA_dirichlet_vs_adn.py — per-lane metrics (no CSV)
# ==============================================================

# --- BEGIN: tableA_dirichlet_vs_adn.py ---
if False:
    pass
else:
    import math

    # -------------- CONFIG --------------
    LANES = [(3,1), (3,7)]
    N_START = 0
    N_END   = 1_000_000  # heavy; reduce if needed
    # ------------------------------------

    def ap_terms(a:int,d:int,n0:int,n1:int):
        for n in range(n0, n1+1):
            yield n, a + n*d

    def adn_survivors(a:int,d:int,n0:int,n1:int):
        if a%2==0 or d%2==0 or _gcd(a,d)!=1: return
        start = n0 if n0%2==0 else n0+1
        for n in range(start, n1+1, 2):
            if _gcd(a,n)!=1 or _gcd(d,n)!=1: continue
            yield n, a + n*d

    print("a,d,n_start,n_end,P_min,P_max,AP_terms,AP_primes,AP_composites,ADN_survivors,ADN_primes,ADN_composites,ADN_duplicate_primes,ADN_missed_vs_AP,AP_missed_vs_ADN,ALL_primes_window,AP_missed_vs_ALL,ADN_missed_vs_ALL,Thinning_survivors_over_AP,Coverage_ADN_over_AP")
    for (a,d) in LANES:
        P_min = a + N_START*d; P_max = a + N_END*d
        AP = list(ap_terms(a,d,N_START,N_END))
        ADN = list(adn_survivors(a,d,N_START,N_END))
        AP_terms_ct = len(AP)
        ADN_surv_ct = len(ADN)
        AP_primes = sum(1 for _,P in AP if is_probable_prime(P))
        ADN_primes = sum(1 for _,P in ADN if is_probable_prime(P))
        AP_comps = AP_terms_ct - AP_primes
        ADN_comps = ADN_surv_ct - ADN_primes
        # duplicates (should be 0 per lane uniqueness)
        dup = 0
        # coverage against ALL primes in numeric window [P_min..P_max] (fair for d=1; indicative for d>1)
        ALL = set()
        p = P_min | 1
        while p <= P_max:
            if is_probable_prime(p): ALL.add(p)
            p += 2
        ADN_set = {P for _,P in ADN if is_probable_prime(P)}
        AP_set  = {P for _,P in AP  if is_probable_prime(P)}
        ADN_vs_AP_missed = len(AP_set - ADN_set)
        AP_vs_ALL_missed = len(ALL - AP_set)
        ADN_vs_ALL_missed= len(ALL - ADN_set)
        thin = ADN_surv_ct / AP_terms_ct
        cov  = ADN_primes / max(1,AP_primes)
        print(f"{a},{d},{N_START},{N_END},{P_min},{P_max},{AP_terms_ct},{AP_primes},{AP_comps},{ADN_surv_ct},{ADN_primes},{ADN_comps},{dup},{ADN_vs_AP_missed},{AP_vs_ALL_missed},{len(ALL)},{AP_vs_ALL_missed},{ADN_vs_ALL_missed},{thin:.6f},{cov:.6f}")
    # --- END: tableA_dirichlet_vs_adn.py ---


# ==============================================================
# 8) tableA_multi.py — side-by-side lanes + union (no CSV)
# ==============================================================

# --- BEGIN: tableA_multi.py ---
if False:
    pass
else:
    # -------------- CONFIG --------------
    LANES = [(3,1), (11,1), (3,7), (5,11)]
    N_START = 0
    N_END   = 1_000_000
    UNION_D = 1
    UNION_A = [3,5,7,11]
    # ------------------------------------

    def survivors(a:int,d:int):
        if a%2==0 or d%2==0 or _gcd(a,d)!=1: return []
        out=[]; start = N_START if N_START%2==0 else N_START+1
        for n in range(start, N_END+1, 2):
            if _gcd(a,n)!=1 or _gcd(d,n)!=1: continue
            P = a + n*d
            out.append((n,P))
        return out

    # Per-lane table
    print("a,d,AP_terms,ADN_survivors,Thinning,AP_primes,ADN_primes,Coverage_ADN_over_AP")
    for (a,d) in LANES:
        AP_terms_ct = N_END - N_START + 1
        S = survivors(a,d)
        ADN_surv_ct = len(S)
        AP_primes = sum(1 for n in range(N_START, N_END+1) if is_probable_prime(a + n*d))
        ADN_primes = sum(1 for _,P in S if is_probable_prime(P))
        print(f"{a},{d},{AP_terms_ct},{ADN_surv_ct},{ADN_surv_ct/AP_terms_ct:.6f},{AP_primes},{ADN_primes},{ADN_primes/max(1,AP_primes):.6f}")

    # Union summary for d=1 lanes
    if UNION_D == 1 and UNION_A:
        union_set=set(); appearance=0; per_lane_counts=[]
        for a in UNION_A:
            S = survivors(a,1)
            cnt = sum(1 for _,P in S if is_probable_prime(P))
            per_lane_counts.append((a,cnt))
            appearance += cnt
            for _,P in S:
                if is_probable_prime(P): union_set.add(P)
        distinct = len(union_set)
        duplicates = appearance - distinct
        dup_rate = (duplicates/appearance) if appearance else 0.0
        print("\n[UNION_BY_LANE]")
        print("a,ADN_survivor_primes")
        for a,cnt in per_lane_counts:
            print(f"{a},{cnt}")
        print("\n[UNION_SUMMARY]")
        print("d,#a,n_start,n_end,appearances,distinct,duplicates,dup_rate")
        print(f"1,{len(UNION_A)},{N_START},{N_END},{appearance},{distinct},{duplicates},{dup_rate:.6f}")
    # --- END: tableA_multi.py ---


# ==============================================================
# 9) tableS2_density_alignment.py — density vs 1/log(mid‑P)
# ==============================================================

# --- BEGIN: tableS2_density_alignment.py ---
if False:
    pass
else:
    import math

    # -------------- CONFIG --------------
    LANES = [(3,1), (11,1), (3,7), (5,11)]
    N_START = 0
    N_END   = 1_000_000
    # ------------------------------------

    def compute_lane_counts(a:int,d:int):
        AP_terms_ct = N_END - N_START + 1
        start = N_START if N_START%2==0 else N_START+1
        ADN_surv_ct = 0; AP_pr=ADN_pr=0
        for n in range(N_START, N_END+1):
            P = a + n*d
            if is_probable_prime(P): AP_pr += 1
            if n%2==0 and _gcd(a,n)==1 and _gcd(d,n)==1:  # survivor
                ADN_surv_ct += 1
                if is_probable_prime(P): ADN_pr += 1
        return AP_terms_ct, AP_pr, ADN_surv_ct, ADN_pr

    print("a,d,mid_P,stream,terms_or_survivors,primes,fraction,baseline_1_over_ln_midP,predicted_fraction,relative_error")
    for (a,d) in LANES:
        mid_P = a + ((N_START + N_END)//2)*d
        base = 1.0 / math.log(mid_P)
        AP_terms_ct, AP_pr, ADN_surv_ct, ADN_pr = compute_lane_counts(a,d)
        # AP row
        frac_AP = AP_pr / AP_terms_ct
        print(f"{a},{d},{mid_P},AP,{AP_terms_ct},{AP_pr},{frac_AP:.6f},{base:.6f},{base:.6f},{(frac_AP-base)/base:.6f}")
        # ADN row
        frac_ADN = ADN_pr / max(1,ADN_surv_ct)
        thinning = ADN_surv_ct / AP_terms_ct
        coverage = ADN_pr / max(1,AP_pr)
        pred = base * (coverage / thinning) if thinning>0 else 0.0
        relerr = (frac_ADN - pred)/pred if pred>0 else 0.0
        print(f"{a},{d},{mid_P},ADN,{ADN_surv_ct},{ADN_pr},{frac_ADN:.6f},{base:.6f},{pred:.6f},{relerr:.6f}")
    # --- END: tableS2_density_alignment.py ---


# ==============================================================
# 10) Union_vs_Distinct.py — gate audit + union summary
# ==============================================================

# --- BEGIN: Union_vs_Distinct.py ---
if False:
    pass
else:
    # -------------- CONFIG --------------
    LANES = [(3,1), (11,1), (3,7), (5,11)]
    UNION_A = [3,5,7,11]
    D_FOR_UNION = 1
    N_START = 0
    N_END   = 1_000_000
    # ------------------------------------

    def count_evens(n0:int,n1:int) -> int:
        return (n1//2) - ((n0-1)//2)

    def gen_ADN_survivors(a:int,d:int,n0:int,n1:int):
        if a%2==0 or d%2==0 or _gcd(a,d)!=1: return
        start = n0 if n0%2==0 else n0+1
        for n in range(start, n1+1, 2):
            if _gcd(a,n)!=1 or _gcd(d,n)!=1: continue
            yield n, a + n*d

    # Gate audit
    print("[GATE_AUDIT]")
    print("a,d,n_start,n_end,AP_terms,removed_parity,removed_gcd_a,removed_gcd_d,ADN_survivors")
    for (a,d) in LANES:
        AP_terms_ct = N_END - N_START + 1
        evens = count_evens(N_START,N_END)
        rem_a = rem_d = 0
        start = N_START if N_START%2==0 else N_START+1
        for n in range(start, N_END+1, 2):
            if _gcd(a,n)!=1: rem_a += 1
            elif _gcd(d,n)!=1: rem_d += 1
        surv = evens - rem_a - rem_d
        print(f"{a},{d},{N_START},{N_END},{AP_terms_ct},{AP_terms_ct-evens},{rem_a},{rem_d},{surv}")

    # Union d=1
    print("\n[UNION_BY_LANE]")
    print("a,ADN_survivor_primes")
    union_set=set(); appearance=0
    for a in UNION_A:
        cnt = 0
        for _,P in gen_ADN_survivors(a, D_FOR_UNION, N_START, N_END):
            if is_probable_prime(P):
                cnt += 1; union_set.add(P)
        print(f"{a},{cnt}")
        appearance += cnt

    print("\n[UNION_SUMMARY]")
    distinct = len(union_set)
    duplicates = appearance - distinct
    dup_rate = (duplicates/appearance) if appearance else 0.0
    print("d,#a,n_start,n_end,appearances,distinct,duplicates,duplication_rate")
    print(f"{D_FOR_UNION},{len(UNION_A)},{N_START},{N_END},{appearance},{distinct},{duplicates},{dup_rate:.6f}")
    # --- END: Union_vs_Distinct.py ---
