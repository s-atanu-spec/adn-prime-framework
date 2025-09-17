#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Strict-ADN Predictor (AP-slice + Lane-aware with congruence alignment)
- G(a)   = 2*rad(a)                  for d=1
- G(a,d) = 2*rad(a)*rad(d)           for d>1
- Enforce P ≡ a (mod d) by aligned k: k ≡ k0 (mod d/g), g = gcd(G,d) (=rad(d))
  <=> step in P by M = lcm(G,d) = G*(d/g) from an aligned base.
- Hybrid k-schedule: quick pass (k/t small) then targeted window around k̂ ~ (φ(G)/G)*log P1.

Deterministic Miller–Rabin is valid for n < 2^64 (covers 1e10+ easily).
Produces two CSVs: per-stream parameters and overall benchmark summary.

Usage (defaults are safe):
  python predictor_hybrid_strict_adn.py \
    --seed 170141183460469231731687303715884105727 \
    --ap 3 5 7 11 \
    --lane 3,7 5,11 \
    --K0 9 --c 2.0 --alpha 1.0 --max-expands 5 \
    --out-prefix content/sample_data/predictor_1e10

Or pin an exact anchor prime:
  python predictor_hybrid_strict_adn.py --P1 170141183460469231731687303715884105727
"""
import argparse, math, time, csv, os, sys
from typing import List, Tuple, Optional, Dict

# ---------------- arithmetic utilities ----------------
def is_probable_prime(n: int) -> bool:
    if n < 2: return False
    smalls = (2,3,5,7,11,13,17,19,23,29)
    for p in smalls:
        if n == p: return True
        if n % p == 0: return n == p
    # write n-1 = d * 2^s
    d = n - 1
    s = (d & -d).bit_length() - 1
    d >>= s
    def check(a):
        x = pow(a, d, n)
        if x == 1 or x == n - 1: return True
        for _ in range(s - 1):
            x = (x * x) % n
            if x == n - 1: return True
        return False
    for a in (2,3,5,7,11,13,17):
        if not check(a): return False
    return True

def prev_prime_at_or_below(n: int) -> int:
    x = n if n % 2 else n-1
    while x >= 3 and not is_probable_prime(x):
        x -= 2
    if x < 2:
        raise ValueError("No prime found below the seed.")
    return x

def next_prime_above(n: int) -> int:
    x = n + 1 if (n % 2 == 0) else n + 2
    while not is_probable_prime(x):
        x += 2
    return x

def rad(n: int) -> int:
    """Square-free kernel"""
    r, m = 1, n
    p = 2
    while p*p <= m:
        if m % p == 0:
            r *= p
            while m % p == 0:
                m //= p
        p = 3 if p == 2 else p + 2
    if m > 1: r *= m
    return r

def phi_from_factor(n: int) -> int:
    """φ(n) via its distinct prime factors (n small here)"""
    res, m, p = n, n, 2
    while p*p <= m:
        if m % p == 0:
            res -= res // p
            while m % p == 0:
                m //= p
        p = 3 if p == 2 else p + 2
    if m > 1: res -= res // m
    return res

def egcd(a: int, b: int) -> Tuple[int,int,int]:
    if a == 0: return (b, 0, 1)
    g, x, y = egcd(b % a, a)
    return (g, y - (b // a) * x, x)

def inv_mod(a: int, m: int) -> int:
    a %= m
    g, x, _ = egcd(a, m)
    if g != 1:
        raise ValueError("no inverse")
    return x % m

# ---------------- predictor streams ----------------
class Stream:
    def __init__(self, kind:str, a:int, d:int, G:int, k_hat:int, H:int, P1:int):
        self.kind, self.a, self.d, self.G = kind, a, d, G
        self.k_hat, self.H = k_hat, H
        self.P1 = P1
        self.tests = 0
        self.name = f"{kind}:{(a if kind=='AP' else (a,d))}"
        # Alignment for d>1
        if d == 1:
            self.aligned = False
            self.base = P1
            self.M = G
        else:
            g = math.gcd(G, d)          # = rad(d)
            m = d // g                  # modulus for k
            # Find base P0 >= P1 with (P0 ≡ a mod d) to ensure lane congruence
            delta = (a - (P1 % d)) % d
            self.base = P1 + (delta if delta != 0 else 0)
            # Solve (G/g) * k ≡ (a - base)/g (mod m) to get phase k0
            rhs = ((a - (self.base % d)) // g) % m
            inv = inv_mod((G // g) % m, m)
            self.k0 = (rhs * inv) % m
            self.M = G * m              # lcm(G,d)
            self.aligned = True
            self.g, self.m = g, m
        # Derived window start for targeted phase (aligned to lane)
        if d == 1:
            self.t0 = max(0, self.k_hat - 0)  # t plays the role of k here
            self.H_aligned = self.H
        else:
            # map k-window [k̂, k̂+H] to aligned t-window using k = k0 + t*m
            if self.k_hat <= self.k0:
                tstart = 0
            else:
                tstart = (self.k_hat - self.k0 + self.m - 1) // self.m
            self.t0 = tstart
            self.H_aligned = max(1, (self.H + self.m - 1) // self.m)

def build_streams(P1:int, ap_set:List[int], lane_set:List[Tuple[int,int]], c_const:float) -> List[Stream]:
    streams = []
    # d=1
    for a in ap_set:
        if a % 2 == 0: continue
        G = 2*rad(a)
        phiG = phi_from_factor(G)
        k_hat = math.ceil((phiG/G) * math.log(P1))
        H = math.ceil(c_const * math.log(P1))
        streams.append(Stream("AP", a, 1, G, k_hat, H, P1))
    # d>1
    for (a,d) in lane_set:
        if (a % 2 == 0) or (d % 2 == 0) or (math.gcd(a,d) != 1): continue
        G = 2*rad(a)*rad(d)
        phiG = phi_from_factor(G)
        k_hat = math.ceil((phiG/G) * math.log(P1))
        H = math.ceil(c_const * math.log(P1))
        streams.append(Stream("LANE", a, d, G, k_hat, H, P1))
    return streams

# ---------------- runner (hybrid schedule) ----------------
def run_predictor(P1:int,
                  P2:Optional[int],
                  ap_set:List[int],
                  lane_set:List[Tuple[int,int]],
                  K0:int=9, c_const:float=2.0, alpha:float=1.0, max_expands:int=5) -> Dict:
    streams = build_streams(P1, ap_set, lane_set, c_const)
    tested_union = set()
    found = None
    phase = None

    # QUICK PASS
    for t in range(1, K0+1):
        for s in streams:
            if s.d == 1:
                cand = s.base + t * s.G
            else:
                # first aligned candidate on lane at or after base: k = k0 + t*m
                k = s.k0 + t*s.m
                cand = s.base + k*s.G
            s.tests += 1
            if cand not in tested_union:
                tested_union.add(cand)
                if is_probable_prime(cand):
                    found, phase = (cand, s, t), "quick"
                    break
        if found: break

    # TARGETED WINDOWS
    expands = 0
    while not found and expands <= max_expands:
        for s in streams:
            if s.d == 1:
                k_start = s.k_hat
                k_end   = s.k_hat + s.H
                for k in range(k_start, k_end + 1):
                    cand = s.base + k * s.G
                    s.tests += 1
                    if cand in tested_union: continue
                    tested_union.add(cand)
                    if is_probable_prime(cand):
                        found, phase = (cand, s, k), "targeted"
                        break
            else:
                t_start = s.t0
                t_end   = s.t0 + s.H_aligned
                for t in range(t_start, t_end + 1):
                    k = s.k0 + t*s.m
                    if k <= 0: continue
                    cand = s.base + k * s.G
                    s.tests += 1
                    if cand in tested_union: continue
                    tested_union.add(cand)
                    if is_probable_prime(cand):
                        found, phase = (cand, s, t), "targeted"
                        break
            if found: break
        if not found:
            expands += 1
            bump = math.ceil(alpha * math.log(P1))
            for s in streams:
                if s.d == 1:
                    s.H += bump
                else:
                    add = max(1, (bump + s.m - 1)//s.m)
                    s.H_aligned += add

    ap_tests   = sum(s.tests for s in streams if s.kind=="AP")
    lane_tests = sum(s.tests for s in streams if s.kind=="LANE")

    return {
        "P1": P1, "P2": P2, "phase": phase or "none",
        "hit": (found[0] if found else None),
        "gap": ((found[0]-P1) if found else None),
        "AP_tests": ap_tests, "LANE_tests": lane_tests,
        "UNION_unique_tested": len(tested_union),
        "streams": [{
            "name": s.name, "a": s.a, "d": s.d, "G": s.G,
            "phi_over_G": phi_from_factor(s.G)/s.G,
            "k_hat": s.k_hat,
            "H_effective": (s.H if s.d==1 else s.H_aligned),
            "tests": s.tests
        } for s in streams]
    }

# ---------------- CLI & CSV ----------------
def main():
    p = argparse.ArgumentParser(description="Strict-ADN predictor (aligned, hybrid schedule)")
    g = p.add_mutually_exclusive_group(required=False)
    g.add_argument("--P1", type=int, help="Exact anchor prime P1 (will be used as-is)")
    g.add_argument("--seed", type=int, default=170141183460469231731687303715884105727,
                   help="Find P1 = prev_prime_at_or_below(seed) (default: 1e10+7)")
    p.add_argument("--ap", type=int, nargs="*", default=[3,5,7,11],
                   help="AP a-values (d=1)")
    p.add_argument("--lane", type=str, nargs="*", default=["3,7","5,11"],
                   help="Lane pairs 'a,d' (both odd; gcd(a,d)=1)")
    p.add_argument("--K0", type=int, default=9, help="Quick-pass count per stream")
    p.add_argument("--c", type=float, default=2.0, help="H = ceil(c * log P1)")
    p.add_argument("--alpha", type=float, default=1.0, help="Expansion increment α·log P1")
    p.add_argument("--max-expands", type=int, default=5, help="Max expansion rounds")
    p.add_argument("--out-prefix", type=str, default="content/sample_data/predictor_output",
                   help="Prefix for CSV files")
    # Add this line to parse known arguments and ignore unknown ones
    args, unknown = p.parse_known_args()

    # Anchor
    if args.P1:
        P1 = args.P1
        if not is_probable_prime(P1):
            print(f"[warn] --P1={P1} is not prime by Miller–Rabin; continue anyway.", file=sys.stderr)
    else:
        P1 = prev_prime_at_or_below(int(args.seed))
    P2 = next_prime_above(P1)

    # Lane set
    lane_set = []
    for s in (args.lane or []):
        a_s, d_s = s.split(",")
        lane_set.append((int(a_s), int(d_s)))

    # Run
    t0 = time.time()
    res = run_predictor(P1, P2, args.ap, lane_set, K0=args.K0,
                        c_const=args.c, alpha=args.alpha, max_expands=args.max_expands)
    elapsed = time.time() - t0

    # Console summary
    print("=== SUMMARY ===")
    print(f"P1={res['P1']}  P2={res['P2']}  phase={res['phase']}")
    print(f"hit={res['hit']}  gap={res['gap']}")
    print(f"AP_tests={res['AP_tests']}  LANE_tests={res['LANE_tests']}  UNION_unique={res['UNION_unique_tested']}")
    print(f"elapsed_sec={elapsed:.3f}")
    print("\nPer-stream:")
    for s in res["streams"]:
        print(f"  {s['name']:>12} | G={s['G']:<4}  φ(G)/G={s['phi_over_G']:.6f}  k̂={s['k_hat']}  H*={s['H_effective']}  tests={s['tests']}")

    # CSVs
    if args.out_prefix:
        # Ensure the directory exists
        output_dir = os.path.dirname(args.out_prefix)
        if output_dir and not os.path.exists(output_dir):
            os.makedirs(output_dir)

        params_path  = f"{args.out_prefix}_stream_params_P1_{P1}.csv"
        overall_path = f"{args.out_prefix}_benchmark_P1_{P1}.csv"
        with open(params_path, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=["stream","a","d","G","phi_over_G","k_hat","H_effective","tests"])
            w.writeheader()
            for s in res["streams"]:
                row = dict(stream=s["name"], a=s["a"], d=s["d"], G=s["G"],
                           phi_over_G=f"{s['phi_over_G']:.6f}",
                           k_hat=s["k_hat"], H_effective=s["H_effective"], tests=s["tests"])
                w.writerow(row)
        with open(overall_path, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=["P1","P2","phase","hit","gap","AP_tests","LANE_tests","UNION_unique_tested","elapsed_sec"])
            w.writeheader()
            w.writerow({
                "P1": res["P1"], "P2": res["P2"], "phase": res["phase"], "hit": res["hit"], "gap": res["gap"],
                "AP_tests": res["AP_tests"], "LANE_tests": res["LANE_tests"],
                "UNION_unique_tested": res["UNION_unique_tested"], "elapsed_sec": f"{elapsed:.3f}"
            })
        print(f"\nCSV written:\n  {params_path}\n  {overall_path}")

if __name__ == "__main__":
    main()