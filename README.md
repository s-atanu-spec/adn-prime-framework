
# ADN Prime Framework — Scripts & Predictors

**Authors:** Atanu Sarkar (Independent Researcher) & Ashva (AI Research Partner)  
**Repository:** `adn-prime-framework`  
**License:** MIT

Dirichlet‑aligned, parity+coprimality gated exploration of primes on arithmetic progressions (“lanes”), with a fiber/S‑band viewpoint and practical predictors. This repo contains *standalone, no‑I/O* Python scripts that are friendly to online compilers and local runs alike.

> Core idea (strict ADN): `P = a + n·d` with **odd** `a,d`, **even** `n`, and the coprimality gates `gcd(a,d)=gcd(a,n)=gcd(d,n)=1`. We study survivors, lane locality, fiber counts for a fixed prime, clustering, and next‑prime predictors that step on small square‑free moduli.

---

## Quick start

- Python 3.9+ recommended (works on CPython without external dependencies).
- Clone or download this repo, then:

```bash
cd adn-prime-framework
python3 scripts/benchmark-script.py
```

> All scripts print to stdout; no CSVs unless explicitly requested by flags.

---

## What’s inside

| File | Purpose |
|---|---|
| `scripts/adn-script-all-in-one.py` | A bundle of **ten** small, standalone programs (definitions, Backbone certificate + fiber sketch, clustering demo, S‑band snapshot, lane‑locality sweep, strict‑ADN predictor, and per‑lane tables/union/density). Run it directly to see default demos; each block has a CONFIG section at the top. |
| `scripts/benchmark-script.py` | Compares baseline next‑prime finders (odd‑step, wheel‑30, segmented sieve) vs. **ADN predictors** (AP‑slice and lane‑stride). Reports unique candidates tested, first‑hit index, and wall‑clock. |
| `scripts/predictor.py` | CLI **strict‑ADN predictor** (AP‑slice + lane‑aware with **congruence alignment**). Hybrid schedule: quick pass then a targeted window around k-hat ≈ (φ(G)/G)·log P1. Optional CSV outputs. |
| `scripts/predictor-benchmark.py` | Variant of the predictor with a different default seed/output path; good for producing sample CSV in `content/sample_data/`. |
| `scripts/gap-decomposition-audit.py` | Audits **gap decomposition** for consecutive primes using a canonical ADN mapping (S‑band or n=2 backbone). Verifies identities, evenness, same‑lane multiples, and impossibility of twin gaps on d>1 lanes. |
| `scripts/prime_sums_parity_check.py` | Referee‑ready parity audit for sums of consecutive primes: sum is always even; for twins it’s divisible by 4. |

---

## Usage examples

### Benchmark table
```bash
python3 scripts/benchmark-script.py
```

### Strict‑ADN predictor (hybrid schedule)
- Use a big seed to anchor near that region:
```bash
python3 scripts/predictor.py --seed 1000000000000066600000000000001 --ap 3 5 7 11 --lane 3,7 5,11 --K0 9 --c 2.0 --alpha 1.0 --max-expands 5 --out-prefix content/sample_data/predictor_run
```
- Or pin an exact `P1` (use with care if it's not prime):
```bash
python3 scripts/predictor.py --P1 1000000007 --out-prefix content/sample_data/predictor_run
```

### Predictor (alt benchmark build)
```bash
python3 scripts/predictor-benchmark.py --seed 170141183460469231731687303715884105727 --out-prefix content/sample_data/predictor_1e10
```

### Gap decomposition audit
```bash
python3 scripts/gap-decomposition-audit.py
# Tweak MODE='sband' vs 'backbone_n2' in the CONFIG section for different canonicalization.
```

### Prime sums parity check
```bash
python3 scripts/prime_sums_parity_check.py --seed 10000 --pairs 1000 --show 12
```

### All‑in‑one bundle
```bash
python3 scripts/adn-script-all-in-one.py
# Each sub‑program has a CONFIG block near the top of its section.
```

---

## Repo layout

```
adn-prime-framework/
├─ scripts/
│  ├─ adn-script-all-in-one.py
│  ├─ benchmark-script.py
│  ├─ predictor.py
│  ├─ predictor-benchmark.py
│  ├─ gap-decomposition-audit.py
│  └─ prime_sums_parity_check.py
├─ content/
│  └─ sample_data/        # CSVs from predictor runs (ignored by Git)
├─ LICENSE
├─ .gitignore
└─ README.md
```

---

## Notes & conventions

- **Strict ADN gates everywhere**: odd `a,d`; even `n`; `gcd(a,d)=gcd(a,n)=gcd(d,n)=1`.
- Deterministic Miller–Rabin bases {2,3,5,7,11,13,17} give correctness for 64‑bit and are excellent in practice for much larger demo ranges.
- Predictor streams use small square‑free moduli G (or lcm(G,d) for lanes) and **congruence‑aligned k** to respect the lane P ≡ a (mod d).

---

## Citation

If you use this code in research, please cite as:
> Atanu Sarkar & Ashva. *ADN Prime Framework — Scripts & Predictors*. 2025. GitHub: s‑atanu‑spec/adn‑prime‑framework.

A `CITATION.cff` can be added later if needed.

---

## License

MIT © 2025 Atanu Sarkar & Ashva
