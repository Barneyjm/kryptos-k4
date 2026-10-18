"""T3: "Change the language base" -> additive cipher over a 27-symbol base instead of Z26.
The 26 letters plus one extra symbol (#: space / '?' / null) are laid out in a keyed sequence of 27 cells and
combined in one of three abelian groups of order 27:
  Z27      (mod-27 Vigenere, e.g. alphabet + '?')
  Z3^3     (trifid cube: layer/row/column digits added mod 3, 'Delastelle trifid coordinates' without fractionation)
  Z9xZ3    (two digit splits: (e//3, e%3) and (e%9, e//9))
Keys: periodic keyword L = 1..26 (vig/beau; variant == vig for equality tests), ciphertext autokey and plaintext
autokey (lag 1..26; vig, beau, variant). Embeddings: every base alphabet (AZ, KA, 116 clue, 842k dictionary) with
# inserted at each of 27 slots. Exact tests. A planted solution is recovered first; shuffled CTs calibrate.
Note: K4 contains no '#'-cell output; for a uniform 27-symbol cipher that happens with p = (26/27)^97 = 2.6%.
"""
import sys, time, numpy as np, k4
import fe_common as F
from multiprocessing import Pool

GROUPS = ["Z27", "Z3^3", "Z9xZ3a", "Z9xZ3b"]
MIN = 7
def digits(e, g):
    if g == "Z27": return [(e, 27)]
    if g == "Z3^3": return [(e // 9, 3), ((e // 3) % 3, 3), (e % 3, 3)]
    if g == "Z9xZ3a": return [(e // 3, 9), (e % 3, 3)]
    if g == "Z9xZ3b": return [(e % 9, 9), (e // 9, 3)]
def combine(c, p, g, op):
    """c, p: int arrays of cell indices. op: 'vig' k=c-p, 'beau' k=c+p, 'var' k=p-c. Returns encoded key 0..26."""
    out = 0; mult = 1
    for (cd, m), (pd, _) in zip(digits(c, g), digits(p, g)):
        k = (cd - pd) % m if op == "vig" else ((cd + pd) % m if op == "beau" else (pd - cd) % m)
        out = out + k * mult; mult *= m
    return out
def encode(e, g):
    """re-encode a cell index into the same mixed-radix code that combine() returns (for autokey comparisons)."""
    out = 0; mult = 1
    for d, m in digits(e, g): out = out + d * mult; mult *= m
    return out
def apply(p, k, g, op):
    """encrypt: return cell index of c from p and encoded key k."""
    ks = []; kk = k
    for _, m in digits(np.int64(0), g): ks.append(kk % m); kk //= m
    cd = []
    for (pd, m), kd in zip(digits(p, g), ks):
        cd.append((pd + kd) % m if op == "vig" else ((kd - pd) % m if op == "beau" else (pd - kd) % m))
    # rebuild cell index from digits
    if g == "Z27": return cd[0]
    if g == "Z3^3": return cd[0] * 9 + cd[1] * 3 + cd[2]
    if g == "Z9xZ3a": return cd[0] * 3 + cd[1]
    if g == "Z9xZ3b": return cd[1] * 9 + cd[0]

ENG = F.Engine(((f"periodic L={L}", [p % L for p in F.POS]) for L in range(1, 27)), MIN)

def embeddings(alphs):
    """alphs: list of (alphabet, name) -> (E (n,26) cell index of each letter, names)."""
    E = []; names = []
    for a, n in alphs:
        base = [ord(c) - 65 for c in a]
        for s in range(27):
            seq = base[:s] + [26] + base[s:]
            pos = np.empty(27, np.int64); pos[seq] = np.arange(27)
            E.append(pos[:26]); names.append(f"{n}#@{s}")
    return np.array(E), names

def test(ct, E):
    C = np.array([ord(ch) - 65 for ch in ct]); Cc = C[F.POS]
    hits = []
    for g in GROUPS:
        pe = E[:, F.PL]; ce = E[:, Cc]
        for op in ("vig", "beau"):
            K = combine(ce, pe, g, op)
            for r, m in ENG.fits(K):
                hits.append((int(r), g, op, ENG.names[m], int(ENG.checks[m])))
        for op in ("vig", "beau", "var"):
            K = combine(ce, pe, g, op)
            for lag in range(1, 27):
                # ciphertext autokey
                cols = [j for j, p in enumerate(F.POS) if p - lag >= 0]
                src = E[:, C[[F.POS[j] - lag for j in cols]]]
                ok = (K[:, cols] == encode(src, g)).all(axis=1)
                for r in np.nonzero(ok)[0]: hits.append((int(r), g, op, f"ct-autokey lag={lag}", len(cols)))
                # plaintext autokey (only where the source letter is a crib letter)
                pcols = [j for j, p in enumerate(F.POS) if (p - lag) in F.POS]
                if len(pcols) >= 6:
                    srcp = E[:, F.PL[[F.POS.index(F.POS[j] - lag) for j in pcols]]]
                    ok = (K[:, pcols] == encode(srcp, g)).all(axis=1)
                    for r in np.nonzero(ok)[0]: hits.append((int(r), g, op, f"pt-autokey lag={lag}", len(pcols)))
    return hits

def job(args):
    ct, chunk = args
    E, names = embeddings(chunk)
    return [(names[r],) + h[1:] for h in test(ct, E) for r in [h[0]]]

def encrypt_periodic(pt, E, kw, g, op):
    cell2let = {int(e): chr(65 + i) for i, e in enumerate(E)}
    out = []
    for i, ch in enumerate(pt):
        c = int(apply(np.int64(E[ord(ch) - 65]), kw[i % len(kw)], g, op))
        if c not in cell2let: return None
        out.append(cell2let[c])
    return "".join(out)

def full(ct, AL, label, pool):
    t = time.time(); hits = []
    for h in pool.imap_unordered(job, [(ct, AL[i:i + 4000]) for i in range(0, len(AL), 4000)]): hits += h
    nemb = len(AL) * 27
    exp = sum(nemb * 2 * len(GROUPS) * 27.0 ** (-c) for c in ENG.checks)
    print(f"{label}: {nemb:,} embeddings x {len(GROUPS)} groups; periodic models {list(ENG.names)[-1]} (checks>={MIN}); "
          f"hits {len(hits)} (periodic chance ~{exp:.2f}) ({time.time()-t:.0f}s)", flush=True)
    return hits

if __name__ == "__main__":
    AL = F.all_alphabets()
    rng = np.random.default_rng(9)
    print(f"periodic models kept (>= {MIN} checks): {ENG.names}", flush=True)
    print("planted-solution validation:", flush=True)
    for g, op, L in (("Z27", "vig", 7), ("Z3^3", "beau", 11), ("Z9xZ3a", "vig", 5), ("Z3^3", "vig", 13)):
        j = int(rng.integers(len(AL))); E, names = embeddings([AL[j]]); s = int(rng.integers(27))
        for attempt in range(20000):
            kw = [int(v) for v in rng.integers(0, 27, L)]
            ct = encrypt_periodic(F.planted_plaintext(rng), E[s], kw, g, op)
            if ct: break
        E2, n2 = embeddings([AL[j]] + [AL[x] for x in rng.choice(len(AL), 300, replace=False)])
        h = [x for x in test(ct, E2) if x[0] == s and x[1] == g and x[2] == op and x[3] == f"periodic L={L}"]
        print(f"  planted {g} {op} L={L} {names[s]}: recovered={bool(h)} (needed {attempt+1} key draws for a #-free CT)", flush=True)
        assert h
    # planted ciphertext autokey
    j = int(rng.integers(len(AL))); E, names = embeddings([AL[j]]); s = 4; e = E[s]; lag = 3
    cell2let = {int(v): chr(65 + i) for i, v in enumerate(e)}
    for attempt in range(100000):
        pt = F.planted_plaintext(rng); seed = [int(v) for v in rng.integers(0, 26, lag)]; out = []
        for i, ch in enumerate(pt):
            k = encode(np.int64(e[seed[i]]) if i < lag else np.int64(e[ord(out[i - lag]) - 65]), "Z3^3")
            c = int(apply(np.int64(e[ord(ch) - 65]), int(k), "Z3^3", "vig"))
            if c not in cell2let: break
            out.append(cell2let[c])
        if len(out) == 97: break
    ct = "".join(out)
    h = [x for x in test(ct, E) if x[0] == s and x[1] == "Z3^3" and x[2] == "vig" and x[3] == f"ct-autokey lag={lag}"]
    print(f"  planted Z3^3 ct-autokey lag={lag}: recovered={bool(h)}", flush=True); assert h
    with Pool(int(sys.argv[1]) if len(sys.argv) > 1 else 3) as pool:
        h = full(k4.CT, AL, "K4", pool)
        for x in sorted(h, key=lambda x: -x[4])[:30]: print("   ", x, flush=True)
        r2 = np.random.default_rng(13)
        for t in range(2):
            full("".join(r2.permutation(list(k4.CT))), AL, f"shuffled{t}", pool)
