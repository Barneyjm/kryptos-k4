"""T2: Wheatstone cryptograph (1867) and its hand-paper equivalent: a 27-cell outer ring (26 letters + blank)
and a 26-cell inner ring geared 1:1 per step, so the inner ring gains one letter per outer revolution
("change the language base": 27 vs 26). Encrypt: move the pointer clockwise from the previous plaintext letter
to the current one (d = 1..27 steps), the inner ring moves d steps too, read the cipher letter on the inner ring.
Within a crib the inner offsets S_k = sum d are fully determined by the outer ring, so for ANY inner alphabet the
crib letters must satisfy  (S_a == S_b mod 26)  <=>  (C_a == C_b), with one unknown shift delta between cribs.
Outer rings: every base alphabet (AZ, KA, 116 clue, 842k dictionary; reversals = counter-clockwise) with the
blank at each of the 26 cyclically distinct slots. Inner ring: arbitrary (solved exactly). Exact test.
"""
import sys, time, numpy as np, k4
import fe_common as F
from multiprocessing import Pool

W1, W2 = "EASTNORTHEAST", "BERLINCLOCK"
def crib_cipher(ct):
    return [ct[21:34], ct[63:74]]

def offsets(posmat, word):
    """posmat (n,27): position on ring of symbol index (0..25 letters, 26 blank). Returns S (n,len) mod 26."""
    w = [ord(c) - 65 for c in word]
    P = posmat[:, w]                                   # (n,len)
    d = ((P[:, 1:] - P[:, :-1] - 1) % 27) + 1
    S = np.concatenate([np.zeros((len(P), 1), np.int64), np.cumsum(d, axis=1)], axis=1)
    return S % 26

def eqmat(s):
    return np.array([[a == b for b in s] for a in s])

def check(posmat, ct):
    c1, c2 = crib_cipher(ct)
    S1 = offsets(posmat, W1); S2 = offsets(posmat, W2)
    ok = np.ones(len(posmat), bool)
    for S, c in ((S1, c1), (S2, c2)):
        E = eqmat(c); n = len(c)
        for a in range(n):
            for b in range(a + 1, n):
                ok &= ((S[:, a] == S[:, b]) == E[a, b])
    idx = np.nonzero(ok)[0]; res = []
    X = np.array([[x == y for y in c2] for x in c1])  # (13,11)
    for r in idx:
        for dl in range(26):
            s2 = (S2[r] + dl) % 26
            if np.array_equal(S1[r][:, None] == s2[None, :], X):
                res.append((int(r), dl))
    return ok.sum(), res

def rings_for(alph):
    """26 cyclically distinct rings: blank inserted after position j (j=0..25)."""
    return [alph[:j + 1] + "#" + alph[j + 1:] for j in range(26)]

def posmat_of(rings):
    M = np.zeros((len(rings), 27), np.int64)
    for i, r in enumerate(rings):
        for p, ch in enumerate(r): M[i, 26 if ch == "#" else ord(ch) - 65] = p
    return M

def job(args):
    ct, chunk = args
    rings = [(r, n) for a, n in chunk for r in rings_for(a)]
    PM = posmat_of([r for r, _ in rings])
    internal, res = check(PM, ct)
    return int(internal), [(rings[r][1], rings[r][0], dl) for r, dl in res]

def encrypt(pt, ring, inner, q0):
    o = ring.index(pt[0]) - 1; q = q0; out = []
    for ch in pt:
        t = ring.index(ch); d = ((t - o - 1) % 27) + 1
        o = t; q = (q + d) % 26; out.append(inner[q])
    return "".join(out)

def full(ct, AL, label, pool):
    t = time.time(); chunks = [(ct, AL[i:i + 5000]) for i in range(0, len(AL), 5000)]
    internal = 0; hits = []
    for a, h in pool.imap_unordered(job, chunks):
        internal += a; hits += h
    print(f"{label}: {len(AL)*26:,} outer rings; pass within-crib checks: {internal:,}; pass both cribs (some delta): {len(hits)} ({time.time()-t:.0f}s)", flush=True)
    return hits

if __name__ == "__main__":
    AL = F.all_alphabets()
    rng = np.random.default_rng(3)
    print("planted-solution validation:", flush=True)
    for trial in range(3):
        a = AL[int(rng.integers(len(AL)))]; ring = rings_for(a[0])[int(rng.integers(26))]
        inner = AL[int(rng.integers(len(AL)))][0]
        pt = F.planted_plaintext(rng)
        ct = encrypt(pt, ring, inner, int(rng.integers(26)))
        internal, res = check(posmat_of([ring] + rings_for(AL[5][0])), ct)
        print(f"  planted ring {ring} ({a[1]}): recovered={any(r == 0 for r, _ in res)}", flush=True)
        assert any(r == 0 for r, _ in res)
    with Pool(int(sys.argv[1]) if len(sys.argv) > 1 else 3) as pool:
        h = full(k4.CT, AL, "K4", pool)
        for x in h[:30]: print("   ", x, flush=True)
        r2 = np.random.default_rng(5)
        for t in range(3):
            sh = "".join(r2.permutation(list(k4.CT)))
            full(sh, AL, f"shuffled{t}", pool)
