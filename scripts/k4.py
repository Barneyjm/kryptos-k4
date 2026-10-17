"""Kryptos K4 crib-consistency harness.
Model A: P -> periodic sub -> transposition -> C
Model B: P -> transposition -> periodic sub -> C
Sub = Vigenere or Beaufort over the standard alphabet or the KRYPTOS-keyed alphabet.
Transposition = irregular columnar (width w, column read order), both directions.
"""
import itertools, sys, numpy as np
CT = "OBKRUOXOGHULBSOLIFBBWFLRVQQPRNGKSSOTWTQSJQSSEKZZWATJKLUDIAWINFBNYPVTTMZFPKWGDKZXTJCDIGKUHUAUEKCAR"
N = len(CT); assert N == 97
CRIBS = {21: "EASTNORTHEAST", 63: "BERLINCLOCK"}   # 0-indexed starts
CPOS = np.array([s + i for s, w in CRIBS.items() for i in range(len(w))])
CPL = "".join(CRIBS.values())
AZ = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
KA = "KRYPTOSABCDEFGHIJLMNQUVWXZ"
ALPHS = {"AZ": AZ, "KA": KA}
PERIODS = range(1, 27)

def order_from(w, colorder):
    return [i for c in colorder for i in range(c, N, w)]

def check_batch(dests, label, out):
    """dests: (n,N) int array, dest[i] = ciphertext position of plaintext position i."""
    dc = dests[:, CPOS]                                   # (n,24) cipher positions of crib letters
    for an, X in ALPHS.items():
        ci = np.array([X.index(c) for c in CT])[dc]       # (n,24)
        pi = np.array([X.index(c) for c in CPL])[None, :]
        for kn, k in (("vig", (ci - pi) % 26), ("beau", (ci + pi) % 26)):
            for model in "AB":
                for L in PERIODS:
                    r = (CPOS[None, :] % L) if model == "A" else (dc % L)
                    if model == "A":
                        r = np.broadcast_to(r, k.shape)
                    same_r = r[:, :, None] == r[:, None, :]
                    diff_k = k[:, :, None] != k[:, None, :]
                    ok = ~(same_r & diff_k).any(axis=(1, 2))
                    # require real constraint: at least 8 crib letters sharing residues
                    for idx in np.nonzero(ok)[0]:
                        rr = r[idx]; cons = len(rr) - len(set(rr.tolist()))
                        if cons >= 8:
                            out.append((label(idx), an, kn, model, L, cons))

def run_columnar(maxw):
    hits = []
    for w in range(2, maxw + 1):
        perms = list(itertools.permutations(range(w)))
        for start in range(0, len(perms), 20000):
            chunk = perms[start:start + 20000]
            orders = np.array([order_from(w, p) for p in chunk])          # C[n] = P[order[n]]
            dest = np.empty_like(orders); rows = np.arange(len(chunk))[:, None]
            dest[rows, orders] = np.arange(N)[None, :]
            for direction, D in (("enc", dest), ("dec", orders)):       # dec: inverse transposition
                check_batch(D, lambda i, w=w, c=chunk, d=direction: f"col w={w} {d} {c[i]}", hits)
        print(f"w={w} done, hits so far {len(hits)}", file=sys.stderr, flush=True)
    return hits

def identity_test():
    hits = []
    check_batch(np.arange(N)[None, :], lambda i: "no transposition", hits)
    return hits

if __name__ == "__main__":
    print("identity:", identity_test())
    hits = run_columnar(int(sys.argv[1]) if len(sys.argv) > 1 else 8)
    for h in hits[:200]: print(h)
    print("total hits", len(hits))
