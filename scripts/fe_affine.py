"""Bonus (idea 5): periodic affine keys / per-class Vigenere-Beaufort mixing.
Position class r = i mod L uses x(C) = a_r * x(P) + b_r (mod 26) with its own unit a_r (12 choices) and shift b_r.
Special case a_r in {+1,-1}: each class is Vigenere or Beaufort (e.g. a 7-letter sign pattern).
x() = alphabet index: Quagmire III over 842k dictionary + 116 clue alphabets, plus P=AZ/KA or C=AZ/KA arrangements.
Exact; planted solution recovered first; shuffled CTs calibrate."""
import sys, time, numpy as np, k4, fe_common as F
UNITS = [a for a in range(26) if np.gcd(a, 26) == 1]
def classes(L):
    cl = {}
    for j, p in enumerate(F.POS): cl.setdefault(p % L, []).append(j)
    return [m for m in cl.values() if len(m) >= 3]
def checks(L, aset):
    return sum(len(m) - 2 for m in classes(L)) if len(aset) > 2 else sum(len(m) - 1 for m in classes(L))
def fits(xP, xC, L, aset):
    ok = np.ones(len(xP), bool)
    for m in classes(L):
        cls_ok = np.zeros(len(xP), bool)
        for a in aset:
            b = (xC[:, m] - a * xP[:, m]) % 26
            cls_ok |= (b == b[:, :1]).all(axis=1)
        ok &= cls_ok
    return ok
def scan(ct, AL, Ls, aset):
    A = F.alph_array([a for a, _ in AL]); INV = F.inv_index(A).astype(np.int64); n = len(A)
    AZi = np.tile(np.arange(26), (n, 1)); KAi = np.tile(np.array([k4.KA.index(chr(65 + x)) for x in range(26)]), (n, 1))
    C = np.array([ord(c) - 65 for c in ct])[F.POS]; hits = []
    for pn, Pi, cn, Ci in (("A", INV, "A", INV), ("AZ", AZi, "A", INV), ("A", INV, "AZ", AZi), ("KA", KAi, "A", INV), ("A", INV, "KA", KAi), ("AZ", AZi[:1], "AZ", AZi[:1]), ("KA", KAi[:1], "KA", KAi[:1])):
        xP = Pi[:, F.PL]; xC = Ci[:, C]
        for L in Ls:
            for r in np.nonzero(fits(xP, xC, L, aset))[0]: hits.append((AL[r][1] if len(xP) > 1 else pn, pn, cn, L, checks(L, aset)))
    return hits
if __name__ == "__main__":
    AL = F.all_alphabets(); rng = np.random.default_rng(21)
    LS_aff = list(range(1, 14)); LS_sgn = [L for L in range(1, 27) if checks(L, [1, 25]) >= 7]
    print("affine periods (classes with >= 3 crib letters checked):", {L: [len(m) for m in classes(L)] for L in LS_aff})
    print("sign-mix periods with >= 7 checks:", {L: checks(L, [1, 25]) for L in LS_sgn})
    # planted
    for L, aset in ((7, UNITS), (9, [1, 25])):
        j = int(rng.integers(len(AL))); X = AL[j][0]
        a = [int(rng.choice(aset)) for _ in range(L)]; b = [int(v) for v in rng.integers(0, 26, L)]
        pt = F.planted_plaintext(rng); ct = "".join(X[(a[i % L] * X.index(c) + b[i % L]) % 26] for i, c in enumerate(pt))
        sub = [AL[x] for x in rng.choice(len(AL), 20000, replace=False)] + [AL[j]]
        h = scan(ct, sub, [L], aset)
        print(f"  planted L={L} {'affine' if len(aset) > 2 else 'sign-mix'} {AL[j][1]}: recovered={any(x[0] == AL[j][1] and x[1] == 'A' and x[2] == 'A' for x in h)}; other hits {len(h) - 1}", flush=True)
    for label, ct in [("K4", k4.CT)] + [(f"shuffled{t}", "".join(rng.permutation(list(k4.CT)))) for t in range(3)]:
        t0 = time.time(); ha = []; hs = []
        for s in range(0, len(AL), 100000):
            ha += scan(ct, AL[s:s + 100000], LS_aff, UNITS); hs += scan(ct, AL[s:s + 100000], LS_sgn, [1, 25])
        nvec = len(AL) * 5
        ea = sum(nvec * np.prod([min(1.0, 12 * 26.0 ** -(len(m) - 1)) for m in classes(L)]) for L in LS_aff)
        es = sum(nvec * np.prod([min(1.0, 2 * 26.0 ** -(len(m) - 1)) for m in classes(L)]) for L in LS_sgn)
        print(f"{label}: affine hits {len(ha)} (chance ~{ea:.2f}), sign-mix hits {len(hs)} (chance ~{es:.2f}) ({time.time()-t0:.0f}s)", flush=True)
        for x in (ha + hs)[:10]: print("    ", x)
