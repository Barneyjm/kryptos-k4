"""Constrained-key decoding: Beaufort/Vigenere with key values restricted to the set seen at crib positions.
Beam search for the most English-like plaintext; compare against null ciphertexts processed identically."""
import numpy as np, random, k4
from collections import Counter
from quad import Q
AZ = k4.AZ; CRIB = dict(zip(k4.CPOS.tolist(), k4.CPL))
F = {"beau": (lambda c, p: (c + p) % 26, lambda c, k: (k - c) % 26), "vig": (lambda c, p: (c - p) % 26, lambda c, k: (c - k) % 26)}
def decode(ct, X, var, beam=3000):
    keyf, ptf = F[var]
    ks = [keyf(X.index(ct[i]), X.index(p)) for i, p in CRIB.items()]
    cnt = Counter(ks); S = list(cnt); logk = {k: np.log10(cnt[k] / len(ks)) for k in S}
    toAZ = [AZ.index(c) for c in X]
    states = [(0.0, ())]                                   # (score, plaintext AZ indices)
    for i in range(97):
        if i in CRIB: cands = [(AZ.index(CRIB[i]), 0.0)]
        else: cands = [(toAZ[ptf(X.index(ct[i]), k)], logk[k]) for k in S]
        new = []
        for sc, pt in states:
            for p, lk in cands:
                s = sc + lk
                if len(pt) >= 3: s += Q[((pt[-3]*26 + pt[-2])*26 + pt[-1])*26 + p]
                new.append((s, pt + (p,)))
        new.sort(key=lambda x: -x[0]); states = new[:beam]
    sc, pt = states[0]; txt = "".join(AZ[p] for p in pt)
    noncrib = [j for j in range(3, 97) if j not in CRIB and j - 1 not in CRIB]
    q = np.mean([Q[((pt[j-3]*26 + pt[j-2])*26 + pt[j-1])*26 + pt[j]] for j in range(3, 97)])
    return q, len(S), txt
if __name__ == "__main__":
    random.seed(5)
    for an, X in k4.ALPHS.items():
        for var in F:
            q, n, txt = decode(k4.CT, X, var)
            print(f"K4  {an} {var:4s} key set {n:2d} values  english score {q:.2f}  {txt}", flush=True)
    nulls = []
    for t in range(12):
        ct = list(k4.CT); random.shuffle(ct); ct = "".join(ct)
        q, n, _ = decode(ct, k4.AZ, "beau"); nulls.append((q, n))
    print("null (shuffled K4, AZ beau): scores", [f"{q:.2f}/{n}" for q, n in nulls])
