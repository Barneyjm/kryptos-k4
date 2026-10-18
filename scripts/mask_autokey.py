"""Arbitrary mask + autokey, solved directly.
(a) plaintext mask x, ciphertext autokey:  x[P_i] = num(C_i) -/+ num(C_{i-lag})  (vig / beau / variant)
(b) ciphertext mask y, plaintext autokey:  y[C_i] = num(P_i) +/- num(P_{i-lag})  (where both plaintext letters known)
Consistent iff the forced value is the same for each repeated letter and distinct across letters."""
import random, k4
CR = dict(zip(k4.CPOS.tolist(), k4.CPL))
def check(forced):
    m = {}
    for letter, v in forced:
        if m.setdefault(letter, v) != v: return None
    return m if len(set(m.values())) == len(m) else None
def run(ct):
    hits = []
    for an, X in (("AZ", k4.AZ), ("KA", k4.KA)):
        n = lambda c: X.index(c)
        for lag in range(1, 97):
            for vn, f in (("C-C'", lambda a, b: (a - b) % 26), ("C+C'", lambda a, b: (a + b) % 26), ("C'-C", lambda a, b: (b - a) % 26)):
                forced = [(p, f(n(ct[i]), n(ct[i - lag]))) for i, p in CR.items() if i - lag >= 0]
                if len(forced) >= 14 and check(forced): hits.append(("plain-mask+ctak", an, lag, vn, len(forced)))
                forced = [(ct[i], f(n(p), n(CR[i - lag]))) for i, p in CR.items() if i - lag in CR]
                if len(forced) >= 8 and check(forced): hits.append(("cipher-mask+ptak", an, lag, vn, len(forced)))
    return hits
if __name__ == "__main__":
    random.seed(3); P = [random.choice(k4.AZ) for _ in range(97)]
    for s, w in k4.CRIBS.items(): P[s:s+len(w)] = w
    mask = random.sample(k4.AZ, 26); X = k4.AZ; c = []
    for i, p in enumerate(P):
        prev = X.index(c[i - 5]) if i >= 5 else 7
        c.append(X[(X.index(mask[X.index(p)]) + prev) % 26])
    print("planted (mask + ciphertext autokey lag 5):", [h for h in run("".join(c)) if h[2] == 5])
    print("K4:", run(k4.CT) or "no consistent mask + autokey at any lag 1-96")
