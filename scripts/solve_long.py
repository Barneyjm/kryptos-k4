import json, sys, numpy as np, k4
from quad import Q
import solve_periodic as SP
def refine(ct, L, X, var, key, rounds=4):
    a = np.array([X.index(c) for c in ct]); conv = np.array([ord(c) - 65 for c in X])
    def dec(k):
        kk = np.array(k)[np.arange(len(a)) % L]
        p = (a - kk) % 26 if var == "vig" else (kk - a) % 26
        return conv[p]
    def sc(pt): return Q[((pt[:-3]*26+pt[1:-2])*26+pt[2:-1])*26+pt[3:]].sum()
    key = list(key); best = sc(dec(key))
    for _ in range(rounds):
        for r in range(L):
            for v in range(26):
                k2 = key[:]; k2[r] = v; s = sc(dec(k2))
                if s > best: best, key = s, k2
    pt = "".join(chr(65 + x) for x in dec(key)); return best / (len(ct) - 3), "".join(X[k] for k in key), pt
if __name__ == "__main__":
    P = {p["slug"]: p["ct"] for p in json.load(open("ctf_puzzles.json"))}
    slug, Ls = sys.argv[1], [int(x) for x in sys.argv[2].split(",")]
    res = []
    for L in Ls:
        for an, X in (("AZ", k4.AZ), ("KA", k4.KA)):
            for var in ("vig", "beau"):
                _, k0, _ = SP.solve(P[slug], L, X, var)
                s, k, pt = refine(P[slug], L, X, var, [X.index(c) for c in k0])
                res.append((s, an, var, L, k, pt))
    for r in sorted(res, reverse=True)[:2]: print(f"{r[0]:.2f} {r[1]} {r[2]} L={r[3]} key={r[4]}\n   {r[5]}")
