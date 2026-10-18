"""Periodic polyalphabetic solver: for period L, try alphabets (AZ, KA, keyword), vig/beau; per-column shift chosen by
chi-square vs English, then score the full decryption with quadgrams."""
import json, numpy as np, sys, k4
from quad import Q
EN = np.array([8.2,1.5,2.8,4.3,12.7,2.2,2.0,6.1,7.0,.15,.77,4.0,2.4,6.7,7.5,1.9,.095,6.0,6.3,9.1,2.8,.98,2.4,.15,2.0,.074]) / 100
def qscore(t): a = np.array([ord(c)-65 for c in t]); return Q[((a[:-3]*26+a[1:-2])*26+a[2:-1])*26+a[3:]].mean()
def solve(ct, L, X, var):
    a = np.array([X.index(c) for c in ct]); key = []
    for r in range(L):
        col = a[r::L]; best = None
        for k in range(26):
            p = (col - k) % 26 if var == "vig" else (k - col) % 26
            letters = np.array([ord(X[v]) - 65 for v in p]); obs = np.bincount(letters, minlength=26)
            chi = ((obs - EN * len(col))**2 / (EN * len(col))).sum()
            if best is None or chi < best[0]: best = (chi, k)
        key.append(best[1])
    pt = "".join(X[((a[i] - key[i % L]) if var == "vig" else (key[i % L] - a[i])) % 26] for i in range(len(ct)))
    return qscore(pt), "".join(X[k] for k in key), pt
if __name__ == "__main__":
    P = {p["slug"]: p["ct"] for p in json.load(open("ctf_puzzles.json"))}
    ct = P[sys.argv[1]]; res = []
    for an, X in (("AZ", k4.AZ), ("KA", k4.KA)):
        for var in ("vig", "beau"):
            for L in range(1, 31): res.append((*solve(ct, L, X, var), an, var, L))
    for r in sorted(res, reverse=True)[:3]: print(f"{r[0]:.2f} {r[3]} {r[4]} L={r[5]} key={r[1]}\n   {r[2][:150]}")
