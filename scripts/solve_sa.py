import json, sys, math, random, numpy as np, k4
from quad import Q
import solve_periodic as SP
def run(ct, L, X, var, iters=60000, seed=0, init=None):
    rnd = random.Random(seed); a = np.array([X.index(c) for c in ct]); conv = np.array([ord(c) - 65 for c in X]); n = len(a)
    idx = np.arange(n) % L
    def dec(k): kk = k[idx]; return conv[((a - kk) if var == "vig" else (kk - a)) % 26]
    def sc(pt): return Q[((pt[:-3]*26+pt[1:-2])*26+pt[2:-1])*26+pt[3:]].sum()
    key = np.array(init if init is not None else [rnd.randrange(26) for _ in range(L)]); cur = sc(dec(key)); best = (cur, key.copy()); T = 20.0
    for it in range(iters):
        r = rnd.randrange(L); old = key[r]; key[r] = rnd.randrange(26); s = sc(dec(key))
        if s >= cur or rnd.random() < math.exp((s - cur) / T): cur = s
        else: key[r] = old
        if cur > best[0]: best = (cur, key.copy())
        T = max(0.5, T * 0.9998)
    k = best[1]; return best[0] / (n - 3), "".join(X[v] for v in k), "".join(chr(65 + x) for x in dec(k))
if __name__ == "__main__":
    P = {p["slug"]: p["ct"] for p in json.load(open("ctf_puzzles.json"))}
    slug, Ls = sys.argv[1], [int(x) for x in sys.argv[2].split(",")]
    res = []
    for L in Ls:
        for an, X in (("KA", k4.KA), ("AZ", k4.AZ)):
            for var in ("vig", "beau"):
                _, k0, _ = SP.solve(P[slug], L, X, var)
                for seed in range(3):
                    res.append((*run(P[slug], L, X, var, seed=seed, init=[X.index(c) for c in k0] if seed == 0 else None), an, var, L))
    for r in sorted(res, reverse=True)[:2]: print(f"{r[0]:.2f} {r[3]} {r[4]} L={r[5]} key={r[1]}\n   {r[2]}")
