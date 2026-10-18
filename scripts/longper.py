"""Fix: periodic keys L=17,18,23,24 were suppressed by the >=8-check threshold. Re-test exactly (any checks) across all
alphabet sets. L=25,26: fill the unknown residues and score the full decryption as English vs shuffled baseline."""
import numpy as np, itertools, time, k4, cokernel, gap1, rk_mixed, gap2
from round3 import ALPHA
from quad import Q
POS = k4.CPOS.tolist()
def models_for(Ls):
    out = []
    for L in Ls:
        M = np.zeros((24, L), int)
        for r, p in enumerate(POS): M[r, p % L] = 1
        N2, N13 = cokernel.left_null(M, 2), cokernel.left_null(M, 13); out.append((f"periodic L={L}", N2, N13, len(N13)))
    return out
D = models_for([17, 18, 23, 24]); print("checks per model:", [(d[0], d[3]) for d in D], flush=True)
gap1.P = D; gap1.names = [d[0] for d in D]
gap1.B13 = np.vstack([d[2] for d in D]).astype(np.float32); gap1.B2 = np.vstack([d[1] for d in D]).astype(np.float32)
gap1.seg13 = np.cumsum([0] + [len(d[2]) for d in D])[:-1]; gap1.seg2 = np.cumsum([0] + [len(d[1]) for d in D])[:-1]
AL = rk_mixed.alph_list() + gap2.build()[0] + [(k4.AZ, "AZ"), (k4.KA, "KA")] + [(v, k) for k, v in ALPHA.items()]
t = time.time(); h = gap1.scan(k4.CT, AL, "K4"); print(f"L=17,18,23,24 exact over {len(AL):,} alphabets: {len(h)} hits ({time.time()-t:.0f}s)", flush=True)
for x in h[:20]: print("  ", x, flush=True)
# L=25, 26 with std/KA pairs: fill unknown residues, score English
def qs(a): return Q[((a[:-3]*26+a[1:-2])*26+a[2:-1])*26+a[3:]].mean()
CT = np.array([ord(c) - 65 for c in k4.CT]); CR = dict(zip(POS, k4.CPL)); rng = np.random.default_rng(0)
def best_fill(ct, L, PA, CA, var):
    ci = np.array([CA.index(chr(65 + c)) for c in ct]); forced = {}
    for i, p in CR.items():
        k = (ci[i] - PA.index(p)) % 26 if var == "vig" else (ci[i] + PA.index(p)) % 26
        if forced.setdefault(i % L, k) != k: return None
    free = [r for r in range(L) if r not in forced]; conv = np.array([ord(c) - 65 for c in PA])
    best = (-9, None)
    for vals in itertools.product(range(26), repeat=len(free)):
        key = np.array([forced.get(r, 0) for r in range(L)]); key[free] = vals
        kk = key[np.arange(97) % L]; pt = conv[((ci - kk) if var == "vig" else (kk - ci)) % 26]
        s = qs(pt)
        if s > best[0]: best = (s, "".join(chr(65 + x) for x in pt))
    return best
for L in (26, 25):
    for (pn, PA), (cn, CA) in itertools.product((("AZ", k4.AZ), ("KA", k4.KA)), repeat=2):
        for var in ("vig", "beau"):
            r = best_fill(CT, L, PA, CA, var)
            if r is None: print(f"L={L} {pn}|{cn} {var}: cribs inconsistent"); continue
            nulls = [best_fill(rng.permutation(CT), L, PA, CA, var) for _ in range(2)]
            nulls = [n[0] for n in nulls if n]
            print(f"L={L} {pn}|{cn} {var}: K4 best {r[0]:.2f} vs shuffled {[round(x, 2) for x in nulls]}  {r[1][:60]}", flush=True)
