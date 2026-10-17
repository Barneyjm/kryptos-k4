"""Fingerprint battery: statistics of K4 vs simulated ciphertexts from candidate hand methods."""
import numpy as np, re, glob, random, k4
from collections import Counter
rng = np.random.default_rng(7); random.seed(7)
ENG = re.sub("[^A-Z]", "", " ".join(open(f, encoding="utf-8", errors="ignore").read().upper() for f in glob.glob("g*.txt")))
KA = k4.KA; KAi = {c: i for i, c in enumerate(KA)}; AZ = k4.AZ
N = 97
def eng(n=N):
    s = random.randrange(10000, len(ENG) - n - 10000); return ENG[s:s+n]
def vig(p, key, X=KA):
    return "".join(X[(X.index(c) + X.index(key[i % len(key)])) % 26] for i, c in enumerate(p))
def beau(p, key, X=KA):
    return "".join(X[(X.index(key[i % len(key)]) - X.index(c)) % 26] for i, c in enumerate(p))
def rkey(L): return "".join(random.choice(AZ) for _ in range(L))
def col(t, w):
    order = list(range(w)); random.shuffle(order); return "".join(t[i] for c in order for i in range(c, len(t), w))
def ptak(p, primer, X=KA):
    k = primer + p; return "".join(X[(X.index(c) + X.index(k[i])) % 26] for i, c in enumerate(p))
MODELS = {
 "random letters":           lambda: "".join(random.choice(AZ) for _ in range(N)),
 "transposition only":       lambda: col(eng(), random.randint(5, 12)),
 "Vig KA period 4-7":        lambda: vig(eng(), rkey(random.randint(4, 7))),
 "Vig KA period 8-14":       lambda: vig(eng(), rkey(random.randint(8, 14))),
 "Beau AZ period 8-14":      lambda: beau(eng(), rkey(random.randint(8, 14)), AZ),
 "Vig KA running key":       lambda: vig(eng(), eng()),
 "Vig KA plaintext autokey": lambda: ptak(eng(), rkey(random.randint(3, 8))),
 "Vig KA p8-14 + columnar":  lambda: col(vig(eng(), rkey(random.randint(8, 14))), random.randint(5, 12)),
 "columnar + Vig KA p8-14":  lambda: vig(col(eng(), random.randint(5, 12)), rkey(random.randint(8, 14))),
 "running key + columnar":   lambda: col(vig(eng(), eng()), random.randint(5, 12)),
}
def stats(c):
    a = np.array([ord(x) - 65 for x in c]); ka = np.array([KAi[x] for x in c]); n = len(a)
    cnt = np.bincount(a, minlength=26)
    s = {}
    s["IC"] = (cnt * (cnt - 1)).sum() / (n * (n - 1))
    s["max letter count"] = cnt.max()
    s["letters absent"] = (cnt == 0).sum()
    s["doubled letters"] = int((a[1:] == a[:-1]).sum())
    s["repeated bigrams"] = sum(v - 1 for v in Counter(c[i:i+2] for i in range(n - 1)).values() if v > 1)
    s["repeated trigrams"] = sum(v - 1 for v in Counter(c[i:i+3] for i in range(n - 2)).values() if v > 1)
    kap = [(a[sh:] == a[:-sh]).mean() for sh in range(1, 49)]
    s["max kappa (shift 1-48)"] = max(kap)
    s["best period IC (2-24)"] = max(np.mean([((lambda cc: (cc*(cc-1)).sum()/max(len(col_)*(len(col_)-1),1))(np.bincount(col_, minlength=26))) for col_ in [a[r::L] for r in range(L)]]) for L in range(2, 25))
    dAZ = (a[1:] - a[:-1]) % 26; dKA = (ka[1:] - ka[:-1]) % 26
    s["IC of AZ differences"] = (lambda cc: (cc*(cc-1)).sum()/(len(dAZ)*(len(dAZ)-1)))(np.bincount(dAZ, minlength=26))
    s["IC of KA differences"] = (lambda cc: (cc*(cc-1)).sum()/(len(dKA)*(len(dKA)-1)))(np.bincount(dKA, minlength=26))
    s["KA-adjacent pairs (diff +-1)"] = int(((dKA == 1) | (dKA == 25)).sum())
    s["chi2 vs uniform"] = ((cnt - n/26)**2 / (n/26)).sum()
    return s
K = stats(k4.CT)
SIMS = 3000
res = {}
for m, f in MODELS.items():
    sims = [stats(f()) for _ in range(SIMS)]
    res[m] = {k: np.array([s[k] for s in sims]) for k in K}
print(f"{'statistic':32s} {'K4':>7s} | " + " | ".join(f"{m[:14]:>14s}" for m in MODELS))
for k, v in K.items():
    row = []
    for m in MODELS:
        x = res[m][k]; lo = (x <= v).mean(); hi = (x >= v).mean(); p = min(1, 2 * min(lo, hi))
        row.append(f"{p:14.3f}")
    print(f"{k:32s} {v:7.3f} | " + " | ".join(row))
print("\n(cells = two-sided p-value of K4's value under that model; < 0.01 means K4 is unusual for it)")
