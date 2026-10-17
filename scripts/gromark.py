import itertools, numpy as np, k4
from texts import KEYWORDS
DEKW = ["WELTZEITUHR","URANIA","ALEXANDERPLATZ","MAUER","BERLINERMAUER","UHR","ZEIT","WELTUHR","OSTBERLIN","NORDOSTEN","SCHATTEN","GEHEIM"]
keyed = lambda w: "".join(dict.fromkeys(w + k4.AZ))
def alphabets():
    A = {"AZ": k4.AZ, "KA": k4.KA}
    for w in set(KEYWORDS + DEKW + ["SHADOW", "TENJ"]): A["kw:" + w] = keyed(w)
    for k in list(A): A[k + "(rev)"] = A[k][::-1]
    return A
BLOCKS = [list(range(21, 34)), list(range(63, 74))]
def digits(ct, cr, PA, CA, sign):
    return {i: (sign * (CA.index(ct[i]) - PA.index(cr[i]))) % 26 for i in cr}
def search(ct, cr, ALPH, general_for=("AZ", "KA", "AZ(rev)", "KA(rev)")):
    hits = []
    for (pn, PA), (cn, CA) in itertools.product(ALPH.items(), repeat=2):
        for sign in (1, -1):
            d = digits(ct, cr, PA, CA, sign); mx = max(d.values())
            for n in range(2, 8):
                eqs = [(d[i], [d[i - j] for j in range(1, n + 1)]) for blk in BLOCKS for i in blk if i - n >= blk[0]]
                tgt = np.array([e[0] for e in eqs]); prev = np.array([e[1] for e in eqs])
                for b in range(max(mx + 1, 2), 27):
                    for p, q in itertools.combinations(range(1, n + 1), 2):
                        if np.all((prev[:, p-1] + prev[:, q-1]) % b == tgt): hits.append((len(eqs), pn, cn, sign, b, n, f"d[i-{p}]+d[i-{q}]"))
                    if n <= 3 and pn in general_for and cn in general_for:
                        for coef in itertools.product(range(b), repeat=n):
                            if np.all((prev @ np.array(coef)) % b == tgt): hits.append((len(eqs), pn, cn, sign, b, n, f"coef{coef}"))
    return hits
if __name__ == "__main__":
    import random, sys
    ALPH = alphabets()
    if "test" in sys.argv:   # plant: base 10 Gromark, primer 26717, plain AZ, cipher keyed KRYPTOS
        random.seed(3); P = [random.choice(k4.AZ) for _ in range(97)]
        for s, w in k4.CRIBS.items(): P[s:s+len(w)] = w
        d = [2, 6, 7, 1, 7]
        while len(d) < 97: d.append((d[-5] + d[-4]) % 10)
        CA = keyed("KRYPTOS"); ct = "".join(CA[(k4.AZ.index(p) + x) % 26] for p, x in zip(P, d))
        cr = dict(zip(k4.CPOS.tolist(), k4.CPL))
        h = search(ct, cr, {"AZ": k4.AZ, "kw:KRYPTOS": CA, "KA": k4.KA}); print("planted found:", [x for x in h if x[4] == 10][:3])
    else:
        h = search(k4.CT, dict(zip(k4.CPOS.tolist(), k4.CPL)), ALPH)
        print(f"alphabets {len(ALPH)}; consistent rules: {len(h)}"); [print(" ", x) for x in sorted(h, reverse=True)[:20]]
