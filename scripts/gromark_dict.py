"""Gromark with keyword-mixed cipher alphabets from a full dictionary (straight + columnar-transposed forms).
Filter: all 24 implied key digits < base; then check two-lag recurrence (any lags up to 7)."""
import numpy as np, itertools, re, sys, k4, multiprocessing as mp
AZ = k4.AZ; CR = dict(zip(k4.CPOS.tolist(), k4.CPL)); POS = k4.CPOS.tolist()
BLOCKS = [list(range(21, 34)), list(range(63, 74))]
def straight(w): return "".join(dict.fromkeys(w + AZ))
def transposed(w):
    k = "".join(dict.fromkeys(w)); s = straight(w); n = len(k)
    cols = sorted(range(n), key=lambda c: k[c])
    return "".join(s[c::n] for c in cols)
def check(word):
    out = []
    for form, CA in (("straight", straight(word)), ("transposed", transposed(word))):
        for PA_name, PA in (("AZ", AZ), ("KA", k4.KA)):
            for sign in (1, -1):
                d = {i: (sign * (CA.index(k4.CT[i]) - PA.index(CR[i]))) % 26 for i in POS}
                mx = max(d.values())
                if mx >= 20: continue                      # only interesting if digits look like a small base
                for b in range(mx + 1, 21):
                    for n in range(2, 8):
                        eqs = [(d[i], [d[i - j] for j in range(1, n + 1)]) for blk in BLOCKS for i in blk if i - n >= blk[0]]
                        for p, q in itertools.combinations(range(1, n + 1), 2):
                            if all((pr[p-1] + pr[q-1]) % b == t for t, pr in eqs):
                                out.append((word, form, PA_name, sign, b, n, p, q, len(eqs)))
                out.append(("SMALLDIGITS", word, form, PA_name, sign, mx))
    return out
if __name__ == "__main__":
    words = sorted({re.sub("[^A-Z]", "", w.upper()) for w in open("/usr/share/dict/words")} - {""})
    words += ["KRYPTOS", "PALIMPSEST", "ABSCISSA", "BERLINCLOCK", "WELTZEITUHR", "SHADOW", "CLOCK", "BERLIN"]
    with mp.Pool(10) as pool: res = pool.map(check, words, chunksize=2000)
    allr = [x for r in res for x in r]
    small = [x for x in allr if x[0] == "SMALLDIGITS"]; rec = [x for x in allr if x[0] != "SMALLDIGITS"]
    print(f"keywords {len(words):,}; alphabets where all 24 digits < 20: {len(small)}")
    for x in sorted(small, key=lambda x: x[-1])[:15]: print("  small:", x)
    print("recurrence matches:", len(rec)); [print("  ", x) for x in rec[:20]]
