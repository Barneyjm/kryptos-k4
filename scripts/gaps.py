import re, numpy as np, multiprocessing as mp, k4, rowcol
from gromark_dict import straight, transposed
def ok(A):
    g = lambda a, b: (A.index(b) - A.index(a)) % 26
    return g("N", "O") == g("S", "T") == g("I", "N")
if __name__ == "__main__":
    words = sorted({re.sub("[^A-Z]", "", w.upper()) for w in open("/usr/share/dict/words")} - {""})
    words += ["KRYPTOS", "PALIMPSEST", "ABSCISSA", "BERLINCLOCK", "WELTZEITUHR", "SHADOW", "BERLIN", "CLOCK"]
    cands = {}
    for w in words:
        for form, A in (("straight", straight(w)), ("transposed", transposed(w)), ("straight-rev", straight(w)[::-1]), ("transposed-rev", transposed(w)[::-1])):
            if ok(A): cands.setdefault(A, f"{w}/{form}")
    print(f"distinct keyword alphabets where N->O, S->T, I->N are equal steps: {len(cands):,}", flush=True)
    rowcol.ALPHA.clear(); rowcol.ALPHA.update({"AZ": k4.AZ, "KA": k4.KA}); rowcol.ALPHA.update({v: k for k, v in cands.items()})
    rowcol.ALPHA = {**{"AZ": k4.AZ, "KA": k4.KA}, **{name: A for A, name in cands.items()}}
    jobs = [(name, c, [7]) for name in list(rowcol.ALPHA)[2:] for c in ("AZ", "KA", name)]
    with mp.Pool(11) as pool: res = pool.map(rowcol.test, jobs, chunksize=200)
    hits = sorted([h for r in res for h in r], reverse=True)
    print("grid-key (width 7) fits using those alphabets:", len(hits)); [print("  ", h) for h in hits[:15]]
