import re, multiprocessing as mp, k4
from linear import rank_consistent
from gromark_dict import straight, transposed
CR = list(zip(k4.CPOS.tolist(), k4.CPL))
def system(W, off):
    ru = sorted({(i + off) // W for i, _ in CR}); RI = {r: j for j, r in enumerate(ru)}; M = []
    for i, _ in CR:
        v = [0] * (W + len(ru)); v[(i + off) % W] = 1; v[W + RI[(i + off) // W]] = 1; M.append(v)
    return M
SYS7 = [system(7, off) for off in range(7)]
def test(args):
    name, PA, CA = args; out = []
    for kn, f in (("vig", lambda c, p: (c - p) % 26), ("beau", lambda c, p: (c + p) % 26)):
        k = [f(CA.index(k4.CT[i]), PA.index(p)) for i, p in CR]
        for off, M in enumerate(SYS7):
            a, r2 = rank_consistent(M, k, 2)
            if not a: continue
            b, r13 = rank_consistent(M, k, 13)
            if b and 24 - max(r2, r13) >= 5: out.append((24 - max(r2, r13), name, kn, off))
    return out
def ok(A):
    g = lambda a, b: (A.index(b) - A.index(a)) % 26
    return g("N", "O") == g("S", "T") == g("I", "N")
if __name__ == "__main__":
    words = sorted({re.sub("[^A-Z]", "", w.upper()) for w in open("/usr/share/dict/words")} - {""}) + ["KRYPTOS", "PALIMPSEST", "ABSCISSA", "BERLINCLOCK", "WELTZEITUHR"]
    cands = {}
    for w in words:
        for form, A in (("straight", straight(w)), ("transposed", transposed(w)), ("straight-rev", straight(w)[::-1]), ("transposed-rev", transposed(w)[::-1])):
            if ok(A): cands.setdefault(A, f"{w}/{form}")
    print(f"keyword alphabets with equal N->O, S->T, I->N steps: {len(cands):,}", flush=True)
    jobs = [(f"{n} | cipher {cn}", A, C) for A, n in cands.items() for cn, C in (("AZ", k4.AZ), ("KA", k4.KA), ("same", A))]
    with mp.Pool(11) as pool: res = pool.map(test, jobs, chunksize=300)
    hits = sorted([h for r in res for h in r], reverse=True)
    print("width-7 grid-key fits:", len(hits)); [print("  ", h) for h in hits[:15]]
