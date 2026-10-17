"""Running key (English or German text) with keyword-mixed alphabets, no transposition.
For alphabet A (dictionary keyword, straight/transposed, fwd/rev) and each Quagmire arrangement (PA, CA, indicator)
in {AZ, A}^3, vig/beau: implied key letters at the 24 crib positions; score as English/German (13 + 11 letter runs)."""
import numpy as np, re, sys, k4, multiprocessing as mp
from quad import Q as QEN
from gromark_dict import straight, transposed
POS = k4.CPOS; PL = np.array([ord(c) - 65 for c in k4.CPL])
def frag_score(Kl, Q):   # Kl (n,24) letter indices 0..25
    def qs(a):
        i = ((a[:, :-3]*26 + a[:, 1:-2])*26 + a[:, 2:-1])*26 + a[:, 3:]; return Q[i].sum(axis=1)
    return (qs(Kl[:, :13]) + qs(Kl[:, 13:])) / 18
def run(ct, alphs, Q, top=5):
    C = np.array([ord(c) - 65 for c in ct])[POS]
    A = np.array([[ord(ch) - 65 for ch in a] for a in alphs])             # (n,26) letter at slot
    INV = np.argsort(A, axis=1)                                            # letter -> slot
    n = len(A); AZ = np.tile(np.arange(26), (n, 1)); AZinv = AZ
    best = []
    for pn, Pinv in (("AZ", AZinv), ("A", INV)):
        for cn, Cinv in (("AZ", AZinv), ("A", INV)):
            if pn == "AZ" and cn == "AZ": continue
            for kn in ("vig", "beau"):
                pc = np.take_along_axis(Pinv, np.tile(PL, (n, 1)), 1); cc = np.take_along_axis(Cinv, np.tile(C, (n, 1)), 1)
                k = (cc - pc) % 26 if kn == "vig" else (cc + pc) % 26
                for iname, IND in (("AZ", AZ), ("A", A)):
                    Kl = np.take_along_axis(IND, k, 1); s = frag_score(Kl, Q)
                    for j in np.argsort(s)[-top:]:
                        best.append((float(s[j]), j, pn, cn, kn, iname, "".join(chr(65 + x) for x in Kl[j])))
    best.sort(reverse=True); return best[:top]
def alph_list():
    words = sorted({re.sub("[^A-Z]", "", w.upper()) for w in open("/usr/share/dict/words")} - {""})
    words += ["KRYPTOS", "PALIMPSEST", "ABSCISSA", "BERLINCLOCK", "WELTZEITUHR", "SHADOW", "URANIA"]
    al = {}
    for w in words:
        for f, a in (("s", straight(w)), ("t", transposed(w))):
            al.setdefault(a, f"{w}/{f}"); al.setdefault(a[::-1], f"{w}/{f}-rev")
    return list(al.items())
def job(args):
    lang, ct, chunk = args
    if lang == "de":
        import german; Q = german.QDE
    else: Q = QEN
    res = run(ct, [a for a, _ in chunk], Q)
    return [(s, chunk[j][1], pn, cn, kn, iname, key) for s, j, pn, cn, kn, iname, key in res]
if __name__ == "__main__":
    AL = alph_list(); print(f"distinct alphabets: {len(AL):,}", flush=True)
    rng = np.random.default_rng(0)
    cts = [("K4", k4.CT)] + [(f"null{t}", "".join(rng.permutation(list(k4.CT)))) for t in range(3)]
    chunks = [AL[i:i + 20000] for i in range(0, len(AL), 20000)]
    with mp.Pool(11) as pool:
        for name, ct in cts:
            for lang in ("en", "de"):
                res = sorted([x for r in pool.map(job, [(lang, ct, ch) for ch in chunks]) for x in r], reverse=True)
                print(f"{name:6s} {lang}: best {res[0][0]:.2f}  {res[0][1:]}", flush=True)
                if name == "K4":
                    for x in res[1:4]: print(f"         {x[0]:.2f}  {x[1:]}", flush=True)
