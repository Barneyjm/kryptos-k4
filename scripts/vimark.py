"""Chain addition (Vimark / VIC-style) from keyword primers: k_i = k_{i-p} + k_{i-q} mod 26, primer = word letters
(A-Z or KRYPTOS numbering). Both alphabets free (Bean's graph check). Survivors ranked by width-21 repeats explained."""
import numpy as np, itertools, re, multiprocessing as mp, k4, gromark2 as G
from gromark_sweep import REP
W = 21
def job(args):
    n, words, numbering = args; X = k4.AZ if numbering == "AZ" else k4.KA
    prim = np.array([[X.index(c) for c in w] for w in words]); out = []
    for p, q in itertools.combinations(range(1, n + 1), 2):
        K = G.expand(prim, 26, (p, q), 97); Kc = K[:, G.POS]
        ok = np.all((Kc @ G.CYC) % 26 == 0, axis=1) & np.all((Kc @ G.DIST) % 26 != 0, axis=1)
        for r in np.nonzero(ok)[0]:
            k = K[r]; ex = sum(1 for i, j in REP if k[i] == k[j] and k[i + W] == k[j + W])
            out.append((ex, words[r], numbering, (p, q)))
    return out
if __name__ == "__main__":
    words = sorted({re.sub("[^A-Z]", "", w.upper()) for w in open("/usr/share/dict/words")} - {""})
    words += ["KRYPTOS", "PALIMPSEST", "ABSCISSA", "BERLINCLOCK", "WELTZEITUHR", "SHADOW", "URANIA", "BERLIN", "CLOCK",
              "EASTNORTHEAST", "LAYERTWO", "IQLUSION", "UNDERGRUUND", "DESPARATLY", "QUAE", "DYAHR", "SCHEIDT", "SANBORN"]
    bylen = {}
    for w in set(words):
        if 3 <= len(w) <= 14: bylen.setdefault(len(w), []).append(w)
    jobs = [(n, ws[i:i + 20000], num) for n, ws in bylen.items() for i in range(0, len(ws), 20000) for num in ("AZ", "KA")]
    with mp.Pool(11) as pool: res = [x for r in pool.map(job, jobs) for x in r]
    res.sort(reverse=True)
    print(f"primer words: {sum(len(v) for v in bylen.values()):,}; survivors of the 24 known letters: {len(res):,}")
    from collections import Counter; print("width-21 pairs explained, distribution:", sorted(Counter(r[0] for r in res).items()))
    for r in res[:12]: print("  ", r)
