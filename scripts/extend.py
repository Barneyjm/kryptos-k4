"""Crib extension: place candidate words at the 4 crib edges, then re-test the masking model for the periods that
were underdetermined. Compare the fit rate against random-letter extensions of the same lengths."""
import re, glob, random, numpy as np, multiprocessing as mp
from collections import Counter
import k4, mask
BASE = dict(zip(k4.CPOS.tolist(), k4.CPL))
UNDET = [13, 16, 18, 19, 20, 23, 24, 25, 26]
EDGES = {"before EAST": ("L", 20), "after NORTHEAST": ("R", 34), "before BERLIN": ("L", 62), "after CLOCK": ("R", 74)}
def candidates():
    corpus = re.findall(r"[A-Z]+", " ".join(open(f, encoding="utf-8", errors="ignore").read().upper() for f in glob.glob("g*.txt")))
    top = [w for w, _ in Counter(corpus).most_common(3000)]
    san = [w for w, c in Counter(re.findall(r"[A-Z]+", open("sanborn_oral_history.txt").read().upper())).most_common(1500)]
    extra = "DEGREES MINUTES SECONDS FEET PACES STEPS MILES METERS YARDS NORTH SOUTH WEST TOWARD TOWARDS FACING POINTS POINTING SHADOW WALL GATE EGYPT CAIRO GIZA PYRAMID SPHINX TOMB DOOR MESSAGE DELIVER DELIVERED TIME HOUR NOON MIDNIGHT ZONE WORLD THEN THERE HERE IS WAS THE OF AT FROM TO LOOK SEE BURIED LOCATION POSITION COMPASS ROSE LODESTONE".split()
    words = sorted({w for w in top + san + extra if 2 <= len(w) <= 10})
    return words
def place(cribs, edge, w):
    side, at = EDGES[edge]; c = dict(cribs)
    rng = range(at - len(w) + 1, at + 1) if side == "L" else range(at, at + len(w))
    for i, ch in zip(rng, w):
        if i in c and c[i] != ch: return None
        c[i] = ch
    return c
def job(args):
    words, tag = args; fits = []; tried = 0
    for w in words:
        for edge in EDGES:
            c = place(BASE, edge, w)
            if c is None: continue
            for side in ("plain", "cipher"):
                for an, X in (("AZ", k4.AZ), ("KA", k4.KA)):
                    for L in UNDET:
                        tried += 1
                        st, sol = mask.test(side, X, L, c)
                        if st == "FITS": fits.append((w, edge, side, an, L))
    return tried, fits
if __name__ == "__main__":
    words = candidates(); print("candidate words:", len(words), flush=True)
    rnd = random.Random(7); fake = ["".join(rnd.choice(k4.AZ) for _ in range(len(w))) for w in words]
    with mp.Pool(11) as pool:
        for tag, ws in (("REAL", words), ("RANDOM", fake)):
            res = pool.map(job, [(ws[i::44], tag) for i in range(44)])
            tried = sum(r[0] for r in res); fits = [f for r in res for f in r[1]]
            print(f"{tag:6s}: {tried:,} tests, {len(fits):,} fits ({len(fits)/max(tried,1):.4%})", flush=True)
            if tag == "REAL":
                from collections import Counter as C
                print("   fits by period:", sorted(C(f[4] for f in fits).items()), flush=True)
                import pickle; pickle.dump(fits, open("extend_fits.pkl", "wb"))
