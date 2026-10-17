import itertools, numpy as np, multiprocessing as mp, k4, tol, round2
from texts import KEYWORDS
from runkey import SOURCES
def dbl():
    ords = {}
    for kw in KEYWORDS:
        ords[kw] = np.array(k4.order_from(len(kw), round2.kw_order(kw, len(kw))))
        for w in (7, 8, 10, 11, 13, 14, 24): ords[f"{kw}@{w}"] = np.array(k4.order_from(w, round2.kw_order(kw, w)))
    names = list(ords); hits = []; B, Lb = [], []
    for a, b in itertools.product(names, repeat=2):
        o = ords[a][ords[b]]
        for d, arr in (("enc", round2.dest_from_order(o)), ("dec", o)): B.append(arr); Lb.append(f"dbl {a}|{b} {d}")
        if len(B) >= 5000:
            tol.tol_check(np.array(B), k4.CT, k4.CPOS, k4.CPL, lambda i, L=list(Lb): L[i], hits, tol.STD); B.clear(); Lb.clear()
    if B: tol.tol_check(np.array(B), k4.CT, k4.CPOS, k4.CPL, lambda i, L=list(Lb): L[i], hits, tol.STD)
    return hits
def runkey_tol():
    hits = []
    for an, X in k4.ALPHS.items():
        ci = np.array([X.index(c) for c in k4.CT])[k4.CPOS]; pi = np.array([X.index(c) for c in k4.CPL])
        for kn, need in (("vig", (ci-pi) % 26), ("beau", (ci+pi) % 26), ("vb", (pi-ci) % 26)):
            for sn, s in SOURCES.items():
                kidx = np.array([X.index(c) for c in s])
                for off in range(-21, len(s) - 74):
                    p = off + k4.CPOS
                    if p.min() < 0 or p.max() >= len(s): continue
                    m = int((kidx[p] == need).sum())
                    if m >= 20: hits.append((an, kn, sn, off, m))
    return hits
if __name__ == "__main__":
    print("running key (K texts), >=20/24 match:", runkey_tol(), flush=True)
    h = dbl(); print("double keyword columnar, 1-2 errors:", len(h), h[:20], flush=True)
    with mp.Pool(8) as pool: res = pool.map(tol.worker, [(9, f) for f in range(9)])
    h = [x for r in res for x in r]; print("columnar w9, 1-2 errors:", len(h), h[:20], flush=True)
