"""Model A running-key: P -> sub with unknown English running key -> transposition -> C.
Score implied key fragments (positions 21-33, 63-73) as English."""
import itertools, sys, numpy as np, multiprocessing as mp, k4
from quad import Q
AZi = {c: i for i, c in enumerate(k4.AZ)}
VARS = []
for an, X in k4.ALPHS.items():
    toAZ = np.array([AZi[c] for c in X]); cta = np.array([X.index(c) for c in k4.CT]); pi = np.array([X.index(c) for c in k4.CPL])
    for fn in ("vig", "beau", "vb"): VARS.append((an, fn, toAZ, cta, pi))
def frag_score(K):  # K (n,24) AZ indices
    def qs(a):
        i = ((a[:, :-3]*26 + a[:, 1:-2])*26 + a[:, 2:-1])*26 + a[:, 3:]; return Q[i].sum(axis=1)
    return (qs(K[:, :13]) + qs(K[:, 13:])) / 18
def scan(D, labels, keep):
    dc = D[:, k4.CPOS]
    for an, fn, toAZ, cta, pi in VARS:
        ci = cta[dc]
        k = {"vig": (ci - pi) % 26, "beau": (ci + pi) % 26, "vb": (pi - ci) % 26}[fn]
        s = frag_score(toAZ[k])
        for i in np.argsort(s)[-3:]:
            keep.append((float(s[i]), labels(i), an, fn, "".join(k4.AZ[x] for x in toAZ[k[i]])))
def worker(args):
    w, first = args; keep = []; rest = [c for c in range(w) if c != first]; batch = []
    def flush():
        orders = np.array([k4.order_from(w, p) for p in batch]); dests = np.empty_like(orders)
        dests[np.arange(len(batch))[:, None], orders] = np.arange(97)
        for d, M in (("enc", dests), ("dec", orders)):
            scan(M, lambda i, b=list(batch), d=d: f"w{w} {d} {b[i]}", keep)
        batch.clear()
    for p in itertools.permutations(rest):
        batch.append((first,) + p)
        if len(batch) == 20000: flush()
    if batch: flush()
    keep.sort(reverse=True); return keep[:10]
if __name__ == "__main__":
    W = int(sys.argv[1]) if len(sys.argv) > 1 else 9
    base = []; scan(np.arange(97)[None, :], lambda i: "none", base); print("no transposition:", sorted(base)[-3:])
    jobs = [(w, f) for w in range(2, W + 1) for f in range(w)]
    with mp.Pool(6) as pool: res = pool.map(worker, jobs)
    allk = sorted([x for r in res for x in r], reverse=True)
    for x in allk[:25]: print("%.3f" % x[0], x[1:])
