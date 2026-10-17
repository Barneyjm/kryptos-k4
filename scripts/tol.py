"""Error-tolerant crib search: accept 1-2 crib letters that disagree (hand-encryption / carving slip)."""
import itertools, sys, numpy as np, multiprocessing as mp, k4, round2
from texts import KEYWORDS
NC = len(k4.CPOS)
def tol_check(D, ct, cpos, cpl, label, out, alphs, maxerr=2, mineff=9):
    D = np.atleast_2d(D); dc = D[:, cpos]; n = len(D); rows = np.arange(n)[:, None]
    for an, (PA, CA) in alphs.items():
        cta = np.array([CA.index(c) for c in ct]); ci = cta[dc]; pi = np.array([PA.index(c) for c in cpl])[None, :]
        for kn, k in (("vig", (ci - pi) % 26), ("beau", (ci + pi) % 26)):
            for model in "AB":
                for L in range(1, 27):
                    r = np.broadcast_to(cpos[None, :] % L, k.shape) if model == "A" else dc % L
                    cnt = np.zeros((n, L * 26), np.int16); np.add.at(cnt, (np.broadcast_to(rows, k.shape), r * 26 + k), 1)
                    cnt = cnt.reshape(n, L, 26)
                    matches = cnt.max(axis=2).sum(axis=1); mism = NC - matches
                    distinct = (cnt.sum(axis=2) > 0).sum(axis=1); cons = NC - distinct
                    ok = (mism >= 1) & (mism <= maxerr) & (cons - mism >= mineff)
                    for i in np.nonzero(ok)[0]:
                        out.append((label(i), an, kn, model, L, int(cons[i]), int(mism[i])))
STD = {"AZ": (k4.AZ, k4.AZ), "KA": (k4.KA, k4.KA)}
def worker(args):
    w, first = args; hits = []; rest = [c for c in range(w) if c != first]; batch = []
    def flush():
        orders = np.array([k4.order_from(w, p) for p in batch]); dests = np.empty_like(orders)
        dests[np.arange(len(batch))[:, None], orders] = np.arange(97)
        for d, M in (("enc", dests), ("dec", orders)):
            tol_check(M, k4.CT, k4.CPOS, k4.CPL, lambda i, b=list(batch), d=d: f"w{w} {d} {b[i]}", hits, STD)
        batch.clear()
    for p in itertools.permutations(rest):
        batch.append((first,) + p)
        if len(batch) == 5000: flush()
    if batch: flush()
    return hits
if __name__ == "__main__":
    # 1. no transposition, all keyword alphabet pairs (Quagmire family)
    from round3 import ALPHA
    pairs = {f"{p}|{c}": (ALPHA[p], ALPHA[c]) for p in ALPHA for c in ALPHA}
    h0 = []; tol_check(np.arange(97), k4.CT, k4.CPOS, k4.CPL, lambda i: "none", h0, pairs)
    print("no-transposition, 1-2 errors:", len(h0)); [print(x) for x in h0[:15]]
    W = int(sys.argv[1]) if len(sys.argv) > 1 else 8
    with mp.Pool(8) as pool: res = pool.map(worker, [(w, f) for w in range(2, W + 1) for f in range(w)])
    h = [x for r in res for x in r]; print(f"columnar w<= {W}, 1-2 errors:", len(h)); [print(x) for x in h[:30]]
