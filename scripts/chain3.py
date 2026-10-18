import numpy as np, itertools, time, multiprocessing as mp, k4, gromark2 as G
from gromark_sweep import REP
W = 21
def expand3(prim, b, lags, upto=97):
    m, n = prim.shape; d = np.zeros((m, upto), np.int64); d[:, :n] = prim
    for i in range(n, upto): d[:, i] = (d[:, i - lags[0]] + d[:, i - lags[1]] + d[:, i - lags[2]]) % b
    return d
def job(args):
    b, n, lags = args; total = b ** n; out = []
    for start in range(0, total, 400_000):
        idx = np.arange(start, min(total, start + 400_000))
        prim = np.stack([(idx // b ** (n - 1 - t)) % b for t in range(n)], axis=1)
        K = expand3(prim, b, lags); Kc = K[:, G.POS]
        ok = np.all((Kc @ G.CYC) % 26 == 0, axis=1) & np.all((Kc @ G.DIST) % 26 != 0, axis=1)
        for r in np.nonzero(ok)[0]:
            k = K[r]; out.append((sum(1 for i, j in REP if k[i] == k[j] and k[i + W] == k[j + W]), b, n, lags, tuple(prim[r].tolist())))
    return out
if __name__ == "__main__":
    # planted check: 3-term rule recovers
    jobs = [(b, n, l) for b in range(3, 27) for n in range(3, 9) if b ** n <= 1.2e7 for l in itertools.combinations(range(1, n + 1), 3)]
    jobs.sort(key=lambda j: j[0] ** j[1]); print("configs", len(jobs), "primers {:,}".format(sum(b ** n for b, n, _ in jobs)), flush=True)
    t = time.time(); allr = []
    with mp.Pool(7) as pool:
        for i, r in enumerate(pool.imap_unordered(job, jobs)):
            allr += r
            if i % 100 == 0: print(f"{i}/{len(jobs)} {len(allr):,} survivors {time.time()-t:.0f}s", flush=True)
    allr.sort(reverse=True); from collections import Counter
    print(f"DONE {time.time()-t:.0f}s survivors {len(allr):,}; width-21 explained distribution {sorted(Counter(x[0] for x in allr).items())}", flush=True)
    for x in allr[:10]: print("  ", x, flush=True)
