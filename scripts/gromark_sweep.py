import numpy as np, itertools, time, multiprocessing as mp, k4, gromark2 as G
CT = np.array([ord(c) - 65 for c in k4.CT]); W = 21
REP = [(i, j) for i, j in itertools.combinations(range(97 - W), 2) if CT[i] == CT[j] and CT[i + W] == CT[j + W]]
def job(args):
    b, n, lags = args; total = b ** n; out = []
    for start in range(0, total, 500_000):
        idx = np.arange(start, min(total, start + 500_000))
        prim = np.stack([(idx // b ** (n - 1 - t)) % b for t in range(n)], axis=1)
        K = G.expand(prim, b, lags, 97)
        Kc = K[:, G.POS]
        ok = np.all((Kc @ G.CYC) % 26 == 0, axis=1) & np.all((Kc @ G.DIST) % 26 != 0, axis=1)
        for r in np.nonzero(ok)[0]:
            k = K[r]; expl = sum(1 for i, j in REP if k[i] == k[j] and k[i + W] == k[j + W])
            out.append((expl, b, n, lags, tuple(prim[r].tolist())))
    return (b, n, lags, total, out)
if __name__ == "__main__":
    print(f"K4 repeated vertical pairs at width 21: {len(REP)} (pairs of positions)", flush=True)
    jobs = [(b, n, (p, q)) for b in range(3, 27) for n in range(2, 9) if b ** n <= 2e7 for p, q in itertools.combinations(range(1, n + 1), 2)]
    jobs.sort(key=lambda j: j[0] ** j[1])
    print("configs:", len(jobs), "primers total: {:,}".format(sum(b ** n for b, n, _ in jobs)), flush=True)
    t = time.time(); allsurv = []; done = 0
    with mp.Pool(11) as pool:
        for b, n, lags, total, out in pool.imap_unordered(job, jobs):
            done += 1; allsurv += out
            if done % 50 == 0: print(f"{done}/{len(jobs)} configs, {len(allsurv):,} survivors, {time.time()-t:.0f}s", flush=True)
    allsurv.sort(reverse=True)
    print(f"DONE {time.time()-t:.0f}s: {len(allsurv):,} survivors. Best at explaining the width-21 repeats:", flush=True)
    for s in allsurv[:25]: print("  ", s, flush=True)
    import pickle; pickle.dump(allsurv, open("gromark_survivors.pkl", "wb"))
