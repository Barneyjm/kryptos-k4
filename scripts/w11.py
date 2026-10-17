import itertools, sys, numpy as np, multiprocessing as mp, k4, round2
W = 11
def worker(prefix):
    hits = []; rest = [c for c in range(W) if c not in prefix]; batch = []
    def flush():
        orders = np.array([k4.order_from(W, p) for p in batch])
        dests = np.empty_like(orders); rows = np.arange(len(batch))[:, None]; dests[rows, orders] = np.arange(97)
        for dname, D in (("enc", dests), ("dec", orders)):
            round2.check_generic(D, k4.CT, k4.CPOS, k4.CPL, lambda i, b=list(batch), d=dname: f"w{W} {d} {b[i]}", hits)
        batch.clear()
    for p in itertools.permutations(rest):
        batch.append(prefix + p)
        if len(batch) == 15000: flush()
    if batch: flush()
    return hits
if __name__ == "__main__":
    prefixes = list(itertools.permutations(range(W), 2))
    with mp.Pool(11) as pool:
        hits = []
        for n, r in enumerate(pool.imap_unordered(worker, prefixes)):
            hits += r
            if n % 10 == 0: print(f"{n+1}/{len(prefixes)} prefixes, hits {len(hits)}", flush=True)
    print("DONE w11 hits", len(hits)); [print(h) for h in hits]
