import sys, time, numpy as np, multiprocessing as mp, prune
def job(args):
    w, d, prefix, lmax = args
    prune.LMAX = lmax; prune.LS = np.arange(1, lmax + 1)
    s = prune.Search(w, d); s.run(prefix); return s.nodes, s.hits
if __name__ == "__main__":
    plan = [(w, 14, ("enc",)) for w in range(12, 26)] + [(w, 14, ("dec",)) for w in (12, 13)] + [(w, 18, ("enc",)) for w in range(8, 15)]
    for w, lmax, dirs in plan:
        t = time.time()
        jobs = [(w, d, (a,), lmax) for d in dirs for a in range(w)]
        with mp.Pool(11) as pool: res = pool.map(job, jobs, chunksize=1)
        nodes = sum(r[0] for r in res); hits = [h for r in res for h in r[1]]
        print(f"{dirs[0]} w={w} periods<= {lmax}: {nodes:,} nodes, hits {len(hits)}, {time.time()-t:.0f}s", flush=True)
        for h in hits[:10]: print("   ", h, flush=True)
