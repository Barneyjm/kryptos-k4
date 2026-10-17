import itertools, multiprocessing as mp, numpy as np, time, rkscan, quad
def job(args):
    lang, w, prefix = args
    if lang == "de":
        import german; rkscan.Q = german.QDE
    else:
        rkscan.Q = quad.Q
    keep = []; rest = [c for c in range(w) if c not in prefix]; batch = []
    def flush():
        orders = np.array([rkscan.k4.order_from(w, p) for p in batch]); dests = np.empty_like(orders)
        dests[np.arange(len(batch))[:, None], orders] = np.arange(97)
        for d, M in (("enc", dests), ("dec", orders)):
            rkscan.scan(M, lambda i, b=list(batch), d=d: f"w{w} {d} {b[i]}", keep)
        batch.clear()
    for p in itertools.permutations(rest):
        batch.append(prefix + p)
        if len(batch) == 20000: flush()
    if batch: flush()
    keep.sort(reverse=True); return [(lang,) + x for x in keep[:5]]
if __name__ == "__main__":
    t = time.time()
    for w in (10, 11):
        jobs = [(lang, w, pre) for lang in ("en", "de") for pre in itertools.permutations(range(w), 2 if w == 11 else 1)]
        with mp.Pool(10) as pool: res = pool.map(job, jobs)
        allk = sorted([x for r in res for x in r], key=lambda x: -x[1])
        for lang in ("en", "de"):
            print(f"w{w} {lang} best:", flush=True)
            for x in [y for y in allk if y[0] == lang][:5]: print("  %.2f" % x[1], x[2:], flush=True)
        print("elapsed %.0fs" % (time.time() - t), flush=True)
