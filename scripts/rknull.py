"""Null model: same scan on shuffled K4 (meaningless ciphertext, same letters). Best scores = noise ceiling."""
import itertools, multiprocessing as mp, numpy as np, random, sys, rkscan, rkwide, k4
def job(args):
    seed, lang, w, prefix = args
    ct = list(k4.CT); random.Random(seed).shuffle(ct); ct = "".join(ct)
    rkscan.VARS[:] = []
    for an, X in k4.ALPHS.items():
        toAZ = np.array([rkscan.AZi[c] for c in X]); cta = np.array([X.index(c) for c in ct]); pi = np.array([X.index(c) for c in k4.CPL])
        for fn in ("vig", "beau", "vb"): rkscan.VARS.append((an, fn, toAZ, cta, pi))
    return rkwide.job((lang, w, prefix))
if __name__ == "__main__":
    for seed in (1, 2, 3):
        jobs = [(seed, lang, 11, pre) for lang in ("en", "de") for pre in itertools.permutations(range(11), 2)]
        with mp.Pool(10) as pool: res = pool.map(job, jobs)
        allk = [x for r in res for x in r]
        for lang in ("en", "de"):
            print(f"null seed{seed} w11 {lang} best %.2f" % max(x[1] for x in allk if x[0] == lang), flush=True)
