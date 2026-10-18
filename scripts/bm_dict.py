import multiprocessing as mp, k4, rk_mixed, bm
from collections import Counter
CR = list(zip(k4.CPOS.tolist(), k4.CPL))
def job(chunk):
    dist = Counter(); low = []
    for a, name in chunk:
        for pn, PA in (("AZ", k4.AZ), ("KA", k4.KA), ("A", a)):
            for cn, CA in (("AZ", k4.AZ), ("KA", k4.KA), ("A", a)):
                if "A" not in (pn, cn): continue
                for kn, sg in (("vig", -1), ("beau", 1)):
                    k = [(CA.index(k4.CT[i]) + sg * PA.index(p)) % 26 for i, p in CR]
                    (l2, j2), (l13, j13) = bm.lc(k); dist[l13] += 1
                    if l13 <= 5 or (j2 and j13): low.append((l13, l2, j2 and j13, name, pn, cn, kn))
    return dist, low
if __name__ == "__main__":
    AL = rk_mixed.alph_list(); chunks = [AL[i:i + 5000] for i in range(0, len(AL), 5000)]
    tot = Counter(); low = []
    with mp.Pool(10) as pool:
        for d, l in pool.imap_unordered(job, chunks): tot.update(d); low += l
    print("sequences:", sum(tot.values()), "mod-13 complexity distribution:", sorted(tot.items()))
    print("lowest / joint-rule cases:", sorted(low)[:15])
