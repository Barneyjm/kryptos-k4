"""Role (iii): PA=A1=D1 (clue alphabet from round3.ALPHA), B1=A2=D (dictionary). Reversed alphabets are exactly
equivalent (k -> -k [+const]) so '-rev' dictionary entries and reversed D1 duplicates are skipped."""
import time, numpy as np, multiprocessing as mp, k4, rk_mixed, ts2_vec as V
from collections import Counter
from round3 import ALPHA
from ts2_run import sweep, report
def main():
    out = open("ts2_results.out", "a")
    def log(s): print(s, flush=True); out.write(s + "\n"); out.flush()
    AL = [(a, n) for a, n in rk_mixed.alph_list() if not n.endswith("-rev")]
    strs = [a for a, _ in AL]; names = [n for _, n in AL]
    seen = {k4.AZ, k4.KA, k4.AZ[::-1], k4.KA[::-1]}; D1s = []
    for n, a in ALPHA.items():
        if len(a) == 26 and len(set(a)) == 26 and a not in seen: D1s.append((n, a)); seen |= {a, a[::-1]}
    log(f"=== role iii {time.ctime()}: {len(D1s)} clue alphabets D1 (excl. AZ/KA, dedup incl. reversal) x {len(strs):,} dictionary alphabets (no -rev)")
    rng = np.random.default_rng(0); sh = "".join(rng.permutation(list(k4.CT)))
    agg = {"K4": [], "shuffled": []}; ncf = Counter(); t0 = time.time()
    with mp.Pool(4) as pool:
        for i, (d1n, d1) in enumerate(D1s):
            for ctn, ct in (("K4", k4.CT), ("shuffled", sh)):
                hits, ncfg = sweep(pool, ct, "i", d1, strs, names); ncf[ctn] += ncfg
                agg[ctn] += [(h, d1n, d1) for h in hits]
                for h in hits:
                    if h[0] >= 5: log(f"   [{ctn}] D1={d1n} net={h[0]} sols={h[1]} D={h[2]} {h[3]} L=({h[4]},{h[5]}) s=({h[6]},{h[7]})")
                    if ctn == "K4" and h[0] >= 6: report(ct, h, "i", d1, log)
            log(f"  [{i+1}/{len(D1s)}] D1={d1n}: K4 hits so far {len(agg['K4'])}, shuffled {len(agg['shuffled'])} ({time.time()-t0:.0f}s)")
    for c in ("K4", "shuffled"):
        log(f"[role iii total {c}] {ncf[c]:,} alphabet-configs; hits {len(agg[c])}; by net {dict(sorted(Counter(h[0] for h, _, _ in agg[c]).items()))}")
if __name__ == "__main__":
    main()
