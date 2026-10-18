import sys, time, itertools, numpy as np, multiprocessing as mp, k4, rk_mixed, twostage, ts2_vec as V, ts2_fast2 as FB
from collections import Counter
CH = 50000; LMAX = 12; NETMIN = 4

def sweep(pool, ct, role, X, strs, names):
    tasks = [(ct, role, X, strs[s:s + CH], names[s:s + CH], LMAX, NETMIN) for s in range(0, len(strs), CH)]
    hits = []; ncfg = 0
    for o, nc in pool.imap_unordered(FB.run_chunk, tasks): hits += o; ncfg += nc
    return hits, ncfg

def roles_of(role, X, D):
    return (X, X, D, D) if role == "i" else (D, D, X, X)

def report(ct, h, role, X, log):
    net, sols, nm, D, L1, L2, s1, s2 = h
    PA, A1, B1, A2 = roles_of(role, X, D)
    comps = V.assignments(ct, L1, L2, PA, A1, B1, A2, s1, s2)
    det = [c for c in comps if len(c) < 26]; und = [c for c in comps if len(c) == 26]
    log(f"    components: {len(comps)} ({len(det)} constrained with {[len(c) for c in det]} solutions, {len(und)} free/26)")
    for combo in itertools.islice(itertools.product(*det), 12):
        k1 = {}; k2 = {}
        for asg in combo:
            for (t, r), val in asg.items(): (k1 if t == "1" else k2)[r] = val
        pt = V.decrypt(ct, L1, L2, PA, A1, B1, A2, s1, s2, k1, k2)
        crib_ok = all(pt[p] == c for p, c in zip(k4.CPOS.tolist(), k4.CPL))
        ks = lambda k, L: "".join(chr(65 + k[r]) if r in k else "?" for r in range(L))
        log(f"    k1={ks(k1, L1)} k2={ks(k2, L2)} cribs_ok={crib_ok}\n      PT {pt}")

def main():
    out = open("ts2_results.out", "a")
    def log(s): print(s, flush=True); out.write(s + "\n"); out.flush()
    AL = rk_mixed.alph_list(); strs = [a for a, _ in AL]; names = [n for _, n in AL]
    log(f"=== ts2 run {time.ctime()}  alphabets={len(strs):,}  L1,L2<=12  netmin={NETMIN}")
    pls = [(L1, L2) for L1 in range(1, LMAX + 1) for L2 in range(1, LMAX + 1) if 24 - V.get_plan(L1, L2)[0] >= NETMIN]
    log(f"(L1,L2) pairs with net>=4: {len(pls)} of {LMAX*LMAX}; x4 sign combos = {4*len(pls)} configs per alphabet")
    rng = np.random.default_rng(0); sh = "".join(rng.permutation(list(k4.CT)))
    with mp.Pool(4) as pool:
        # ---- validation ----
        r2 = np.random.default_rng(42); tgt = names.index("HOROLOGE/t"); D = strs[tgt]
        pt = [chr(65 + x) for x in r2.integers(0, 26, 97)]
        for p, c in zip(k4.CPOS.tolist(), k4.CPL): pt[p] = c
        pt = "".join(pt)
        k1 = {i: ord(c) - 65 for i, c in enumerate("PALIMPSEST")}; k2 = {i: ord(c) - 65 for i, c in enumerate("ABSCISSA")}
        pct = V.encrypt(pt, 10, 8, k4.KA, k4.KA, D, D, 1, 1, k1, k2)
        t = time.time(); hits, ncfg = sweep(pool, pct, "i", k4.KA, strs, names)
        log(f"[validate] planted KA/PALIMPSEST(+) -> {names[tgt]}/ABSCISSA(+): {ncfg:,} alphabet-configs, {len(hits)} hits ({time.time()-t:.0f}s)")
        for h in sorted(hits, reverse=True)[:8]: log(f"   net={h[0]} sols={h[1]} {h[2]} L=({h[4]},{h[5]}) s=({h[6]},{h[7]})")
        agree = all(twostage.count(pct, h[4], h[5], k4.KA, k4.KA, h[3], h[3], h[6], h[7]) == (h[1], 24 - h[0]) for h in hits)
        planted = [h for h in hits if h[3] == D and h[4:] == (10, 8, 1, 1)]
        log(f"[validate] planted flagged at (10,8,+,+): {bool(planted)} net={planted[0][0] if planted else None}; twostage.count agrees on all hits: {agree}")
        if planted: report(pct, planted[0], "i", k4.KA, log); log(f"    true PT {pt}")
        # ---- real + shuffled ----
        summary = {}
        for ctn, ct in (("K4", k4.CT), ("shuffled", sh)):
            for role in ("i", "ii"):
                for Xn, X in (("AZ", k4.AZ), ("KA", k4.KA)):
                    t = time.time(); hits, ncfg = sweep(pool, ct, role, X, strs, names)
                    byn = Counter(h[0] for h in hits)
                    summary[(ctn, role, Xn)] = (ncfg, hits)
                    log(f"[{ctn}] role {role} X={Xn}: {ncfg:,} alphabet-configs, significant hits {len(hits)} by net {dict(sorted(byn.items()))} ({time.time()-t:.0f}s)")
                    for h in sorted(hits, reverse=True)[:6]: log(f"   net={h[0]} sols={h[1]} {h[2]} {h[3]} L=({h[4]},{h[5]}) s=({h[6]},{h[7]})")
                    if ctn == "K4":
                        for h in sorted(hits, reverse=True):
                            if h[0] >= 6:
                                log(f"  ** K4 net>=6 hit: {h[:3]} {h[3]} L=({h[4]},{h[5]}) s=({h[6]},{h[7]})"); report(ct, h, role, X, log)
        tot = {c: sum(v[0] for k, v in summary.items() if k[0] == c) for c in ("K4", "shuffled")}
        th = {c: [h for k, v in summary.items() if k[0] == c for h in v[1]] for c in ("K4", "shuffled")}
        for c in ("K4", "shuffled"):
            log(f"[total {c}] {tot[c]:,} alphabet-configs; hits {len(th[c])}; net>=5: {sum(h[0]>=5 for h in th[c])}; net>=6: {sum(h[0]>=6 for h in th[c])}; by net {dict(sorted(Counter(h[0] for h in th[c]).items()))}")

if __name__ == "__main__":
    main()
