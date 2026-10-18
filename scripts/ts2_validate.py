import numpy as np, time, k4, rk_mixed, twostage, ts2_vec as V
AL = rk_mixed.alph_list(); strs = [a for a, _ in AL]; names = [n for _, n in AL]
rng = np.random.default_rng(42)
tgt = names.index("HOROLOGE/t") if "HOROLOGE/t" in names else 123457
D = strs[tgt]; print("planted alphabet", names[tgt], D)
pt = [chr(65 + x) for x in rng.integers(0, 26, 97)]
for p, c in zip(k4.CPOS.tolist(), k4.CPL): pt[p] = c
pt = "".join(pt)
k1 = {i: ord(c) - 65 for i, c in enumerate("PALIMPSEST")}; k2 = {i: ord(c) - 65 for i, c in enumerate("ABSCISSA")}
ct = V.encrypt(pt, 10, 8, k4.KA, k4.KA, D, D, 1, 1, k1, k2)
assert V.decrypt(ct, 10, 8, k4.KA, k4.KA, D, D, 1, 1, k1, k2) == pt
print("planted CT", ct)
t = time.time(); hits = []; ncfg = 0
for s in range(0, len(strs), 50000):
    o, nc = V.run_chunk((ct, "i", k4.KA, strs[s:s + 50000], names[s:s + 50000], 12, 4)); hits += o; ncfg += nc
print(f"role i X=KA: {ncfg:,} alphabet-configs, {len(hits)} significant hits, {time.time()-t:.0f}s")
for h in sorted(hits, reverse=True)[:20]: print("  ", h[:3], h[4:])
ok = [h for h in hits if h[3] == D and (h[4], h[5], h[6], h[7]) == (10, 8, 1, 1)]
print("planted flagged at (10,8,+,+):", bool(ok), ok[:1])
# cross-check counts with twostage.count for every hit + a sample of non-hits
for h in hits:
    net, sols, nm, A, L1, L2, s1, s2 = h
    c, nodes = twostage.count(ct, L1, L2, k4.KA, k4.KA, A, A, s1, s2)
    assert c == sols and 24 - nodes == net, (h, c, nodes)
print("twostage.count agrees on all", len(hits), "hits")
# sample of random alphabets/configs: vector says 0 -> scalar says 0 or net<4
idx = rng.choice(len(strs), 300, replace=False)
sub_s = [strs[i] for i in idx]; sub_n = [names[i] for i in idx]
P = V.prep(np.broadcast_to(V.to_arr([k4.KA]), (300, 26)), np.broadcast_to(V.to_arr([k4.KA]), (300, 26)), V.to_arr(sub_s), V.to_arr(sub_s), ct)
bad = 0; checked = 0
for L1, L2 in [(10, 8), (5, 4), (3, 7), (12, 3), (2, 2), (1, 9)]:
    for s1 in (1, -1):
        for s2 in (1, -1):
            sols = V.solve_cfg(P, 300, L1, L2, s1, s2, V.get_plan(L1, L2))
            for a in range(0, 300, 7):
                c, _ = twostage.count(ct, L1, L2, k4.KA, k4.KA, sub_s[a], sub_s[a], s1, s2); checked += 1
                bad += c != sols[a]
print(f"random sample cross-check: {checked} checks, {bad} mismatches")
# also role ii sanity: planted with D in stage 1 and KA in stage 2
ct2 = V.encrypt(pt, 10, 8, D, D, k4.KA, k4.KA, 1, 1, k1, k2)
o, _ = V.run_chunk((ct2, "ii", k4.KA, strs[tgt - 100:tgt + 100], names[tgt - 100:tgt + 100], 12, 4))
print("role ii planted:", [(h[0], h[1], h[2], h[4:]) for h in o])
