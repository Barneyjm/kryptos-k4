"""T1 alphabet-free extension: every interrupted-key model of fe_interrupt, but each key position (class) may use an
ARBITRARY independent alphabet (covers unknown tableaux / strip ciphers / any keyed alphabets at once).
Necessary condition: plain->cipher is a partial bijection inside every class. Exact. Planted + shuffled calibration."""
import numpy as np, collections, k4, fe_common as F, fe_interrupt as T
def ok(ct, labels):
    f, g = {}, {}
    for j, (p, ch) in enumerate(zip(F.POS, k4.CPL)):
        l = labels[j]
        if l is None: continue
        c = ct[p]
        if f.setdefault((l, ch), c) != c or g.setdefault((l, c), ch) != ch: return False
    return True
def strength(labels):   # number of same-class crib pairs (rough constraint count)
    c = collections.Counter(l for l in labels if l is not None); return sum(v * (v - 1) // 2 for v in c.values())
rng = np.random.default_rng(4)
# planted: ct_restart_after T=K, L=6, random independent alphabets per key position
L = 6; al = [rng.permutation(26) for _ in range(L)]; pt = F.planted_plaintext(rng); out = []; last = -1
for i, ch in enumerate(pt):
    c = chr(65 + al[(i - last - 1) % L][ord(ch) - 65]); out.append(c)
    if c == "K": last = i
pct = "".join(out)
names = {n: l for n, l in T.gen_models(pct)}
print("planted ct_restart_after T=K L=6 (independent alphabets) recovered:", ok(pct, names["ct_restart_after T=K L=6"]))
def run(ct):
    res = collections.Counter(); tot = collections.Counter(); strong = []
    for n, l in T.gen_models(ct):
        fam = n.split()[0]; s = strength(l)
        if s < 20: continue
        tot[fam] += 1
        if ok(ct, l): res[fam] += 1; strong.append((s, n))
    return res, tot, sorted(strong, reverse=True)
r, t, s = run(k4.CT)
print("K4: consistent / tested (models with >= 20 same-class pairs):", {f: (r[f], t[f]) for f in t})
print("K4 strongest consistent models:", s[:15])
for k in range(3):
    sh = "".join(rng.permutation(list(k4.CT))); r2, t2, s2 = run(sh)
    print(f"shuffled{k}:", {f: (r2[f], t2[f]) for f in t2}, " strongest:", s2[:3])
