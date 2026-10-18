"""Calibrate the one notable alphabet-free fit: line_restart L=7 (key restarts at each carved line; 7 arbitrary alphabets)."""
import numpy as np, k4, fe_common as F
LS = [0, 4, 35, 66]
def lab(L): return [(p - max(b for b in LS if b <= p)) % L for p in F.POS]
def ok(ct, labels):
    f, g = {}, {}
    for l, p, ch in zip(labels, F.POS, k4.CPL):
        c = ct[p]
        if f.setdefault((l, ch), c) != c or g.setdefault((l, c), ch) != ch: return False
    return True
rng = np.random.default_rng(99); N = 20000
for L in range(2, 27):
    k = ok(k4.CT, lab(L))
    fr = np.mean([ok("".join(rng.permutation(list(k4.CT))), lab(L)) for _ in range(N if L == 7 else 2000)])
    print(f"line_restart L={L:2d}: K4 consistent={k}  shuffled fraction consistent={fr:.4f}")
