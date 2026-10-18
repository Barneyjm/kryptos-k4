"""Per-period breakdown of fe_affine (periodic affine key) hits: K4 vs chance vs one shuffled CT."""
import numpy as np, collections, k4, fe_common as F, fe_affine as A
AL = F.all_alphabets(); rng = np.random.default_rng(21)
for label, ct in (("K4", k4.CT), ("shuffled", "".join(rng.permutation(list(k4.CT))))):
    h = []
    for s in range(0, len(AL), 100000): h += A.scan(ct, AL[s:s + 100000], range(1, 14), A.UNITS)
    c = collections.Counter(x[3] for x in h)
    print(label, {L: (c[L], round(len(AL) * 5 * float(np.prod([min(1.0, 12 * 26.0 ** -(len(m) - 1)) for m in A.classes(L)])), 3)) for L in range(1, 14)}, flush=True)
    if label == "K4": print("  K4 hits with L<=10:", [x for x in h if x[3] <= 10][:20])
