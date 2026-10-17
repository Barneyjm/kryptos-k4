import numpy as np, itertools, k4, round2
from texts import KEYWORDS
def colorder(kw, w=None):
    w = w or len(kw); return np.array(k4.order_from(w, round2.kw_order(kw, w)))
ords = {}
for kw in KEYWORDS:
    ords[kw] = colorder(kw)
    for w in (7, 8, 10, 11, 13, 14, 24): ords[f"{kw}@{w}"] = colorder(kw, w)
names = list(ords); hits = []
batch, labels = [], []
def flush():
    D = np.array(batch)
    round2.check_generic(D, k4.CT, k4.CPOS, k4.CPL, lambda i, L=list(labels): L[i], hits)
    batch.clear(); labels.clear()
for a, b in itertools.product(names, repeat=2):
    o = ords[a][ords[b]]                     # C[n] = P[o[n]] (apply a then b)
    for d, arr in (("enc", round2.dest_from_order(o)), ("dec", o)):
        batch.append(arr); labels.append(f"dbl {a} | {b} {d}")
    if len(batch) >= 20000: flush()
flush()
print("double-columnar pairs:", len(names)**2, "hits:", len(hits)); [print(h) for h in hits[:20]]
