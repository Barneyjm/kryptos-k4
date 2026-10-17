import numpy as np, itertools, k4, dumb
rows = open("tableau.txt").read().split()
print("tableau letters:", sum(map(len, rows)))
body = rows[1:-1]
SRC = {"rows": "".join(rows), "rows_nolabel": "".join(r[1:] for r in body), "rows_rev": "".join(rows)[::-1],
       "boustro": "".join(r if i % 2 == 0 else r[::-1] for i, r in enumerate(rows))}
W = max(map(len, rows))
SRC["cols"] = "".join(r[c] for c in range(W) for r in rows if c < len(r))
SRC["cols_up"] = "".join(r[c] for c in range(W) for r in rows[::-1] if c < len(r))
SRC["diag"] = "".join(body[i][(i + j) % len(body[i])] for j in range(31) for i in range(26))
# ---- (a) tableau as running key: exact >= 20/24, with and without simple transpositions
lib = dumb.library(); names = list(lib); O = np.array([lib[n] for n in names])
D = np.empty_like(O); D[np.arange(len(names))[:, None], O] = np.arange(97)
DEST = np.vstack([D, O]); LAB = [n + " enc" for n in names] + [n + " dec" for n in names]
hits = []
for an, X in k4.ALPHS.items():
    cta = np.array([X.index(c) for c in k4.CT]); pi = np.array([X.index(c) for c in k4.CPL]); ci = cta[DEST[:, k4.CPOS]]
    for kn, need in (("vig", (ci - pi) % 26), ("beau", (ci + pi) % 26), ("vb", (pi - ci) % 26)):
        for sn, s in SRC.items():
            kk = np.array([X.index(c) for c in s + s[:200]]); offs = np.arange(len(s))
            for model in "AB":
                pos = np.broadcast_to(k4.CPOS, need.shape) if model == "A" else DEST[:, k4.CPOS]
                m = (kk[offs[:, None, None] + pos[None]] == need[None]).sum(axis=2)
                for o, t in zip(*np.nonzero(m >= 19)): hits.append((int(m[o, t]), LAB[t], an, kn, sn, model, int(o)))
print("(a) tableau running key, best:", sorted(hits, reverse=True)[:5] or "none >= 19/24")
# ---- (b) LAYER TWO: K4 = m interleaved layers (i mod m), each its own periodic key/alphabet/variant
F = {"vig": lambda c, p: (c - p) % 26, "beau": lambda c, p: (c + p) % 26}
cr = list(zip(k4.CPOS.tolist(), k4.CPL))
res = []
for m in range(2, 7):
    layer_ok = []
    for cls in range(m):
        pts = [(i // m, i, p) for i, p in cr if i % m == cls]
        opts = []
        for an, X in k4.ALPHS.items():
            for kn, f in F.items():
                ks = [(j, f(X.index(k4.CT[i]), X.index(p))) for j, i, p in pts]
                for L in range(1, 27):
                    d = {}; ok = all(d.setdefault(j % L, v) == v for j, v in ks)
                    if ok: opts.append((len(ks) - len(d), an, kn, L))
        opts = [o for o in opts if o[0] >= 2]
        layer_ok.append(max(opts) if opts else None)
    if all(layer_ok):
        tot = sum(o[0] for o in layer_ok)
        res.append((tot, m, layer_ok))
print("(b) interleaved layers (every layer must fit with >= 2 checks):", sorted(res, reverse=True)[:5] or "none")
