import itertools, sys, numpy as np, multiprocessing as mp, k4
from texts import KEYWORDS

def check_generic(dests, ct, cpos, cpl, label, out, min_cons=8):
    dests = np.atleast_2d(dests); dc = dests[:, cpos]
    for an, X in k4.ALPHS.items():
        cta = np.array([X.index(c) if c in X else -100 for c in ct])
        ci = cta[dc]; pi = np.array([X.index(c) for c in cpl])[None, :]
        bad = (ci < 0).any(axis=1)
        for kn, k in (("vig", (ci - pi) % 26), ("beau", (ci + pi) % 26)):
            for model in "AB":
                for L in range(1, 27):
                    r = np.broadcast_to(cpos[None, :] % L, k.shape) if model == "A" else dc % L
                    ok = ~((r[:, :, None] == r[:, None, :]) & (k[:, :, None] != k[:, None, :])).any(axis=(1, 2)) & ~bad
                    for i in np.nonzero(ok)[0]:
                        cons = len(r[i]) - len(set(r[i].tolist()))
                        if cons >= min_cons: out.append((label(i), an, kn, model, L, cons))

def dest_from_order(order):
    order = np.asarray(order); d = np.empty_like(order); d[order] = np.arange(len(order)); return d

def kw_order(word, w):
    word = (word * (w // len(word) + 1))[:w]
    return [i for _, i in sorted((c, i) for i, c in enumerate(word))]

# ---- part 1: width 10 exhaustive, parallel by first column
def w10_worker(first):
    hits = []; w = 10; rest = [c for c in range(w) if c != first]
    batch = []
    def flush():
        orders = np.array([k4.order_from(w, p) for p in batch])
        dests = np.empty_like(orders); rows = np.arange(len(batch))[:, None]; dests[rows, orders] = np.arange(97)
        for dname, D in (("enc", dests), ("dec", orders)):
            check_generic(D, k4.CT, k4.CPOS, k4.CPL, lambda i, b=list(batch), d=dname: f"w10 {d} {b[i]}", hits)
        batch.clear()
    for p in itertools.permutations(rest):
        batch.append((first,) + p)
        if len(batch) == 15000: flush()
    if batch: flush()
    return hits

# ---- part 2/3: keyword columnar + K3-style rotations
def rotations(n, rows, cols):
    """Yield (name, order) for writing n chars row-wise in rows x cols grid and reading with a rotation."""
    g = np.arange(rows * cols).reshape(rows, cols)
    for k in range(4):
        for flip in (False, True):
            h = np.rot90(g, k); h = np.fliplr(h) if flip else h
            yield f"rot{k}{'f' if flip else ''}", h.flatten()

def structured():
    hits = []
    ct, cpos, cpl = k4.CT, k4.CPOS, k4.CPL
    for w in range(2, 31):
        for kw in KEYWORDS:
            order = np.array(k4.order_from(w, kw_order(kw, w)))
            for dname, D in (("enc", dest_from_order(order)), ("dec", order)):
                check_generic(D, ct, cpos, cpl, lambda i, w=w, kw=kw, d=dname: f"kwcol {kw} w={w} {d}", hits)
    # rotations: 97 + pads, and '?'+K4 (98)
    variants = [("K4", k4.CT, 0)] + [("?K4", "?" + k4.CT, 1)] + [(f"K4+{p}pad", k4.CT + "#" * p, 0) for p in range(1, 15)]
    for vname, s, shift in variants:
        n = len(s); cp = cpos + shift
        dims = [(r, n // r) for r in range(2, n) if n % r == 0]
        for r1, c1 in dims:
            for n1, o1 in rotations(n, r1, c1):
                for dname, D in (("enc", dest_from_order(o1)), ("dec", o1)):
                    check_generic(D, s, cp, cpl, lambda i, a=(vname, r1, c1, n1, dname): f"rot {a}", hits)
                for r2, c2 in dims:  # double rotation (K3 used two)
                    for n2, o2 in rotations(n, r2, c2):
                        o = o1[o2]
                        for dname, D in (("enc", dest_from_order(o)), ("dec", o)):
                            check_generic(D, s, cp, cpl, lambda i, a=(vname, r1, c1, n1, r2, c2, n2, dname): f"rot2 {a}", hits)
    return hits

if __name__ == "__main__":
    h = structured(); print("structured hits:", len(h)); [print(x) for x in h[:30]]
    with mp.Pool(12) as pool:
        res = pool.map(w10_worker, range(10))
    h10 = [x for r in res for x in r]; print("w10 hits:", len(h10)); [print(x) for x in h10[:30]]
