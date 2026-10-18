"""Gaps 3+4: (a) keyword expanded by writing it into a w-wide grid and reading columns (text stays in place);
(b) block-constant key (new arbitrary setting every n letters, any phase); (c) keyword + per-letter drift."""
import numpy as np, itertools, time, pickle, k4, cokernel, gap1, rk_mixed, gap2
from round3 import ALPHA
POS = k4.CPOS.tolist()
def onehot(vals):
    u = sorted(set(vals)); I = {v: j for j, v in enumerate(u)}; M = np.zeros((24, len(u)), int)
    for r, v in enumerate(vals): M[r, I[v]] = 1
    return M
def models():
    out = {}
    for w in range(2, 31):
        rows = -(-97 // w); cells = [(r, c) for r in range(rows) for c in range(w) if r * w + c < 97]
        colorder = sorted(cells, key=lambda rc: (rc[1], rc[0]))          # read columns top-down, left to right
        for name, seq in (("rows->cols", colorder), ("rows->cols R2L", sorted(cells, key=lambda rc: (-rc[1], rc[0])))):
            cell_of_n = {n: r * w + c for n, (r, c) in enumerate(seq)}  # n-th key letter came from cell index
            n_of_cell = {v: k for k, v in cell_of_n.items()}
            for L in range(1, 27):
                out[f"keygrid w={w} {name} L={L}"] = onehot([cell_of_n[p] % L for p in POS])      # key written row-wise, read col-wise
                out[f"keygrid w={w} {name} inv L={L}"] = onehot([n_of_cell[p] % L for p in POS])  # key written col-wise, read row-wise
    for n in range(2, 31):
        for off in range(n): out[f"block n={n} off={off}"] = onehot([(p + off) // n for p in POS])
    for L in range(1, 27): out[f"drift L={L}"] = np.hstack([onehot([p % L for p in POS]), np.array([[p] for p in POS])])
    return out
if __name__ == "__main__":
    D = []
    for name, M in models().items():
        N2, N13 = cokernel.left_null(M, 2), cokernel.left_null(M, 13)
        if min(len(N2), len(N13)) >= 5: D.append((name, N2, N13, len(N13)))
    print("models with >= 5 checks:", len(D), flush=True)
    # std pairs directly
    hits = []
    for (pn, PA), (cn, CA) in itertools.product((("AZ", k4.AZ), ("KA", k4.KA)), repeat=2):
        for kn, f in (("vig", lambda c, p: (c - p) % 26), ("beau", lambda c, p: (c + p) % 26)):
            k = np.array([f(CA.index(k4.CT[i]), PA.index(p)) for i, p in zip(POS, k4.CPL)])
            for name, N2, N13, c in D:
                if not (N2 @ k % 2).any() and not (N13 @ k % 13).any(): hits.append((name, pn, cn, kn, c))
    print("std/KA pairs:", len(hits), hits[:10], flush=True)
    gap1.P = D; gap1.names = [d[0] for d in D]
    gap1.B13 = np.vstack([d[2] for d in D]).astype(np.float32); gap1.B2 = np.vstack([d[1] for d in D]).astype(np.float32)
    gap1.seg13 = np.cumsum([0] + [len(d[2]) for d in D])[:-1]; gap1.seg2 = np.cumsum([0] + [len(d[1]) for d in D])[:-1]
    AL = [(v, k) for k, v in ALPHA.items()]
    t = time.time(); h = gap1.scan(k4.CT, AL, "K4"); print(f"clue alphabets: {len(h)} hits ({time.time()-t:.0f}s)", flush=True); [print("  ", x) for x in h[:10]]
    AL = rk_mixed.alph_list() + gap2.build()[0]
    t = time.time(); h = gap1.scan(k4.CT, AL, "K4"); print(f"dictionary+constructed ({len(AL):,}): {len(h)} hits ({time.time()-t:.0f}s)", flush=True)
    for x in h[:20]: print("  ", x, flush=True)
    pickle.dump(h, open("keyperm_hits.pkl", "wb"))
