"""Library of 'dumb' hand transpositions; each tested with every substitution family."""
import itertools, numpy as np, k4, round2, tol
from sympy import Matrix
N = 97
def grid_orders(w, lead=0):
    """Positions laid out in rows of width w; first row starts at column `lead` (physical layout)."""
    cells = [(r, c) for r, c in ((divmod(i + lead, w)) for i in range(N))]
    rows = {}
    for i, (r, c) in enumerate(cells): rows.setdefault(r, {})[c] = i
    R = sorted(rows); out = {}
    def rd(seq_rows, colwise=False, snake=False, rev_cols=False, rev_rows=False):
        o = []
        if not colwise:
            for k, r in enumerate(seq_rows):
                cs = sorted(rows[r]); cs = cs[::-1] if (rev_rows ^ (snake and k % 2)) else cs
                o += [rows[r][c] for c in cs]
        else:
            cols = sorted({c for r in R for c in rows[r]}); cols = cols[::-1] if rev_cols else cols
            for k, c in enumerate(cols):
                rr = [r for r in seq_rows if c in rows[r]]; rr = rr[::-1] if (snake and k % 2) else rr
                o += [rows[r][c] for r in rr]
        return o
    for top in (False, True):
        seq = R[::-1] if top else R
        for snake in (False, True):
            out[f"rows{'Up' if top else ''}{' snake' if snake else ''}"] = rd(seq, snake=snake)
            out[f"rows{'Up' if top else ''} rev{' snake' if snake else ''}"] = rd(seq, snake=snake, rev_rows=True)
            for rc in (False, True):
                out[f"cols{'Up' if top else ''}{' R2L' if rc else ''}{' snake' if snake else ''}"] = rd(seq, colwise=True, snake=snake, rev_cols=rc)
    return {f"w{w}+{lead} {k}": v for k, v in out.items() if sorted(v) == list(range(N))}
def library():
    lib = {"identity": list(range(N)), "reverse": list(range(N))[::-1]}
    for w in range(2, 49):
        lib.update(grid_orders(w))
    for lead in range(31): lib.update(grid_orders(31, lead))      # carved layout: 31 wide, K4 starts mid-row
    for k in range(2, N):                                          # decimation: every k-th letter (97 prime)
        for s in (0,):
            lib[f"every{k}th"] = [(s + k * i) % N for i in range(N)]
    for rails in range(2, 12):                                     # rail fence (zigzag) and simple even/odd split
        pat = list(range(rails)) + list(range(rails - 2, 0, -1))
        rail = [pat[i % len(pat)] for i in range(N)]
        lib[f"railfence{rails}"] = [i for r in range(rails) for i in range(N) if rail[i] == r]
        lib[f"split{rails}"] = [i for r in range(rails) for i in range(r, N, rails)]
    # also each of these reversed afterwards
    for k, v in list(lib.items()): lib[k + " +rev"] = v[::-1]
    uniq = {}
    for k, v in lib.items(): uniq.setdefault(tuple(v), k)
    return {k: list(v) for v, k in uniq.items()}
def hill_test(Xin, Y, n, X):
    blocks = []
    for a in range(n):
        for s in range(a, N - n + 1, n):
            if all(Xin[j] >= 0 for j in range(s, s + n)): blocks.append((a, [Xin[j] for j in range(s, s + n)], [Y[j] for j in range(s, s + n)]))
    hits = []
    for a in range(n):
        bl = [(p, c) for aa, p, c in blocks if aa == a]
        if len(bl) <= n: continue
        for combo in itertools.combinations(range(len(bl)), n):
            A = Matrix([bl[c][0] for c in combo]).T; B = Matrix([bl[c][1] for c in combo]).T
            try: K = (B * A.inv_mod(26)).applyfunc(lambda x: x % 26)
            except Exception: continue
            if all(list((K * Matrix(p)).applyfunc(lambda x: x % 26)) == c for p, c in bl): hits.append((n, a, K.tolist()))
            break
    return hits
if __name__ == "__main__":
    lib = library(); print("distinct dumb transpositions:", len(lib))
    names = list(lib); orders = np.array([lib[k] for k in names])
    dests = np.empty_like(orders); dests[np.arange(len(names))[:, None], orders] = np.arange(N)
    exact, tolh = [], []
    for D in (dests, orders):
        round2.check_generic(D, k4.CT, k4.CPOS, k4.CPL, lambda i: names[i], exact, min_cons=6)
        tol.tol_check(D, k4.CT, k4.CPOS, k4.CPL, lambda i: names[i], tolh, tol.STD, mineff=7)
    from round3 import ALPHA
    quag = {f"{p}|{c}": (ALPHA[p], ALPHA[c]) for p in ("AZ", "KA", "kw:PALIMPSEST", "kw:ABSCISSA", "kw:BERLINCLOCK", "kw:KRYPTOS")
            for c in ("AZ", "KA", "kw:PALIMPSEST", "kw:ABSCISSA", "kw:BERLINCLOCK", "kw:KRYPTOS")}
    qh = []; tol.tol_check(np.vstack([dests, orders]), k4.CT, k4.CPOS, k4.CPL, lambda i: names[i % len(names)], qh, quag, maxerr=0 + 2, mineff=7)
    print("exact periodic hits:", exact[:20]); print("1-2 error hits:", tolh[:20]); print("quagmire tolerant hits:", qh[:20])
    # Hill 2/3/4 after/before each dumb transposition
    hh = []
    for an, X in (("AZ", k4.AZ), ("KA", k4.KA)):
        Cn = [X.index(c) for c in k4.CT]; Pk = {int(i): X.index(c) for i, c in zip(k4.CPOS, k4.CPL)}
        for k, o in lib.items():
            d = [0] * N
            for nn, i in enumerate(o): d[i] = nn
            # Model A: Hill on plaintext, then transposition: H[i] = C[d[i]]
            XinA = [Pk.get(i, -1) for i in range(N)]; YA = [Cn[d[i]] for i in range(N)]
            # Model B: transposition, then Hill: T[n] = P[o[n]]
            XinB = [Pk.get(o[nn], -1) for nn in range(N)]
            for n in (2, 3, 4):
                for model, (xi, y) in (("A", (XinA, YA)), ("B", (XinB, Cn))):
                    for h in hill_test(xi, y, n, X): hh.append((k, an, model, h))
    print("Hill hits:", hh[:20])
