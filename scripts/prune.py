"""Pruned columnar search: assign columns slot by slot, check crib consistency as soon as letters land.
Configs: (alphabet AZ/KA) x (vig/beau) x (model A: key by plaintext pos | B: key by ciphertext pos) x period 1..26."""
import numpy as np, sys, time, multiprocessing as mp, itertools, k4
N = 97; LMAX = 26; MINCONS = 8
CT = k4.CT; CRIB = dict(zip(k4.CPOS.tolist(), k4.CPL))
ALPHS = [k4.AZ, k4.KA]
LS = np.arange(1, LMAX + 1)                                   # periods
def vals(dest, p):                                             # key value per (alphabet, variant) -> (4,)
    out = []
    for X in ALPHS:
        c, q = X.index(CT[dest]), X.index(p); out += [(c - q) % 26, (c + q) % 26]
    return np.array(out)
class Search:
    def __init__(self, w, direction):
        self.w, self.dir = w, direction
        q, r = divmod(N, w); self.len = [q + 1 if c < r else q for c in range(w)]
        self.nodes = 0; self.hits = []
    def letters(self, c, off):
        """crib letters placed when column c is assigned and its block starts at ciphertext/plaintext offset `off`."""
        out = []
        if self.dir == "enc":                                  # plaintext rows-> columns read in order: dest = off + row
            for row in range(self.len[c]):
                i = c + row * self.w
                if i in CRIB: out.append((i, off + row))
        else:                                                  # plaintext read by columns: P[off+row] sits at grid cell c+row*w
            for row in range(self.len[c]):
                i = off + row
                if i in CRIB: out.append((i, c + row * self.w))
        return out
    def run(self, prefix=()):
        S = -np.ones((4, 2, LMAX, LMAX), np.int16)            # key value at residue
        cons = np.zeros((4, 2, LMAX), np.int16); alive = np.ones((4, 2, LMAX), bool)
        for L in LS: alive[:, :, L - 1] &= True
        used = [False] * self.w; order = []
        state = (S, cons, alive)
        for c in prefix:
            state = self.place(state, c, sum(self.len[x] for x in order))
            if state is None: return
            used[c] = True; order.append(c)
        self.dfs(state, used, order)
    def place(self, state, c, off):
        S, cons, alive = state; S = S.copy(); cons = cons.copy(); alive = alive.copy()
        for i, d in self.letters(c, off):
            v = vals(d, CRIB[i])                               # (4,)
            for m, pos in ((0, i), (1, d)):
                r = pos % LS                                    # (26,)
                cur = S[:, m, LS - 1, r]                        # (4,26)
                conflict = (cur >= 0) & (cur != v[:, None])
                agree = (cur >= 0) & (cur == v[:, None])
                alive[:, m, :] &= ~conflict
                cons[:, m, :] += agree
                S[:, m, LS - 1, r] = np.where(cur < 0, v[:, None], cur)
        if not alive.any(): return None
        return S, cons, alive
    def dfs(self, state, used, order):
        self.nodes += 1
        if len(order) == self.w:
            S, cons, alive = state
            ok = alive & (cons >= MINCONS)
            for a, m, L in zip(*np.nonzero(ok)):
                self.hits.append((self.w, self.dir, tuple(order), ["AZvig", "AZbeau", "KAvig", "KAbeau"][a], "AB"[m], int(L) + 1, int(cons[a, m, L])))
            return
        off = sum(self.len[x] for x in order)
        for c in range(self.w):
            if used[c]: continue
            ns = self.place(state, c, off)
            if ns is None: continue
            used[c] = True; order.append(c); self.dfs(ns, used, order); order.pop(); used[c] = False
def job(args):
    w, d, prefix = args; s = Search(w, d); s.run(prefix); return s.nodes, s.hits
if __name__ == "__main__":
    for w in map(int, sys.argv[1:]):
        t = time.time(); jobs = [(w, d, (a, b)) for d in ("enc", "dec") for a in range(w) for b in range(w) if a != b]
        with mp.Pool(11) as pool: res = pool.map(job, jobs, chunksize=1)
        nodes = sum(r[0] for r in res); hits = [h for r in res for h in r[1]]
        print(f"w={w}: {nodes:,} nodes visited (vs {2*np.prod(np.arange(1, w+1), dtype=float):,.0f} full orders), hits {len(hits)}, {time.time()-t:.0f}s", flush=True)
        for h in hits[:10]: print("  ", h, flush=True)
