"""Transposition solver: (1) grid readings for every factorization (rows/cols, snake, rotations);
(2) keyed columnar (both directions) hill-climbing on column order with quadgram score."""
import json, sys, random, numpy as np
from quad import Q
def qs(t): a = np.frombuffer(t.encode(), np.uint8).astype(np.int64) - 65; return Q[((a[:-3]*26+a[1:-2])*26+a[2:-1])*26+a[3:]].mean()
def col_decrypt(ct, order, enc=True):
    n = len(ct); w = len(order); q, r = divmod(n, w); lens = [q + (1 if c < r else 0) for c in range(w)]
    if enc:   # ciphertext = columns read in `order`; rebuild rows
        cols = {}; i = 0
        for c in order: cols[c] = ct[i:i + lens[c]]; i += lens[c]
        return "".join(cols[c][row] for row in range(q + 1) for c in range(w) if row < lens[c])
    else:     # ciphertext written in rows, plaintext read down columns in `order`
        return "".join(ct[c + w * row] for c in order for row in range(lens[c]))
def climb(ct, w, enc, iters=4000, restarts=8, seed=0):
    rnd = random.Random(seed); best = (-99, None)
    for _ in range(restarts):
        o = list(range(w)); rnd.shuffle(o); cur = qs(col_decrypt(ct, o, enc))
        for _ in range(iters):
            n = o[:]; i, j = rnd.sample(range(w), 2)
            if rnd.random() < 0.5: n[i], n[j] = n[j], n[i]
            else: i, j = min(i, j), max(i, j); n[i:j+1] = n[i:j+1][::-1]
            s = qs(col_decrypt(ct, n, enc))
            if s >= cur: o, cur = n, s
        if cur > best[0]: best = (cur, o)
    return best
def grids(ct):
    n = len(ct); out = []
    for r in range(2, n):
        if n % r: continue
        g = np.array(list(ct)).reshape(r, n // r)
        for k in range(4):
            for flip in (0, 1):
                h = np.rot90(g, k); h = np.fliplr(h) if flip else h
                out.append((qs("".join(h.flatten())), f"{r}x{n//r} rot{k} flip{flip}", "".join(h.flatten())))
                sn = "".join("".join(row if i % 2 == 0 else row[::-1]) for i, row in enumerate(h))
                out.append((qs(sn), f"{r}x{n//r} rot{k} flip{flip} snake", sn))
    return sorted(out, reverse=True)
if __name__ == "__main__":
    ct = {p["slug"]: p["ct"] for p in json.load(open("ctf_puzzles.json"))}[sys.argv[1]]
    g = grids(ct)[0]; print(f"best simple grid reading: {g[0]:.2f} {g[1]}\n   {g[2][:120]}")
    best = []
    for w in range(2, int(sys.argv[2]) if len(sys.argv) > 2 else 16):
        for enc in (True, False):
            s, o = climb(ct, w, enc); best.append((s, w, enc, o))
    for s, w, enc, o in sorted(best, reverse=True)[:3]:
        print(f"columnar w={w} {'enc' if enc else 'dec'} score {s:.2f} order {o}\n   {col_decrypt(ct, o, enc)[:160]}")
