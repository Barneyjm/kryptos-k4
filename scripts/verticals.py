"""Vertical keys: each K4 letter (row r, col c in the carved layout) keyed by the letter d rows above/below in the same
column, taken from the cipher panel or the tableau (plain or mirrored), with small column shifts. Exact + best partial."""
import numpy as np, k4
CL = open("panel.txt").read().split(); TL = open("tableau.txt").read().split()
cells = [(r, c, ch) for r, line in enumerate(CL) for c, ch in enumerate(line) if ch != "?"]
k4cells = cells[-97:]; CR = dict(zip(k4.CPOS.tolist(), k4.CPL))
def get(src, r, c, mirror):
    if not 0 <= r < 28: return None
    row = (CL if src == "cipher" else TL)[r]
    cc = len(row) - 1 - c if mirror else c
    ch = row[cc] if 0 <= cc < len(row) else None
    return ch if ch and ch.isalpha() else None
res = []
for src in ("cipher", "tableau"):
    for mirror in (False, True):
        for d in range(-27, 28):
            if src == "cipher" and d == 0: continue
            for cs in range(-3, 4):
                keys = {i: get(src, r + d, c + cs, mirror) for i, (r, c, _) in enumerate(k4cells)}
                usable = [i for i in CR if keys[i]]
                if len(usable) < 16: continue
                for an, X in (("AZ", k4.AZ), ("KA", k4.KA)):
                    for vn, f in (("vig", lambda c, k: (c - k) % 26), ("beau", lambda c, k: (k - c) % 26), ("vb", lambda c, k: (c + k) % 26)):
                        m = sum(X[f(X.index(k4.CT[i]), X.index(keys[i]))] == CR[i] for i in usable)
                        res.append((m, len(usable), src, mirror, d, cs, an, vn))
res.sort(reverse=True)
print(f"configurations: {len(res):,}; exact fits: {sum(1 for x in res if x[0] == x[1])}")
print("best partial matches (matched, checkable, source, mirrored, rows offset, col shift, alphabet, type):")
for x in res[:5]: print("  ", x)
rng = np.random.default_rng(0); print("expected best by chance across this many tries: ~4-5 of 24")
