"""Physical overlay: cipher panel (28 lines) vs tableau (28 lines). Key for each K4 letter = tableau letter at the same
(row, column), optionally mirrored (seen from behind), with row shifts and column offsets. Exact crib test, vig/beau/vb, AZ/KA.
Also: coincidence rate between the two halves at aligned positions across the whole sculpture."""
import numpy as np, k4
CL = open("panel.txt").read().split(); TL = open("tableau.txt").read().split()
assert len(CL) == 28 and len(TL) == 28
# K4 cells: (row, col) of each K4 letter in the carved layout ('?' occupies a cell)
cells = []
for r, line in enumerate(CL):
    for c, ch in enumerate(line): cells.append((r, c, ch))
letters = [(r, c, ch) for r, c, ch in cells if ch != "?"]
k4cells = letters[-97:]; assert "".join(ch for _, _, ch in k4cells) == k4.CT
CR = dict(zip(k4.CPOS.tolist(), k4.CPL)); hits = []; best = []
def tab(r, c, mirror, rs, cs):
    rr = r + rs
    if not 0 <= rr < 28: return None
    row = TL[rr]; cc = (len(row) - 1 - c if mirror else c) + cs
    return row[cc] if 0 <= cc < len(row) else None
for mirror in (False, True):
    for rs in range(-3, 4):
        for cs in range(-6, 7):
            keys = {i: tab(r, c, mirror, rs, cs) for i, (r, c, _) in enumerate(k4cells)}
            if any(keys[i] is None for i in CR): continue
            for an, X in (("AZ", k4.AZ), ("KA", k4.KA)):
                for vn, f in (("vig", lambda c, k: (c - k) % 26), ("beau", lambda c, k: (k - c) % 26), ("vb", lambda c, k: (c + k) % 26)):
                    m = sum(X[f(X.index(k4.CT[i]), X.index(keys[i]))] == p for i, p in CR.items())
                    best.append((m, mirror, rs, cs, an, vn))
                    if m == 24: hits.append((mirror, rs, cs, an, vn))
print("overlay-key exact fits:", hits or "none"); print("best partial:", sorted(best, reverse=True)[:3], "(chance ~1/24)")
# coincidence between halves at aligned positions (whole sculpture)
for mirror in (False, True):
    same = tot = 0
    for r in range(28):
        a, b = CL[r], TL[r]
        for c in range(min(len(a), len(b))):
            bc = b[len(b) - 1 - c] if mirror else b[c]
            if a[c] != "?" and bc != " ": tot += 1; same += a[c] == bc
    print(f"{'mirrored' if mirror else 'direct  '} overlay: {same}/{tot} aligned letters identical ({same/tot:.3%}; chance ~3.8%)")
