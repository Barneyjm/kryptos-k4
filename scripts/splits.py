"""Scheidt-style key splits: key_i = a[i mod L1] + b[i mod L2] (+ c[i mod L3]) mod 26.
(1) unknown splits: linear algebra over Z2 x Z13 for every (L1, L2); (2) known clue words as splits."""
import itertools, numpy as np, k4
from linear import rank_consistent
from texts import KEYWORDS
F = {"vig": lambda c, p: (c - p) % 26, "beau": lambda c, p: (c + p) % 26}
CR = list(zip(k4.CPOS.tolist(), k4.CPL))
res = []
for an, X in k4.ALPHS.items():
    for kn, f in F.items():
        k = [f(X.index(k4.CT[i]), X.index(p)) for i, p in CR]
        for L1 in range(2, 27):
            for L2 in range(L1 + 1, 27):
                rows = []
                for i, _ in CR:
                    r = [0] * (L1 + L2); r[i % L1] = 1; r[L1 + i % L2] = 1; rows.append(r)
                ok2, r2 = rank_consistent(rows, k, 2); ok13, r13 = rank_consistent(rows, k, 13)
                slack = 24 - max(r2, r13)
                if ok2 and ok13 and slack >= 4: res.append((slack, an, kn, L1, L2))
print("(1) two unknown splits, consistent with >= 4 spare checks:", sorted(res, reverse=True)[:10] or "none")
words = sorted(set(KEYWORDS + ["SHADOW", "LUCID", "MEMORY", "POSITION", "QUAE", "DYAHR", "YAR", "TENJ", "WELTZEITUHR", "URANIA"]))
hits = []
for an, X in k4.ALPHS.items():
    for kn, f in F.items():
        need = [f(X.index(k4.CT[i]), X.index(p)) for i, p in CR]
        for r in (2, 3):
            for combo in itertools.combinations(words, r):
                for sgn in itertools.product((1, -1), repeat=r - 1):
                    for off in itertools.product(range(26), repeat=0):
                        m = sum(((X.index(combo[0][i % len(combo[0])]) + sum(s * X.index(w[i % len(w)]) for s, w in zip(sgn, combo[1:]))) % 26) == v
                                for (i, _), v in zip(CR, need))
                        if m >= 18: hits.append((m, an, kn, combo, sgn))
print("(2) clue-word splits (pairs and triples, added or subtracted), >= 18/24:", sorted(hits, reverse=True)[:10] or "none")
