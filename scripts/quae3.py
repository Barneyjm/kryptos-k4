import itertools, k4
from linear import rank_consistent
CR = list(zip(k4.CPOS.tolist(), k4.CPL))
pairs = [(L1, L2) for L1 in range(2, 27) for L2 in range(L1 + 1, 27)]
ROWS = {}
for L1, L2 in pairs:
    rows = []
    for i, _ in CR:
        r = [0] * (L1 + L2); r[i % L1] = 1; r[L1 + i % L2] = 1; rows.append(r)
    ROWS[(L1, L2)] = rows
hits = []
for q in sorted({"".join(p) for p in itertools.permutations("QUAE")}):
    for an, X in k4.ALPHS.items():
        for kn, f in (("vig", lambda c, p: (c - p) % 26), ("beau", lambda c, p: (c + p) % 26)):
            for sgn in (1, -1):
                k = [(f(X.index(k4.CT[i]), X.index(p)) - sgn * X.index(q[i % 4])) % 26 for i, p in CR]
                for (L1, L2), rows in ROWS.items():
                    a, r2 = rank_consistent(rows, k, 2)
                    if not a: continue
                    b, r13 = rank_consistent(rows, k, 13)
                    if b and 24 - max(r2, r13) >= 4: hits.append((24 - max(r2, r13), q, an, kn, sgn, L1, L2))
print("QUAE (any order) + two unknown splits, consistent with >= 4 spare checks:", sorted(hits, reverse=True)[:10] or "none")
