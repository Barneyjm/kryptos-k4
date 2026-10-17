"""Model-free pattern search: is there ANY linear rule mod 26 linking each letter to its neighbours/position?
target = sum(a_f * feature_f) mod 26, features drawn from: const, i, i^2, P[i], P[i-1], P[i-2], P[i+1], C[i-1], C[i-2], C[i+1], C[i+2].
Solved exactly over Z2 and Z13 (Z26 = Z2 x Z13). Report rules consistent with every usable equation and with spare equations."""
import itertools, numpy as np, k4
P = dict(zip(k4.CPOS.tolist(), k4.CPL)); CT = k4.CT
def rank_consistent(A, b, p):
    A = [list(r) + [bb] for r, bb in zip(A, b)]; m, n = len(A), len(A[0]) - 1; r = 0
    for c in range(n):
        pr = next((i for i in range(r, m) if A[i][c] % p), None)
        if pr is None: continue
        A[r], A[pr] = A[pr], A[r]; inv = pow(A[r][c], -1, p); A[r] = [x * inv % p for x in A[r]]
        for i in range(m):
            if i != r and A[i][c] % p: f = A[i][c]; A[i] = [(x - f * y) % p for x, y in zip(A[i], A[r])]
        r += 1
    consistent = not any(all(x % p == 0 for x in row[:n]) and row[n] % p for row in A)
    return consistent, r
results = []
for an, X in (("AZ", k4.AZ), ("KA", k4.KA)):
    pv = {i: X.index(c) for i, c in P.items()}; cv = {i: X.index(c) for i, c in enumerate(CT)}
    feats = {"1": lambda i: 1, "i": lambda i: i, "i2": lambda i: i * i % 26,
             "P0": lambda i: pv.get(i), "P-1": lambda i: pv.get(i - 1), "P-2": lambda i: pv.get(i - 2), "P+1": lambda i: pv.get(i + 1),
             "C0": lambda i: cv.get(i), "C-1": lambda i: cv.get(i - 1), "C-2": lambda i: cv.get(i - 2), "C+1": lambda i: cv.get(i + 1), "C+2": lambda i: cv.get(i + 2)}
    for target in ("C0", "P0"):
        others = [f for f in feats if f != target]
        for k in range(1, 6):
            for S in itertools.combinations(others, k):
                if target == "C0" and "P0" not in S: continue
                if target == "P0" and "C0" not in S: continue
                rows, rhs = [], []
                for i in P:
                    vals = [feats[f](i) for f in S]; t = feats[target](i)
                    if None in vals or t is None: continue
                    rows.append(vals); rhs.append(t)
                ok2, r2 = rank_consistent(rows, rhs, 2); ok13, r13 = rank_consistent(rows, rhs, 13)
                slack = len(rows) - max(r2, r13)
                if ok2 and ok13 and slack >= 6:
                    results.append((slack, an, target, S, len(rows)))
print("linear rules tested over both alphabets; consistent rules with >= 6 spare equations:", len(results))
for r in sorted(results, reverse=True)[:15]: print(" ", r)
