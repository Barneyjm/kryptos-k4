import itertools, numpy as np, k4, hill
from quad import score
X = k4.AZ
def solve_mod_p(A, b, p):
    """All x with A x = b (mod p prime). A: m x n list. Returns list of solutions."""
    A = [list(r) + [bb] for r, bb in zip(A, b)]; m, n = len(A), len(A[0]) - 1; piv = []; r = 0
    for c in range(n):
        pr = next((i for i in range(r, m) if A[i][c] % p), None)
        if pr is None: continue
        A[r], A[pr] = A[pr], A[r]; inv = pow(A[r][c], -1, p); A[r] = [x * inv % p for x in A[r]]
        for i in range(m):
            if i != r and A[i][c] % p: f = A[i][c]; A[i] = [(x - f * y) % p for x, y in zip(A[i], A[r])]
        piv.append(c); r += 1
    if any(all(x % p == 0 for x in row[:n]) and row[n] % p for row in A): return []
    free = [c for c in range(n) if c not in piv]; sols = []
    for fv in itertools.product(range(p), repeat=len(free)):
        x = [0] * n
        for c, v in zip(free, fv): x[c] = v
        for i, c in enumerate(piv): x[c] = (A[i][n] - sum(A[i][f] * x[f] for f in free)) % p
        sols.append(x)
    return sols
def crt(a2, a13): return [next(v for v in range(26) if v % 2 == x and v % 13 == y) for x, y in zip(a2, a13)]
for shift, align in ((0, 0), (1, 1)):
    ct = ("?" + k4.CT) if shift else k4.CT
    blocks = []
    for s in range(align, 97 + shift, 4):
        idx = [j - shift for j in range(s, s + 4)]
        if all(j in hill.P for j in idx): blocks.append(([X.index(hill.P[j]) for j in idx], [X.index(k4.CT[j]) for j in idx]))
    Cm = np.array([b[1] for b in blocks]); Pm = np.array([b[0] for b in blocks])   # rows = blocks
    # decryption matrix D: D c = p for each block -> for row r of D: Cm @ d_r = Pm[:, r]
    rows = []
    for r in range(4):
        s2 = solve_mod_p(Cm.tolist(), Pm[:, r].tolist(), 2); s13 = solve_mod_p(Cm.tolist(), Pm[:, r].tolist(), 13)
        rows.append([crt(a, b) for a in s2 for b in s13])
    print(f"shift{shift} align{align}: solutions per row", [len(x) for x in rows])
    best = []
    for D in itertools.product(*rows):
        D = np.array(D)
        det = int(round(np.linalg.det(D))) % 26
        if det % 2 == 0 or det % 13 == 0: continue
        pt = "".join("".join(X[v] for v in (D @ np.array([X.index(c) for c in ct[s:s+4]])) % 26)
                     for s in range(align, len(ct) - 3, 4) if "?" not in ct[s:s+4])
        best.append((score(pt), pt))
    best.sort(reverse=True); print(" valid invertible keys:", len(best))
    for sc, pt in best[:5]: print("  %.2f %s" % (sc, pt))
