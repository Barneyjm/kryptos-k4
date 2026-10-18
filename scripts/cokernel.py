"""Linear key models: key k (24 crib values) must lie in the column space of a fixed design matrix M (24 x u) over Z26.
Consistency <=> N2 k = 0 (mod 2) and N13 k = 0 (mod 13), N = left null space of M. Precompute N; apply to millions of k."""
import numpy as np, itertools, k4
POS = k4.CPOS.tolist()
def left_null(M, p):
    """Rows N with N M = 0 mod p (basis of null space of M^T)."""
    A = (np.array(M).T % p).astype(np.int64); m, n = A.shape   # A: u x 24 ; want x with A x = 0
    A = A.copy(); piv = []; r = 0
    for c in range(n):
        pr = next((i for i in range(r, m) if A[i, c] % p), None)
        if pr is None: continue
        A[[r, pr]] = A[[pr, r]]; A[r] = (A[r] * pow(int(A[r, c]), -1, p)) % p
        for i in range(m):
            if i != r and A[i, c] % p: A[i] = (A[i] - A[i, c] * A[r]) % p
        piv.append(c); r += 1
    free = [c for c in range(n) if c not in piv]; basis = []
    for f in free:
        x = np.zeros(n, np.int64); x[f] = 1
        for i, c in enumerate(piv): x[c] = (-A[i, f]) % p
        basis.append(x)
    return np.array(basis, np.int64).reshape(-1, n)
def models():
    out = {}
    def onehot(vals):
        u = sorted(set(vals)); I = {v: j for j, v in enumerate(u)}; M = np.zeros((24, len(u)), int)
        for r, v in enumerate(vals): M[r, I[v]] = 1
        return M
    for L in range(1, 27): out[f"periodic L={L}"] = onehot([p % L for p in POS])
    for L in range(1, 27): out[f"progressive L={L}"] = np.hstack([onehot([p % L for p in POS]), np.array([[p // L] for p in POS])])
    for W in range(2, 31):
        for off in range(W):
            out[f"grid W={W} off={off}"] = np.hstack([onehot([(p + off) % W for p in POS]), onehot([(p + off) // W for p in POS])])
    for L1, L2 in itertools.combinations(range(2, 27), 2):
        out[f"splits {L1}+{L2}"] = np.hstack([onehot([p % L1 for p in POS]), onehot([p % L2 for p in POS])])
    out["linear a*i+b"] = np.array([[1, p] for p in POS]); out["quadratic"] = np.array([[1, p, p * p] for p in POS])
    return out
def prepared(min_checks=8):
    res = []
    for name, M in models().items():
        N2, N13 = left_null(M, 2), left_null(M, 13)
        if min(len(N2), len(N13)) >= min_checks: res.append((name, N2, N13, len(N13)))
    return res
if __name__ == "__main__":
    P = prepared(); print(len(P), "models with >= 8 independent checks")
    # self-test on a planted periodic + grid + split key with std alphabet
    import random; random.seed(0)
    for name, gen in (("periodic L=9", lambda i: [3,1,4,1,5,9,2,6,5][i % 9]), ("grid W=7 off=2", lambda i: [4,0,19,2,7,11,8][(i+2) % 7] + [3,14,1,5,9,2,6,5,3,5,8,9,7,9][(i+2)//7]),
                      ("splits 5+11", lambda i: [1,7,2,8,3][i % 5] + [9,4,6,0,13,2,21,5,17,11,8][i % 11])):
        k = np.array([gen(p) % 26 for p in POS])
        hits = [m for m, N2, N13, c in P if not (N2 @ k % 2).any() and not (N13 @ k % 13).any()]
        print("planted", name, "->", name in hits, f"({len(hits)} models consistent)")
