"""Masking model: x[P_i] + z[i mod L] = num(C_i) (mod 26), x = ARBITRARY injective map on plaintext letters (the mask),
z = periodic key. Also the cipher-side version y[C_i] = num(P_i) + z[i mod L]. Solved exactly over Z2 and Z13."""
import itertools, numpy as np, k4
def solve_mod(A, b, p):
    """Return (consistent, particular solution or None, nullity) of A x = b mod prime p."""
    A = [list(r) + [bb % p] for r, bb in zip(A, b)]; m, n = len(A), len(A[0]) - 1; piv = []; r = 0
    for c in range(n):
        pr = next((i for i in range(r, m) if A[i][c] % p), None)
        if pr is None: continue
        A[r], A[pr] = A[pr], A[r]; inv = pow(A[r][c], -1, p); A[r] = [x * inv % p for x in A[r]]
        for i in range(m):
            if i != r and A[i][c] % p: f = A[i][c]; A[i] = [(x - f * y) % p for x, y in zip(A[i], A[r])]
        piv.append(c); r += 1
    if any(all(v == 0 for v in row[:n]) and row[n] for row in A): return False, None, None
    x = [0] * n
    for i, c in enumerate(piv): x[c] = A[i][n]
    return True, x, n - len(piv)
def crt(a2, a13): return [next(v for v in range(26) if v % 2 == x and v % 13 == y) for x, y in zip(a2, a13)]
def test(side, X, L, cribs):
    """side 'plain': unknown mask on plaintext letters; 'cipher': unknown mask on ciphertext letters."""
    pos = sorted(cribs); letters = sorted({(cribs[i] if side == "plain" else k4.CT[i]) for i in pos})
    LI = {c: j for j, c in enumerate(letters)}; nx = len(letters)
    rows, rhs = [], []
    for i in pos:
        row = [0] * (nx + L)
        if side == "plain": row[LI[cribs[i]]] = 1; row[nx + i % L] = 1; rhs.append(X.index(k4.CT[i]))
        else: row[LI[k4.CT[i]]] = 1; row[nx + i % L] = -1 % 26; rhs.append(X.index(cribs[i]))
        rows.append(row)
    # fix the gauge: x[letters[0]] = 0
    g = [0] * (nx + L); g[0] = 1; rows.append(g); rhs.append(0)
    ok2, s2, n2 = solve_mod(rows, rhs, 2); ok13, s13, n13 = solve_mod(rows, rhs, 13)
    if not (ok2 and ok13): return "inconsistent", None
    if n2 or n13: return "underdetermined", None
    sol = crt(s2, s13); x = sol[:nx]
    if len(set(x)) < nx: return "not one-to-one", None
    return "FITS", dict(zip(letters, x))
if __name__ == "__main__":
    CR = dict(zip(k4.CPOS.tolist(), k4.CPL))
    for side in ("plain", "cipher"):
        for an, X in (("AZ", k4.AZ), ("KA", k4.KA)):
            summary = {}
            for L in range(1, 27):
                st, sol = test(side, X, L, CR); summary.setdefault(st, []).append(L)
                if st == "FITS": print("   FIT:", side, an, L, sol)
            print(f"mask on {side:6s} side, numbering {an}: " + "; ".join(f"{k}: L={v}" for k, v in summary.items()))
    # planted check
    import random; random.seed(1); P = [random.choice(k4.AZ) for _ in range(97)]
    for s, w in k4.CRIBS.items(): P[s:s+len(w)] = w
    mask = list(k4.AZ); random.shuffle(mask); key = [3, 14, 7, 22, 9, 0, 17]
    ct = "".join(k4.AZ[(k4.AZ.index(mask[k4.AZ.index(p)]) + key[i % 7]) % 26] for i, p in enumerate(P))
    real = k4.CT; k4.CT = ct; print("planted (mask + period 7):", test("plain", k4.AZ, 7, CR)[0]); k4.CT = real
