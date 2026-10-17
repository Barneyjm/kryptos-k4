import itertools, numpy as np, k4
from sympy import Matrix
P = dict(zip(k4.CPOS.tolist(), k4.CPL))
NUMS = {"AZ": k4.AZ, "KA": k4.KA, "AZ1": None}
def solve(n, X, align, ct, shift):
    blocks = []
    for s in range(align, 97 + shift, n):
        idx = [j - shift for j in range(s, s + n)]
        if all(j in P for j in idx):
            pv = [X.index(P[j]) for j in idx]; cv = [X.index(ct[j]) for j in idx]
            blocks.append((pv, cv))
    res = []
    if len(blocks) <= n: return len(blocks), res
    for direction in ("enc", "dec"):           # C = K P  or  P = K C
        for combo in itertools.combinations(range(len(blocks)), n):
            A = Matrix([blocks[c][0 if direction == "enc" else 1] for c in combo]).T
            B = Matrix([blocks[c][1 if direction == "enc" else 0] for c in combo]).T
            try: Ai = A.inv_mod(26)
            except Exception: continue
            K = (B * Ai).applyfunc(lambda x: x % 26)
            good = sum(all(((K * Matrix(b[0 if direction == "enc" else 1])).applyfunc(lambda x: x % 26)) == Matrix(b[1 if direction == "enc" else 0]) for _ in [0]) for b in blocks)
            if good >= len(blocks) - 0: res.append((direction, K.tolist())); break
            if good == len(blocks) - 1 and len(blocks) - n >= 2: res.append(("1-error", direction, K.tolist()))
    return len(blocks), res
for n in (2, 3, 4, 5):
    for an, X in (("AZ", k4.AZ), ("KA", k4.KA)):
        for shift in (0, 1):                       # 1 = count the '?' before K4 as part of the block stream
            for align in range(n):
                nb, r = solve(n, X, align, k4.CT, shift)
                if r: print("HIT", n, an, "shift", shift, "align", align, r)
                elif n == 4: print(f"n=4 {an} shift{shift} align{align}: {nb} known blocks, no consistent key")
print("done")
