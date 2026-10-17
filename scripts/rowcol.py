"""Grid key: key_i = col[(i+off) % W] + row[(i+off) // W] (mod 26). Quagmire-style: CA(C) - PA(P) = key (vig)
or CA(C) + PA(P) = key (beau). Linear consistency over Z2 x Z13 with only the unknowns the cribs touch."""
import itertools, numpy as np, sys, multiprocessing as mp, k4
from linear import rank_consistent
from round3 import ALPHA
CR = list(zip(k4.CPOS.tolist(), k4.CPL))
def system(W, off):
    rows_used = sorted({(i + off) // W for i, _ in CR}); RI = {r: j for j, r in enumerate(rows_used)}
    M = []
    for i, _ in CR:
        v = [0] * (W + len(rows_used)); v[(i + off) % W] = 1; v[W + RI[(i + off) // W]] = 1; M.append(v)
    return M
SYS = {(W, off): system(W, off) for W in range(2, 31) for off in range(W)}
def test(args):
    pn, cn, Ws = args; PA, CA = ALPHA[pn], ALPHA[cn]; out = []
    for kn, f in (("vig", lambda c, p: (c - p) % 26), ("beau", lambda c, p: (c + p) % 26)):
        k = [f(CA.index(k4.CT[i]), PA.index(p)) for i, p in CR]
        for W in Ws:
            for off in range(W):
                M = SYS[(W, off)]
                a, r2 = rank_consistent(M, k, 2)
                if not a: continue
                b, r13 = rank_consistent(M, k, 13)
                slack = 24 - max(r2, r13)
                if b and slack >= 5: out.append((slack, pn, cn, kn, W, off))
    return out
if __name__ == "__main__":
    names = list(ALPHA)
    jobs = [(p, c, range(2, 31)) for p in ("AZ", "KA") for c in ("AZ", "KA")]
    jobs += [(p, c, [7, 14, 21]) for p in names for c in names]
    with mp.Pool(11) as pool: res = pool.map(test, jobs, chunksize=50)
    hits = sorted([h for r in res for h in r], reverse=True)
    print(f"alphabet pairs: {len(jobs):,}; grid keys consistent with >= 5 spare checks: {len(hits)}")
    for h in hits[:20]: print("  ", h)
