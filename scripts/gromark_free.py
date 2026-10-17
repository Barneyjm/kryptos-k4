"""Gromark with an UNKNOWN cipher alphabet (any permutation). For every primer, key digits at the 24 crib
positions force slot(C_i) = PA.index(P_i) + d_i (mod 26). Survive only if consistent as a permutation."""
import numpy as np, itertools, sys, time, k4
POS = k4.CPOS; PL = k4.CPL; CL = [k4.CT[i] for i in POS]
same = [(a, b) for a, b in itertools.combinations(range(24), 2) if CL[a] == CL[b]]
diff = [(a, b) for a, b in itertools.combinations(range(24), 2) if CL[a] != CL[b]]
SA, SB = np.array([p[0] for p in same]), np.array([p[1] for p in same])
DA, DB = np.array([p[0] for p in diff]), np.array([p[1] for p in diff])
def expand(primers, b, n, lags):
    """primers: (m, n) digits -> (m, 74) digits (only need up to index 73)."""
    m = len(primers); d = np.zeros((m, 74), np.int16); d[:, :n] = primers
    p, q = lags
    for i in range(n, 74): d[:, i] = (d[:, i - p] + d[:, i - q]) % b
    return d
def survivors(b, n, lags, PA, sign):
    total = b ** n; out = []
    pidx = np.array([PA.index(c) for c in PL])
    for start in range(0, total, 2_000_000):
        idx = np.arange(start, min(total, start + 2_000_000))
        prim = np.stack([(idx // b ** (n - 1 - k)) % b for k in range(n)], axis=1)
        d = expand(prim, b, n, lags)[:, POS]
        slot = (pidx[None, :] + sign * d) % 26
        ok = np.all(slot[:, SA] == slot[:, SB], axis=1) & np.all(slot[:, DA] != slot[:, DB], axis=1)
        out += [tuple(int(x) for x in prim[i]) for i in np.nonzero(ok)[0]]
    return out
if __name__ == "__main__":
    t = time.time(); results = []
    for b in range(3, 17):
        for n in range(2, 8):
            if b ** n > 3e7: continue
            for lags in sorted({(n - 1, n), (1, n), (1, 2)} if n > 2 else {(1, 2)}):
                for pa_name, PA in (("AZ", k4.AZ), ("KA", k4.KA)):
                    for sign in (1, -1):
                        s = survivors(b, n, lags, PA, sign)
                        results.append((b, n, lags, pa_name, sign, len(s), s[:5]))
                        if s: print(f"base {b} primer len {n} lags {lags} plain {pa_name} sign {sign:+d}: {len(s)} survivors e.g. {s[:3]}", flush=True)
    print(f"done in {time.time()-t:.0f}s; configs {len(results)}; total survivors {sum(r[5] for r in results)}")
