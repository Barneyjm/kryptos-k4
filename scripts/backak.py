"""Gap 5: backward autokey. key_i = value(C_{i+lag}) or value(P_{i+lag}), lag 1..96, all alphabets. Exact + near (>=18/24)."""
import numpy as np, time, k4, rk_mixed, gap2
from round3 import ALPHA
from collections import Counter
POS = k4.CPOS; PL = np.array([ord(c) - 65 for c in k4.CPL]); CTall = np.array([ord(c) - 65 for c in k4.CT]); CR = dict(zip(POS.tolist(), PL.tolist()))
def scan(AL, ct=CTall):
    hits = []; dist = Counter()
    for s in range(0, len(AL), 40000):
        chunk = AL[s:s + 40000]; A = np.array([[ord(ch) - 65 for ch in a] for a, _ in chunk]); INV = np.argsort(A, axis=1); n = len(A)
        AZ = np.tile(np.arange(26), (n, 1)); KAi = np.tile(np.array([k4.KA.index(chr(65 + x)) for x in range(26)]), (n, 1))
        for pn, Pi in (("AZ", AZ), ("KA", KAi), ("A", INV)):
            for cn, Ci in (("AZ", AZ), ("KA", KAi), ("A", INV)):
                if "A" not in (pn, cn) and chunk[0][1] not in ("AZ", "KA"): continue
                pc = np.take_along_axis(Pi, np.tile(PL, (n, 1)), 1); cc = np.take_along_axis(Ci, np.tile(ct[POS], (n, 1)), 1)
                for kn, k in (("vig", (cc - pc) % 26), ("beau", (cc + pc) % 26)):
                    for vn, V in (("in-plain", Pi), ("in-cipher", Ci), ("AZ", AZ), ("KA", KAi)):
                        Vct = np.take_along_axis(V, np.tile(ct, (n, 1)), 1)
                        for lag in range(1, 97):
                            idx = [j for j, p in enumerate(POS) if p + lag < 97]
                            if len(idx) < 12: continue
                            m = (k[:, idx] == Vct[:, POS[idx] + lag]).sum(axis=1); dist.update((m * 24 // len(idx)).tolist())
                            for r in np.nonzero(m >= len(idx) - (0 if len(idx) < 18 else 2))[0]: hits.append((int(m[r]), len(idx), "ctak-back", chunk[r][1], pn, cn, kn, vn, lag))
                            pidx = [j for j, p in enumerate(POS) if p + lag in CR]
                            if len(pidx) >= 8:
                                nxt = np.array([CR[POS[j] + lag] for j in pidx])
                                m2 = (k[:, pidx] == np.take_along_axis(V, np.tile(nxt, (n, 1)), 1)).sum(axis=1)
                                for r in np.nonzero(m2 == len(pidx))[0]: hits.append((int(m2[r]), len(pidx), "ptak-back", chunk[r][1], pn, cn, kn, vn, lag))
    return hits, dist
if __name__ == "__main__":
    AL = [(k4.AZ, "AZ"), (k4.KA, "KA")] + [(v, k) for k, v in ALPHA.items()]
    t = time.time(); h, d = scan(AL); print(f"std+clue alphabets: hits {len(h)} ({time.time()-t:.0f}s)", flush=True); [print("  ", x) for x in sorted(h, reverse=True)[:10]]
    AL = rk_mixed.alph_list(); t = time.time(); h, d = scan(AL); print(f"dictionary: hits {len(h)}; scaled-match tail {dict((k, v) for k, v in sorted(d.items()) if k >= 10)} ({time.time()-t:.0f}s)", flush=True)
    for x in sorted(h, reverse=True)[:10]: print("  ", x, flush=True)
