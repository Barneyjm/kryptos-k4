"""Autokey with keyword-mixed alphabets. key_i = value of C[i-lag] (ciphertext autokey) or P[i-lag] (plaintext autokey,
only where both letters are known), value measured in the plain, cipher, A-Z or KRYPTOS alphabet. Exact on all checkable cribs."""
import numpy as np, time, pickle, sys, k4, rk_mixed
POS = k4.CPOS; PL = np.array([ord(c) - 65 for c in k4.CPL]); CTall = np.array([ord(c) - 65 for c in k4.CT])
CR = dict(zip(POS.tolist(), PL.tolist()))
def scan(AL):
    hits = []
    for s in range(0, len(AL), 40000):
        chunk = AL[s:s + 40000]; A = np.array([[ord(ch) - 65 for ch in a] for a, _ in chunk]); INV = np.argsort(A, axis=1); n = len(A)
        AZ = np.tile(np.arange(26), (n, 1)); KAi = np.tile(np.array([k4.KA.index(chr(65 + x)) for x in range(26)]), (n, 1))
        for pn, Pi in (("AZ", AZ), ("KA", KAi), ("A", INV)):
            for cn, Ci in (("AZ", AZ), ("KA", KAi), ("A", INV)):
                if "A" not in (pn, cn): continue
                pc = np.take_along_axis(Pi, np.tile(PL, (n, 1)), 1); cc = np.take_along_axis(Ci, np.tile(CTall[POS], (n, 1)), 1)
                for kn, k in (("vig", (cc - pc) % 26), ("beau", (cc + pc) % 26)):
                    for vn, V in (("in-plain", Pi), ("in-cipher", Ci), ("AZ", AZ), ("KA", KAi)):
                        Vct = np.take_along_axis(V, np.tile(CTall, (n, 1)), 1)          # value of every ciphertext letter
                        for lag in range(1, 97):
                            idx = [j for j, p in enumerate(POS) if p - lag >= 0]
                            if len(idx) < 12: continue
                            ok = np.all(k[:, idx] == Vct[:, POS[idx] - lag], axis=1)
                            for r in np.nonzero(ok)[0]: hits.append(("ctak", chunk[r][1], pn, cn, kn, vn, lag))
                            pidx = [j for j, p in enumerate(POS) if p - lag in CR]
                            if len(pidx) >= 8:
                                prev = np.array([CR[POS[j] - lag] for j in pidx])
                                ok = np.all(k[:, pidx] == np.take_along_axis(V, np.tile(prev, (n, 1)), 1), axis=1)
                                for r in np.nonzero(ok)[0]: hits.append(("ptak", chunk[r][1], pn, cn, kn, vn, lag))
    return hits
if __name__ == "__main__":
    AL = rk_mixed.alph_list(); t = time.time(); h = scan(AL)
    print(f"autokey x {len(AL):,} dictionary alphabets: {len(h)} exact fits ({time.time()-t:.0f}s)"); [print("  ", x) for x in h[:20]]
