"""Near-miss version of the dictionary-alphabet autokey search: record the match count (out of all checkable cribs)
for every configuration; compare the top of the distribution with the same search on shuffled ciphertext."""
import numpy as np, time, sys, k4, rk_mixed
from collections import Counter
POS = k4.CPOS; PL = np.array([ord(c) - 65 for c in k4.CPL])
def scan(ct, AL):
    CTall = np.array([ord(c) - 65 for c in ct]); dist = Counter(); best = []
    for s in range(0, len(AL), 40000):
        chunk = AL[s:s + 40000]; A = np.array([[ord(ch) - 65 for ch in a] for a, _ in chunk]); INV = np.argsort(A, axis=1); n = len(A)
        AZ = np.tile(np.arange(26), (n, 1)); KAi = np.tile(np.array([k4.KA.index(chr(65 + x)) for x in range(26)]), (n, 1))
        for pn, Pi in (("AZ", AZ), ("KA", KAi), ("A", INV)):
            for cn, Ci in (("AZ", AZ), ("KA", KAi), ("A", INV)):
                if "A" not in (pn, cn): continue
                pc = np.take_along_axis(Pi, np.tile(PL, (n, 1)), 1); cc = np.take_along_axis(Ci, np.tile(CTall[POS], (n, 1)), 1)
                for kn, k in (("vig", (cc - pc) % 26), ("beau", (cc + pc) % 26)):
                    for vn, V in (("in-plain", Pi), ("in-cipher", Ci), ("AZ", AZ), ("KA", KAi)):
                        Vct = np.take_along_axis(V, np.tile(CTall, (n, 1)), 1)
                        for lag in range(1, 22):             # all 24 cribs checkable
                            m = (k == Vct[:, POS - lag]).sum(axis=1)
                            dist.update(m.tolist())
                            top = np.nonzero(m >= 13)[0]
                            for r in top: best.append((int(m[r]), chunk[r][1], pn, cn, kn, vn, lag))
    return dist, sorted(best, reverse=True)
if __name__ == "__main__":
    AL = rk_mixed.alph_list()
    rng = np.random.default_rng(5)
    for label, ct in (("K4", k4.CT), ("shuffled-1", "".join(rng.permutation(list(k4.CT)))), ("shuffled-2", "".join(rng.permutation(list(k4.CT))))):
        t = time.time(); dist, best = scan(ct, AL)
        tail = {m: c for m, c in sorted(dist.items()) if m >= 9}
        print(f"{label}: configs {sum(dist.values()):,}; matches>=9 tail {tail}; best {best[:3]} ({time.time()-t:.0f}s)", flush=True)
