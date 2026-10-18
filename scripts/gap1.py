import numpy as np, time, sys, k4, cokernel, rk_mixed
MIN = int(sys.argv[1]) if len(sys.argv) > 1 else 8
P = cokernel.prepared(MIN)
names = [p[0] for p in P]
B13 = np.vstack([p[2] for p in P]).astype(np.float32); B2 = np.vstack([p[1] for p in P]).astype(np.float32)
seg13 = np.cumsum([0] + [len(p[2]) for p in P])[:-1]; seg2 = np.cumsum([0] + [len(p[1]) for p in P])[:-1]
POS = k4.CPOS; PL = np.array([ord(c) - 65 for c in k4.CPL])
def scan(ct, AL, label):
    C = np.array([ord(c) - 65 for c in ct])[POS]; hits = []
    for s in range(0, len(AL), 40000):
        chunk = AL[s:s + 40000]; A = np.array([[ord(ch) - 65 for ch in a] for a, _ in chunk]); INV = np.argsort(A, axis=1); n = len(A)
        AZ = np.tile(np.arange(26), (n, 1)); KAi = np.tile(np.array([k4.KA.index(chr(65 + x)) for x in range(26)]), (n, 1))
        for pn, Pi in (("AZ", AZ), ("KA", KAi), ("A", INV)):
            for cn, Ci in (("AZ", AZ), ("KA", KAi), ("A", INV)):
                if "A" not in (pn, cn): continue
                pc = np.take_along_axis(Pi, np.tile(PL, (n, 1)), 1); cc = np.take_along_axis(Ci, np.tile(C, (n, 1)), 1)
                for kn, k in (("vig", (cc - pc) % 26), ("beau", (cc + pc) % 26)):
                    kf = k.astype(np.float32)
                    z13 = (np.mod(kf @ B13.T, 13) == 0); z2 = (np.mod(kf @ B2.T, 2) == 0)
                    ok = np.logical_and.reduceat(z13, seg13, axis=1) & np.logical_and.reduceat(z2, seg2, axis=1)
                    for r, m in zip(*np.nonzero(ok)): hits.append((chunk[r][1], pn, cn, kn, names[m]))
    return hits
if __name__ == "__main__":
    AL = rk_mixed.alph_list(); print(f"models: {len(P)} (>= {MIN} checks); alphabets: {len(AL):,}", flush=True)
    t = time.time(); h = scan(k4.CT, AL, "K4"); print(f"K4: {len(h)} hits in {time.time()-t:.0f}s", flush=True)
    for x in h[:20]: print("  ", x)
    rng = np.random.default_rng(1)
    for i in range(2):
        hn = scan("".join(rng.permutation(list(k4.CT))), AL, "null"); print(f"shuffled baseline {i}: {len(hn)} hits", flush=True)
    import pickle; pickle.dump(h, open(f"gap1_hits_{MIN}.pkl", "wb"))
