"""One-slip tolerance for linear key models: for each model and each crib j, build the checksums on the other 23
cribs (column j zeroed). A hit means the model fits with crib j treated as an encryption/carving slip."""
import numpy as np, time, pickle, sys, k4, cokernel, rk_mixed, gap2
POS = k4.CPOS.tolist(); MIN = 8
def drop_models():
    out = []
    for name, M in cokernel.models().items():
        if not (name.startswith("periodic") or name.startswith("progressive")): continue
        for j in range(24):
            keep = [r for r in range(24) if r != j]; Mj = M[keep]
            N2, N13 = cokernel.left_null(Mj, 2), cokernel.left_null(Mj, 13)
            if min(len(N2), len(N13)) < MIN: continue
            Z2 = np.zeros((len(N2), 24)); Z2[:, keep] = N2; Z13 = np.zeros((len(N13), 24)); Z13[:, keep] = N13
            out.append((f"{name} slip@{POS[j]+1}", Z2, Z13, len(N13)))
    return out
if __name__ == "__main__":
    D = drop_models(); print("model x dropped-crib variants:", len(D), flush=True)
    import gap1
    gap1.P = D; gap1.names = [d[0] for d in D]
    gap1.B13 = np.vstack([d[2] for d in D]).astype(np.float32); gap1.B2 = np.vstack([d[1] for d in D]).astype(np.float32)
    gap1.seg13 = np.cumsum([0] + [len(d[2]) for d in D])[:-1]; gap1.seg2 = np.cumsum([0] + [len(d[1]) for d in D])[:-1]
    # planted check: periodic L=9 key with one corrupted crib
    k = np.array([[3,1,4,1,5,9,2,6,5][p % 9] for p in POS]); k[5] = (k[5] + 7) % 26
    ok = [D[i][0] for i in range(len(D)) if not (D[i][1] @ k % 2).any() and not (D[i][2] @ k % 13).any()]
    print("planted periodic L=9 with slip at crib 6 ->", [o for o in ok if o.startswith("periodic L=9")], flush=True)
    AL = rk_mixed.alph_list() + gap2.build()[0] + [("".join(k4.AZ), "AZ"), (k4.KA, "KA")]
    t = time.time(); h = gap1.scan(k4.CT, AL, "K4"); print(f"K4 one-slip: {len(h)} hits over {len(AL):,} alphabets ({time.time()-t:.0f}s)", flush=True)
    for x in h[:25]: print("  ", x, flush=True)
    pickle.dump(h, open("slip1_hits.pkl", "wb"))
